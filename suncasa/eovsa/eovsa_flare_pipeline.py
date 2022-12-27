import os, sys, glob
import numpy as np
from astropy.time import Time
import time
from astropy.io import fits
import datetime as dt
import subprocess
from casatasks import *
from casatools import image, ms, msmetadata, table
from suncasa.utils import qlookplot
from suncasa.utils import helioimage2fits as hf
import timeit
from sunpy.time import parse_time
import matplotlib

matplotlib.use('Agg')

start = timeit.default_timer()
ia = image()
ms = ms()
msmd = msmetadata()
tb = table()


class FlareSelfCalib():
    def __init__(self, vis=None, workpath='./', logfile=None):
        ##========================= initial setups =================================
        self.vis = vis
        if not logfile:
            logfile = workpath + "selfcal_{0:s}.log".format(Time.now().isot[:-4].replace(':', '').replace('-', ''))
        self.logfile = logfile
        self.identify_data_gap = True
        self.slfcal_spws = [3, 5, 8, 10, 15, 20, 24, 30, 35, 40, 42, 45]
        # self.slfcal_spws = [15, 24]  ## spw used for flare finding and final imaging
        self.maximum_spw = max(self.slfcal_spws)
        self.minimum_spw = min(self.slfcal_spws)
        self.slfcal_spwstr = ','.join(
            [str(s) for s in self.slfcal_spws])  ## string format of self.spws for, e.g., split
        self.selfcal_spw = '0~' + str(self.maximum_spw)  ## actual spw range used for self-calibration

        ##===================== final imaging parameters ===========================
        self.total_duration = 480.  ### seconds to image; time will be centred on detected flare peak
        self.final_image_cadence = 10  ### Cadence of the images after final self-calibration
        self.final_image_int = 10  ### Integration time of final images
        self.min_restoring_beam = 6.  ### minimum size of the restoring beam in arcsec
        self.beam_1GHz = '89.7arcsec'  ### beam size at 1 GHz (and scales with 1/nu)
        self.cell_size = 2.  ### cell size of the final images
        self.final_imsize = 512  ### number of pixels for the final images

        # ============ declaring the working directories ============
        ### remember / is necessary in all the folder names
        # self.workpath = '/data1/bchen/flare_pipeline/tmp/'  ## / is needed at the end of all paths
        self.workpath = workpath
        self.slfcaldir = self.workpath + 'slfcal/'  # place to put all selfcalibration products
        self.imagedir = self.slfcaldir + 'images/'  # place to put all selfcalibration images
        self.caltbdir = self.slfcaldir + 'caltables/'  # place to put calibration tables
        self.maskdir = self.slfcaldir + 'masks/'  # place to put masks
        self.imagedir_slfcaled = self.slfcaldir + 'images_slfcal'  # place to put all selfcalibrated images

        if os.path.exists(self.slfcaldir):
            print("{0:s} already exists. Re-initialize it.".format(self.slfcaldir))
            os.system('rm -rf {0:s}'.format(self.slfcaldir))

        os.makedirs(self.slfcaldir)
        os.makedirs(self.imagedir)
        os.makedirs(self.caltbdir)
        os.makedirs(self.maskdir)
        os.makedirs(self.imagedir_slfcaled)

        # ============ selfcal parameters ===============
        self.refantenna = '0'  ### reference antenna
        self.calc_cell = True  ### If set to False use the value in beam given below
        # TODO what is this?
        self.cell = [10]  ### size needs to be same as the number of spws
        self.calc_imsize = True  ### is False uses the value given below
        self.max_frac_freq_avg = 0.5  ### Frequency averaging factor
        self.maxiter = 10  ### maximum selfcal iteration rounds
        self.uvlim = 25  ### maximum uv distance in kilo-lambda
        self.avg_spw_max = 5  ### maximum number of spws for averaging
        self.flag_antennas = ''  ###anything except 13~15. Those antennas are always flagged.

        # ============ the following could be populated using the flare data base =======
        self.phasecenter = ''  ### phasecenter of the flare in RA and DEC
        self.flare_time_available = False  ### Flag on whether the flare time (peak and duration) is available
        self.flare_loc_available = False  ### Flag on whether the flare location is available

    @property
    def vis(self):
        """
        Getting the input visibility for self-calibration
        """
        return self._vis

    @vis.setter
    def vis(self, value):
        """
        Setting the input visibility for self-calibration
        Parameters
        ----------
        value: full path for the input visibility data
        -------
        """
        if value:
            if not os.path.exists(value):
                raise ValueError('input visibility {0:s} does not exist! Please supply a valid one.'.format(value))
        else:
            print("No visibility provided. Please supply a valid one by FlareSelfCalib.vis = 'your visibility'")
        self._vis = value

    def vis_info(self):
        if os.path.exists(self.vis):
            listobs(self.vis, listfile=self.vis + '.lisobs')
            print('visibility info stored in {0:s}'.format(self.vis + '.listobs'))
        else:
            print('input visibility {0:s} do not exist! Abort...'.format(self.vis))

    @staticmethod
    def get_img_center_heliocoords(images):
        """
        Provide a set of images in helioprojective coordinates (at different frequencies), find the peak location
        Parameters
        ----------
        images: list of fits image files

        Returns
        -------
        xycen: solar x and y coordinates, in arcsec
        """
        num_img = len(images)
        x = np.zeros(num_img)
        y = np.zeros(num_img)
        for j, img in enumerate(images):
            head = fits.getheader(img)
            data = fits.getdata(img)[0, 0, :, :]
            x0 = head['CRVAL1']
            y0 = head['CRVAL2']
            cell_x = head['CDELT1']
            cell_y = head['CDELT2']
            pos = np.where(data == np.nanmax(data))
            xmax = pos[1][0]
            ymax = pos[0][0]
            xcen = head['CRPIX1']
            ycen = head['CRPIX2']
            dx = xmax - xcen
            dy = ymax - ycen
            dx_asec = dx * cell_x + x0
            dy_asec = dy * cell_y + y0
            x[j] = dx_asec
            y[j] = dy_asec
        xycen = [np.median(x), np.median(y)]
        return xycen

    @staticmethod
    def find_sidelobe_level(image):
        head = fits.getheader(image)
        imsize = head['naxis1']
        cellsize = float(abs(head['CDELT1'])) * 3600

        psf_image = image[:-5] + ".psf"
        ia.open(psf_image)
        psf_data = ia.getchunk()[:, :, 0, 0]
        ia.close()

        x0 = head['CRPIX2'] - 1
        y0 = head['CRPIX1'] - 1

        posang = head['BPA'] * np.pi / 180  # in degrees --> radian
        bmaj = head['BMAJ'] * 3600 / cellsize  # in degrees --> pix Full maj axis extent
        bmin = head['BMIN'] * 3600 / cellsize  # ''
        R = np.array([[np.cos(posang), -np.sin(posang)], [np.sin(posang), np.cos(posang)]])
        a = bmaj / 2.  # Chosing nsig region for calculating flux density
        b = bmin / 2.

        low_limit_x = int(x0 - a)
        low_limit_y = int(y0 - a)
        upper_limit_x = int(x0 + a)
        upper_limit_y = int(y0 + a)
        psf_region_values = []

        for k in range(low_limit_y, upper_limit_y):
            for j in range(low_limit_x, upper_limit_x):
                x_v = np.matrix([[k - x0], [j - y0]])
                R = np.matrix(R)
                xp_v = R * x_v
                xp_v = np.array(xp_v.transpose()).reshape(2)
                if xp_v[0] ** 2 / a ** 2 + xp_v[1] ** 2 / b ** 2 < 1:
                    psf_region_values.append(psf_data[k, j])
                    psf_data[k, j] = np.nan
        max_val = np.nanmax(psf_data)
        min_lev = np.nanmin(np.array(psf_region_values))
        return max_val, min_lev

    @staticmethod
    def check_shift(image, shift, cell):
        if type(cell) == str:
            # input cell is in the format of '5.0000arcsec'
            cell = float(cell.split('a')[0])
        header = imhead(image)
        major = header['restoringbeam']['major']['value']
        minor = header['restoringbeam']['minor']['value']
        pa = header['restoringbeam']['positionangle']['value'] * np.pi / 180
        v2 = np.array([np.sin(pa), np.cos(pa)])
        s = np.dot(v2, shift) / (np.sqrt(v2[0] ** 2 + v2[1] ** 2) * np.sqrt(shift[0] ** 2 + shift[1] ** 2))
        if shift[0] ** 2 + shift[1] ** 2 < 1e-3:
            return True
        theta = np.arccos(s)
        abs_shift = np.sqrt(shift[0] ** 2 + shift[1] ** 2)
        major_shift = abs(abs_shift * np.cos(theta))
        minor_shift = abs(abs_shift * np.cos(theta))
        print(major_shift, minor_shift, major / cell, minor / cell)
        if major_shift > 0.75 * major / cell:
            return False
        if minor_shift > 0.75 * minor / cell:
            return False
        return True

    @staticmethod
    def grow_mask(image, mask, thres):
        image_data = fits.getdata(image)[0, 0, :, :]
        ia.open(mask)
        mask_data = ia.getchunk()[:, :, 0, 0].T
        ia.close()

        shape = np.shape(mask_data)
        rows = shape[0]
        cols = shape[1]

        pos = np.where(mask_data == 1)
        if len(pos) == 0:
            print("Mask blank. First run gen_mask. Exiting")
            return False

        xpos1 = pos[1]
        ypos1 = pos[0]

        while np.size(xpos1) != 0:
            for x, y in zip(xpos1, ypos1):
                j = x
                i = y + 1
                while i < rows:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        i += 1
                i = y
                j = x + 1
                while j < cols:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        j += 1

            del xpos1
            del ypos1
            del pos
            pos = np.where(mask_data == 2)
            xpos1 = pos[1]
            ypos1 = pos[0]
            mask_data[pos] = 1

        pos = np.where(mask_data == 1)
        xpos1 = pos[1]
        ypos1 = pos[0]

        while np.size(xpos1) != 0:
            for x, y in zip(xpos1, ypos1):
                j = x
                i = y - 1
                while i < rows:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        i -= 1
                i = y
                j = x - 1
                while j < cols:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        j -= 1

            del xpos1
            del ypos1
            del pos
            pos = np.where(mask_data == 2)
            xpos1 = pos[1]
            ypos1 = pos[0]
            mask_data[pos] = 1

        pos = np.where(mask_data == 1)
        xpos1 = pos[1]
        ypos1 = pos[0]

        while np.size(xpos1) != 0:
            for x, y in zip(xpos1, ypos1):
                j = x
                i = y - 1
                while i < rows:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        i -= 1
                i = y
                j = x - 1
                while j < cols:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        j += 1

            del xpos1
            del ypos1
            del pos
            pos = np.where(mask_data == 2)
            xpos1 = pos[1]
            ypos1 = pos[0]
            mask_data[pos] = 1

        pos = np.where(mask_data == 1)
        xpos1 = pos[1]
        ypos1 = pos[0]

        while np.size(xpos1) != 0:
            for x, y in zip(xpos1, ypos1):
                j = x
                i = y + 1
                while i < rows:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        i += 1
                i = y
                j = x - 1
                while j < cols:
                    if mask_data[i, j] == 1 or mask_data[i, j] == 2 or image_data[i, j] < thres:
                        break
                    else:
                        mask_data[i, j] = 2
                        j -= 1

            del xpos1
            del ypos1
            del pos
            pos = np.where(mask_data == 2)
            xpos1 = pos[1]
            ypos1 = pos[0]
            mask_data[pos] = 1

        pos = np.where(mask_data == 1)
        xpos1 = pos[1]
        ypos1 = pos[0]

        del mask_data
        ia.open(mask)
        mask_data = ia.getchunk()

        for x, y in zip(pos[1], pos[0]):
            mask_data[x, y, 0, 0] = 1
        ia.putchunk(mask_data)
        ia.close()

        return True

    @staticmethod
    def get_spw_num(visibility):
        msmd.open(visibility)
        nspws = msmd.nspw()
        msmd.done()
        return nspws

    @staticmethod
    def calc_cellsize(visibility):
        nspws = FlareSelfCalib.get_spw_num(visibility)
        beams = []
        ms.open(visibility)
        for i in range(nspws):
            ms.selectinit(datadescid=0, reset=True)
            ms.selectinit(datadescid=i)
            data = ms.getdata(['u', 'v', 'axis_info'], ifraxis=True)
            fghz = data['axis_info']['freq_axis']['chan_freq'][
                       0, 0] / 1e9  ### taking the first frequency. EOVSA spws bandwidths are small. Hence this is fine
            uvdist = np.sqrt(data['u'] ** 2 + data['v'] ** 2)
            max_uvlambda = np.max(uvdist / 299792458.0 * fghz * 1e9)
            beam_val = 1.0 / max_uvlambda * 180 / 3.14159 * 3600 / 4.0  ### we take that the beam is resolved by 4 pixels by default
            beams.append(beam_val)
        ms.close()
        return beams

    @staticmethod
    def get_ref_freqlist(visibility):
        tb.open(visibility + '/SPECTRAL_WINDOW')
        reffreqs = tb.getcol('REF_FREQUENCY')
        bdwds = tb.getcol('TOTAL_BANDWIDTH')
        cfreqs = reffreqs + bdwds / 2.
        tb.close()
        return cfreqs

    @staticmethod
    def read_bandpass(bptable, nant=16):
        bptb = tb.open(bptable, nomodify=True)  # Open the bandpass table
        bpass = tb.getcol('CPARAM')  # Extract the CPARAM=GAIN column from the bandpass table
        flag = tb.getcol('FLAG')  # Extract FLAG column from bandpass table
        snr = tb.getcol('SNR');
        dims = bpass.shape
        npol = dims[0]
        nchan = dims[1]
        ntime = int(dims[2] / nant)
        flag = np.reshape(flag, [npol, nchan, ntime, nant])
        bpass = np.reshape(bpass, [npol, nchan, ntime, nant])
        success = True
        tb.close()
        return bpass, flag, success

    @staticmethod
    def combine_groups(group, pos):
        groups_com = []
        deleted_groups = []

        for g1 in range(len(group)):
            if g1 in deleted_groups:
                continue
            for g2 in range(g1 + 1, len(group)):
                combine = False
                for mem1 in group[g1]:
                    for mem2 in group[g2]:
                        x1 = pos[1][mem1]
                        x2 = pos[1][mem2]
                        y1 = pos[0][mem1]
                        y2 = pos[0][mem2]
                        # if x1==x2 and y1==y2:
                        if np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) < 3:
                            combine = True
                            break
                    if combine == True:
                        break
                if combine == True:
                    for mem2 in group[g2]:
                        group[g1].append(mem2)
                    deleted_groups.append(g2)
                    continue
            temp = []
            for mem in group[g1]:
                temp.append(mem)
            groups_com.append(temp)
        return groups_com

    @staticmethod
    def gen_fof_groups(data3, thres):
        pos = np.where(data3 > thres)
        y = pos[0]
        x = pos[1]
        size = len(x)

        groups = [[0]]
        for i in range(1, size):
            x0 = pos[1][i]
            y0 = pos[0][i]
            neighbour = False
            for g in groups:
                for mem in g:
                    x1 = pos[1][mem]
                    y1 = pos[0][mem]

                    if np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2) < 3:
                        g.append(i)
                        neighbour = True
                        break
                if neighbour == True:
                    break
            if neighbour == False:
                groups.append([i])
        # print("Combining groups")
        # print(groups)
        group_com = FlareSelfCalib.combine_groups(groups, pos)
        # print("Groups combined")
        final_groups_x = []
        final_groups_y = []
        for g1 in group_com:
            # print (len(g1))
            if len(g1) > 3:
                for mem in g1:
                    final_groups_x.append(pos[1][mem])
                    final_groups_y.append(pos[0][mem])
        return [final_groups_y, final_groups_x]

    def gen_mask(self, image1, image2, mask1, mask2, threshold, imsize, s,
                 make_shifted_mask=False, grow_threshold=0.5):  ### threshold corresponds to image1
        '''
        We will allow for small shifts and small change of size ere
        '''
        data2 = fits.getdata(image2)[0, 0, :, :]
        sidelobe_level, min_lev = FlareSelfCalib.find_sidelobe_level(image2)
        if image1 != '':
            data1 = fits.getdata(image1)[0, 0, :, :]
            head1 = fits.getheader(image1)
            data2 = fits.getdata(image2)[0, 0, :, :]
            pos1 = np.where(data2 == np.nanmax(data2))
            head2 = fits.getheader(image2)
            # print (mask1)
            ia.open(mask1)
            mask1_data = ia.getchunk()[:, :, 0, 0].T
            ia.close()
            pos = np.where(mask1_data == 1)
            # print (pos)
            if np.min(pos[0]) == 0 and np.min(pos[1]) == 0 and np.max(pos[1]) == imsize - 1 and np.max(
                    pos[1]) == imsize - 1:
                del pos
                pos = np.where(data1 > threshold)
            if np.size(pos) == 0:
                return False
            x_min = np.min(pos[1])
            x_max = np.max(pos[1])
            y_min = np.min(pos[0])
            y_max = np.max(pos[0])
            box_x = x_max - x_min
            box_y = y_max - y_min
            # print(x_min,x_max,y_min,y_max)
            self.logf.write("spw=" + str(s).zfill(2) + " Previous image supplied\n")
            self.logf.write("spw=" + str(s).zfill(2) + " Mask boundaries in previous image:" + \
                            str(x_min) + "," + str(x_max) + "," + str(y_min) + "," + str(y_max))

            cell1 = head1['CDELT1']
            cell2 = head2['CDELT1']
            # print (cell1,cell2)
            ref_ra1 = head1['CRVAL1']
            ref_dec1 = head1['CRVAL2']
            ref_ra2 = head2['CRVAL1']
            ref_dec2 = head2['CRVAL2']
            ref_xpix1 = head1['CRPIX1']
            ref_ypix1 = head1['CRPIX2']
            ref_xpix2 = head2['CRPIX1']
            ref_ypix2 = head2['CRPIX2']

            ra_min = (x_min - ref_xpix1) * cell1 + ref_ra1
            ra_max = (x_max - ref_xpix1) * cell1 + ref_ra1
            dec_min = (y_min - ref_ypix1) * cell1 + ref_dec1
            dec_max = (y_max - ref_ypix1) * cell1 + ref_dec1

            x2_min = max(int((ra_min - ref_ra2) / cell2 + ref_xpix2), 0)
            x2_max = min(int((ra_max - ref_ra2) / cell2 + ref_xpix2), imsize - 1)
            y2_min = max(int((dec_min - ref_dec2) / cell2 + ref_ypix2), 0)
            y2_max = min(int((dec_max - ref_dec2) / cell2 + ref_ypix2), imsize - 1)
            # print (x2_min,x2_max,y2_min,y2_max)
            if x2_min > x2_max or y2_min > y2_max:  #### I do not know what to do here.
                #### Using the full image. Hopefull
                x2_min = 0
                x2_max = imsize - 1
                y2_min = 0
                y2_max = imsize - 1
            self.logf.write("spw=" + str(s).zfill(2) + " Somehow lower boundary is greater than " + \
                            "the outer boundary. Masing the full image\n")
            data3 = data2[y2_min:y2_max, x2_min:x2_max]
            pos2 = np.where(data3 == np.nanmax(data3))
            if pos2[0][0] + y2_min != pos1[0][0] or pos2[1][0] + x2_min != pos1[1][0]:
                self.logf.write("spw=" + str(s).zfill(2) + "Peak pixel outside the mask." \
                                                           " Check for position shifts\n")
                # filename.write("Trying to phase up at the location of the previous mask.\n")
                self.logf.write("spw=" + str(s).zfill(2) + "Difference in x:" + \
                                str(int(pos1[1][0] - pos2[1][0])) + "pixels\n")
                self.logf.write("spw=" + str(s).zfill(2) + "Difference in y:" + \
                                str(int((pos1[0][0] - pos2[0][0]))) + "pixels\n")
                ##### for testing #####
                if make_shifted_mask == False:
                    self.logf.write("spw=" + str(s).zfill(2) + "Mske_shifted_mask is False. Will " + \
                                    "mask the entire region within the region of previous image\n")
                    mask_y = np.arange(y2_min, y2_max + 1)
                    mask_x = np.arange(x2_min, x2_max + 1)
                    ia.open(mask2)
                    data = ia.getchunk()
                    data *= 0
                    # for i,j in zip(mask_y,mask_x):
                    for j in mask_x:
                        for i in mask_y:
                            data[j, i, 0, 0] = 1
                    ia.putchunk(data)
                    ia.close()
                    del data3
                    return False

        else:
            y2_min = 0
            y2_max = imsize
            x2_min = 0
            x2_max = imsize
        data3 = data2[y2_min:y2_max, x2_min:x2_max]
        # print (y2_min,y2_max,x2_min,x2_max)
        self.logf.write("spw=" + str(s).zfill(2) + " Mask boundaries in previous image:" + \
                        "{}\n".format(x2_min, x2_max, y2_min, y2_max))
        # print ("shift successfull")
        max_data3 = np.nanmax(data3)
        chosen_thres = sidelobe_level + 0.01  ### some more caution
        mask_y, mask_x = FlareSelfCalib.gen_fof_groups(data3, chosen_thres * max_data3)
        # print (mask_y, mask_x)
        if len(mask_x) == 0:
            self.logf.write("spw=" + str(s).zfill(2) + " Mask not found. Masking the region" + \
                            " covered by the previous image\n")
            mask_y = np.arange(y2_min, y2_max + 1)
            mask_x = np.arange(x2_min, x2_max + 1)
            ia.open(mask2)
            data = ia.getchunk()
            data *= 0
            # for i,j in zip(mask_y,mask_x):
            for j in mask_x:
                for i in mask_y:
                    data[j, i, 0, 0] = 1
            ia.putchunk(data)
            ia.close()
            return False
        for i in range(len(mask_y)):
            mask_y[i] += y2_min
            mask_x[i] += x2_min
        ia.open(mask2)
        data = ia.getchunk()
        data *= 0
        for i, j in zip(mask_y, mask_x):
            data[j, i, 0, 0] = 1
        # print (i,j)

        ia.putchunk(data)
        ia.close()
        self.logf.write("spw=" + str(s).zfill(2) + " Mask generated. Growing the mask" + \
                        " with threshold {}\n".format(grow_threshold))
        success = FlareSelfCalib.grow_mask(image2, mask2,
                                           max(grow_threshold, min_lev - 0.1) * max_data3)  ### again some safety here.
        self.logf.write("spw=" + str(s).zfill(2) + " Mask generated\n")
        return True

    def confirm_maximum_pixel(self, imagename, mask, spwran, msname, uvrange,
                              imsize, cell, s):
        try:
            first_freq = int(spwran.split('~')[0])
            last_freq = int(spwran.split('~')[1])
        except IndexError:
            first_freq = int(spwran)
            last_freq = first_freq
        # print (first_freq, last_freq)

        self.logf.write("spw=" + str(s).zfill(2) + " confirm_maximum_pixel called with {}\n".format(spwran))
        maxpix = imstat(imagename, mask="mask(\'" + mask + "\')")["maxpos"]

        first_freq = max(self.minimum_spw, first_freq - 1)
        last_freq = min(last_freq + 1, self.maximum_spw)

        spwran = str(first_freq) + "~" + str(last_freq)
        self.logf.write("spw=" + str(s).zfill(2) + " Checking against " + str(first_freq) + "~" + str(last_freq) + "\n")

        tclean(vis=msname, imagename="check_maxpix_robust", uvrange=uvrange,
               spw=spwran, imsize=imsize, cell=cell, niter=1, gain=0.05,
               cyclefactor=10, weighting='briggs', robust=0.0, savemodel='none',
               pbcor=False, pblimit=0.001, phasecenter=self.phasecenter, stokes='XX')

        maxpix1 = imstat("check_maxpix_robust.image", mask="mask(\'" + mask + "\')")["maxpos"]

        shift = np.array([maxpix[0] - maxpix1[0], maxpix[1] - maxpix1[1]])
        self.logf.write("spw=" + str(s).zfill(2) + " Shift is {}\n".format(shift))
        cellsize = float(cell.split('arcsec')[0])
        success = FlareSelfCalib.check_shift("check_maxpix_robust.image", shift, cellsize)
        self.logf.write("spw=" + str(s).zfill(2) + " Shift ok? {}\n".format(success))
        os.system("rm -rf check_maxpix_robust.*")
        return success

    @staticmethod
    def restore_previous_condition(imagename):
        temp = imagename.split('_')
        prefix = '_'.join(temp[:-1])
        iteration_num = int(temp[-1])
        psf = glob.glob(prefix + "*.psf")
        num_psf = len(psf)
        for i in range(iteration_num + 1, num_psf):
            image = prefix + '_' + str(i).zfill(2)
            for str1 in ['.image', '.mask', '.residual', '.psf',
                         '.model', '.sumwt', '.pb', '.fits']:
                os.system("rm -rf " + image + str1)
            os.system("rm -rf " + image + ".gcal")
        return

    @staticmethod
    def flag_data_gap(visibility, sp):
        # TODO would it be easier to just use CASA's flagdata task? Or is this way more efficient?
        ms.open(visibility)
        ms.selectinit(datadescid=0, reset=True)
        ms.selectinit(datadescid=sp)
        data = ms.getdata('amplitude')['amplitude']
        ms.close()
        ms.open(visibility, nomodify=False)
        ms.selectinit(datadescid=0, reset=True)
        ms.selectinit(datadescid=sp)
        flag = ms.getdata('flag')
        pos = np.where(data < 1e-3)
        flag['flag'][pos] = True
        ms.putdata(flag)
        ms.close()
        return

    def get_img_stat(self, imagename):
        ia.open(imagename + ".image")
        data = ia.getchunk()[:, :, 0, 0]
        (ny, nx) = data.shape
        max_pix = np.nanmax(data)
        min_pix = np.nanmin(data)
        data[int(ny * 0.25):int(ny * 0.75), int(nx * 0.5):int(nx * 0.75)] = np.nan
        rms = np.nanstd(data)
        ia.close()
        return max_pix, min_pix, rms

    def flare_finder(self):
        """
        Provide input visibility, find out the flare peak time, flare duration, and
        suitable time ranges for performing self-calibration
        Parameters
        ----------

        Returns
        -------

        """
        ###TODO: currently it iterates over ALL the supplied spw list to find the flare peak,
        ### which is likely an overkill.
        ### Also self-calbration time range needs to be a shorter time around the peak.
        ### Now it is as long as 1-min, which may include significant flare evolution
        flare_times = []
        found_flares = []
        flare_peak = []
        start_time = []
        end_time = []
        start_datetime = []
        end_datetime = []

        if os.path.isfile(self.vis[:-3] + ".listobs") == False:
            listobs(self.vis, listfile=self.vis[:-3] + ".listobs")

        t = \
            subprocess.check_output("grep Observed " + self.vis[:-3] + ".listobs", shell=True).decode(
                'utf8').strip().split(
                ' ')[4]
        t = parse_time(t.replace('/', ' ')).value
        hmy = t.split('T')[1].split(':')
        hour = int(hmy[0])
        minute = int(hmy[1])
        second = int(float(hmy[2]))

        # ymd = starttime.split(' ')[0].split('-')
        # year = int(ymd[0])
        # month = int(ymd[1])
        # day = int(ymd[2])

        # start = dt.datetime(year, month, day, hour, minute, second)

        peak = 0
        peak_time = 0

        ms.open(self.vis)
        # for sp in range(self.minimum_spw, self.maximum_spw + 1):
        for sp in self.slfcal_spws:
            print('processing spw ', sp)
            ms.selectinit(datadescid=0, reset=True)
            ms.selectinit(datadescid=sp)
            data = ms.getdata(['amplitude', 'axis_info'], ifraxis=True)
            flag = ms.getdata('flag', ifraxis=True)['flag'][0, :, :, :]
            ms_startmjd = data['axis_info']['time_axis']['MJDseconds'][0]
            # print('startmjd: ', startmjd)
            if sp == self.minimum_spw:
                start = Time(ms_startmjd / 3600. / 24., format='mjd').datetime
            ms_endmjd = data['axis_info']['time_axis']['MJDseconds'][-1]
            tot_data = data['amplitude'][0, :, :, :]
            pos = np.where(flag == True)
            tot_data[pos] = np.nan
            #### first taking median over time
            median = np.nanmedian(tot_data, axis=2)
            median = np.expand_dims(median, axis=2)
            med_subtracted_power = (tot_data - median) / median
            power = np.nanmedian(med_subtracted_power, axis=(0, 1))
            # power=data['amplitude'][0,3,0,:]
            # print (data['axis_info']['time_axis'].keys())
            median_power = np.nanmedian(power)
            smoothed = np.convolve(power, np.ones(5), mode='same') * 1. / 5
            peak_pow = np.nanmax(smoothed)
            mad = np.nanmedian(abs(power - median_power))
            pos = np.where(np.isnan(smoothed) == True)
            smoothed[pos] = 0.0
            power[pos] = 0.0
            thres = 5.0
            y = (smoothed - median_power) / mad
            pos = np.argmax(y)
            tot_times = np.size(y)
            if y[pos] > peak:
                peak = y[pos]
                peak_time = data['axis_info']['time_axis']['MJDseconds'][pos]
            del pos
            pos = np.where(y > thres)[0]
            duration = 60
            if len(pos) == 0:
                found_flares.append(False)
            else:
                found_flares.append(True)
                temp = np.convolve(pos, np.array([1, -1]), mode='same')
                #### Couting the number of continuous zeros in temp. This will give the event duration
                counts = []
                count = 1
                for i in temp:
                    if i == 1:
                        count += 1
                    else:
                        if count != 1:
                            counts.append(count)
                            count = 0
                counts.append(count)

                counts.sort()
                duration = counts[-1]
                print('flare duration is {0:d}s'.format(duration))

                #### max integration time=60s#####
                if duration / 3 >= 60:
                    duration = 60
                elif duration / 3 > 20 and duration / 3 < 60:
                    duration = 40
                else:
                    duration = 20

            peak_pos = np.argmax(smoothed - median_power)
            flare_peak.append(smoothed[peak_pos] - median_power)
            # flare_peak.append(np.max(smoothed-median_power))
            ### assume that the duration is maximum when the peak is highest
            # peak_pos=np.where(abs(smoothed-flare_peak[-1]-median_power)<0.1)[0][0]

            ### check if peak_pos can be kept at middle of integration window.
            start_fine = True
            end_fine = True

            try:
                if smoothed[int(peak_pos - duration / 2)] <= median_power + thres * mad:
                    start_fine = False
            except IndexError:
                start_fine = False
            try:
                if smoothed[int(peak_pos + duration / 2)] <= median_power + thres * mad:
                    end_fine = False
            except IndexError:
                end_fine = False

            if (start_fine == True and end_fine == True) or (start_fine == False and end_fine == False):
                start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(max(0, peak_pos - duration / 2))])
                end_time.append(
                    data['axis_info']['time_axis']['MJDseconds'][int(min(peak_pos + duration / 2, np.size(power) - 1))])
            elif start_fine == False:
                diff = peak_pos - duration / 2 + 1
                end_found = False
                while end_found == False:
                    try:
                        if smoothed[int(diff)] > median_power + thres * mad:
                            end_found = True
                        else:
                            diff += 1
                    except IndexError:
                        diff += 1

                distance = peak_pos - diff + 1
                distance_left = duration - distance
                start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(diff)])
                end_time.append(data['axis_info']['time_axis']['MJDseconds'][
                                    int(min(peak_pos + distance_left, np.size(power) - 1))])

            else:
                diff = peak_pos + duration / 2 - 1
                end_found = False
                while end_found == False:
                    try:
                        if smoothed[int(diff)] > median_power + thres * mad:
                            end_found = True
                        else:
                            diff -= 1
                    except IndexError:
                        diff -= 1
                distance = diff - peak_pos + 1
                distance_left = duration - distance
                start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(max(0, peak_pos - distance_left))])
                end_time.append(data['axis_info']['time_axis']['MJDseconds'][int(diff)])
            time_diff = start_time[-1] - ms_startmjd
            startstr = (start + dt.timedelta(seconds=time_diff)).strftime('%Y/%m/%d/%H:%M:%S')
            start_datetime.append(start + dt.timedelta(seconds=time_diff))
            time_diff = end_time[-1] - ms_startmjd
            endstr = (start + dt.timedelta(seconds=time_diff)).strftime('%Y/%m/%d/%H:%M:%S')
            flare_times.append(startstr + '~' + endstr)
            print(sp, flare_times[-1])
            end_datetime.append(start + dt.timedelta(seconds=time_diff))
            del tot_data, median, med_subtracted_power, power, smoothed
        ms.close()
        # flare_peak_time = start + dt.timedelta(seconds=peak_time - startmjd)
        self.flare_time_available = True
        self.flare_times = flare_times
        self.found_flares = found_flares
        self.flare_start_datetime = start_datetime
        self.flare_end_datetime = end_datetime
        self.flare_peak_mjd = peak_time  # in mjd seconds
        self.flare_peak_datetime = Time(peak_time / 3600. / 24., format='mjd').datetime
        self.ms_startmjd = ms_startmjd
        self.ms_endmjd = ms_endmjd

    def find_previous_image(self, spw):
        imgprefix = self.imagedir + "selfcal"
        imagename = imgprefix + '_spw_*.image'
        images = glob.glob(imagename)
        num_img = len(images)
        if num_img == 0:
            return False, ''
        s = np.zeros(num_img)
        for i in range(num_img):
            chunks = images[i].split('_')
            for n, c in enumerate(chunks):
                if c == 'spw':
                    break

            if '~' in chunks[n + 1]:
                temp = chunks[n + 1].split('~')
                freq1 = int(temp[0])
                freq2 = int(temp[1])
                s[i] = 0.5 * (freq1 + freq2)
            else:
                s[i] = int(chunks[n + 1])  ## the way I have named the images, the spw number will
            ## always be just after the str "spw"
        diff = np.abs(-s + spw)
        ind = np.argsort(diff)
        non_zero = True
        for index in ind:
            if diff[index] != 0:
                return True, images[index]
        return False, ''

    def gen_blank_cal(self, spw):
        calprefix = self.caltbdir + 'selfcal'
        iteration_num = 0
        caltable = calprefix + '_spw_' + spw.zfill(2) + '_' + str(iteration_num).zfill(2) + '.gcal'
        os.system("rm -rf " + calprefix + '*')
        gencal(vis=self.vis, caltable=caltable, caltype='amp', spw=spw, antenna='', pol='', parameter=[1.0])
        return

    def find_phasecenter(self):
        """
        The purpose of this module is to find a new phasecenter at the flare location for imaging
        Returns
        -------
        Updates self.phasecenter (in J2000 RA and DEC) to be the flare location
        """
        imagename = self.workpath + "tmpimg"
        ra = []
        dec = []
        ###TODO: spw range is hard coded to use up to the first 3 spws of the supplied spw list,
        ### It could fail if they do not have a clear flare response. Need to be more adaptive.
        for j, s in enumerate(self.slfcal_spws):
            if j == 4:
                break
            current_trange = self.flare_times[j]

            tclean(vis=self.vis, imagename=imagename, timerange=current_trange, spw=str(s), uvrange=self.uvranges[j],
                   imsize=4096, cell=self.cell_vals[j], niter=0, gain=0.05, cyclefactor=10,
                   weighting='briggs', robust=0.0, savemodel='none', pbcor=False,
                   pblimit=0.01, stokes='XX')

            max_pix, min_pix, rms = self.get_img_stat(imagename)
            if max_pix / rms > 10:
                j += 1
                flare_loc_available = True
                pos = imstat(imagename + ".image")['maxposf']
                temp = pos.split(',')
                ra_str = temp[0].split(':')
                dec_str = temp[1].split('.')
                ra.append((abs(int(ra_str[0])) + int(ra_str[1]) / 60.0 + float(ra_str[2]) / 3600.) * 15)
                sign = 1
                if '-' in ra_str[0]:
                    sign = -1
                ra[-1] = ra[-1] * sign
                sign = 1
                if '-' in dec_str[0]:
                    sign = -1
                try:
                    dec.append(abs(int(dec_str[0])) + int(dec_str[1]) / 60.0 + float(
                        dec_str[2] + '0.' + dec_str[3]) / 3600.)
                except IndexError:
                    dec.append(abs(int(dec_str[0])) + int(dec_str[1]) / 60.0 + float(dec_str[2]) / 3600.)
                dec[-1] = dec[-1] * sign
                os.system("rm -rf " + imagename + ".*")
            else:
                os.system("rm -rf " + imagename + ".*")

        if flare_loc_available == True:
            ra_final = np.median(np.array(ra))
            dec_final = np.median(np.array(dec))
            phasecenter = 'J2000 ' + str(ra_final) + "deg " + str(dec_final) + "deg"
            self.phasecenter = phasecenter
            self.flare_loc_available = True
            print(phasecenter)
            # logf.write("Phasecenter:{}\n".format(phasecenter))
            os.system("rm -rf " + imagename + ".*")
        else:
            # logf.write("Appropriate phase center could not be found." + \
            #           "Please restart after providing an appropriate one.\n")
            print("Appropriate phase center could not be found." +
                  "Please restart after providing one manually by." +
                  "setting FlareSelfCalib.phasecenter='J2000 00h00m00s +00d00m00s'. ")
            self.phasecenter = ''
            self.flare_loc_available = False

    def do_selfcal(self, slfcalms, sp, spwran, uvrange='', cell_val='2arcsec', imsize=2048, ref_image='',
                   make_shifted_mask=False, combine_spws=False):
        clearcal(slfcalms)
        calprefix = self.caltbdir + 'selfcal'
        imgprefix = self.imagedir + 'selfcal'
        os.system('rm -rf ' + self.caltbdir + '*_spw_' + sp.zfill(2) + '*')
        os.system('rm -rf ' + self.imagedir + '*spw_' + sp.zfill(2) + '*')
        iteration_num = 0
        found_mask = False
        change_grow_mask_threshold = False
        grow_mask_threshold = 0.7
        image_to_go_back_to = ''
        self.logf.write("spw=" + str(sp).zfill(2) + " grow_mask_threshold={}\n".format(grow_mask_threshold))
        while iteration_num < self.maxiter:
            # print("starting round "+str(iteration_num))
            caltable = calprefix + '_spw_' + sp.zfill(2) + '_' + str(iteration_num).zfill(2) + '.gcal'
            imagename = imgprefix + '_spw_' + sp.zfill(2) + '_' + str(iteration_num).zfill(2)

            self.logf.write("spw=" + str(sp).zfill(2) + " Creating the first diry map\n")
            if iteration_num == 0:
                tclean(vis=slfcalms, imagename=imagename, uvrange=uvrange, spw=spwran, imsize=imsize, \
                       cell=cell_val, niter=1, gain=0.05, cyclefactor=10, weighting='briggs', stokes='XX', \
                       robust=0.0, savemodel='none', pbcor=False, pblimit=0.001, phasecenter=self.phasecenter)

                # print ("Check if image produced or not")
                if os.path.isdir(imagename + ".image") == False:
                    # print ("image not produced")
                    self.logf.write("spw=" + str(sp).zfill(2) + " Image not produced\n")
                    os.system('rm -rf ' + self.caltbdir + '*_spw_' + sp.zfill(2) + '*')
                    os.system('rm -rf ' + self.imagedir + '*spw_' + sp.zfill(2) + '*')
                    return False, False
                # print ("Getting image statistics")
                self.logf.write("spw=" + str(sp).zfill(2) + " Getting image statistics\n")
                max_pix, min_pix, rms = self.get_img_stat(imagename)
                # print (max_pix)
                threshold = 10 * rms
                if max_pix < 15 * rms and len(ref_image) == 0:
                    self.logf.write("spw=" + str(sp).zfill(2) + " Image SNR too low for selfcal with no mask\n")
                    self.logf.write("spw=" + str(sp).zfill(2) + " Trying to find past image\n")
                    # print("Trying to find past image")
                    past_image_found, past_image = self.find_previous_image(int(sp))
                    if past_image_found == True:
                        os.system('rm -rf ' + self.caltbdir + '*_spw_' + sp.zfill(2) + '*')
                        os.system('rm -rf ' + self.imagedir + '*spw_' + sp.zfill(2) + '*')
                        self.logf.write("spw=" + str(sp).zfill(2) + " Previous mask files found." + \
                                        "Trying to use the previous mask\n")
                        return False, True
                    self.logf.write("spw=" + str(sp).zfill(2) + " Previous mask does not exist. Autogenerating mask\n")
                    threshold = 0.5 * max_pix
                    exportfits(imagename=imagename + ".image", fitsimage=imagename + ".fits")
                    found_mask = self.gen_mask('', imagename + ".fits", '', imagename + ".mask", 10 * rms,
                                               imsize, sp, make_shifted_mask, grow_mask_threshold)
                    if found_mask == False:
                        self.logf.write("spw=" + str(sp).zfill(2) + " Mask not produced. Exiting\n")
                        # print("Mask not produced. Exitting.")
                        os.system('rm -rf ' + self.caltbdir + '*_spw_' + sp.zfill(2) + '*')
                        os.system('rm -rf ' + self.imagedir + '*spw_' + sp.zfill(2) + '*')
                        return False, True

            # elif max_pix<10*rms and len(ref_image)==0:
            #   print ("Source not detected")
            #  return False, False
            else:
                threshold = 10 * rms
                if combine_spws == False:
                    spwran = sp

            # print ("Starting actual clean")
            mask = ''
            if len(ref_image) != 0 and iteration_num == 0:
                # print (len(ref_image),ref_image)
                self.logf.write("spw=" + str(sp).zfill(2) + " Generating mask using the previous mask\n")
                # print ("Generating mask using the previous mask\n")
                threshold = 0.5 * max_pix
                exportfits(imagename=imagename + ".image", fitsimage=imagename + ".fits")
                max_pix, min_pix, rms = self.get_img_stat(ref_image[:-5])
                found_mask = self.gen_mask(ref_image, imagename + ".fits", ref_image[:-5] + ".mask",
                                           imagename + ".mask", 10 * rms, imsize, sp,
                                           make_shifted_mask, grow_mask_threshold)

                if found_mask == False:
                    self.logf.write("spw=" + str(sp).zfill(2) + " Mask not generated." + \
                                    " Trying to do more averaging in frequency\n ")
                    # print("Mask  not generated. Trying to do more averaging in frequency")
                    return False, True

                # print(imagename+".image")
                self.logf.write("spw=" + str(sp).zfill(2) + " Checking if the maximum pixel inside the " + \
                                "mask is in the correct position\n")
                success = self.confirm_maximum_pixel(imagename + ".image", imagename + ".mask",
                                                     spwran, slfcalms, uvrange, imsize, cell_val,
                                                     sp)

                if success == False:
                    self.logf.write("spw=" + str(sp).zfill(2) + " Probable shift in image plane." + \
                                    " Do more averaging in frequency.\n ")
                    print("Probable shift in image plane. Do more averaging in frequency.")
                    return False, True

            elif len(ref_image) != 0:
                # max_pix,min_pix,rms=get_img_stat(imagename+'.image')
                threshold = 0.2 * max_pix
                max_pix, min_pix, rms = self.get_img_stat(ref_image[:-5])
                exportfits(imagename=imagename[:-2] + str(iteration_num - 1).zfill(2) + ".image", fitsimage=
                imagename[:-2] + str(iteration_num - 1).zfill(2) + ".fits", overwrite=True)
                # print("Fits file generated")
                # print (imagename, iteration_num)
                found_mask = self.gen_mask(ref_image,
                                           imagename[:-2] + str(iteration_num - 1).zfill(2) + ".fits",
                                           ref_image[:-5] + '.mask', imagename[:-2] + "00.mask", 10 * rms,
                                           imsize, sp, make_shifted_mask, grow_mask_threshold)
                ### if found_mask=False, \
                #### the earlier mask was \
                #### not updated and hence can be used.
                # print (found_mask)

                mask = imagename[:-2] + "00.mask"
                # print (max_pix)

            elif iteration_num != 0 and len(ref_image) == 0 and max_pix < 15 * rms:
                threshold = 0.2 * max_pix
                exportfits(imagename=imagename[:-2] + str(iteration_num - 1).zfill(2) + ".image", fitsimage=
                imagename[:-2] + str(iteration_num - 1).zfill(2) + ".fits", overwrite=True)
                found_mask = self.gen_mask('', imagename[:-2] + str(iteration_num - 1).zfill(2) + ".fits",
                                           '', imagename[:-2] + "00.mask", 10 * rms, imsize, sp,
                                           make_shifted_mask, grow_mask_threshold)
                mask = imagename[:-2] + "00.mask"
            # print (str(threshold)+"Jy")
            self.logf.write("spw=" + str(sp).zfill(2) + " Cleaning threshold:{}\n".format(threshold))
            # if no_shift_mask==True and iteration_num>3:
            #   raise RuntimeError("shift not in image")
            if found_mask == True or max_pix > 15 * rms:
                tclean(vis=slfcalms, imagename=imagename, spw=spwran, uvrange=uvrange,
                       imsize=imsize, threshold=str(threshold) + "Jy", cell=cell_val, niter=10000,
                       gain=0.05, cyclefactor=10, weighting='briggs', robust=0.0,
                       savemodel='modelcolumn', pbcor=False, pblimit=0.01, mask=mask,
                       phasecenter=self.phasecenter, stokes='XX')

                ### sometimes it happens that due to some error, a blank model is generated.
                ### But this should not be the case, as at least the maximum pixel should be
                ### picked. I saw this happening in the case of 20210507 spw 45. CASA said
                ### and I am copying it " Caught Exception: NonLinearFitLM: error in loop
                ### solution" and then it restored a blank image. In these situations, gaincal
                ### will essentially use a point source at the phase center as a model. But, this
                ### should not be the case. Hence, I will check if the model is blank or not.
                ### If yes, I will assume that some problem has happened and I will not proceed

                ia.open(imagename + ".model")
                modeldata = ia.getchunk()
                ia.close()
                if np.nansum(modeldata) < 1e-6:
                    # print("Blank model. Checking if this is because of very high threshold.")
                    self.logf.write("spw=" + str(sp).zfill(2) + " Blank model. Checking if this" + \
                                    " is because of very high threshold.\n")
                    max_pix = imstat(imagename=imagename + ".image", mask="mask(\'" + mask + "\')")["max"]
                    if max_pix < threshold:
                        threshold = max_pix * 0.5
                        self.logf.write("spw=" + str(sp).zfill(2) + " Threshold was too high." + \
                                        "new threshold={}\n".format(threshold))
                        tclean(vis=slfcalms, imagename=imagename, spw=spwran, uvrange=uvrange, \
                               imsize=imsize, threshold=str(threshold) + "Jy", cell=cell_val, niter=10000, \
                               gain=0.05, cyclefactor=10, weighting='briggs', robust=0.0, \
                               savemodel='modelcolumn', pbcor=False, pblimit=0.01, mask=mask, \
                               phasecenter=self.phasecenter, stokes='XX')
                        ia.open(imagename + ".model")
                        modeldata = ia.getchunk()
                        ia.close()
                        if np.nansum(modeldata) < 1e-6:
                            self.logf.write("spw=" + str(sp).zfill(2) + " Model still blank. Exiting\n")
                            return False, True
                    else:
                        return False, True  ### I am returning True here, as I am not really sure
                    ### if it cannot recover ever. It is better to let it
                    ### try, rather than leaving it False. However, my gut
                    ### feeling is that, it would not be able to recover.

            else:
                # print("Mask not found or rms is very high")
                return False, True

            max_pix, min_pix, rms = self.get_img_stat(imagename)
            dyn_range = max_pix / rms
            print("dyn_range=", dyn_range, iteration_num)
            self.logf.write("spw=" + str(sp).zfill(2) + " dynamic range:" + str(dyn_range) + ", Iteration" + \
                            " number:" + str(iteration_num) + "\n")
            if iteration_num <= 1:
                dyn_range1 = dyn_range
                min_pix1 = min_pix
                max_pix1 = max_pix
            else:
                if dyn_range < 10 and dyn_range < dyn_range1 * 0.95:  ##len(ref_image)==0 and mask=='':
                    self.logf.write("spw=" + str(sp).zfill(2) + " SNR too low for selfcal and SNR decreasing\n")
                    print("SNR too low for selfcal")
                    if change_grow_mask_threshold == True:
                        print("spw=" + str(sp).zfill(2) + " Reverting back and exiting with success")
                        FlareSelfCalib.restore_previous_condition(image_to_go_back_to)
                        return True, True
                    return False, True
                elif dyn_range / dyn_range1 > 0.95 and dyn_range / dyn_range1 < 1.05 and iteration_num > 2:  ### otherwise, this means no more improvement
                    print(dyn_range, max_pix / abs(min_pix))

                    if change_grow_mask_threshold == False:
                        change_grow_mask_threshold = True
                        grow_mask_threshold = 0.2
                        self.logf.write("spw=" + str(sp).zfill(2) + " Changing grow_mask_threshold to" + \
                                        "{}\n".format(grow_mask_threshold))
                        image_to_go_back_to = imagename
                        print(imagename)
                    else:
                        self.logf.write("spw=" + str(sp).zfill(2) + " Converged\n")
                        return True, True
                elif dyn_range / dyn_range1 < 1.1 and iteration_num == 5 and spwran == sp:
                    self.logf.write("spw=" + str(sp).zfill(2) + " Not enough improvement.\n")
                    print("Not enough improvement")
                    if change_grow_mask_threshold == False:
                        change_grow_mask_threshold = True
                        grow_mask_threshold = 0.5
                        image_to_go_back_to = imagename
                        self.logf.write("spw=" + str(sp).zfill(2) + " Changing grow_mask_threshold to" + \
                                        "{}\n".format(grow_mask_threshold))
                    else:
                        return True, True

            # if iteration_num>=1:
            #	flagmanager(vis=slfcalms,mode='restore',versionname='applycal_1')
            #	flagmanager(vis=slfcalms,mode='delete',versionname='applycal_1')
            if os.path.isdir(caltable):
                os.system("rm -rf " + caltable)
            gaincal(vis=slfcalms, refant=self.refantenna, spw=sp, caltable=caltable, uvrange=uvrange, \
                    calmode='p', combine='scan', minblperant=4, minsnr=3, append=False, \
                    solnorm=True, solmode='L1R', rmsthresh=[10, 7], normtype='median')

            if os.path.isdir(caltable) == False:
                self.logf.write("spw=" + str(sp).zfill(2) + " Caltable not produced\n")
                print("Solution not found")
                if change_grow_mask_threshold == True:
                    FlareSelfCalib.restore_previous_condition(image_to_go_back_to)
                    print("spw=" + str(sp).zfill(2) + " Reverting back and exiting with success")
                    return True, True
                return False, True

            bpass, flag, success = FlareSelfCalib.read_bandpass(caltable)
            pos = np.where(flag[0, 0, 0, :] == True)[0]
            num_flagged_ant = np.size(pos)
            if iteration_num == 0:
                num_flagged_ant1 = num_flagged_ant
            self.logf.write("Antennas flagged in caltable:" + str(pos) + '\n')
            if num_flagged_ant > num_flagged_ant1 + 3:
                self.logf.write("spw=" + str(sp).zfill(2) + " Flagged antennas have increased with iteration. " + \
                                "Possible divergence. Exiting\n")
                print("Flagged anttenae number increased a lot")
                if change_grow_mask_threshold == True:
                    print("spw=" + str(sp).zfill(2) + " Reverting back and exiting with success")
                    FlareSelfCalib.restore_previous_condition(image_to_go_back_to)
                    return True, True
                return False, True

            if 15 - len(pos) < 5:  ### will not proceed with this selfcal if 5 or more antennas are flaggged
                self.logf.write("spw=" + str(sp).zfill(2) + " Too many flagged antennas. Exiting\n")
                print("Too many flagged antennas")
                if change_grow_mask_threshold == True:
                    print("spw=" + str(sp).zfill(2) + " Reverting back and exiting with success")
                    FlareSelfCalib.restore_previous_condition(image_to_go_back_to)
                    return True, True
                return False, True

            ### TODO Need to check if this helps the process or makes it worse. Initial guess is it makes it worse.
            applycal(vis=slfcalms, spw=sp, gaintable=caltable, applymode='calonly')
            if iteration_num != self.maxiter - 1:
                delmod(slfcalms, scr=True)
            iteration_num += 1
            dyn_range1 = dyn_range
            min_pix1 = min_pix
            max_pix1 = max_pix
            num_flagged_ant1 = num_flagged_ant
        self.logf.write("spw=" + str(sp).zfill(2) + " Successfull exit\n")
        return True, True

    def calling_do_selfcal(self, slfcalms, s, uvrange='', cell_val='2arcsec'):
        sp = str(s)
        ref_image = ''
        print('processing spw: ' + sp)
        self.logf.write("spw=" + str(s).zfill(2) + " Cell size=" + str(cell_val) + "\n")
        # f.write("Imsize="+str(imsize)+'\n')
        self.logf.write("spw=" + str(s).zfill(2) + " Uvrange= " + uvrange + '\n')

        success, signal = self.do_selfcal(slfcalms, sp, sp, uvrange=uvrange, cell_val=cell_val)

        if success == False and signal == False:
            return success

        if success == False:
            self.logf.write("spw=" + str(s).zfill(2) + " Trying with previous mask\n")
            print("Trying with previous mask \n \n \n \n")
            ### here I will try to find the image at nearest spw. Assumption is
            ### that if it is at the nearest frequency, it is very likely that
            ### the sources at both of them will be located close by.
            past_image_found, past_image = self.find_previous_image(s)
            if past_image_found == True:
                self.logf.write("spw=" + str(s).zfill(2) + " Past image found. Using already used masks\n")
                fitsimage = past_image[:-6] + ".fits"
                if os.path.isfile(fitsimage) == False:
                    exportfits(imagename=past_image, fitsimage=fitsimage)

                ref_image = fitsimage
                success, signal = self.do_selfcal(slfcalms, sp, sp, uvrange=uvrange, cell_val=cell_val,
                                                  ref_image=ref_image)

        if success == False:  ### this means that we would need to do
            ### multi-frequency synthesis. To do that
            ### first I would make sure that the maximum
            ### pixel within the previous mask in the
            ### dirty map is a real feature and not an
            ### artifact.
            print("Trying to do mfs \n \n \n \n")
            self.logf.write("spw=" + str(s).zfill(2) + " Trying to do mfs\n")
            imgprefix = self.imagedir + 'selfcal'
            imagename = imgprefix + '_spw_' + sp.zfill(2) + '_00'
            mask = imagename + ".mask"
            maxpix = imstat(imagename + ".image", mask="mask(\'" + mask + "\')")["maxpos"]

            ###TODO A big assumption is that the source is at the same location at all
            ### frequencies, in the sense that they lie in the box created by the lower
            ### frequency. If the box is very small and the true source location moves,
            ###  then there is a problem. Hence, we need to do another check. We shall
            ### take a box of size lets say 5 arcminutes centred on the original mask and
            ### also caluclate the maximum pixel there. If the maximum pixel inside this
            ### box remains same even when we do mfs, then probably the true source is
            ### also shifted.
            freq_int_found = False
            for avg_spw in range(1, self.avg_spw_max):
                min_spw = max(self.minimum_spw, s - avg_spw)
                max_spw = min(self.maximum_spw, s + avg_spw)
                #### implementing a check where I ensure that I will never average more
                #### more than 0.5 times the central frequency
                if abs(self.freqs_ms[max_spw] - self.freqs_ms[min_spw]) > self.max_frac_freq_avg * self.freqs_ms[s]:
                    self.logf.write("Did not a suitable averaging range in frequency.",
                                    "Trying with maximum possible frequency bandwidth.\n")
                    break
                if s - avg_spw == min_spw - 1 and s + avg_spw == max_spw:
                    self.logf.write("Did not a suitable averaging range in frequency. Exiting.")
                    break
                if min_spw != max_spw:
                    spwran = str(min_spw) + '~' + str(max_spw)
                else:
                    spwran = sp
                print('processling spw {0:s} using spw {1:s} as model'.format(sp, spwran))
                self.logf.write('Processing spw {0:s} using spw {1:s} as model\n'.format(sp, spwran))

                if os.path.isdir("check_maxpix.psf"):
                    os.system("rm -rf check_maxpix*")
                tclean(vis=slfcalms, imagename='check_maxpix', spw=spwran, uvrange=uvrange,
                       imsize=2048, cell=cell_val, niter=0, gain=0.05, cyclefactor=10,
                       weighting='briggs', stokes='XX',
                       robust=0.0, savemodel='none', pbcor=False, pblimit=0.01,
                       mask='user', phasecenter=self.phasecenter)
                maxpix1 = imstat("check_maxpix.image", mask="mask(\'" + mask + "\')")["maxpos"]
                print(maxpix, maxpix1)
                if maxpix[0] != maxpix1[0] and maxpix[1] != maxpix1[1]:
                    shift_ok = FlareSelfCalib.check_shift("check_maxpix.image",
                                                          np.array([maxpix[0] - maxpix1[0], maxpix[1] - maxpix1[1]]),
                                                          cell_val)  ## 2 is the cellsize
                else:
                    shift_ok = True
                print("not an artefact", shift_ok)
                self.logf.write("spw=" + str(s).zfill(2) + " Not an artefact? {}\n".format(shift_ok))
                if shift_ok == True:
                    freq_int_found = True
                    break
                else:
                    maxpix = maxpix1
                    os.system("rm -rf check_maxpix.*")
            if freq_int_found == False:
                self.logf.write("spw=" + str(s).zfill(2) + " Doing the best we can get and trying it out\n")
                s -= 1
                min_spw = max(self.minimum_spw, s - avg_spw)
                max_spw = min(self.maximum_spw, s + avg_spw)
                freq_int_found = True
            print(min_spw, max_spw)
            self.logf.write("spw=" + str(s).zfill(2) + ":" + str(min_spw) + "," + str(max_spw))
            os.system("rm -rf check_maxpix.*")
            if freq_int_found == True:
                for avg_spw1 in range(avg_spw, self.avg_spw_max):
                    min_spw = max(self.minimum_spw, s - avg_spw1)
                    max_spw = min(self.maximum_spw, s + avg_spw1)
                    if s - avg_spw == min_spw - 1 and s + avg_spw == max_spw:
                        success = False
                        break
                    if min_spw != max_spw:
                        spwran = str(min_spw) + '~' + str(max_spw)
                    else:
                        spwran = sp

                    print("Calling do_selfcal with mfs")
                    self.logf.write("spw=" + str(s).zfill(2) + ":" + str(min_spw) + "," + str(max_spw) + "\n")
                    success, signal = self.do_selfcal(slfcalms, sp, spwran, uvrange=uvrange, cell_val=cell_val,
                                                      ref_image=ref_image, make_shifted_mask=True)
                    print("success=", success)
                    self.logf.write("spw=" + str(s).zfill(2) + " success={}".format(success))

                    if success == True and signal == True:
                        break
            else:
                success = False

        if success == False and freq_int_found == True:
            print("Calling do_selfcal with mfs and solving gains at all frequencies")
            self.logf.write("spw=" + str(s).zfill(2) + ":" + str(min_spw) + "," + str(max_spw) + "\n")
            success, signal = self.do_selfcal(slfcalms, sp, spwran, uvrange=uvrange, cell_val=cell_val,
                                              ref_image=ref_image, make_shifted_mask=True, combine_spws=True)
            print("success=", success)
            self.logf.write("spw=" + str(s).zfill(2) + " success={}".format(success))

        return success

    def slfcal_init(self):
        ##=========== Obtain information from the input visibility =======
        nspw_ms = FlareSelfCalib.get_spw_num(self.vis)
        freqs_ms = FlareSelfCalib.get_ref_freqlist(self.vis)
        if len(self.flag_antennas) > 0:
            flagdata(vis=self.vis, mode='manual', antenna=self.flag_antennas)

        ##=========== Determine the flare peak and selfcal timerange ########
        # TODO: if flare_peak_time and duration is already given, this (time-consuming) step can be skipped
        flare_finder_timer = timeit.default_timer()
        if not self.flare_time_available:
            print('No flare time provided. Trying to find flare times from the visibility')
            self.flare_finder()
            time1 = timeit.default_timer()
            self.logf.write("Finding flare time took {0:d} seconds:" + str(time1 - flare_finder_timer))

        tmpms = self.workpath + 'temp_ms.ms'
        if os.path.exists(tmpms):
            os.system('rm -rf {0:s}'.format(tmpms))

        split(vis=self.vis, outputvis=tmpms, spw=self.slfcal_spwstr,
              datacolumn='data', timerange=self.flare_times[0])

        nspw_slfcal = FlareSelfCalib.get_spw_num(tmpms)
        freqs_slfcal = FlareSelfCalib.get_ref_freqlist(tmpms)
        cell_vals = []
        uvranges = []
        if self.calc_cell:
            cell = FlareSelfCalib.calc_cellsize(tmpms)
            for i in range(len(cell)):
                cell_vals.append(str(cell[i]) + 'arcsec')
                uvranges.append('>' + str(self.uvlim * freqs_slfcal[i] / freqs_slfcal[0]) + 'lambda')
        elif len(self.cell) < nspw_slfcal:
            print("Number of beams provided does not match number of spw." +
                  " Using the first beam value and then I will scale with frequency")

            for i in range(len(freqs_slfcal)):
                cell_vals.append(str(self.cell[0] * freqs_slfcal[0] / freqs_slfcal[i]) + 'arcsec')
                uvranges.append('>' + str(self.uvlim * freqs_slfcal[i] / freqs_slfcal[0]) + 'lambda')
        ### scaled the beam here with frequency
        else:
            for i in range(nspw_slfcal):
                cell_vals.append(str(self.cell[i]) + 'arcsec')
                uvranges.append('>' + str(self.uvlim * freqs_slfcal[i] / freqs_slfcal[0]) + 'lambda')

        self.uvranges = uvranges
        self.cell_vals = cell_vals
        self.freqs_slfcal = freqs_slfcal
        self.nspw_ms = nspw_ms
        self.freqs_ms = freqs_ms
        self.nspw_slfcal = nspw_slfcal
        os.system('rm -rf {0:s}'.format(tmpms))

        ## ======= find flare location and change phase center ========
        if not self.flare_loc_available:
            flare_loc_timer = timeit.default_timer()
            print('No flare location provided. Trying to find flare location from the visibility')
            self.find_phasecenter()
            time1 = timeit.default_timer()
            self.logf.write("Finding flare location took {0:d} seconds:" + str(time1 - flare_loc_timer))

    def slfcal_pipeline(self, doselfcal=True, doimaging=False):
        print("Starting the self-calibration process")
        logf = open(self.logfile, "a")
        self.logf = logf

        ## ====== perform flare detection, initialize the pipeline ============
        self.slfcal_init()

        if doselfcal:
            ## ====== perform self-calibration =========
            for j, s in enumerate(self.slfcal_spws):
                success = False
                uvrange = self.uvranges[j]
                cell_val = self.cell_vals[j]
                start_time = self.flare_start_datetime[j]
                end_time = self.flare_end_datetime[j]
                time_delta = (end_time - start_time).seconds
                if self.found_flares[j]:
                    max_time_delta = 180  ## maximum is set to 3 minutes
                    if time_delta < 20:
                        time_increase = 5
                    elif 20 < time_delta < 50:
                        time_increase = 10
                    else:
                        time_increase = 15
                else:
                    max_time_delta = 300  ## 5 minutes
                    time_increase = 30
                while not success:
                    time_delta = (end_time - start_time).seconds
                    if time_delta > max_time_delta:
                        break
                    startstr = start_time.strftime('%Y/%m/%d/%H:%M:%S')
                    endstr = end_time.strftime('%Y/%m/%d/%H:%M:%S')
                    current_trange = startstr + "~" + endstr
                    self.logf.write("spw=" + str(s).zfill(2) + " Self-calibrating spw " + str(s) + \
                                    " with time range " + current_trange + '\n')

                    ms_slfcaled = self.vis.replace('.ms', '.slfcaled_v2.ms')  # output, selfcaled, visibility
                    slfcalms = self.slfcaldir + 'slfcalms.XX.slfcal'  # intermediate small visibility for selfcalbration
                    slfcaledms = self.slfcaldir + 'slfcalms.XX.slfcaled'
                    print("Spliting from " + self.vis)
                    self.logf.write("spw=" + str(s).zfill(2) + " Splitting..\n")
                    ##TODO: is this only performed for the first spw?
                    if not os.path.exists(slfcalms):
                        split(vis=self.vis, outputvis=slfcalms, datacolumn='data', timerange=current_trange,
                              correlation='XX', spw=self.selfcal_spw)
                    if self.identify_data_gap:
                        self.flag_data_gap(slfcalms, s)

                    print('Processing ' + current_trange)
                    success = self.calling_do_selfcal(slfcalms, s, uvrange=uvrange, cell_val=cell_val)

                    self.logf.write("success=" + str(success))
                    if success != True:
                        start_time = start_time - dt.timedelta(seconds=time_increase)
                        end_time = end_time + dt.timedelta(seconds=time_increase)
                    os.system("rm -rf " + slfcalms + "*")
                    break

                sp = str(s)
                self.logf.write("spw=" + str(s).zfill(2) + " success={}\n".format(success))
                calprefix = self.caltbdir + 'selfcal'
                caltable = calprefix + '_spw_' + sp.zfill(2)
                imgprefix = self.imagedir + "selfcal"
                image = imgprefix + '_spw_' + sp.zfill(2)

                if success == False:
                    for str1 in ['.image', '.mask', '.residual', '.psf', \
                                 '.model', '.sumwt', '.pb', '.fits']:
                        os.system("rm -rf " + image + "*" + str1)
                    os.system("rm -rf " + caltable + "*.gcal")
                    self.gen_blank_cal(sp)

                caltables = glob.glob(caltable + "*.gcal")
                num_caltable = len(caltables)
                num_img = len(glob.glob(image + "*.image"))
                if num_caltable != 1:
                    for i in range(0, num_caltable - 1):
                        os.system("rm -rf " + caltable + "_" + str(i).zfill(2) + ".gcal")

                for i in range(0, num_img - 1):
                    for str1 in ['.image', '.mask', '.residual', '.psf', \
                                 '.model', '.sumwt', '.pb', '.fits']:
                        os.system("rm -rf " + image + "_" + str(i).zfill(2) + str1)
                final_img = image + "_" + str(num_img - 1).zfill(2) + ".image"
                print(self.vis, current_trange)
                hf.imreg(vis=self.vis, imagefile=final_img, timerange=current_trange,
                         fitsfile=final_img[:-6] + "_helio.fits", \
                         usephacenter=False, verbose=False)  ### converting final image to heliocentric coordinates
                final_cal = self.caltbdir + "final_cal_spw_" + str(s).zfill(2) + ".gcal"
                if not os.path.isdir(final_cal):
                    os.system("mv " + caltable + "_" + str(num_caltable - 1).zfill(2) + ".gcal " + final_cal)
                else:
                    gaincal(vis=slfcalms, refant=self.refantenna, spw=sp, caltable=caltable, uvrange=uvrange, \
                            calmode='p', combine='scan', minblperant=4, minsnr=3, append=True, solnorm=True)
                # print('Calibration done in {0:s}'.format(slfcalms))
                os.system("rm -rf " + slfcalms + "*")

            selfcal_timer = timeit.default_timer()
            self.logf.write("Time taken for selfcal in seconds is " + str(selfcal_timer - start))

            final_ms = self.vis.replace('.ms', '_selfcaled.ms')
            print(final_ms, self.vis)
            print("Applying caltables")
            if os.path.isdir(final_ms) == False:
                for s in self.slfcal_spws:
                    caltable = self.caltbdir + 'final_cal_spw_' + str(s).zfill(2) + '.gcal'
                    if os.path.isdir(caltable) == False:
                        continue
                    applycal(vis=self.vis, gaintable=caltable, spw=str(s), applymode='calonly', interp='nearest')
                ##TODO: perhaps we only need to split out the spws that are self-calibrated
                # split(vis=self.vis, outputvis=final_ms, spw=self.spwstr, correlation='XX', datacolumn='corrected')
                split(vis=self.vis, outputvis=final_ms, correlation='XX', datacolumn='corrected')
                self.logf.write("Calibrated MS split\n")
                time1 = timeit.default_timer()
                self.logf.write("Calibrated data split in seconds:" + str(time1 - selfcal_timer))

            self.final_ms = final_ms

        if doimaging:
            imaging_timer = timeit.default_timer()
            ##TODO: the following relies on the completion of the prior self-calibration produced images.
            ## Need to remove such dependence if these images are not available.
            image_list = glob.glob(imgprefix + "*_helio.fits")
            if hasattr(self, 'final_ms'):
                if os.path.exists(self.final_ms):
                    msname = self.final_ms
                else:
                    print('self-calibrated ms does not exist. Using original visibility for imaging.')
                    msname = self.vis
            xycen = self.get_img_center_heliocoords(image_list)

            imaging_start_mjd = self.flare_peak_mjd - self.total_duration / 2
            if imaging_start_mjd < self.ms_startmjd:
                diff = imaging_start_mjd - self.ms_startmjd
                imaging_start_time = (self.flare_peak_datetime - dt.timedelta(seconds=diff)).strftime(
                    "%Y/%m/%d/%H:%M:%S")
            else:
                imaging_start_time = (
                        self.flare_peak_datetime - dt.timedelta(seconds=self.total_duration / 2)).strftime(
                    "%Y/%m/%d/%H:%M:%S")
            imaging_end_mjd = self.flare_peak_mjd + self.total_duration / 2
            if imaging_end_mjd > self.ms_endmjd:
                diff = self.ms_endmjd - imaging_start_mjd
                imaging_end_time = (self.flare_peak_datetime + dt.timedelta(seconds=diff)).strftime("%Y/%m/%d/%H:%M:%S")
            else:
                imaging_end_time = (self.flare_peak_datetime + dt.timedelta(seconds=self.total_duration / 2)).strftime(
                    "%Y/%m/%d/%H:%M:%S")
            time_str = imaging_start_time + "~" + imaging_end_time

            specfile = msname[:-3] + "_dspec.npz"
            spws = []
            for s in self.slfcal_spws:
                spws.append(str(s))
            cell_size1 = [str(self.cell_size) + "arcsec"]
            stokes = 'XX'
            mkmovie = True
            docompress = True
            opencontour = False
            clevels = [0.6, 1.0]
            plotaia = True
            aiawave = 1600
            movieformat = 'mp4'
            overwrite = False
            subregion = ''  # box[[128pix,128pix],[284pix,384pix]]'
            qlookplot.qlookplot(specfile=specfile, vis=msname, timerange=time_str, spw=spws,
                                ncpu=1, imsize=[self.final_imsize], cell=cell_size1, restoringbeam=[self.beam_1GHz],
                                robust=0.5, opencontour=opencontour, clevels=clevels, plotaia=plotaia,
                                aiawave=aiawave, mkmovie=mkmovie, twidth=int(self.final_image_int),
                                docompress=docompress,
                                stokes=stokes, movieformat=movieformat, uvrange='',
                                niter=300, overwrite=overwrite, xycen=xycen, fov=[256, 256])
            final_clean_timer = timeit.default_timer()
            self.logf.write("Final clean done in seconds:" + str(final_clean_timer - imaging_timer))

        self.logf.close()
