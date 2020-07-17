import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sunpy import map as smap
import sunpy.cm.cm as cm_smap
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar as colorbar
import matplotlib.patches as patches
from datetime import timedelta
from datetime import datetime
from glob import glob
import numpy as np
from astropy.time import Time
import calendar
from suncasa.utils import plot_mapX as pmX
import urllib.request
from sunpy.physics.differential_rotation import diffrot_map
from suncasa.utils import DButil

# imgfitsdir = '/Users/fisher/myworkspace/'
# imgfitstmpdir = '/Users/fisher/myworkspace/fitstmp/'
# pltfigdir = '/Users/fisher/myworkspace/SynopticImg/eovsamedia/eovsa-browser/'

imgfitsdir = '/data1/eovsa/fits/synoptic/'
imgfitstmpdir = '/data1/workdir/fitstmp/'
pltfigdir = '/common/webplots/SynopticImg/eovsamedia/eovsa-browser/'


def clearImage():
    for (dirpath, dirnames, filenames) in os.walk(pltfigdir):
        for filename in filenames:
            for k in ['0094', '0193', '0335', '4500', '0171', '0304', '0131', '1700', '0211', '1600', '_HMIcont',
                      '_HMImag']:
                # for k in ['_Halph_fr']:
                if k in filename:
                    print(os.path.join(dirpath, filename))
                    os.system('rm -rf ' + os.path.join(dirpath, filename))


def pltEovsaQlookImageSeries(timobjs, spw, vmax, vmin, aiawave, fig=None, axs=None, imgoutdir=None, overwrite=False,
                             verbose=False):
    plt.ioff()
    imgfiles = []
    dpi = 512. / 4
    aiaDataSource = {"0094": 8,
                     "0193": 11,
                     "0335": 14,
                     # "4500": 17,
                     "0171": 10,
                     "0304": 13,
                     "0131": 9,
                     "1700": 16,
                     "0211": 12,
                     # "1600": 15,
                     "_HMIcont": 18,
                     "_HMImag": 19}

    tmjd = timobjs.mjd
    tmjd_base = np.floor(tmjd)
    tmjd_hr = (tmjd - tmjd_base) * 24
    for tidx, timobj in enumerate(timobjs):
        dateobj = timobj.to_datetime()
        timestr = dateobj.strftime("%Y-%m-%dT%H:%M:%SZ")
        tstrname = dateobj.strftime("%Y%m%dT%H%M%SZ")
        datestrdir = dateobj.strftime("%Y/%m/%d/")
        imgindir = imgfitsdir + datestrdir
        timobj_prevday = Time(tmjd[tidx] - 1, format='mjd')
        dateobj_prevday = timobj_prevday.to_datetime()
        datestrdir_prevday = dateobj_prevday.strftime("%Y/%m/%d/")
        imgindir_prevday = imgfitsdir + datestrdir_prevday

        # if not os.path.exists(imgindir): os.makedirs(imgindir)
        cmap = cm_smap.get_cmap('sdoaia304')

        if verbose: print('Processing EOVSA images for date {}'.format(dateobj.strftime('%Y-%m-%d')))
        key, sourceid = aiawave, aiaDataSource[aiawave]
        s, sp = 0, spw[0]

        figoutname = os.path.join(imgoutdir, 'eovsa_bd{:02d}_aia{}_{}.jpg'.format(s + 1, key, tstrname))

        if overwrite or (not os.path.exists(figoutname)):
            ax = axs[0]
            ax.cla()
            spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
            eofile = imgindir + 'eovsa_{}.spw{}.tb.disk.fits'.format(dateobj.strftime('%Y%m%d'), spwstr)
            if not os.path.exists(eofile):
                continue

            norm = colors.Normalize(vmin=vmin[s], vmax=vmax[s])
            eomap = smap.Map(eofile)
            eomap = eomap.resample(u.Quantity(eomap.dimensions) / 2)
            eomap_rot = diffrot_map(eomap, time=timobj)
            eomap_rot.data[np.where(eomap_rot.data < 0)] = np.nan

            t_hr = tmjd_hr[tidx]
            t_hr_st_blend = 0.0
            t_hr_ed_blend = 12.0
            if t_hr_st_blend <= t_hr < t_hr_ed_blend:
                eofile_prevday = imgindir_prevday + 'eovsa_{}.spw{}.tb.disk.fits'.format(
                    dateobj_prevday.strftime('%Y%m%d'), spwstr)
                if not os.path.exists(eofile_prevday):
                    continue
                eomap_prevd = smap.Map(eofile_prevday)
                eomap_prevd = eomap_prevd.resample(u.Quantity(eomap_prevd.dimensions) / 2)
                eomap_rot_prevd = diffrot_map(eomap_prevd, time=timobj)
                eomap_rot_prevd.data[np.where(eomap_rot_prevd.data < 0)] = np.nan
                alpha = (t_hr - t_hr_st_blend) / (t_hr_ed_blend-t_hr_st_blend)
                alpha_prevd = 1.0 - alpha
                eomap_plt = smap.Map((eomap.data*alpha+eomap_prevd.data*alpha_prevd),eomap.meta)
                eomap_rot_plt = smap.Map((eomap_rot.data*alpha+eomap_rot_prevd.data*alpha_prevd), eomap_rot.meta)
            else:
                eomap_plt = eomap
                eomap_rot_plt = eomap_rot

            eomap_plt.plot(axes=ax, cmap=cmap, norm=norm)
            eomap_rot_plt.plot(axes=ax, cmap=cmap, norm=norm)

            # eomap_ = pmX.Sunmap(eomap)
            # eomap_.imshow(axes=ax, cmap=cmap, norm=norm)
            # eomap_.draw_limb(axes=ax, lw=0.25, alpha=0.5)
            # eomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.25)

            # eomap_rot_ = pmX.Sunmap(eomap_rot)
            # eomap_rot_.imshow(axes=ax, cmap=cmap, norm=norm)
            # eomap_rot.plot(axes=ax, cmap=cmap, norm=norm, alpha=alpha)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.text(0.02, 0.02,
                    'EOVSA {:.1f} GHz  {}'.format(eomap.meta['CRVAL3'] / 1e9, dateobj.strftime('%d-%b-%Y %H:%M UT')),
                    transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=9)
            ax.text(0.98, 0.02, 'Max Tb {:.0f} K'.format(np.nanmax(eomap.data)),
                    transform=ax.transAxes, color='w', ha='right', va='bottom', fontsize=9)
            ax.set_xlim(-1227, 1227)
            ax.set_ylim(-1227, 1227)

            ax = axs[1]
            # for key, sourceid in aiaDataSource.items():
            ax.cla()

            if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
            sdourl = 'https://api.helioviewer.org/v2/getJP2Image/?date={}&sourceId={}'.format(timestr, sourceid)
            sdofile = os.path.join(imgoutdir, 'AIA' + key + '.{}.jp2'.format(tstrname))
            if not os.path.exists(sdofile):
                urllib.request.urlretrieve(sdourl, sdofile)

            if not os.path.exists(sdofile): continue
            sdomap = smap.Map(sdofile)
            norm = colors.Normalize()
            # sdomap_ = pmX.Sunmap(sdomap)
            if "HMI" in key:
                cmap = plt.get_cmap('gray')
            else:
                cmap = cm_smap.get_cmap('sdoaia' + key.lstrip('0'))
            sdomap.plot(axes=ax, cmap=cmap, norm=norm)
            # sdomap_.imshow(axes=ax, cmap=cmap, norm=norm)
            # sdomap_.draw_limb(axes=ax, lw=0.25, alpha=0.5)
            # sdomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.25)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.text(0.02, 0.02,
                    '{}/{} {}  {}'.format(sdomap.observatory, sdomap.instrument.split(' ')[0], sdomap.measurement,
                                          sdomap.date.strftime('%d-%b-%Y %H:%M UT')),
                    transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=9)
            ax.set_xlim(-1227, 1227)
            ax.set_ylim(-1227, 1227)
            fig.savefig(figoutname, dpi=np.int(dpi), quality=85)
        imgfiles.append(figoutname)

    return imgfiles


def main(year, month, day=None, dayspan=30):
    '''
    By default, the subroutine create EOVSA monthly movie
    '''
    # tst = datetime.strptime("2017-04-01", "%Y-%m-%d")
    # ted = datetime.strptime("2019-12-31", "%Y-%m-%d")
    if day is None:
        dst, ded = calendar.monthrange(year, month)
        tst = datetime(year, month, 1)
        ted = datetime(year, month, ded)
    else:
        if year:
            ted = datetime(year, month, day)
        else:
            ted = datetime.now() - timedelta(days=2)
        tst = Time(np.fix(Time(ted).mjd) - dayspan, format='mjd').datetime
    tsep = datetime.strptime('2019-02-22', "%Y-%m-%d")

    vmaxs = [22.0e4, 8.0e4, 6.0e4, 3.5e4, 2.3e4, 1.8e4, 1.5e4]
    vmins = [-9.0e3, -5.5e3, -3.4e3, -2.5e3, -2.5e3, -2.5e3, -2.5e3]
    plt.ioff()
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0, hspace=0.0, wspace=0.0)

    dateobs = tst
    imgfileslist = []
    while dateobs < ted:
        monthstrdir = dateobs.strftime("%Y/%m/")
        imgoutdir = imgfitstmpdir + monthstrdir
        movieoutdir = pltfigdir + monthstrdir
        for odir in [imgoutdir, movieoutdir]:
            if not os.path.exists(odir):
                os.makedirs(odir)

        # dateobs = datetime.strptime("2019-12-21", "%Y-%m-%d")

        if dateobs > tsep:
            spws = ['0~1', '2~5', '6~10', '11~20', '21~30', '31~43', '44~49']
        else:
            spws = ['1~3', '4~9', '10~16', '17~24', '25~30']
        spw = spws[2:3]
        vmax = vmaxs[2:3]
        vmin = vmins[2:3]
        aiawave = '0304'

        tdateobs = Time(Time(dateobs).mjd + np.arange(0, 24, 1) / 24, format='mjd')
        imgfiles = pltEovsaQlookImageSeries(tdateobs, spw, vmax, vmin, aiawave, fig=fig, axs=axs, overwrite=False,
                                            imgoutdir=imgoutdir)
        imgfileslist = imgfileslist + imgfiles
        dateobs = dateobs + timedelta(days=1)
    moviename = 'eovsa_bd01_aia304_{}'.format(dateobs.strftime("%Y%m"))
    DButil.img2movie(imgprefix=imgfileslist, outname=moviename, img_ext='jpg')
    os.system('mv {} {}'.format(os.path.join(imgoutdir + '../', moviename + '.mp4'),
                                os.path.join(movieoutdir, moviename + '.mp4')))
    plt.close('all')


# if __name__ == '__main__':
#     import sys
#     import numpy as np
#
#     # import subprocess
#     # shell = subprocess.check_output('echo $0', shell=True).decode().replace('\n', '').split('/')[-1]
#     # print("shell " + shell + " is using")
#
#     print(sys.argv)
#     try:
#         argv = sys.argv[1:]
#         if '--clearcache' in argv:
#             clearcache = True
#             argv.remove('--clearcache')  # Allows --clearcache to be either before or after date items
#         else:
#             clearcache = False
#
#         year = np.int(argv[0])
#         month = np.int(argv[1])
#         day = np.int(argv[2])
#         if len(argv) == 3:
#             dayspan = 30
#         else:
#             dayspan = np.int(argv[3])
#     except:
#         print('Error interpreting command line argument')
#         year = None
#         month = None
#         day = None
#         dayspan = 30
#         clearcache = True
#     print("Running pipeline_plt for date {}-{}-{}".format(year, month, day))
#     main(year, month, None, dayspan)
