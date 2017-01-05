import os
# import sys
import gc
import numpy as np
# import matplotlib.pyplot as plt
import scipy.constants as constants
import time
# import glob
import aipy
import eovsapy.chan_util_bc as chan_util_bc
import eovsapy.read_idb as ri
from eovsapy.util import Time
from taskinit import *
from split_cli import split_cli as split
from concat_cli import concat_cli as concat


# from parallel.parallel_data_helper import ParallelDataHelper

def bl_list2(nant=16):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = -1.
    '''
    bl2ord = np.ones((nant, nant), dtype='int') * (-1)
    k = 0
    for i in range(nant):
        for j in range(i, nant):
            bl2ord[i, j] = k
            # bl2ord[j,i] = k
            k += 1
    return bl2ord


def get_band_edge(nband=34):
    # Input the frequencies from UV, returen the indices frequency edges of all bands
    idx_start_freq = [0]
    ntmp = 0
    for i in range(1, nband + 1):
        ntmp += len(chan_util_bc.start_freq(i))
        idx_start_freq.append(ntmp)
    return np.asarray(idx_start_freq)


def creatms(idbfile, outpath, timebin=None, width=None):
    uv = aipy.miriad.UV(idbfile)
    uv.rewind()

    start_time = 0  # The start and stop times are referenced to ref_time_jd in second
    end_time = 600
    time0 = time.time()

    if 'antlist' in uv.vartable:
        ants = uv['antlist']
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    good_idx = np.arange(len(uv['sfreq']))

    ref_time_jd = uv['time']
    freq = uv['sfreq'][good_idx]
    sdf = uv['sdf']
    project = uv['proj']
    source_id = uv['source']
    bandedge = get_band_edge(nband=34)
    msname = list(idbfile.split('/')[-1])
    msname = outpath + ''.join(msname) + '-10m.ms'

    if os.path.exists(msname):
        os.system("rm -fr %s" % msname)

    """ Creates an empty measurement set using CASA simulate (sm) tool. """
    sm.open(msname)

    enu = np.reshape(uv['antpos'], (16, 3)) * constants.speed_of_light / 1e9
    refpos_wgs84 = me.position('wgs84',
                               '-118.286952892965deg',
                               '37.2331698901026deg',
                               '1207.1339m')
    lon, lat, rad = [me.measure(refpos_wgs84, 'itrf')[x]['value'] for x in 'm0', 'm1', 'm2']
    # 3x3 transform matrix. Each row is a normal vector, i.e. the rows are (dE,dN,dU)
    # ----------- local xyz ------------
    xform = np.array([
        [0, -np.sin(lat), np.cos(lat)],
        [1, 0, 0],
        [0, np.cos(lat), np.sin(lat)]])
    xyz = enu.dot(xform)  # + xyz0[np.newaxis,:]

    # ----------- global xyz ------------
    # xyz0 = rad*np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    # # 3x3 transform matrix. Each row is a normal vector, i.e. the rows are (dE,dN,dU)
    # xform = np.array([
    #     [-np.sin(lon),np.cos(lon),0],
    #     [-np.cos(lon)*np.sin(lat),-np.sin(lon)*np.sin(lat),np.cos(lat)],
    #     [np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
    # ])
    # xyz = xyz0[np.newaxis,:] + enu.dot(xform)

    dishdiam = np.asarray([2.1] * uv['nants'])
    dishdiam[-3:-1] = 27
    dishdiam[-1] = 2.1
    station = uv['telescop']
    mount = ['ALT-AZ'] * uv['nants']
    for l in [8, 9, 10, 12, 13, 14]:
        mount[l] = 'EQUATORIAL'
    sm.setconfig(telescopename=station,
                 x=np.asarray(xyz)[:, 0],
                 y=np.asarray(xyz)[:, 1],
                 z=np.asarray(xyz)[:, 2],
                 dishdiameter=dishdiam,
                 mount=mount,
                 antname=['eo' + "{0:02d}".format(l) for l in antlist],
                 padname=station,
                 coordsystem='local', referencelocation=refpos_wgs84)

    sm.setfield(sourcename=source_id,
                sourcedirection=me.direction('J2000',
                                             '{:22.19f}'.format(uv['obsra']) + 'rad',
                                             '{:22.19f}'.format(uv['obsdec']) + 'rad'))
    sm.setfeed(mode='perfect X Y')

    ref_time = me.epoch('tai',
                        '{:20.13f}'.format(ref_time_jd - 2400000.5) + 'd')

    sm.settimes(integrationtime='1s',
                usehourangle=False,
                referencetime=ref_time)

    for l, bdedge in enumerate(bandedge[:-1]):
        nchannels = (bandedge[l + 1] - bandedge[l])
        stokes = 'XX YY XY YX'
        df = sdf[bandedge[l]]
        st_freq = freq[bandedge[l]]

        sm.setspwindow(spwname='band%02d' % (l + 1),
                       freq='{:22.19f}'.format(st_freq) + 'GHz',
                       deltafreq='{:22.19f}'.format(df) + 'GHz',
                       freqresolution='{:22.19f}'.format(df) + 'GHz',
                       nchannels=nchannels,
                       stokes=stokes)

    nband = len(bandedge) - 1
    for bdid in range(nband):
        sm.observe(source_id, 'band%02d' % (bdid + 1),
                   starttime=start_time, stoptime=end_time,
                   project=project,
                   state_obs_mode='')

    if sm.done():
        casalog.post('Empty MS {0} created in --- {1:10.2f} seconds ---'.format(msname, (time.time() - time0)))
    else:
        raise RuntimeError('Failed to create MS. Look at the log file. '
                           'Double check you settings.')

    modelms = msname + '.MSmodel'
    os.system('mv {} {}'.format(msname, modelms))

    # if timebin != '0s' or width != 1:
    #     modelms = msname + '.MSmodel'
    #     split(vis=msname, outputvis=modelms, datacolumn='data', timebin=timebin, width=width)
    #     os.system('rm -rf {}'.format(msname))

    return modelms


def importeovsa(idbfiles, timebin=None, width=None, visprefix=None, nocreatms=True, doconcat=False, modelms=None):
    casalog.origin('importeovsa')

    # # Initialize the helper class
    # pdh = ParallelDataHelper("importeovsa", locals())
    #
    # # Validate input and output parameters
    # try:
    #     pdh.setupIO()
    # except Exception, instance:
    #     casalog.post('%s' % instance, 'ERROR')
    #     return False


    if type(idbfiles) == Time:
        filelist = ri.get_trange_files(idbfiles)
    else:
        # If input type is not Time, assume that it is the list of files to read
        filelist = idbfiles

    if type(filelist) == str:
        filelist = [filelist]

    for f in filelist:
        if not os.path.exists(f):
            casalog.post("Some files in filelist are invalid. Aborting...")
            return False
    if not visprefix:
        visprefix = './'
    if not timebin:
        timebin = '0s'
    if not width:
        width = 1

    if not modelms:
        if nocreatms:
            filename = filelist[0]
            modelms = creatms(filename, visprefix)
    else:
        if not os.path.exists(modelms):
            if nocreatms:
                filename = filelist[0]
                modelms = creatms(filename, visprefix)

    msfile = []
    for filename in filelist:
        uv = aipy.miriad.UV(filename)
        # if uv['source'].lower() == 'sun':
        #     visprefix = visprefix + 'sun/'
        #     if not os.path.exists(visprefix):
        #         os.mkdir(visprefix)
        # else:
        #     visprefix = visprefix + 'calibrator/'
        #     if not os.path.exists(visprefix):
        #         os.mkdir(visprefix)
        # uv.rewind()

        start_time = 0  # The start and stop times are referenced to ref_time_jd in second
        end_time = 600
        delta_time = 1
        time_steps = (end_time - start_time) / delta_time
        time0 = time.time()

        if 'antlist' in uv.vartable:
            ants = uv['antlist']
            antlist = map(int, ants.split())
        else:
            antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

        good_idx = np.arange(len(uv['sfreq']))

        ref_time_jd = uv['time']
        ref_time_mjd = (ref_time_jd - 2400000.5) * 24. * 3600. + 0.5 * delta_time
        nf = len(good_idx)
        npol = uv['npol']
        nants = uv['nants']
        source_id = uv['source']
        ra, dec = uv['ra'], uv['dec']
        nbl = nants * (nants - 1) / 2
        bl2ord = bl_list2(nants)
        npairs = nbl + nants
        flag = np.ones((npol, nf, time_steps, npairs), dtype=bool)
        out = np.zeros((npol, nf, time_steps, npairs), dtype=np.complex64)  # Cross-correlations
        uvwarray = np.zeros((3, time_steps, npairs), dtype=np.float)
        bandedge = get_band_edge(nband=34)
        nband = len(bandedge) - 1

        uv.rewind()
        l = -1
        for preamble, data in uv.all():
            uvw, t, (i0, j0) = preamble
            i = antlist.index(i0 + 1)
            j = antlist.index(j0 + 1)
            if i > j:
                # Reverse order of indices
                j = antlist.index(i0 + 1)
                i = antlist.index(j0 + 1)
            # Assumes uv['pol'] is one of -5, -6, -7, -8
            k = -5 - uv['pol']
            l += 1
            out[k, :, l / (npairs * npol), bl2ord[i0, j0]] = data.data
            flag[k, :, l / (npairs * npol), bl2ord[i0, j0]] = data.mask
            if i != j:
                if k == 3:
                    uvwarray[:, l / (npairs * npol), bl2ord[i0, j0]] = -uvw * constants.speed_of_light / 1e9

        nrows = time_steps * npairs
        out = out.reshape(npol, nf, nrows)
        flag = flag.reshape(npol, nf, nrows)
        uvwarray = uvwarray.reshape(3, nrows)
        uvwarray = np.tile(uvwarray, (1, nband))
        sigma = np.ones((4, nrows), dtype=np.float) + 1
        sigma = np.tile(sigma, (1, nband))

        casalog.post('IDB File {0} is readed in --- {1:10.2f} seconds ---'.format(filename, (time.time() - time0)))

        msname = list(filename.split('/')[-1])
        msname = visprefix + ''.join(msname) + '-10m.ms'

        if not nocreatms:
            modelms = creatms(filename, visprefix)
            os.system('mv {} {}'.format(modelms, msname))
        else:
            casalog.post('----------------------------------------')
            casalog.post('copying standard MS to {0}'.format(msname, (time.time() - time0)))
            casalog.post('----------------------------------------')
            os.system("rm -fr %s" % msname)
            os.system("cp -r " + " %s" % modelms + " %s" % msname)
            casalog.post(
                'Standard MS is copied to {0} in --- {1:10.2f} seconds ---'.format(msname, (time.time() - time0)))

        tb.open(msname, nomodify=False)
        casalog.post('----------------------------------------')
        casalog.post("Updating the main table of" '%s' % msname)
        casalog.post('----------------------------------------')
        for l, bdedge in enumerate(bandedge[:-1]):
            time1 = time.time()
            nchannels = (bandedge[l + 1] - bandedge[l])
            for row in range(nrows):
                tb.putcell('DATA', (row + l * nrows), out[:, bandedge[l]:bandedge[l + 1], row])
                tb.putcell('FLAG', (row + l * nrows), flag[:, bandedge[l]:bandedge[l + 1], row])
            casalog.post('---spw {0:02d} is updated in --- {1:10.2f} seconds ---'.format((l + 1), time.time() - time1))
        tb.putcol('UVW', uvwarray)
        tb.putcol('SIGMA', sigma)
        tb.putcol('WEIGHT', 1.0 / sigma ** 2)
        timearr = np.arange((time_steps), dtype=np.float)
        timearr = timearr.reshape(1, time_steps, 1)
        timearr = np.tile(timearr, (nband, 1, npairs))
        timearr = timearr.reshape(nband * npairs * time_steps) + ref_time_mjd
        tb.putcol('TIME', timearr)
        tb.putcol('TIME_CENTROID', timearr)
        scan_id = tb.getcol('SCAN_NUMBER')
        scan_id *= 0
        tb.putcol('SCAN_NUMBER', scan_id)
        colnames = tb.colnames()
        cols2rm = ["MODEL_DATA", "CORRECTED_DATA"]
        for l in range(len(cols2rm)):
            if cols2rm[l] in colnames:
                tb.removecols(cols2rm[l])
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the OBSERVATION table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/OBSERVATION', nomodify=False)
        tb.putcol('TIME_RANGE',
                  np.asarray([ref_time_mjd - 0.5 * delta_time, ref_time_mjd + end_time - 0.5 * delta_time]).reshape(
                      2, 1))
        tb.putcol('OBSERVER', ['EOVSA team'])
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the POINTING table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/POINTING', nomodify=False)
        timearr = np.arange((time_steps), dtype=np.float).reshape(1, time_steps, 1)
        timearr = np.tile(timearr, (nband, 1, nants))
        timearr = timearr.reshape(nband * time_steps * nants) + ref_time_mjd
        tb.putcol('TIME', timearr)
        tb.putcol('TIME_ORIGIN', timearr - 0.5 * delta_time)
        direction = tb.getcol('DIRECTION')
        direction[0, 0, :] = ra
        direction[1, 0, :] = dec
        tb.putcol('DIRECTION', direction)
        target = tb.getcol('TARGET')
        target[0, 0, :] = ra
        target[1, 0, :] = dec
        tb.putcol('TARGET', target)
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the SOURCE table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/SOURCE', nomodify=False)
        radec = tb.getcol('DIRECTION')
        radec[0], radec[1] = ra, dec
        tb.putcol('DIRECTION', radec)
        name = np.array([source_id], dtype='|S{0}'.format(len(source_id) + 1))
        tb.putcol('NAME', name)
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the DATA_DESCRIPTION table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/DATA_DESCRIPTION/', nomodify=False)
        pol_id = tb.getcol('POLARIZATION_ID')
        pol_id *= 0
        tb.putcol('POLARIZATION_ID', pol_id)
        # spw_id = tb.getcol('SPECTRAL_WINDOW_ID')
        # spw_id *= 0
        # tb.putcol('SPECTRAL_WINDOW_ID', spw_id)
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the POLARIZATION table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/POLARIZATION/', nomodify=False)
        tb.removerows(rownrs=np.arange(1, nband, dtype=int))
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the FIELD table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/FIELD/', nomodify=False)
        delay_dir = tb.getcol('DELAY_DIR')
        delay_dir[0], delay_dir[1] = ra, dec
        tb.putcol('DELAY_DIR', delay_dir)
        phase_dir = tb.getcol('PHASE_DIR')
        phase_dir[0], phase_dir[1] = ra, dec
        tb.putcol('PHASE_DIR', phase_dir)
        reference_dir = tb.getcol('REFERENCE_DIR')
        reference_dir[0], reference_dir[1] = ra, dec
        tb.putcol('REFERENCE_DIR', reference_dir)
        name = np.array([source_id], dtype='|S{0}'.format(len(source_id) + 1))
        tb.putcol('NAME', name)
        tb.close()

        # FIELD: DELAY_DIR, PHASE_DIR, REFERENCE_DIR, NAME


        del out, flag, uvwarray, uv, timearr, sigma
        gc.collect()  #

        if not (timebin == '0s' and width == 1):
            split(vis=msname, outputvis=msname + '.split', datacolumn='data', timebin=timebin, width=width)
            os.system('rm -rf {}'.format(msname))
            msfile.append(msname + '.split')
        else:
            msfile.append(msname)
        casalog.post("finished in --- %s seconds ---" % (time.time() - time0))

    if doconcat:
        msname = list(filelist[0].split('/')[-1])
        concatvis = visprefix + ''.join(msname) + '-{:d}m.ms'.format(10 * len(msfile))
        concat(vis=msfile, concatvis=concatvis)
        for ll in msfile:
            os.system('rm -rf {}'.format(ll))

        return True
