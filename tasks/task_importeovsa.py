import os
import gc
import numpy as np
import scipy.constants as constants
import time
import aipy
import eovsapy.read_idb as ri
from eovsapy.util import Time
from taskinit import tb, casalog
from split_cli import split_cli as split
from concat_cli import concat_cli as concat
from suncasa.utils import impteovsa as ipe



def importeovsa(idbfiles, timebin=None, width=None, visprefix=None, nocreatms=False, doconcat=False, modelms=''):
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
            modelms = ipe.creatms(filename, visprefix)
    else:
        if not os.path.exists(modelms):
            if nocreatms:
                filename = filelist[0]
                modelms = ipe.creatms(filename, visprefix)

    msfile = []
    time_concat = []
    for filename in filelist:
        uv = aipy.miriad.UV(filename)
        # if filename.split('/')[-1][0:3] == 'UDB':
        #     uv_str = uv_hex_rm(uv)
        # else:
        #     uv_str = uv
        uv.select('antennae', 0, 1, include=True)
        uv.select('polarization', -5, -5, include=True)
        times = []
        uv.rewind()
        for preamble, data in uv.all():
            uvw, t, (i, j) = preamble
            times.append(t)

        uv.select('clear', -1, -1, include=True)
        times = ipe.jd2mjds(np.asarray(times))
        inttime = np.median((times - np.roll(times, 1))[1:]) / 60

        time_steps = len(times)
        time_concat.append(int((times[-1] - times[0]) / 60 + inttime))
        time0 = time.time()

        if 'antlist' in uv.vartable:
            ants = uv['antlist'].replace('\x00','')
            antlist = map(int, ants.split())
        else:
            antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

        good_idx = np.where(uv['sfreq'] > 0)[0]

        nf = len(good_idx)
        npol = uv['npol']
        nants = uv['nants']
        source_id = uv['source'].replace('\x00','')
        sfreq = uv['sfreq'][good_idx]
        sdf = uv['sdf'][good_idx]
        ra, dec = uv['ra'], uv['dec']
        nbl = nants * (nants - 1) / 2
        bl2ord = ipe.bl_list2(nants)
        npairs = nbl + nants
        flag = np.ones((npol, nf, time_steps, npairs), dtype=bool)
        out = np.zeros((npol, nf, time_steps, npairs), dtype=np.complex64)  # Cross-correlations
        uvwarray = np.zeros((3, time_steps, npairs), dtype=np.float)
        chan_band = ipe.get_band(sfreq=sfreq, sdf=sdf)
        nband = len(chan_band)

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
            # if i != j:
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
        msname = visprefix + ''.join(msname) + '.ms'

        if not nocreatms:
            modelms = ipe.creatms(filename, visprefix)
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
        for l, cband in enumerate(chan_band):
            time1 = time.time()
            nchannels = len(cband['cidx'])
            for row in range(nrows):
                tb.putcell('DATA', (row + l * nrows), out[:, cband['cidx'][0]:cband['cidx'][-1] + 1, row])
                tb.putcell('FLAG', (row + l * nrows), flag[:, cband['cidx'][0]:cband['cidx'][-1] + 1, row])
            casalog.post('---spw {0:02d} is updated in --- {1:10.2f} seconds ---'.format((l + 1), time.time() - time1))
        tb.putcol('UVW', uvwarray)
        tb.putcol('SIGMA', sigma)
        tb.putcol('WEIGHT', 1.0 / sigma ** 2)
        timearr = times
        timearr = timearr.reshape(1, time_steps, 1)
        timearr = np.tile(timearr, (nband, 1, npairs))
        timearr = timearr.reshape(nband * npairs * time_steps)
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
                  np.asarray([times[0] - 0.5 * inttime, times[-1] + 0.5 * inttime]).reshape(
                      2, 1))
        tb.putcol('OBSERVER', ['EOVSA team'])
        tb.close()

        casalog.post('----------------------------------------')
        casalog.post("Updating the POINTING table of" '%s' % msname)
        casalog.post('----------------------------------------')
        tb.open(msname + '/POINTING', nomodify=False)
        timearr = times.reshape(1, time_steps, 1)
        timearr = np.tile(timearr, (nband, 1, nants))
        timearr = timearr.reshape(nband * time_steps * nants)
        tb.putcol('TIME', timearr)
        tb.putcol('TIME_ORIGIN', timearr)  # - 0.5 * delta_time)
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
            split(vis=msname, outputvis=msname + '.split', datacolumn='data', timebin=timebin, width=width,
                  keepflags=False)
            os.system('rm -rf {}'.format(msname))
            msfile.append(msname + '.split')
        else:
            msfile.append(msname)
        casalog.post("finished in --- %s seconds ---" % (time.time() - time0))

    if doconcat:
        msname = list(filelist[0].split('/')[-1])
        concatvis = visprefix + ''.join(msname) + '-{:d}m.ms'.format(int(sum(time_concat)))
        concat(vis=msfile, concatvis=concatvis, timesort=True)
        # Change all observation ids to be the same (zero)
        tb.open(concatvis+'/OBSERVATION',nomodify=False)
        nobs=tb.nrows()
        tb.removerows([i+1 for i in range(nobs-1)])
        tb.close()
        tb.open(concatvis,nomodify=False)
        obsid=tb.getcol('OBSERVATION_ID')
        newobsid=np.zeros(len(obsid),dtype='int')
        tb.putcol('OBSERVATION_ID',newobsid)
        tb.close()
        for ll in msfile:
            os.system('rm -rf {}'.format(ll))

        return True
