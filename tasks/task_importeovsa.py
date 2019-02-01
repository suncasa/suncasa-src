import os
import numpy as np
import numpy.ma as ma
import scipy.constants as constants
import time
import aipy
from taskinit import tb, casalog
from split_cli import split_cli as split
from suncasa.eovsa import impteovsa as ipe


def importeovsa_iter(filelist, timebin, width, visprefix, nocreatms, modelms, doscaling, keep_nsclms, fileidx):
    from taskinit import tb, casalog
    filename = filelist[fileidx]
    uv = aipy.miriad.UV(filename)
    # try:
    msname0 = list(filename.split('/')[-1])
    msname = visprefix + ''.join(msname0) + '.ms'
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
    durtim = int((times[-1] - times[0]) / 60 + inttime)
    time0 = time.time()

    if 'antlist' in uv.vartable:
        ants = uv['antlist'].replace('\x00', '')
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    good_idx = np.where(uv['sfreq'] > 0)[0]

    nf = len(good_idx)
    npol = uv['npol']
    nants = uv['nants']
    source_id = uv['source'].replace('\x00', '')
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
        data = ma.masked_array(ma.masked_invalid(data), fill_value=0.0)
        out[k, :, l / (npairs * npol), bl2ord[i0, j0]] = data.data
        flag[k, :, l / (npairs * npol), bl2ord[i0, j0]] = data.mask
        # if i != j:
        if k == 3:
            uvwarray[:, l / (npairs * npol), bl2ord[i0, j0]] = -uvw * constants.speed_of_light / 1e9

    nrows = time_steps * npairs
    if doscaling:
        out2 = out.copy()
        for i0 in antlist:
            for j0 in antlist:
                if i0 < j0:
                    i, j = i0 - 1, j0 - 1
                    out2[:, :, :, bl2ord[i, j]] = out[:, :, :, bl2ord[i, j]] / np.sqrt(
                        np.abs(out[:, :, :, bl2ord[i, i]]) * np.abs(out[:, :, :, bl2ord[j, j]]))
        out2 = out2.reshape(npol, nf, nrows)
        out2[np.isnan(out2)] = 0
        out2[np.isinf(out2)] = 0
    # out2 = ma.masked_array(ma.masked_invalid(out2), fill_value=0.0)
    out = out.reshape(npol, nf, nrows) * 1e4
    flag = flag.reshape(npol, nf, nrows)
    uvwarray = uvwarray.reshape(3, nrows)
    uvwarray = np.tile(uvwarray, (1, nband))
    sigma = np.ones((4, nrows), dtype=np.float) + 1
    sigma = np.tile(sigma, (1, nband))

    casalog.post('IDB File {0} is readed in --- {1:10.2f} seconds ---'.format(filename, (time.time() - time0)))

    if not nocreatms:
        modelms = ipe.creatms(filename, visprefix)
        os.system('mv {} {}'.format(modelms, msname))
    else:
        casalog.post('----------------------------------------')
        casalog.post('copying standard MS to {0}'.format(msname, (time.time() - time0)))
        casalog.post('----------------------------------------')
        os.system("rm -fr %s" % msname)
        os.system("cp -r " + " %s" % modelms + " %s" % msname)
        casalog.post('Standard MS is copied to {0} in --- {1:10.2f} seconds ---'.format(msname, (time.time() - time0)))

    tb.open(msname, nomodify=False)
    casalog.post('----------------------------------------')
    casalog.post("Updating the main table of" '%s' % msname)
    casalog.post('----------------------------------------')
    for l, cband in enumerate(chan_band):
        time1 = time.time()
        nchannels = len(cband['cidx'])
        for row in range(nrows):
            if not doscaling or keep_nsclms:
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
    tb.putcol('TIME_RANGE', np.asarray([times[0] - 0.5 * inttime, times[-1] + 0.5 * inttime]).reshape(2, 1))
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

    # casalog.post('----------------------------------------')
    # casalog.post("Updating the POLARIZATION table of" '%s' % msname)
    # casalog.post('----------------------------------------')
    # tb.open(msname + '/POLARIZATION/', nomodify=False)
    # tb.removerows(rownrs=np.arange(1, nband, dtype=int))
    # tb.close()

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

    # del out, flag, uvwarray, uv, timearr, sigma
    # gc.collect()  #
    if doscaling:
        if keep_nsclms:
            msname_scl = visprefix + ''.join(msname0) + '_scl.ms'
            os.system('cp -r {} {}'.format(msname, msname_scl))
        else:
            msname_scl = msname
        tb.open(msname_scl, nomodify=False)
        casalog.post('----------------------------------------')
        casalog.post("Updating the main table of" '%s' % msname_scl)
        casalog.post('----------------------------------------')
        for l, cband in enumerate(chan_band):
            time1 = time.time()
            for row in range(nrows):
                tb.putcell('DATA', (row + l * nrows), out2[:, cband['cidx'][0]:cband['cidx'][-1] + 1, row])
            casalog.post('---spw {0:02d} is updated in --- {1:10.2f} seconds ---'.format((l + 1), time.time() - time1))
        tb.close()

    if not (timebin == '0s' and width == 1):
        msfile = msname + '.split'
        if doscaling:
            split(vis=msname_scl, outputvis=msname_scl + '.split', datacolumn='data', timebin=timebin, width=width, keepflags=False)
            os.system('rm -rf {}'.format(msname_scl))
            msfile_scl = msname_scl + '.split'
        if not (doscaling and not keep_nsclms):
            split(vis=msname, outputvis=msname + '.split', datacolumn='data', timebin=timebin, width=width, keepflags=False)
            os.system('rm -rf {}'.format(msname))
    else:
        msfile = msname
        if doscaling:
            msfile_scl = msname_scl
    casalog.post("finished in --- %s seconds ---" % (time.time() - time0))
    if doscaling:
        return [True, msfile, msfile_scl, durtim]
    else:
        return [True, msfile, durtim]


def importeovsa(idbfiles=None, ncpu=None, timebin=None, width=None, visprefix=None, udb_corr=True, nocreatms=None, doconcat=None, modelms=None,
                doscaling=False, keep_nsclms=False):
    casalog.origin('importeovsa')

    # if type(idbfiles) == Time:
    #     filelist = ri.get_trange_files(idbfiles)
    # else:
    #     # If input type is not Time, assume that it is the list of files to read
    filelist = idbfiles

    if type(filelist) == str:
        filelist = [filelist]

    filelist_tmp = []
    for ll in filelist:
        if not os.path.exists(ll):
            casalog.post("Warning: {} not exist.".format(ll))
        else:
            filelist_tmp.append(ll)

    filelist = filelist_tmp
    if not filelist:
        casalog.post("No file in idbfiles list exists. Abort.")
        return False

    for idx, ll in enumerate(filelist):
        if ll[-1] == '/':
            filelist[idx] = ll[:-1]

    if not visprefix:
        visprefix = './'
    else:
        if os.path.exists(visprefix):
            pass
        else:
            casalog.post("The output path {} does not exist. Abort.".format(visprefix))
            return False
    if not timebin:
        timebin = '0s'
    if not width:
        width = 1
    import sys
    sys.stdout.flush()
    if udb_corr:
        udbcorr_path = visprefix + '/tmp_UDBcorr/'

        sys.stdout.flush()
        if not os.path.exists(udbcorr_path):
            os.makedirs(udbcorr_path)

        sys.stdout.flush()
        from eovsapy import pipeline_cal as pc

        sys.stdout.flush()
        filelist_tmp = []
        for ll in filelist:
            filelist_tmp.append(pc.udb_corr(ll, outpath=udbcorr_path, calibrate=True))
        filelist = filelist_tmp

        sys.stdout.flush()

    sys.stdout.flush()
    if not modelms:
        if nocreatms:
            filename = filelist[0]
            modelms = ipe.creatms(filename, visprefix)
    else:
        if not os.path.exists(modelms):
            if nocreatms:
                filename = filelist[0]
                modelms = ipe.creatms(filename, visprefix)

    iterable = range(len(filelist))

    t0 = time.time()
    casalog.post('Perform importeovsa in parallel with {} CPUs...'.format(ncpu))

    if ncpu == 1:
        res = []
        for fidx, ll in enumerate(filelist):
            res.append(importeovsa_iter(filelist, timebin, width, visprefix, nocreatms, modelms, doscaling, keep_nsclms, fidx))
    if ncpu > 1:
        import multiprocessing as mprocs
        from functools import partial
        imppart = partial(importeovsa_iter, filelist, timebin, width, visprefix, nocreatms, modelms, doscaling, keep_nsclms)
        pool = mprocs.Pool(ncpu)
        res = pool.map(imppart, iterable)
        pool.close()
        pool.join()

    # print res
    t1 = time.time()
    timelapse = t1 - t0
    print 'It took %f secs to complete' % timelapse

    # results = pd.DataFrame({'succeeded': [], 'msfile': [], 'durtim': []})
    # for r in res:
    #     results = results.append(pd.DataFrame({'succeeded': [r[0]], 'msfile': [r[1]], 'durtim': [r[2]]}))
    # try:
    succeeded = []
    msfile = []
    durtim = []
    if doscaling:
        msfile_scl = []
        for r in res:
            succeeded.append(r[0])
            msfile.append(r[1])
            msfile_scl.append(r[2])
            durtim.append(r[3])
        results = {'succeeded': succeeded, 'msfile': msfile, 'msfile_scl': msfile_scl, 'durtim': durtim}
    else:
        for r in res:
            succeeded.append(r[0])
            msfile.append(r[1])
            durtim.append(r[2])
        results = {'succeeded': succeeded, 'msfile': msfile, 'durtim': durtim}
    # except:
    #     print 'errors occurred when creating the output summary.'

    if doconcat:
        from suncasa.tasks import concateovsa_cli as ce
        msname = os.path.basename(filelist[0])
        durtim = int(np.array(results['durtim']).sum())
        if doscaling:
            msfiles = list(np.array(results['msfile_scl'])[np.where(np.array(results['succeeded']) == True)])
            if keep_nsclms:
                concatvis = visprefix + msname + '-{:d}m{}.ms'.format(durtim, '_scl')
            else:
                concatvis = visprefix + msname + '-{:d}m{}.ms'.format(durtim, '')
        else:
            msfiles = list(np.array(results['msfile'])[np.where(np.array(results['succeeded']) == True)])
            concatvis = visprefix + msname + '-{:d}m{}.ms'.format(durtim, '')
        ce.concateovsa(msfiles, concatvis, datacolumn='data', keep_orig_ms=True, cols2rm="model,corrected")
        return True

    if udb_corr:
        os.system('rm -rf {}'.format(udbcorr_path))
