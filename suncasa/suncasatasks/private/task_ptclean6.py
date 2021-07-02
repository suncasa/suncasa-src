import os
import numpy as np
import shutil
from functools import partial
from time import time
import glob
import sys
from suncasa.utils import helioimage2fits as hf

pversion = sys.version_info.major
if pversion < 3:
    ## CASA version < 6
    from taskinit import ms, tb, qa, casalog
    from tclean_cli import tclean_cli as tclean

    is_casa6 = False
else:
    ## CASA version >= 6
    from casatasks import casalog, tclean
    from casatools import table as tbtool
    from casatools import ms as mstool
    from casatools import quanta as qatool

    tb = tbtool()
    ms = mstool()
    qa = qatool()
    is_casa6 = True


def clean_iter(tim, vis, imageprefix, imagesuffix,
               twidth, doreg, docompress, usephacenter, reftime, ephem, msinfo, toTb, sclfactor, subregion, overwrite,
               selectdata, field, spw, timerange, uvrange, antenna, scan, observation, intent, datacolumn,
               imagename, imsize, cell, phasecenter, stokes, projection, startmodel, specmode, reffreq, nchan,
               start,
               width, outframe, veltype, restfreq, interpolation, perchanweightdensity, gridder, facets,
               psfphasecenter,
               wprojplanes, vptable, mosweight, aterm, psterm, wbawp, conjbeams, cfcache, usepointing,
               computepastep,
               rotatepastep, pointingoffsetsigdev, pblimit, normtype, deconvolver, scales, nterms,
               smallscalebias,
               restoration, restoringbeam, pbcor, outlierfile, weighting, robust, noise, npixels, uvtaper, niter,
               gain,
               threshold, nsigma, cycleniter, cyclefactor, minpsffraction, maxpsffraction, interactive, usemask,
               mask,
               pbmask, sidelobethreshold, noisethreshold, lownoisethreshold, negativethreshold, smoothfactor,
               minbeamfrac,
               cutthreshold, growiterations, dogrowprune, minpercentchange, verbose, fastnoise, restart,
               savemodel,
               calcres, calcpsf, psfcutoff, parallel, btidx):
    bt = btidx  # 0
    if bt + twidth < len(tim) - 1:
        et = btidx + twidth - 1
    else:
        et = len(tim) - 1

    if bt == 0:
        bt_d = tim[bt] - ((tim[bt + 1] - tim[bt]) / 2)
    else:
        bt_d = tim[bt] - ((tim[bt] - tim[bt - 1]) / 2)
    if et == (len(tim) - 1) or et == -1:
        et_d = tim[et] + ((tim[et] - tim[et - 1]) / 2)
    else:
        et_d = tim[et] + ((tim[et + 1] - tim[et]) / 2)

    timerange = qa.time(qa.quantity(bt_d, 's'), prec=9, form='ymd')[0] + '~' + \
                qa.time(qa.quantity(et_d, 's'), prec=9, form='ymd')[0]
    btstr = qa.time(qa.quantity(bt_d, 's'), prec=9, form='fits')[0]
    etstr = qa.time(qa.quantity(et_d, 's'), prec=9, form='fits')[0]
    print('cleaning timerange: ' + timerange)

    image0 = btstr.replace(':', '').replace('-', '')
    imname = imageprefix + image0 + imagesuffix

    if overwrite or (len(glob.glob(imname + '*')) == 0):
        os.system('rm -rf {}*'.format(imname))
        # try:
        tclean(vis=vis, selectdata=selectdata, field=field, spw=spw, timerange=timerange, uvrange=uvrange,
               antenna=antenna, scan=scan, observation=observation, intent=intent, datacolumn=datacolumn,
               imagename=imname, imsize=imsize, cell=cell, phasecenter=phasecenter, stokes=stokes,
               projection=projection, startmodel=startmodel, specmode=specmode, reffreq=reffreq, nchan=nchan,
               start=start, width=width, outframe=outframe, veltype=veltype, restfreq=restfreq,
               interpolation=interpolation, perchanweightdensity=perchanweightdensity, gridder=gridder, facets=facets,
               psfphasecenter=psfphasecenter, wprojplanes=wprojplanes, vptable=vptable, mosweight=mosweight,
               aterm=aterm, psterm=psterm, wbawp=wbawp, conjbeams=conjbeams, cfcache=cfcache, usepointing=usepointing,
               computepastep=computepastep, rotatepastep=rotatepastep, pointingoffsetsigdev=pointingoffsetsigdev,
               pblimit=pblimit, normtype=normtype, deconvolver=deconvolver, scales=scales, nterms=nterms,
               smallscalebias=smallscalebias, restoration=restoration, restoringbeam=restoringbeam, pbcor=pbcor,
               outlierfile=outlierfile, weighting=weighting, robust=robust, noise=noise, npixels=npixels,
               uvtaper=uvtaper, niter=niter, gain=gain, threshold=threshold, nsigma=nsigma, cycleniter=cycleniter,
               cyclefactor=cyclefactor, minpsffraction=minpsffraction, maxpsffraction=maxpsffraction,
               interactive=interactive, usemask=usemask, mask=mask, pbmask=pbmask, sidelobethreshold=sidelobethreshold,
               noisethreshold=noisethreshold, lownoisethreshold=lownoisethreshold, negativethreshold=negativethreshold,
               smoothfactor=smoothfactor, minbeamfrac=minbeamfrac, cutthreshold=cutthreshold,
               growiterations=growiterations, dogrowprune=dogrowprune, minpercentchange=minpercentchange,
               verbose=verbose, fastnoise=fastnoise, restart=restart, savemodel=savemodel, calcres=calcres,
               calcpsf=calcpsf, psfcutoff=psfcutoff, parallel=parallel)
        if pbcor:
            clnjunks = ['.flux', '.mask', '.model', '.psf', '.residual', '.pb', '.sumwt', '.image']
        else:
            clnjunks = ['.flux', '.mask', '.model', '.psf', '.residual', '.pb', '.sumwt', '.image.pbcor']
        for clnjunk in clnjunks:
            if os.path.exists(imname + clnjunk):
                shutil.rmtree(imname + clnjunk)
        if pbcor:
            os.system('mv {} {}'.format(imname + '.image.pbcor', imname + '.image'))
        # except:
        #     print('error in cleaning image: ' + btstr)
        #     return [False, btstr, etstr, '']
    else:
        print(imname + ' exists. Clean task aborted.')

    if doreg:
        # ephem.keys()
        # msinfo.keys()
        if os.path.isfile(imname + '.fits'):
            return [True, btstr, etstr, imname + '.fits']
        else:
            try:
                # check if ephemfile and msinfofile exist
                if not ephem:
                    print("ephemeris info does not exist, querying from JPL Horizons on the fly")
                    ephem = hf.read_horizons(vis=vis)
                if not msinfo:
                    print("ms info not provided, generating one on the fly")
                    msinfo = hf.read_msinfo(vis)
                hf.imreg(vis=vis, ephem=ephem, msinfo=msinfo, timerange=timerange, reftime=reftime,
                         imagefile=imname + '.image', fitsfile=imname + '.fits', overwrite=True,
                         toTb=toTb, sclfactor=sclfactor, usephacenter=usephacenter, subregion=subregion,
                         docompress=docompress)
                if os.path.exists(imname + '.fits'):
                    shutil.rmtree(imname + '.image')
                    return [True, btstr, etstr, imname + '.fits']
                else:
                    return [False, btstr, etstr, imname + '.fits']
            except Exception as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                print('error in registering image: ' + btstr)
                return [False, btstr, etstr, imname + '.image']
    else:
        if os.path.exists(imname + '.image'):
            return [True, btstr, etstr, imname + '.image']
        else:
            return [False, btstr, etstr, imname + '.image']


def ptclean6(vis, imageprefix, imagesuffix, ncpu, twidth, doreg, usephacenter, reftime, toTb, sclfactor, subregion,
             docompress, overwrite,
             selectdata, field, spw, timerange, uvrange, antenna, scan, observation, intent, datacolumn,
             imagename, imsize, cell, phasecenter, stokes, projection, startmodel, specmode, reffreq, nchan, start,
             width, outframe, veltype, restfreq, interpolation, perchanweightdensity, gridder, facets, psfphasecenter,
             wprojplanes, vptable, mosweight, aterm, psterm, wbawp, conjbeams, cfcache, usepointing, computepastep,
             rotatepastep, pointingoffsetsigdev, pblimit, normtype, deconvolver, scales, nterms, smallscalebias,
             restoration, restoringbeam, pbcor, outlierfile, weighting, robust, noise, npixels, uvtaper, niter, gain,
             threshold, nsigma, cycleniter, cyclefactor, minpsffraction, maxpsffraction, interactive, usemask, mask,
             pbmask, sidelobethreshold, noisethreshold, lownoisethreshold, negativethreshold, smoothfactor, minbeamfrac,
             cutthreshold, growiterations, dogrowprune, minpercentchange, verbose, fastnoise, restart, savemodel,
             calcres, calcpsf, psfcutoff, parallel):
    if not (type(ncpu) is int):
        casalog.post('ncpu should be an integer')
        ncpu = 1

    if doreg:
        # check if ephem and msinfo exist. If not, generate one on the fly
        try:
            ephem = hf.read_horizons(vis=vis)
        except ValueError:
            print("error in obtaining ephemeris")
        try:
            msinfo = hf.read_msinfo(vis)
        except ValueError:
            print("error in getting ms info")
    else:
        ephem = None
        msinfo = None

    if imageprefix:
        workdir = os.path.dirname(imageprefix)
    else:
        workdir = './'
    # get number of time pixels
    ms.open(vis)
    ms.selectinit()
    timfreq = ms.getdata(['time', 'axis_info'], ifraxis=True)
    tim = timfreq['time']
    ms.close()

    if twidth < 1:
        casalog.post('twidth less than 1. Change to 1')
        twidth = 1

    if twidth > len(tim):
        casalog.post(
            'twidth greater than # of time pixels in the dataset. Change to the timerange of the entire dateset')
        twidth = len(tim)
    # find out the start and end time index according to the parameter timerange
    # if not defined (empty string), use start and end from the entire time of the ms

    if not timerange:
        btidx = 0
        etidx = len(tim) - 1
    else:
        try:
            (tstart, tend) = timerange.split('~')
            bt_s = qa.convert(qa.quantity(tstart, 's'), 's')['value']
            et_s = qa.convert(qa.quantity(tend, 's'), 's')['value']
            # only time is given but not date, add the date (at 0 UT) from the first record
            if bt_s < 86400. or et_s < 86400.:
                bt_s += np.fix(qa.convert(qa.quantity(tim[0], 's'), 'd')['value']) * 86400.
                et_s += np.fix(qa.convert(qa.quantity(tim[0], 's'), 'd')['value']) * 86400.
            btidx = np.argmin(np.abs(tim - bt_s))
            etidx = np.argmin(np.abs(tim - et_s))
            # make the indice back to those bracket by the timerange
            if tim[btidx] < bt_s:
                btidx += 1
            if tim[etidx] > et_s:
                etidx -= 1
            if etidx <= btidx:
                print("ending time must be greater than starting time")
                print("reinitiating to the entire time range")
                btidx = 0
                etidx = len(tim) - 1
        except ValueError:
            print("keyword 'timerange' has a wrong format")


    btstr = qa.time(qa.quantity(tim[btidx], 's'), prec=9, form='fits')[0]
    etstr = qa.time(qa.quantity(tim[etidx], 's'), prec=9, form='fits')[0]

    iterable = range(btidx, etidx + 1, twidth)
    print('First time pixel: ' + btstr)
    print('Last time pixel: ' + etstr)
    print(str(len(iterable)) + ' images to clean...')

    res = []
    # partition
    clnpart = partial(clean_iter, tim, vis, imageprefix, imagesuffix,
                      twidth, doreg, docompress, usephacenter, reftime, ephem, msinfo, toTb, sclfactor, subregion,
                      overwrite,
                      selectdata, field, spw, timerange, uvrange, antenna, scan, observation, intent, datacolumn,
                      imagename, imsize, cell, phasecenter, stokes, projection, startmodel, specmode, reffreq, nchan,
                      start,
                      width, outframe, veltype, restfreq, interpolation, perchanweightdensity, gridder, facets,
                      psfphasecenter,
                      wprojplanes, vptable, mosweight, aterm, psterm, wbawp, conjbeams, cfcache, usepointing,
                      computepastep,
                      rotatepastep, pointingoffsetsigdev, pblimit, normtype, deconvolver, scales, nterms,
                      smallscalebias,
                      restoration, restoringbeam, pbcor, outlierfile, weighting, robust, noise, npixels, uvtaper, niter,
                      gain,
                      threshold, nsigma, cycleniter, cyclefactor, minpsffraction, maxpsffraction, interactive, usemask,
                      mask,
                      pbmask, sidelobethreshold, noisethreshold, lownoisethreshold, negativethreshold, smoothfactor,
                      minbeamfrac,
                      cutthreshold, growiterations, dogrowprune, minpercentchange, verbose, fastnoise, restart,
                      savemodel,
                      calcres, calcpsf, psfcutoff, parallel)
    timelapse = 0
    t0 = time()
    # parallelization
    if ncpu > 1:
        import multiprocessing as mprocs
        casalog.post('Perform clean in parallel ...')
        print('Perform clean in parallel ...')
        pool = mprocs.Pool(ncpu)
        res = pool.map(clnpart, iterable)
        pool.close()
        pool.join()
    else:
        casalog.post('Perform clean in single process ...')
        print('Perform clean in single process ...')
        for i in iterable:
            res.append(clnpart(i))

    t1 = time()
    timelapse = t1 - t0
    print('It took %f secs to complete' % timelapse)
    # repackage this into a single dictionary
    results = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': []}
    for r in res:
        results['Succeeded'].append(r[0])
        results['BeginTime'].append(r[1])
        results['EndTime'].append(r[2])
        results['ImageName'].append(r[3])

    return results
