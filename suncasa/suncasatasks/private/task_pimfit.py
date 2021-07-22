import numpy as np
import numpy.ma as ma
import shutil
import os, struct
import bisect
from math import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.gridspec as gs
from fnmatch import fnmatch
from time import sleep
from time import time
from taskinit import casalog
import multiprocessing as mprocs
from suncasa.utils import DButil


def imfit_iter(imgfiles, doreg, tims, msinfofile, ephem, box, region, chans, stokes, mask, includepix, excludepix,
               residual, model, estimates, logfile, append, newestimates, complist,
               overwrite, dooff, offset, fixoffset, stretch, rms, noisefwhm, summary,
               imidx):
    from taskinit import iatool,rgtool
    import pdb
    try:
        from astropy.io import fits as pyfits
    except:
        try:
            import pyfits
        except ImportError:
            raise ImportError('Neither astropy nor pyfits exists in this CASA installation')
    img = imgfiles[imidx]

    if doreg:
        # check if ephemfile and msinfofile exist
        if not ephem:
            print("ephemeris info does not exist!")
            return
        if not tims:
            print("timestamp of the image does not exist!")
            return
        try:
            tim = tims[imidx]
            helio = vla_prep.ephem_to_helio(msinfo=msinfofile, ephem=ephem, reftime=tim)
            fitsfile = [img.replace('.image', '.fits')]
            vla_prep.imreg(imagefile=[img], fitsfile=fitsfile, helio=helio, toTb=False, scl100=True)
            img = img.replace('.image', '.fits')
        except:
            print 'Failure in vla_prep. Skipping this image file: ' + img

    myia = iatool()
    myrg = rgtool()
    try:
        if (not myia.open(img)):
            raise Exception, "Cannot create image analysis tool using " + img
        print('Processing image: ' + img)
        hdr = pyfits.getheader(img)
        pols = DButil.polsfromfitsheader(hdr)
        ndx, ndy, nchans, npols = myia.shape()
        results = {}
        for itpp in pols:
            results[itpp] = {}
        for pp, itpp in enumerate(pols):
            r = myrg.box(blc=[0, 0, 0, pp], trc=[ndx, ndy, nchans, pp])
            myiapol = myia.subimage(region=r, dropdeg=True)
            results[itpp] = myiapol.fitcomponents(
                box=box, region=region, chans=chans, stokes=stokes,
                mask=mask, includepix=includepix,
                excludepix=excludepix, residual=residual,
                model=model, estimates=estimates, logfile=logfile,
                append=append, newestimates=newestimates,
                complist=complist, overwrite=overwrite, dooff=dooff,
                offset=offset, fixoffset=fixoffset, stretch=stretch,
                rms=rms, noisefwhm=noisefwhm, summary=summary
            )
        # update timestamp

        timstr = hdr['date-obs']
        return [True, timstr, img, results]
    except Exception, instance:
        casalog.post(str('*** Error in imfit ***') + str(instance))
        # raise instance
        return [False, timstr, img, {}]
    finally:
        myia.done()


def pimfit(imagefiles, ncpu, doreg, timestamps, msinfofile, ephemfile, box, region, chans, stokes, mask, includepix,
           excludepix,
           residual, model, estimates, logfile, append, newestimates, complist,
           overwrite, dooff, offset, fixoffset, stretch, rms, noisefwhm, summary):
    from functools import partial
    import pdb
    # check if imagefiles is a single file or a list of files

    if isinstance(imagefiles, str):
        imagefiles = [imagefiles]
    if (not isinstance(imagefiles, list)):
        print 'input "imagefiles" is not a list. Abort...'

    if doreg:
        # check if ephemfile and msinfofile exist
        try:
            ephem = vla_prep.read_horizons(ephemfile=ephemfile)
        except ValueError:
            print("error in reading ephemeris file")
        if not os.path.isfile(msinfofile):
            print("msinfofile does not exist!")
        # check if timestamps exist 
        if (not isinstance(timestamps, list)):
            print('input "timestamps" is not a list. Abort...')
        if (len(timestamps) != len(imagefiles)):
            print('imagefiles and timestamps should have the same length! Abort...')
    else:
        ephem = None

    # check file existence
    imgfiles = []
    tims = []
    for i, img in enumerate(imagefiles):
        if os.path.exists(img):
            imgfiles.append(img)
            if doreg:
                tims.append(timestamps[i])
        else:
            casalog.post(img + 'does not exist. Skipping this one...')

    iterable = range(len(imgfiles))
    res = []

    if not (type(ncpu) is int):
        casalog.post('ncpu should be an integer')
        ncpu = 8

    # partition
    imfit_part = partial(imfit_iter, imgfiles, doreg, tims, msinfofile, ephem, box, region, chans, stokes, mask,
                         includepix, excludepix, \
                         residual, model, estimates, logfile, append, newestimates, complist, \
                         overwrite, dooff, offset, fixoffset, stretch, rms, noisefwhm, summary)

    # parallelization
    para = 1
    timelapse = 0
    t0 = time()
    if para:
        casalog.post('Perform imfit in parallel ...')
        pool = mprocs.Pool(ncpu)
        # res = pool.map_async(imfit_part, iterable)
        res = pool.map(imfit_part, iterable)
        pool.close()
        pool.join()
    else:
        for i in iterable:
            res.append(imfit_part(i))

    t1 = time()
    timelapse = t1 - t0
    print 'It took %f secs to complete' % timelapse
    # repackage this into a single dictionary
    results = {'succeeded': [], 'timestamps': [], 'imagenames': [], 'outputs': []}
    for r in res:
        results['succeeded'].append(r[0])
        results['timestamps'].append(r[1])
        results['imagenames'].append(r[2])
        results['outputs'].append(r[3])

    return results
