#
# This file was generated using xslt from its XML file
#
# Copyright 2009, Associated Universities Inc., Washington DC
#
import sys
import os
from  casac import *
import string
from taskinit import casalog
from taskinit import xmlpath
#from taskmanager import tm
import task_ptclean
def ptclean(vis='', imagedir='./', ncpu=8, twidth=1, doreg=False, ephemfile='', msinfofile='', outlierfile='', field='', spw='', selectdata=True, timerange='', uvrange='', antenna='', scan='', observation='', intent='', mode='mfs', resmooth=False, gridmode='', wprojplanes=-1, facets=1, cfcache='cfcache.dir', rotpainc=5.0, painc=360.0, aterm=True, psterm=False, mterm=True, wbawp=False, conjbeams=True, epjtable='', interpolation='linear', niter=500, gain=0.1, threshold='0.0mJy', psfmode='clark', imagermode='csclean', ftmachine='mosaic', mosweight=False, scaletype='SAULT', multiscale=[0], negcomponent=-1, smallscalebias=0.6, interactive=False, mask=[], nchan=-1, start=0, width=1, outframe='', veltype='radio', imsize=[256, 256], cell=['1.0arcsec'], phasecenter='', restfreq='', stokes='I', weighting='natural', robust=0.0, uvtaper=False, outertaper=[''], innertaper=['1.0'], modelimage='', restoringbeam=[''], pbcor=False, minpb=0.2, usescratch=False, noise='1.0Jy', npixels=0, npercycle=100, cyclefactor=1.5, cyclespeedup=-1, nterms=1, reffreq='', chaniter=False, flatnoise=True, allowchunk=False):

        """Parallelized clean in consecutive time steps
  
        """
        if type(multiscale)==int: multiscale=[multiscale]
        if type(imsize)==int: imsize=[imsize]
        if type(cell)==float: cell=[cell]
        if type(outertaper)==str: outertaper=[outertaper]
        if type(innertaper)==str: innertaper=[innertaper]
        if type(restoringbeam)==str: restoringbeam=[restoringbeam]

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['imagedir'] = imagedir
        mytmp['ncpu'] = ncpu
        mytmp['twidth'] = twidth
        mytmp['doreg'] = doreg
        mytmp['ephemfile'] = ephemfile
        mytmp['msinfofile'] = msinfofile
        mytmp['outlierfile'] = outlierfile
        mytmp['field'] = field
        mytmp['spw'] = spw
        mytmp['selectdata'] = selectdata
        mytmp['timerange'] = timerange
        mytmp['uvrange'] = uvrange
        mytmp['antenna'] = antenna
        mytmp['scan'] = scan
        mytmp['observation'] = observation
        mytmp['intent'] = intent
        mytmp['mode'] = mode
        mytmp['resmooth'] = resmooth
        mytmp['gridmode'] = gridmode
        mytmp['wprojplanes'] = wprojplanes
        mytmp['facets'] = facets
        mytmp['cfcache'] = cfcache
        mytmp['rotpainc'] = rotpainc
        mytmp['painc'] = painc
        mytmp['aterm'] = aterm
        mytmp['psterm'] = psterm
        mytmp['mterm'] = mterm
        mytmp['wbawp'] = wbawp
        mytmp['conjbeams'] = conjbeams
        mytmp['epjtable'] = epjtable
        mytmp['interpolation'] = interpolation
        mytmp['niter'] = niter
        mytmp['gain'] = gain
        if type(threshold) == str :
           mytmp['threshold'] = casac.quanta().quantity(threshold)
        else :
           mytmp['threshold'] = threshold
        mytmp['psfmode'] = psfmode
        mytmp['imagermode'] = imagermode
        mytmp['ftmachine'] = ftmachine
        mytmp['mosweight'] = mosweight
        mytmp['scaletype'] = scaletype
        mytmp['multiscale'] = multiscale
        mytmp['negcomponent'] = negcomponent
        mytmp['smallscalebias'] = smallscalebias
        mytmp['interactive'] = interactive
        mytmp['mask'] = mask
        mytmp['nchan'] = nchan
        mytmp['start'] = start
        mytmp['width'] = width
        mytmp['outframe'] = outframe
        mytmp['veltype'] = veltype
        mytmp['imsize'] = imsize
        if type(cell) == str :
           mytmp['cell'] = casac.quanta().quantity(cell)
        else :
           mytmp['cell'] = cell
        mytmp['phasecenter'] = phasecenter
        mytmp['restfreq'] = restfreq
        mytmp['stokes'] = stokes
        mytmp['weighting'] = weighting
        mytmp['robust'] = robust
        mytmp['uvtaper'] = uvtaper
        mytmp['outertaper'] = outertaper
        mytmp['innertaper'] = innertaper
        mytmp['modelimage'] = modelimage
        mytmp['restoringbeam'] = restoringbeam
        mytmp['pbcor'] = pbcor
        mytmp['minpb'] = minpb
        mytmp['usescratch'] = usescratch
        mytmp['noise'] = noise
        mytmp['npixels'] = npixels
        mytmp['npercycle'] = npercycle
        mytmp['cyclefactor'] = cyclefactor
        mytmp['cyclespeedup'] = cyclespeedup
        mytmp['nterms'] = nterms
        mytmp['reffreq'] = reffreq
        mytmp['chaniter'] = chaniter
        mytmp['flatnoise'] = flatnoise
        mytmp['allowchunk'] = allowchunk
	pathname="file:///Users/binchen/Dropbox/bc_python/casa_task/"
	trec = casac.utils().torecord(pathname+'ptclean.xml')

        casalog.origin('ptclean')
        if trec.has_key('ptclean') and casac.utils().verify(mytmp, trec['ptclean']) :
	    result = task_ptclean.ptclean(vis, imagedir, ncpu, twidth, doreg, ephemfile, msinfofile, outlierfile, field, spw, selectdata, timerange, uvrange, antenna, scan, observation, intent, mode, resmooth, gridmode, wprojplanes, facets, cfcache, rotpainc, painc, aterm, psterm, mterm, wbawp, conjbeams, epjtable, interpolation, niter, gain, threshold, psfmode, imagermode, ftmachine, mosweight, scaletype, multiscale, negcomponent, smallscalebias, interactive, mask, nchan, start, width, outframe, veltype, imsize, cell, phasecenter, restfreq, stokes, weighting, robust, uvtaper, outertaper, innertaper, modelimage, restoringbeam, pbcor, minpb, usescratch, noise, npixels, npercycle, cyclefactor, cyclespeedup, nterms, reffreq, chaniter, flatnoise, allowchunk)

	else :
	  result = False
        return result
