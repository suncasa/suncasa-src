#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import copy, struct
import shutil
import os
from math import *
import dspec 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from fnmatch import fnmatch
import pdb
from taskinit import *
import optparse
import pyfits
from glob import glob

stokesdict={'I':1,'Q':2,'U':3,'V':4,'RR':-1,'LL':-2,'RL':-3,'LR':-4,'XX':-5,'YY':-6,'XY':-7,'YX':-8}

def get_vspec(fitspath=None,vspecfile=None,tinc=1,rgns=None):
    #setup initial parametesr
    #fitspath = '/home/bchen/work/EVLA/20120303/S16-20/fits/rgn1/'
    #rgns=[[[95,112],[120,140]],
    #      [[129,109],[163,141]],
    #      [[111,66],[135,98]]]
    #      [[182,135],[213,168]]]
    # check if the path exist
    if not os.path.exists(fitspath):
        print 'specified path '+fitspath+' does not exist. aborting...'
        return
    fitsfiles_= glob(fitspath+'/*.fits')
    if len(fitsfiles_) < 1:
        print 'no fits files found. aborting...'
        return
    fitsfiles_.sort()
    fitsfiles=fitsfiles_[::tinc]
    ntim=len(fitsfiles)
    tims=np.zeros(ntim)
    timstrs=[]

    #determine the frequency axis
    fitsfile0=fitsfiles[0]
    ia.open(fitsfile0)
    info0=ia.summary()
    dim=info0['shape']
    nx0=dim[0]
    ny0=dim[1]
    nf0=dim[2]
    npol0=dim[3]
    stokes=[]
    svals=(np.arange(npol0)*info0['incr'][3]+info0['refval'][3])
    for sval in svals:
        stokes.append(stokesdict.keys()[stokesdict.values().index(sval)])
    ia.close()
          
    nrgn=len(rgns)
    fmaxs=np.zeros((nf0,ntim,npol0,nrgn))
    fsums=np.zeros((nf0,ntim,npol0,nrgn))
    fmeans=np.zeros((nf0,ntim,npol0,nrgn))
    print str(ntim)+' times to process...'
    timstrs=[]
    tims = np.zeros(ntim)
    for t in range(ntim):
        fitsfile=fitsfiles[t]
        #get fits information 
        hdu=pyfits.open(fitsfile)
        hdr = hdu[0].header
        timstr = hdr['date-obs']
        timstrs.append(timstr)
        tims[t] = qa.convert(timstr,'s')['value']
        exptime=hdr['exptime']
        ndim = hdr['naxis']
        if ndim != 4:
            print 'input dimension is not 4. skipping this fits file'
            break
        nx = hdr['naxis1']
        ny = hdr['naxis2']
        nf = hdr['naxis3']
        npol = hdr['naxis4']
        crpixs = [hdr['crpix1'],hdr['crpix2'],hdr['crpix3'],hdr['crpix4']]
        crvals = [hdr['crval1'],hdr['crval2'],hdr['crval3'],hdr['crval4']]
        cdelts = [hdr['cdelt1'],hdr['cdelt2'],hdr['cdelt3'],hdr['cdelt4']]
        freqs=(np.arange(nf)*cdelts[2]+crvals[2])/1e6
        hdu.close()

        if t % 50 == 0 and t > 0:
            print 'Processed '+str(t)+' times'
            print 'Current time: ', timstr
        ia.open(fitsfile)
        for p in range(npol):
            for k in range(nrgn):
                rgn=rgns[k]
                r1=rg.box(blc=rgn[0],trc=rgn[1])
                try:
                    mystat=ia.statistics(axes=[0,1],region=r1)
                    fmax=mystat['max']
                    fsum=mystat['sum']
                    fmean=mystat['mean']
                    npt=mystat['npts']
                    fmaxs[:,t,p,k]=fmax
                    fsums[:,t,p,k]=fsum
                    fmeans[:,t,p,k]=fmean
                except:
                    print 'failure for this time: ',timstrs[t]
        ia.close()
    if not vspecfile:
        vspecfile='./vspec.t'+timstrs[0].replace(':','')+'-'+\
                  timstrs[-1].replace(':','')+'.npz'
    np.savez(vspecfile,rgns=rgns,freqs=freqs,tims=tims,timstrs=timstrs,
             stokes=stokes,fmaxs=fmaxs,fsums=fsums,fmeans=fmeans)

def wrt_vspec(vspecfile=None,vspecdat=None):
    try:
        vspecfile
    except NameError:
        print 'No input centfile specified for reading. Aborting...'
    if not vspecdat:
        print 'Output file name is not specified, use the default convention'
        vspecdat=vspecfile.replace('npz','dat')
    vspecdata=np.load(vspecfile)
    rgns=vspecdata['rgns']
    nrgn=rgns.shape[0]
    stokes=vspecdata['stokes']
    nstokes=len(stokes)
    fmx=copy.deepcopy(vspecdata['fmaxs'][:,:,:,:])
    fmean=copy.deepcopy(vspecdata['fmeans'][:,:,:,:])
    nf,nt,npl,nr=fmx.shape
    print 'Dimension of the data cube -- # of freq, # of time, # of pol, # of regions:'
    print nf, nt, npl, nr
    # Turn fmx into a 1d array for writing
    fmx.shape = (nf*nt*npl*nr)
    fmean.shape = (nf*nt*npl*nr)
    # Pack the data into a byte buffer.  Note the asterisks in the next three lines!
    # If your data includes integers or other data types, use 'i' or other instead of 'f'
    buf = struct.pack(str(nf)+'f',*vspecdata['freqs'])
    buf += struct.pack(str(nt)+'f',*vspecdata['tims'])
    buf += struct.pack(str(nf*nt*npl*nr)+'f',*fmx)
    buf += struct.pack(str(nf*nt*npl*nr)+'f',*fmean)
    with open(vspecdat,'wb') as f:
        f.write(buf)
    f.close()

def plt_vspec(vspecfile=None,imgfile=None,moment=None, pltmax=True, pltmean=False, bkgsub=False, \
              dmin=None,dmax=None,savepng=True,savepdf=False,\
              timerange=None, freqrange=None,timestr=False):
    # vspecfile: saved npz file containing the vector dynamic spectrum data (from get_max)
    # imagefile: image files to overplot selected regions, the number of images 
    #               must equal to len(stokes)
    # pltmax: plot the maximum value of each selected region
    # pltmin: plot the mean value of each selected region
    # bkgsub: if True, subtract max/mean of the first region from all the other regions
    # dmin: minimum flux to show, one- or len(stokes) element list 
    # dmax: maximum flux to show, one- or len(stokes) element list
    # savepng: if True, save the resulted figure as a png file
    # savepdf: if True, save the resulted figure as a pdf file
    # timerange: format: ['18:00:00','19:00:00']
    # freqrange: format: [1000.,1500.]
    try:
        vspecfile
    except NameError:
        print 'No input centfile specified for reading. Aborting...'
    vspecdata=np.load(vspecfile)
    rgns=vspecdata['rgns']
    nrgn=rgns.shape[0]
    fmaxs=vspecdata['fmaxs']
    fmeans=vspecdata['fmeans']
    tims=vspecdata['tims']
    ntim=tims.shape
    freqs=vspecdata['freqs']
    nfreq=freqs.shape
    stokes=vspecdata['stokes']
    nstokes=len(stokes)
    if pltmax:
        fluxs=fmaxs
    if pltmean:
        fluxs=fmeans
    if not moment:
        moment = 0
    if bkgsub:
        if nrgn < 2:
            print 'Background subtraction requires more than one regions'
            return
        for r in range(nrgn-1):
            fluxs[:,:,:,r+1]-=fluxs[:,:,:,0] #the first region is the background
    if (not dmin) or (not dmax):
        print 'will use automatic brightness scaling'
    elif len(dmin) == 1 and len(dmax) == 1:
        dmin=nstokes*dmin
        dmax=nstokes*dmax
    elif len(dmin) != nstokes or len(dmax) != nstokes:
        raise ValueError('Number of dmin or dmax does not agree with number of stokes')
    if timerange:
        if type(timerange[0]) is str:
            timerange=[qa.convert(qa.quantity(t),'s')['value'] for t in timerange]
        tidx=np.where((tims >= timerange[0]) & (tims <= timerange[1]))[0]
    else:
        tidx=range(ntim[0])
    if freqrange:
        fidx=np.where((freqs >= freqrange[0]) & (freqs <= freqrange[1]))[0]
    else:
        fidx=range(nfreq[0])
    
    #set up time axis
    #1 sec per tick
    tstart=np.fix(tims[0])
    xticks=np.arange(ntim[0]/20)+tstart
    nticks=xticks.shape
    xticknames=[]
    for i in range(nticks[0]):
        tim=xticks[i]
        #timstr=qa.time(qa.quantity(tim,'s'),prec=6)[0]
        timstr="%d" % (tim-tstart)
        xticknames.append(timstr)
    xticks=list(xticks)

    #do the plot
    #plot moment image
    if imgfile:
        try:
            ia.open(imgfile)
        except NameError:
            print 'Specified image file can not be opened'
        f,ax = plt.subplots(figsize=(5,5))
        im2=ia.moments(moments=[moment],excludepix=[-1000,0])
        pix=im2.getchunk()
        im2.done()
        ia.close()
        imdim=pix.shape
        immom=np.transpose(pix[:,:,0])
        ax.imshow(immom,aspect='auto',origin='lower',extent=[0,imdim[0]-1,0,imdim[1]-1])
        ax.set_autoscale_on(False)
        ax.text(0.03,0.93,'Momment '+str(moment)+' image',color='w')
        for rgn in rgns:
            xbg=rgn[0][0]
            xend=rgn[1][0]
            ybg=rgn[0][1]
            yend=rgn[1][1]
            ax.axvspan(xbg,xend,ymin=float(ybg)/float(imdim[1]),ymax=float(yend)/float(imdim[1]),
                        alpha=1.0,fill=False,lw=1,color='w')
    figsize=(8,10)
    nrgns=rgns.shape[0]
    if bkgsub:
        nplot=nrgns-1
    else:
        nplot=nrgns
    f2,axarr = plt.subplots(nplot,nstokes,figsize=figsize)
    for p in range(nstokes):
        for n in range(nplot):
            if nplot == 1 and nstokes == 1:
                ax=axarr
            if nplot > 1 and nstokes == 1:
                ax=axarr[n]
            if nplot > 1 and nstokes > 1:
                ax=axarr[n,p]
            if bkgsub:
                flux=fluxs[:,:,p,n+1]
            else:
                flux=fluxs[:,:,p,n]
            mflux=ma.masked_outside(flux,10,1e6)
            flux_med=ma.median(mflux)
            if not dmin or not dmax:
                vmin=flux_med/3.
                vmax=flux_med*8.
            else:
                vmin=dmin[p]
                vmax=dmax[p]
            ax.pcolormesh(tims,freqs,flux,vmin=vmin,vmax=vmax)
            ax.set_xlim(tims[tidx[0]],tims[tidx[-1]])
            ax.set_ylim(freqs[fidx[0]],freqs[fidx[-1]])
            ax.set_ylabel('Frequency (MHz)')
            if timestr:
                labels=ax.get_xticks().tolist()
                newlabels=[qa.time(qa.quantity(lb,'s'))[0] for lb in labels] 
                ax.set_xticklabels(newlabels)
            #ax.text(0.05,0.04,'1',color='k',transform=ax.transAxes)
    if savepng:
        figfile=vspecfile[:(vspecfile.find('npz'))]+'png'
        f2.savefig(figfile)
    if savepdf:
        figfile=vspecfile[:(vspecfile.find('npz'))]+'pdf'
        f2.savefig(figfile)

def plt_max(vspecfile='',fitsfile='',pltmax=True, pltmean=False, bkgsub=True,rgns='',dmin='',dmax='',figfile=''):
    #vspecfile = vspecpath+'maxfit.t183100.025-183109.974.I.npz'
    #fitspath = '/lustre/bchen/CASA/MAR3/S16-20/fits/rgn1/'
    #fitsfile = 'SUN_66660.025.fits'
    vspecdata=np.load(vspecfile)
    rgns=vspecdata['rgns']
    nrgn=rgns.shape[0]
    fmaxs=vspecdata['fmaxs']
    fmeans=vspecdata['fmeans']
    tims=vspecdata['tims']
    freqs=vspecdata['freqs']
    if pltmax:
        fluxs=fmaxs
    if pltmean:
        fluxs=fmeans
    if bkgsub:
        if nrgn < 2:
            print 'Background subtraction requires more than one regions'
            return
        for r in range(nrgn-1):
            fluxs[:,:,:,r+1]-=fluxs[:,:,:,0] #the first region is the background
    
    #setup plot parameters
    fluxs_med=np.median(fluxs)
    if not dmin:
        dmin=fluxs_med/20.
    if not dmax:
        dmax=fluxs_med*2.
    #fill missing times with zeros
    #find out the time range
    timeran=tims[-1]-tims[0]
    #find out the integration interval
    ntim=len(tims)
    t_ints=[]
    iter=0
    for i in range(ntim-1):
        t_ints.append(tims[i+1]-tims[i])
    t_int=np.median(t_ints)
    ntim_plt=round(timeran/t_int)

    xi=np.linspace(tims[0],tims[-1],900)
    yi=np.linspace(freqs[0],freqs[-1],128)
    fluxs_plt=griddata((tims,freqs),fluxs,(xi,yi),method='cubic')
    #1 sec per tick
    tstart=np.fix(tims[0])
    xticks=np.arange(ntim[0]/20)+tstart
    nticks=xticks.shape
    xticknames=[]
    for i in range(nticks[0]):
        tim=xticks[i]
        #timstr=qa.time(qa.quantity(tim,'s'),prec=6)[0]
        timstr="%d" % (tim-tstart)
        xticknames.append(timstr)
    xticks=list(xticks)

    #do the plot
    figsize=(7.36,8)
    f=plt.figure(figsize=figsize,dpi=100)
    #rect1=[0.12,0.52,0.27,0.46]
    #rect2=[0.45,0.52,0.50,0.46]
    #rect3=[0.12,0.06,0.27,0.42]
    #rect4=[0.40,0.06,0.27,0.42]
    #rect5=[0.68,0.06,0.27,0.42]
    rect1=[0.12,0.52,0.27,0.46]
    rect2=[0.45,0.52,0.50,0.46]
    rect3=[0.12,0.06,0.20,0.42]
    rect4=[0.33,0.06,0.20,0.42]
    rect5=[0.54,0.06,0.20,0.42]
    rect6=[0.75,0.06,0.20,0.42]
    ax1=f.add_axes(rect1)
    ax2=f.add_axes(rect2)
    ax3=f.add_axes(rect3)
    ax4=f.add_axes(rect4)
    ax5=f.add_axes(rect5)
    ax6=f.add_axes(rect6)
    plt.figtext(0.3,0.01,'Time in seconds, starting from 18:30:02')
    #plot cross-power dspec
    specpath='/lustre/bchen/CASA/MAR3/S16-20/'
    specfile=specpath+'sun.mar3.S16-20.cal.grd.ms.bl19-22.spec.npz'
    specdata=np.load(specfile)
    spec=specdata['spec']
    xtim=specdata['tim']
    xtimsecs=np.zeros(xtim.shape)
    #convert to seconds from 2012-03-03 00:00
    for i in range(xtim.shape[0]):
        xtimsecs[i]=qa.sub(qa.quantity(xtim[i],'s'),qa.quantity('2012-03-03','s'))['value']
    xfreq=specdata['freq']
    tbgind=np.where(xtimsecs>tims[0])[0][0]
    tendind=np.where(xtimsecs<tims[-1])[0][-1]
    fbgind=np.where(xfreq/1e6>=freqs[0])[0][0]
    fendind=np.where(xfreq/1e6<=freqs[-1])[0][-1]
    specsel=np.absolute(spec[:,fbgind:fendind+1,tbgind:tendind+1])
    #pdb.set_trace()
    ax1.imshow((specsel[0,:,:]+specsel[1,:,:])/2,aspect='auto',origin='lower',extent=[tims[0],tims[-1],freqs[0],freqs[-1]])
    ax1.set_autoscale_on(False)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticknames)
    ax1.set_ylabel('Frequency (MHz)')
    ax1.text(0.05,0.93,'x-power spec',color='k',transform=ax1.transAxes)
    #plot moment image
    ia.open(fitspath_+fitsfile)
    im2=ia.moments(moments=[0],excludepix=[-1000,0])
    pix=im2.getchunk()
    im2.done()
    ia.close()
    imdim=pix.shape
    immom=np.transpose(pix[:,:,0])
    ax2.imshow(immom,aspect='auto',origin='lower',extent=[0,imdim[0]-1,0,imdim[1]-1])
    ax2.set_autoscale_on(False)
    ax2.text(100,150,'1',color='w')
    ax2.text(150,150,'2',color='w')
    ax2.text(120,50,'3',color='w')
    ax2.text(170,160,'4',color='w')
    ax2.text(0.03,0.93,'Image at 18:31:00 (moment '+str(moment)+')',color='w',transform=ax2.transAxes)
    #rgns=[[[95,112,0,0],[120,140,0,0]],
    #      [[129,109,0,0],[163,141,0,0]],
    #      [[111,66,0,0],[135,98,0,0]]]
    for rgn in rgns:
        xbg=rgn[0][0]
        xend=rgn[1][0]
        ybg=rgn[0][1]
        yend=rgn[1][1]
        ax2.axvspan(xbg,xend,ymin=float(ybg)/float(imdim[1]),ymax=float(yend)/float(imdim[1]),
                    alpha=1.0,fill=False,lw=1,color='w')
    ax3.imshow(fluxs[:,:,0],aspect='auto',origin='lower',extent=[tims[0],tims[-1],freqs[0],freqs[-1]])
    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticknames)
    #ax2.set_xlabel('Universal Time (50 ms/pixel)')
    ax3.set_ylabel('Frequency (MHz)')
    #ax2.set_title('region 3')
    ax3.text(0.05,0.04,'1',color='k',transform=ax3.transAxes)
    ax4.imshow(fluxs[:,:,1],aspect='auto',origin='lower',extent=[tims[0],tims[-1],freqs[0],freqs[-1]])
    ax4.set_xticks(xticks)
    ax4.set_xticklabels(xticknames)
    #ax3.set_xlabel('Universal Time (50 ms/pixel)')
    ax4.set_yticks([])
    ax4.set_ylabel('')
    #ax3.set_title('region 3')
    ax4.text(0.05,0.04,'2',color='w',transform=ax4.transAxes)
    ax5.imshow(fluxs[:,:,2],aspect='auto',origin='lower',extent=[tims[0],tims[-1],freqs[0],freqs[-1]])
    ax5.set_xticks(xticks)
    ax5.set_xticklabels(xticknames)
    #ax4.set_xlabel('Universal Time (50 ms/pixel)')
    ax5.set_yticks([])
    ax5.set_ylabel('')
    #ax4.set_title('region 3')
    ax5.text(0.05,0.04,'3',color='k',transform=ax5.transAxes)
    ax6.imshow(fluxs[:,:,3],aspect='auto',origin='lower',extent=[tims[0],tims[-1],freqs[0],freqs[-1]])
    ax6.set_xticks(xticks)
    ax6.set_xticklabels(xticknames)
    #ax4.set_xlabel('Universal Time (50 ms/pixel)')
    ax6.set_yticks([])
    ax6.set_ylabel('')
    #ax4.set_title('region 3')
    ax6.text(0.05,0.04,'4',color='k',transform=ax6.transAxes)
    if not figfile:
        figfile=centfile[:(centfile.find('npz'))]+'pdf'
    else:
        figfile=figfile
    f.savefig(figfile,dpi=300)

if __name__ == "__main__":
    import sys
    get_max(sys.argv[1],sys.argv[2])
