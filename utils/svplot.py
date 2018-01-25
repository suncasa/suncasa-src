import numpy as np
import matplotlib.pyplot as plt
import os
#from config import get_and_create_download_dir
import shutil
from astropy.io import fits
import urllib2
from suncasa.utils import helioimage2fits as hf
from suncasa.utils import dspec2 as ds
import sunpy.map as smap
from sunpy.net import vso
from astropy import units as u
from astropy.time import Time
from taskinit import ms, tb, qa, iatool
from clean_cli import clean_cli as clean
from sunpy import lightcurve
from sunpy.time import TimeRange
from matplotlib.dates import DateFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord
from sunpy import lightcurve as lc
from sunpy.time import TimeRange,parse_time
import pickle
import datetime
import matplotlib
import matplotlib.cm as cm
import matplotlib.patches as patches
import glob

def svplot(vis, specfile=None, timerange=None, spw=None, pol='RRLL', dmin=None, dmax=None,
           goestime=None, reftime=None, fov=None, aiawave=171,
           imagefile=None, savefig=False, aiafits=None,
           changeheader=True, redoclean=False, fitsfile=None):
	'''
	Required inputs:
            vis: calibrated CASA measurement set 
            timerange: timerange for clean. Standard CASA time selection format
            spw: frequency selection for clean. Standard CASA spectral window selection format
	Optional inputs:
            specfile: supply dynamic spectrum save file (from suncasa.utils.dspec2.get_dspec()). Otherwise 
                      generate on the fly
            pol: pol of the dynamic spectrum,can be 'RR','LL','I','V','IV','RRLL', default is 'RRLL'
            dmin,dmax: color bar parameter
            goestime: goes plot time, example ['2016/02/18 18:00:00','2016/02/18 23:00:00']
            rhessisav: rhessi savefile
            reftime: reftime for the image
            fov: field of view in aia image, in unit of arcsec, example:[[-400,-200],[100,300]]
            aiawave: wave length of aia file in a
            imagefile: cleaned image file
            fitsfile: exist vla fitsfile
            savefig: whether to save the figure
	Example:
	
	'''
	#first part of dyn spec
	if pol !='RR' and pol !='LL' and pol !='I' and pol !='V' and pol !='RRLL' and pol !='IV':
		print 'wrong pol(only LL,RR,RRLL,I,V and IV)'
		return 0
	
	if not os.path.exists(vis):
		print 'input measurement not exist'
		return -1

	
	#split the data
        #generating dynamic spectrum
	if specfile:
            try:
                specdata=np.load(specfile)
            except:
                print('Provided dynamic spectrum file not numpy npz. Generating one from the visibility data')
                specdata=ds.get_dspec(vis,domedian=True,verbose=True)
        else:
            print('Dynamic spectrum file not provided; Generating one from the visibility data')
            specdata=ds.get_dspec(vis,domedian=True,verbose=True)

	tb.open(vis)
	starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
	tb.close()
	datstr = starttim.iso[:10]

        # I do not really like the following section. All we need here is a timerange and frequency range, which does
        # not require splitting and gridding the data at all.  --Bin
	ms.open(vis, nomodify=True)
	vis_sp=vis+'_tmp'
        if os.path.exists(vis_sp):
            os.system('rm -rf '+vis_sp)
	ms.split(outputms=vis_sp,time=timerange,spw=spw,whichcol='DATA')
	ms.close()
	tb.open(vis_sp)
	starttim1 = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
	endtim1 = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
	tb.close()
	ms.open(vis_sp, nomodify=False)
	ms.cvel(outframe='LSRK', mode='frequency', interp='nearest')
	ms.selectinit(datadescid=0, reset=True)
	data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
	nfreq = data['amplitude'].shape[1]
	freq = data['axis_info']['freq_axis']['chan_freq'].reshape(nfreq)
	ms.close()
	ret1=Time(Time(starttim1),format='jd').plot_date
	ret2=Time(Time(endtim1),format='jd').plot_date
	req1=freq[0]/1e9
	req2=freq[-1]/1e9
        #########

	spec = specdata['spec']
        (npol,nbl,nfreq,ntim)=spec.shape
	tidx=range(ntim)
	fidx=range(nfreq)
	tim = specdata['tim']
	freq = specdata['freq']
        try:
            bl = specdata['bl'].item()
        except:
            bl = specdata['bl']
	
	spec_med = np.median(np.absolute(spec))

        midtime_mjd = (starttim1.mjd + endtim1.mjd) / 2.
	
        f=plt.figure(figsize=(9,6),dpi=100)
	if pol !='RRLL' and pol !='IV':
		if pol =='RR':
			spec_plt=spec[0, 0, :, :]
		elif pol =='LL':
			spec_plt=spec[1, 0, :, :]
		elif pol =='I':
			spec_plt= (spec[0, 0, :, :] + spec[1, 0, :, :]) / 2.
		elif pol == 'V':
			spec_plt = (spec[0, 0, :, :] - spec[1, 0, :, :]) / 2.
	
		print 'plot the dynamic spectrum in pol '+pol
		f1=f.add_subplot(221)
		freqg= freq /1e9
		timstrr=range(tidx[-1]+1)
		for i in tidx:
			timstr = qa.time(qa.quantity(tim[i], 's'), form='clean', prec=9)
			timstrr[i]=Time(Time(datstr+' '+timstr[0]),format='jd').plot_date
		timstrr=np.array(timstrr)
		#f1.pcolormesh(timstrr, freqg, spec_plt, cmap='jet', vmin=dmin, vmax=dmax)
		f1.pcolormesh(timstrr, freqg, spec_plt, cmap='jet')
		#f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
		f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,alpha=0.4))
		f1.xaxis_date()
		f1.xaxis.set_major_formatter(DateFormatter("%H:%M"))
		f1.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
		f1.set_ylim(freqg[fidx[0]], freqg[fidx[-1]])
		
		f1.set_ylabel('Frequency (GHz)',fontsize=10)
		f1.set_title('Dynamic spectrum @ pol '+pol,fontsize=12)
		for tick in f1.get_xticklabels():
			#tick.set_fontname('Comic Sans MS')
			tick.set_fontsize(8)
		for tick in f1.get_yticklabels():
			tick.set_fontsize(8)
		f1.set_autoscale_on(False)
	else:
		R_plot=np.absolute(spec[0,0,:,:])
		L_plot=np.absolute(spec[1,0,:,:])
		I_plot=(R_plot+L_plot)/2.
		V_plot=(R_plot-L_plot)/2.
		if pol == 'RRLL':
			spec_plt_1=R_plot
			spec_plt_2=L_plot
			polstr=['RR','LL']
		if pol == 'IV':
			spec_plt_1=I_plot
			spec_plt_2=V_plot
			polstr=['I','V']
		
		print 'plot the dynamic spectrum in pol '+pol
		f1=f.add_subplot(321)
		freqg = freq/1e9
		timstrr=range(tidx[-1]+1)
		for i in tidx:
			timstr = qa.time(qa.quantity(tim[i], 's'), form='clean', prec=9)
			timstrr[i]=Time(Time(datstr+' '+timstr[0]),format='jd').plot_date
		timstrr=np.array(timstrr)
		f1.pcolormesh(timstrr,freqg,spec_plt_1,cmap='jet',vmin=dmin,vmax=dmax)
		f1.set_xlim(timstrr[tidx[0]],timstrr[tidx[-1]])
		f1.xaxis_date()
                f1.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
		f1.set_ylim(freqg[fidx[0]],freqg[fidx[-1]])
		f1.set_ylabel('Frequency (GHz)',fontsize=10)
		f1.set_title('Dynamic spectrum @pol '+polstr[0]+' in upper and '+polstr[1]+' in bottom',fontsize=12)
		f1.set_autoscale_on(False)
		#f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
		f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,alpha=0.4))
		for tick in f1.get_xticklabels():
			tick.set_fontsize(8)
		for tick in f1.get_yticklabels():
			tick.set_fontsize(8)
		f2=f.add_subplot(323)
		f2.pcolormesh(timstrr,freqg,spec_plt_2,cmap='jet',vmin=dmin,vmax=dmax)
		f2.set_xlim(timstrr[tidx[0]],timstrr[tidx[-1]])
		f2.xaxis_date()
                f2.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
		f2.set_ylim(freqg[fidx[0]],freqg[fidx[-1]])
		f2.set_ylabel('Frequency (GHz)',fontsize=10)
		for tick in f2.get_xticklabels():
			tick.set_fontsize(8)
		for tick in f2.get_yticklabels():
			tick.set_fontsize(8)
		f2.set_autoscale_on(False)
		#f2.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
		f2.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,alpha=0.4))
	
	#second part of Goes plot
	print 'plot the Goes soft X-ray curve and derivate'
	if goestime:
		goesplottim=TimeRange(goestime[0],goestime[1])
	else:
		datstrg=datstr.replace('-','/')
		goestime0=datstrg+' '+qa.time(qa.quantity(tim[0]-1800, 's'), form='clean', prec=9)[0]
		goestime1=datstrg+' '+qa.time(qa.quantity(tim[tidx[-1]-1]+1800, 's'), form='clean', prec=9)[0]
		goesplottim=TimeRange(goestime0,goestime1)
		
	if pol == 'RRLL' or pol == 'IV' :
		f3=f.add_subplot(427)
	else:
		f3=f.add_subplot(223)
	
	os.system('rm -rf goes.py')
	fi=open('goes.py','wb')
	fi.write('import os \n')
	fi.write('from sunpy.time import TimeRange \n')
	fi.write('from sunpy import lightcurve as lc \n')
	fi.write('import pickle \n')
	fi.write('goes = lc.GOESLightCurve.create(goesplottim) \n')
	fi.write('fi2 = open("goes.dat", "wb") \n')
	fi.write('pickle.dump(goes, fi2) \n')
	fi.write('fi2.close()')
	fi.close()
	
        try:
            execfile('goes.py')
	
            fi1 = file('goes.dat', 'rb')
            goest = pickle.load(fi1)
            fi1.close()
            
            dates = matplotlib.dates.date2num(parse_time(goest.data.index))
            goesdif=np.diff(goest.data['xrsb'])
            gmax=np.nanmax(goesdif)
            gmin=np.nanmin(goesdif)
            ra=gmax-gmin
            db=3e-4/ra
            goesdifp=np.diff(goest.data['xrsb'])*db+gmin+1e-4
            f3.plot_date(dates[0:-1], goesdifp, '-',label='derivate', color='blue', lw=0.4)
            f3.plot_date(dates, goest.data['xrsb'], '-',label='1.0--8.0 $\AA$', color='red', lw=2)
            
            f3.set_yscale("log")
            f3.set_ylim(1e-7, 1e-3)
            f3.set_title('Goes Soft X-ray and derivate')
            f3.set_ylabel('Watts m$^{-2}$')
            f3.set_xlabel(datetime.datetime.isoformat(goest.data.index[0])[0:10])
            f3.axvspan(dates[899],dates[dates.size-899],alpha=0.2)
            
            ax2 = f3.twinx()
            ax2.set_yscale("log")
            ax2.set_ylim(1e-7, 1e-3)
            ax2.set_yticks(( 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))
            ax2.set_yticklabels(( 'B', 'C', 'M', 'X'))

            f3.yaxis.grid(True, 'major')
            f3.xaxis.grid(False, 'major')
            f3.legend(prop={'size': 6})
            
            formatter = matplotlib.dates.DateFormatter('%H:%M')
            f3.xaxis.set_major_formatter(formatter)
            
            f3.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        except:
            print 'error in plotting GOES. Continue without GOES'
	
	#third part
        # start to download the fits files
	if not aiafits:
            newlist=[]
            items=glob.glob('*.fits')
            for names in items:
                str1=starttim1.iso[:4]+'_'+starttim1.iso[5:7]+'_'+starttim1.iso[8:10]+'t'+starttim1.iso[11:13]+'_'+starttim1.iso[14:16]
                str2=str(aiawave)
                if names.endswith(".fits"):
                    if names.find(str1) !=-1 and names.find(str2) !=-1:
                        newlist.append(names)
		newlist.append('0')
		if os.path.exists(newlist[0]):
                    aiafits=newlist[0]
		else:
                    print 'downloading the aiafits file'
                    client=vso.VSOClient()
                    wave1=aiawave-3
                    wave2=aiawave+3
                    t1 = Time(starttim1.mjd - 0.02/24.,format='mjd')
                    t2 = Time(endtim1.mjd + 0.02/24.,format='mjd')
                    qr=client.query(vso.attrs.Time(t1.iso,t2.iso),vso.attrs.Instrument('aia'),vso.attrs.Wave(wave1*u.AA, wave2*u.AA))
                    res=client.get(qr,path='{file}')
        # Here something is needed to check whether it has finished downloading the fits files or not
	
	if fitsfile:
		fitsfile=fitsfile
	else:
		if not imagefile:
                    eph=hf.read_horizons(t0=Time(midtime_mjd,format='mjd'))
                    phasecenter='J2000 '+str(eph['ra'][0])[:15]+'rad '+str(eph['dec'][0])[:15]+'rad'
                    imagename=vis+'__'
                    if os.path.exists(imagename+'.image') or os.path.exists(imagename+'.flux'):
                        os.system('rm -rf '+imagename+'*')
                    print 'do clean for '+timerange+' in spw '+spw
                    print 'use phasecenter: '+phasecenter
                    clean(vis=vis,imagename=imagename,selectdata=True,spw=spw,timerange=timerange,niter=500,
                          interactive=False,npercycle=50,imsize=[512,512],cell=['5.0arcsec'],phasecenter=phasecenter)
                    os.system('rm -rf '+imagename+'.psf')
                    os.system('rm -rf '+imagename+'.flux')
                    os.system('rm -rf '+imagename+'.model')
                    os.system('rm -rf '+imagename+'.mask')
                    os.system('rm -rf '+imagename+'.residual')
                    imagefile=imagename+'.image'
		fitsfile=imagefile+'.fits'
		hf.imreg(vis=vis,ephem=eph,imagefile=imagefile,timerange=timerange,reftime=reftime,fitsfile=fitsfile,verbose=True)
		if changeheader:
			data, header = fits.getdata(fitsfile, header=True)
			header['date-obs']=starttim1.iso[:10]+'T'+starttim1.iso[12:]
			outfits='n'+fitsfile
			fits.writeto(outfits, data, header, clobber=True)
			os.system('rm -rf '+fitsfile)
			fitsfile=outfits
	print 'fits file '+fitsfile+' selected'
	
	if not aiafits:
            newlist=[]
            items=glob.glob('*.fits')
            for nm in items:
                str1=starttim1.iso[:4]+'_'+starttim1.iso[5:7]+'_'+starttim1.iso[8:10]+'t'+starttim1.iso[11:13]+'_'+starttim1.iso[14:16]
                str2=str(aiawave)
                if nm.find(str1) !=-1 and nm.find(str2) !=-1:
                    newlist.append(nm)
            if newlist:
                aiafits=newlist[0]
                print 'AIA fits '+aiafits+' selected'
            else:
                print 'no AIA fits files found. Proceed without AIA'

        try:
            aiamap=smap.Map(aiafits)
        except:
            print 'error in reading aiafits. Proceed without AIA'

	f4=f.add_subplot(222)
	vlafits=fitsfile
	vlamap=smap.Map(vlafits)
	vlamap.data = vlamap.data.reshape(vlamap.meta['naxis1'],vlamap.meta['naxis2'])
	vlamap=vlamap.submap([-1200,1200]*u.arcsec,[-1200,1200]*u.arcsec)
	cmap=smap.CompositeMap(aiamap)
	cmap.add_map(vlamap,levels=np.array([0.5,0.7,0.9])*np.nanmax(vlamap.data))
	cmap.set_colors(1,cm=cm.rainbow)
	f4=cmap.plot(title='overview image')
	#f4[2].set_autoscale_on(False)  gca
	
	if fov:
		fov=fov
	else:
		raw, col = vlamap.data.shape
		positon = np.nanargmax(vlamap.data)
		m, n = divmod(positon, col)
		length = 200 * u.arcsec
		x0=vlamap.xrange[0]+vlamap.scale[1]*(n+0.5)*u.pix
		y0=vlamap.yrange[0]+vlamap.scale[0]*(m+0.5)*u.pix
		x1=x0-length
		x2=x0+length
		y1=y0-length
		y2=y0+length
		fov=[[x1.value,x2.value],[y1.value,y2.value]]
	
	aiamap.draw_rectangle((fov[0][0],fov[1][0])*u.arcsec, 400*u.arcsec, 400*u.arcsec)
	f5=f.add_subplot(224)
	subaiamap=aiamap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec)
	subvlamap=vlamap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec)
	scmap=smap.CompositeMap(subaiamap)
	scmap.add_map(subvlamap,levels=np.array([0.5,0.7,0.9])*np.nanmax(subvlamap.data))
	scmap.set_colors(1,cm=cm.rainbow)
	f5=scmap.plot(title='zoomed view')
	
		
	f.show()
	#os.system('rm -rf '+vis_sp)
	os.system('rm -rf goes.py')
	os.system('rm -rf goes.dat')
