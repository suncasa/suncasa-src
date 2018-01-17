import numpy as np
import matplotlib.pyplot as plt
import os
from config import get_and_create_download_dir
import shutil
from astropy.io import fits
import urllib2
from suncasa.utils import helioimage2fits
import sunpy.map as smap
from sunpy.net import vso
from astropy import units as u
from astropy.time import Time
from taskinit import ms, tb, qa, iatool
import clean
from clean import clean
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

def svplot(specfile,timerange,spw,vis,pol='RRLL',dmin=None,dmax=None,goestime=None,reftime=None,fov=None,aiawave=171,imagefile=None,savefig=False,aiafits=None,changeheader=True,redoclean=False,fitsfile=None):
	'''
	Required input:
		specfile:name of the dynamic spectrum file
		vis:name of the visbility measurement file
		timerange:timerange for the clean
		spw:frequency range for the clean
	Optional input:
		pol:pol of the dynamic spectrum,can be 'RR','LL','I','V','IV','RRLL', default is 'RRLL'
		dmin,dmax:color bar parameter
		goestime:goes plot time,example ['2016/02/18 18:00:00','2016/02/18 23:00:00'](goes plot time, example ['2016-02-18 18:00:00','2016-02-18 23:00:00']-old)
		rhessisav:rhessi savefile
		reftime:reftime for the image
		fov:field of view in aia image, in unit of arcsec, example:[[-400,-200],[100,300]]
		aiawave:wave length of aia file in a
		imagefile: cleaned image file
		fitsfile: exist vla fitsfile
		savefig:weather to save the figure
	Example:
	
	'''
	#first part of dyn spec
	if pol !='RR' and pol !='LL' and pol !='I' and pol !='V' and pol !='RRLL' and pol !='IV':
		print 'wrong pol(only LL,RR,RRLL,I,V and IV)'
		return 0
	
	if not os.path.exists(specfile):
		print 'input specfile not exist'
		return -1

	if not os.path.exists(vis):
		print 'input measurement not exist'
		return -1

	
	#get prepared work
	tb.open(vis)
	starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
	tb.close()
	datstr = starttim.iso[:10]
	ms.open(vis, nomodify=True)
	vis_sp=timerange+'-'+vis
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
	specdata = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
	nfreq = specdata['amplitude'].shape[1]
	freq = specdata['axis_info']['freq_axis']['chan_freq'].reshape(nfreq)
	ms.close()
	ret1=Time(Time(starttim1),format='jd').plot_date
	ret2=Time(Time(endtim1),format='jd').plot_date
	req1=freq[0]/1e9
	req2=freq[-1]/1e9


	specdata=np.load(specfile)
	spec = specdata['spec']
	npol = specdata['npol']
	ntim = specdata['ntim']
	tidx=range(ntim)
	nfreq = specdata['nfreq']
	fidx=range(nfreq)
	tim = specdata['tim']
	freq = specdata['freq']
	bl = specdata['bl'].item()
	
	spec_med = np.median(np.absolute(spec))
	

	if not dmin:
		dmin = spec_med / 20.
	if not dmax:
		dmax = spec_med * 5.
	
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
		f=plt.figure(figsize=(24,16),dpi=120)
		f1=f.add_subplot(221)
		freqg= freq /1e9
		timstrr=range(tidx[-1]+1)
		for i in tidx:
			timstr = qa.time(qa.quantity(tim[i], 's'), form='clean', prec=9)
			timstrr[i]=Time(Time(datstr+' '+timstr[0]),format='jd').plot_date
		timstrr=np.array(timstrr)
		f1.pcolormesh(timstrr, freqg, spec_plt, cmap='jet', vmin=dmin, vmax=dmax)
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
		f=plt.figure(figsize=(24,16),dpi=120)
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
		f1.xaxis.set_major_formatter(DateFormatter("%H:%M"))
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
		f2.xaxis.set_major_formatter(DateFormatter("%H:%M"))
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
	
	#third part
	#download aia
	if aiafits:
		aiafits=aiafits
	else:
		newlist=[]
		items=os.listdir('.')
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
	
	if fitsfile:
		fitsfile=fitsfile
	else:
		if imagefile:
			imagefile=imagefile
		else:
			eph=helioimage2fits.read_horizons(vis_sp)
			phasecenter='J2000 '+str(eph['ra'][30])[:5]+'rad '+str(eph['dec'][30])[:5]+'rad'
			imagename=vis+timerange+spw
			os.system('rm -rf '+imagename+'*')
			print 'do the clean at '+timerange+' '+spw
			clean(vis=vis,imagename=imagename,selectdata=True,spw=spw,timerange=timerange,interactive=True,npercycle=50,imsize=[512,512],cell=['4.0arcsec'],phasecenter=phasecenter)
			os.system('rm -rf '+imagename+'.psf')
			os.system('rm -rf '+imagename+'.flux')
			os.system('rm -rf '+imagename+'.model')
			os.system('rm -rf '+imagename+'.mask')
			os.system('rm -rf '+imagename+'.residual')
			imagefile=imagename+'.image'
		fitsfile=imagefile+'.fits'
		helioimage2fits.imreg(vis=vis,imagefile=imagefile,timerange=timerange,reftime=reftime,fitsfile=fitsfile)
		if changeheader:
			data, header = fits.getdata(fitsfile, header=True)
			header['date-obs']=starttim1.iso[:10]+'T'+starttim1.iso[12:]
			outfits='n'+fitsfile
			fits.writeto(outfits, data, header, clobber=True)
			os.system('rm -rf '+fitsfile)
			fitsfile=outfits
	print 'vlafits file '+fitsfile+' selected'
	
	if aiafits:
		aiafits=aiafits
	else:
		newlist=[]
		items=os.listdir('.')
		for names in items:
			str1=starttim1.iso[:4]+'_'+starttim1.iso[5:7]+'_'+starttim1.iso[8:10]+'t'+starttim1.iso[11:13]+'_'+starttim1.iso[14:16]
			str2=str(aiawave)
			if names.endswith(".fits"):
				if names.find(str1) !=-1 and names.find(str2) !=-1:
					newlist.append(names)
		aiafits=newlist[0]
	print 'aiafits '+aiafits+' selected'

	f4=f.add_subplot(222)
	aiamap=smap.Map(aiafits)
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
	os.system('rm -rf '+vis_sp)
	os.system('rm -rf goes.py')
	os.system('rm -rf goes.dat')