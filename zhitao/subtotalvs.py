import sys
import os
from pylab import *
from numpy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from  casac import *
import string
from taskinit import casalog
from taskmanager import tm
import task_subvs_lv1
import task_smooth
import numpy as np
from taskinit import ms, cb, tb, qa, casalog, write_history

def pltvs(vis='',outputvis='',plttimerange='',pltspw1='',pltspw2='',pltspw3='',timoffset=4,windowlen=10,windowtype='hamming',pol='LL',bl='19&22'):
	ms.open(vis,nomodify=True)
	""" This function can do two steps: (1) plot the uv amplitude vs. time 
		on three selected channel and compare the original and smoothed signals.
		You can select appropriate window length and type, and the offset time
		to match the peaks and valleys for the three channel curves. (2) Confirm
		the use of specified window length, type and offset time and smooth the
		data channel by channel. The background-subtracted measurement set is 
		then generated.
    """
	timfreq=ms.getdata(['time','axis_info'],ifraxis=True)
	tim = timfreq['time']
	timerange=str(qa.time(qa.quantity(tim[0],'s'),prec=8)[0])+'~'+str(qa.time(qa.quantity(tim[-1],'s'),prec=8)[0])
	# check plotting timerange
	if plttimerange and (type(plttimerange)== str):
		print 'plotting the specified timerange: ', plttimerange
	else:	
		plttimerange=timerange
		print 'plotting the entire timerange: ', plttimerange

	if pltspw1 and (type(pltspw1)== str):
		print 'Using the specified channel 1:', pltspw1
	else:
		pltspw1 = '0:400'

	if pltspw2 and (type(pltspw2)== str):
		print 'Using the specified channel 2:', pltspw2
	else:
		pltspw2 = '0:500'

	if pltspw3 and (type(pltspw3)== str):
		print 'Using the specified channel 3:', pltspw3
	else:
		pltspw3 = '0:600'

	[spwid1,chan1]=pltspw1.split(':')
	[spwid2,chan2]=pltspw2.split(':')
	[spwid3,chan3]=pltspw3.split(':')
	chanid=[chan1,chan2,chan3]

	if not (spwid1 == spwid2 and spwid1 == spwid3):
		print 'Please use the same spectral window'
		exit()

	if not (timoffset and (type(timoffset)==int)):
		timoffset=int(4)
	# initialize the loop-out status
	status='n'
	while True:
		timoffset = int(raw_input("Please specify the offset after smoothing:" ))
		windowlen = int(raw_input("Please specify window width for smoothing:" ))
		windowtype = str(raw_input("Please specify window type for smoothing: (e.g. 'flat', 'hanning', 'hamming', 'bartlett', 'blackman')" ))
		pol = str(raw_input("Please specify polarization: (e.g. RR/LL)" ))
		bl = str(raw_input("Please specify baseline: (e.g. '19&22')" ))
		j=331
		for i in range(len(chanid)):
			ms.selectinit(datadescid=int(spwid1))
			if bl and (type(bl)==str):
				ms.msselect({'baseline':bl})
			if timerange and (type(timerange)==str):
				ms.msselect({'time':timerange})
			if chanid[i] and (type(chanid[i])==str):
				ms.selectchannel(1,int(chanid[i]),1,1)
			specdata=ms.getdata(['data','time','axis_info'],ifraxis=True)
			if pol=='RR':
				spec=specdata['data'][0,0,0,:]
			if pol=='LL':
				spec=specdata['data'][1,0,0,:]

			ms.selectinit(datadescid=int(spwid1))
			if bl and (type(bl)==str):
				ms.msselect({'baseline':bl})
			if plttimerange and (type(plttimerange)==str):
				ms.msselect({'time':plttimerange})
			if chanid[i] and (type(chanid[i])==str):
				ms.selectchannel(1,int(chanid[i]),1,1)
			specdata_plt=ms.getdata(['data','time','axis_info'],ifraxis=True)
			if pol=='RR':
				spec_plt=specdata_plt['data'][0,0,0,:]
			if pol=='LL':
				spec_plt=specdata_plt['data'][1,0,0,:]

			spec_plt_smooth=task_smooth.smooth(spec_plt,timoffset,windowlen,windowtype)
			spec_smooth=task_smooth.smooth(spec,timoffset,windowlen,windowtype)
			spec_plt_amp=np.absolute(spec_plt)
			spec_plt_smooth_amp=np.absolute(spec_plt_smooth)
			#spec_plt_amp=sqrt(spec_plt.real**2+spec_plt.imag**2)
			#spec_plt_smooth_amp=sqrt(spec_plt_smooth.real**2+spec_plt_smooth.imag**2)
			#print len(spec)
			#print len(spec_smooth)
			#print type(spec)
			#print spec[0]
			sp1 = fft(spec)
			sp2 = fft(spec_smooth)
			sp3 = sp1-sp2

			freq1 = fftfreq(len(sp1), d=0.001)
			freq2 = fftfreq(len(sp2), d=0.001)
			freq3 = fftfreq(len(sp3), d=0.001)
			freq1_index = np.argsort(freq1)
			freq2_index = np.argsort(freq2)
			freq3_index = np.argsort(freq3)
			#print min(freq1),max(freq1)

			subplot(j)
			plot(spec_plt_amp)
			plot(spec_plt_smooth_amp)
			title("Signal vs Time")
			j=j+1
			subplot(j)
			plot(freq1[freq1_index],log10(sp1[freq1_index]))
			plot(freq2[freq2_index],log10(sp2[freq2_index]))
			ylim([0,6])
			title("FFT signal vs Frequency")
			j=j+1
			subplot(j)
			#plot(subspec)
			plot(freq3[freq3_index],log10(sp3[freq3_index]))
			ylim([0,6])
			title("FFT smoothed signal vs Frequency")
			j=j+1
			#print "number of original data points: ",len(spec)
			#print "number of smoothed data points: ",len(spec_smooth)
		
		status = str(raw_input("Confirm to use current parameters? (y/n/abort) "))
		if status == 'y':
				flag1 = str(raw_input("Smooth all the channels and time range? (y/n) " ))
				if flag1 == 'y':
					smtimerange=''
					smspw=''
					splitsel=False
				else:
					print 'confirm using window width: ',windowlen
					print 'confirm using window type: ', windowtype
					smtimerange = str(raw_input("Please specify the time range for smoothing (HH:MM:SS) :" ))
					smspw = str(raw_input("Please specify spectral window and channel (e.g. 0:0~1032) :" ))
					splitsel = True
				break
		elif status =='abort':
				print 'Abort background subtraction.'
				sys.exit()

	if not outputvis:		
		outputvis=str(timoffset)+'_'+str(windowlen)+'_'+str(windowtype)+'.ms'
		print "Generating output: ",outputvis


	ms.close()
	result1 = subvs(vis,outputvis,smtimerange,smspw,timoffset,windowlen,windowtype,'','','')

def subvs(vis='', outputvis='', timerange='', spw='',timoffset=4, windowlen=5, windowtype='hamming', splitsel=True, reverse=False, overwrite=False):
	"""Perform vector subtraction for visibilities
    Keyword arguments:
	vis -- Name of input visibility file (MS)
				default: none; example: vis='ngc5921.ms'
	outputvis -- Name of output uv-subtracted visibility file (MS)
				default: none; example: outputvis='ngc5921_src.ms'
	timerange -- Time range of performing the UV subtraction:
				default='' means all times.  examples:
				timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
				timerange = 'hh:mm:ss~hh:mm:ss'
	spw -- Select spectral window/channel.
	windowlen -- Specify the width of window for smoothing
	windowtype --The type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
				flat window will produce a moving average smoothing.
	splitsel -- True or False. default = False. If splitsel = False, then the entire input
				measurement set is copied as the output measurement set (outputvis), with 
				background subtracted at selected timerange and spectral channels. 
				If splitsel = True,then only the selected timerange and spectral channels 
				are copied into the output measurement set (outputvis).
	reverse -- True or False. default = False. If reverse = False, then the times indicated
				by subtime1 and/or subtime2 are treated as background and subtracted; If reverse
				= True, then reverse the sign of the background-subtracted data. The option can 
				be used for mapping absorptive structure.
	overwrite -- True or False. default = False. If overwrite = True and
				outputvis already exists, the selected subtime and spw in the 
				output measurment set will be replaced with background subtracted 
				visibilities

    """
	# Get the time and frequency axis of the input ms
	# Open the ms and plot dynamic spectrum
	print 'using window length: ',windowlen
	print 'using window type: ', windowtype
	ms.open(vis,nomodify=True)
	# ms.selectinit(datadescid=0)
	timfreq=ms.getdata(['time','axis_info'],ifraxis=True)
	tim = timfreq['time']
	# check timerange input; default: entire timerange
	if timerange and (type(timerange)==str):
		[btimeo,etimeo]=timerange.split('~')
		btimeosec=qa.getvalue(qa.convert(qa.totime(btimeo),'s'))
		etimeosec=qa.getvalue(qa.convert(qa.totime(etimeo),'s'))
		timebinosec=etimeosec-btimeosec
		if timebinosec < 0:
			raise Exception, 'Negative timebin! Please check the "timerange" parameter.'
		else:
			casalog.post('Selected timerange: '+timerange+' as the time for UV subtraction.')
	else:
		casalog.post('Output timerange not specified, using the entire timerange')
		timerange=str(qa.time(qa.quantity(tim[0],'s'),prec=8)[0])+'~'+str(qa.time(qa.quantity(tim[-1],'s'),prec=8)[0])
		print 'Output timerange not specified, using the entire timerange',timerange
	# check spectral window input; default: entire channels of spectral window 0
	if spw and (type(spw)==str):
		[spwid,chanran]=spw.split(':')
		[bchan,echan]=chanran.split('~')
		nchan=int(echan)-int(bchan)+1
	else:
		casalog.post('spw not specified, use all frequency channels')
		freq = timfreq['axis_info']['freq_axis']['chan_freq'].flatten()
		nchan=len(freq)
		spwid='0'
		bchan='0'
		echan=str(nchan-1)
		print 'spw not specified, use all frequency channels', spwid+':'+bchan+'~'+str(nchan-1)
	
	ntimergn=len(timerange)
	# To avoid memory error, split the channel into smaller segements for smoothing
	cellstep=2
	chancell=int(nchan/cellstep)
	l=range(nchan)
	chunks=[l[x:x+cellstep] for x in xrange(0, len(l), cellstep)]
	#spwrange='0:0~'+str(chancell)
	ms.close()

	if not (timoffset and (type(timoffset)==int)):
		timoffset=int(4)

	for i in range(len(chunks)):
		spwrange=spwid+':'+str(int(bchan)+min(chunks[i]))+'~'+str(int(bchan)+max(chunks[i]))
		print 'Subtracting visibility from spectral range: ', spwrange
		result2 = task_subvs_lv1.subvs(vis,outputvis,timerange,spwrange,timoffset,windowlen,windowtype,splitsel,False,True)

