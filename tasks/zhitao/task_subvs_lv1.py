import os
import shutil
import time
import stat
import numpy as np
import task_smooth
from taskinit import ms, cb, tb, qa, casalog, write_history
#import pdb

def subvs(vis=None,outputvis=None,timerange=None,spw=None,timoffset=4,
          windowlen=10,windowtype='flat',splitsel=True,reverse=False,overwrite=False):
    """Vector-subtraction in UV using selected time ranges and spectral channels as background
    subvs is a function to do UV vector-subtraction. By selecting gliding averaging window,
    only the low-frequency signals corresponding to the background continuum emission remain.
    As a result, a uv subtraction of original dynamic spectrum from the background can improve
    fine stucture such as fibers. Subvs can be used to subtract the background
    continuum emission to separate the time-dependent emission, e.g. solar 
    coherent radio bursts.  
    
        Keyword arguments:
        vis -- Name of input visibility file
                default: none; example: vis='sun_type3.ms'
        outputvis -- Name of output visibility file
                default: none; example: outputvis='sun_type3.sub.ms'
    timerange -- Select the time range in the data to be subtracted from.
               timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
               Note: if YYYY/MM/DD is missing date, timerange defaults to the
               first day in the dataset
               timerange='09:14:0~09:54:0' picks 40 min on first day
               timerange='25:00:00~27:30:00' picks 1 hr to 3 hr 30min
               on next day
    spw -- Select spectral window/channel.
            default = '' all the spectral channels. Example: spw='0:1~20'
    timoffset -- After the convolution, each single channel of smoothed uv vectors 
                in time series are misaligned from original data. Setting the 
                timoffset allow you to shift by the specifed amount of data points. 
    windowlen -- Choose the width the gliding window for smoothing.             
    windowtype -- Choose the type o gliding window. Options available are  'flat', 
                'hanning', 'hamming', 'bartlett', 'blackman'
    splitsel -- True of False. default = False. If splitsel = False, then the entire input 
           measurement set is copied as the output measurement set (outputvis), with 
           background subtracted at selected timerange and spectral channels. 
           If splitsel = True,then only the selected timerange and spectral channels 
           are copied into the output measurement set (outputvis).
    reverse -- True or False. default = False. If reverse = False, then the times indicated
            by subtime1 and/or subtime2 are treated as background and subtracted; If reverse
            = True, then reverse the sign of the background-subtracted data. The option can 
            be used for mapping absorptive structure.
    overwrite -- True or False. default = False. If overwrite = True and outputvis
           already exists, the selected subtime and spw in the already existing 
           output measurement set will be replaced with subtracted visibilities

        """
    #check the visbility ms
    if not outputvis or outputvis.isspace():
        raise ValueError, 'Please specify outputvis.'

    if os.path.exists(outputvis):
        if overwrite:
            print "The already existing output measurement set will be updated."
        else:
            raise ValueError, "Output MS %s already exists - will not overwrite." % outputvis
    else:
        if not splitsel:
            shutil.copytree(vis,outputvis)
        else:
            ms.open(vis,nomodify=True)
            ms.split(outputvis,spw=spw,time=timerange,whichcol='DATA')
            ms.close()
    # check and specify time range and channel
    if timerange and (type(timerange)==str):
        [btimeo,etimeo]=timerange.split('~')
        btimeosec=qa.getvalue(qa.convert(qa.totime(btimeo),'s'))
        etimeosec=qa.getvalue(qa.convert(qa.totime(etimeo),'s'))
        timebinosec=etimeosec-btimeosec
        if timebinosec < 0:
            raise Exception, 'Negative timebin! Please check the "timerange" parameter.'
        casalog.post('Selected timerange: '+timerange+' as the time for UV subtraction.')
    else:
        casalog.post('Output timerange not specified, using the entire timerange')

    if spw and (type(spw)==str):
        [spwid,chanran]=spw.split(':')
        [bchan,echan]=chanran.split('~')
        nchan=int(echan)-int(bchan)+1
    else:
        casalog.post('spw not specified, use all frequency channels')

    #select data range to be smoothed
    ms.open(vis,nomodify=True)
    #Select the spw id
    ms.msselect({'time':timerange})
    ms.selectinit(datadescid=int(spwid))
    ms.selectchannel(nchan,int(bchan),1,1)
    rec=ms.getdata(['data','time','axis_info'],ifraxis=True)
    #print 'shape of the frequency matrix ',rec1['axis_info']['freq_axis']['chan_freq'].shape
    sz=rec['data'].shape
    print 'dimension of selected background for smoothing', rec['data'].shape
    #the data shape is (n_pol,n_channel,n_baseline,n_time), no need to reshape
    #rec1['data']=rec1['data'].reshape(sz1[0],sz1[1],sz1[2],nspw,sz1[3]/nspw,order='F')
    #print 'reshaped rec1 ', rec1['data'].shape
    if not (timoffset and (type(timoffset)==int)):
        timoffset=int(4)

    for i in range(rec['data'].shape[0]):
     for j in range(rec['data'].shape[1]):
     	for k in range (rec['data'].shape[2]):
    		rec['data'][i,j,k,:]=task_smooth.smooth(rec['data'][i,j,k,:],timoffset,windowlen,windowtype)
    
    casalog.post('Smoothing the visibilities in timerange: '+timerange)
    ms.close()


    #do UV subtraction, according to timerange and spw    
    ms.open(outputvis,nomodify=False)
    if not splitsel:
        #outputvis is identical to input visibility, do the selection
        if timerange and (type(timerange==str)):
            ms.msselect({'time':timerange})
        if spw and (type(spw)==str):
            ms.selectinit(datadescid=int(spwid))
            nchan=int(echan)-int(bchan)+1
            ms.selectchannel(nchan,int(bchan),1,1)
    else:
        #outputvis is splitted, selections have already applied, select all the data
        ms.selectinit(datadescid=0)
    orec=ms.getdata(['data','time','axis_info'],ifraxis=True)
    b_rows=orec['data'].shape[2]
    nchan=orec['data'].shape[1]
    #szo=orec['data'].shape
    print 'dimension of output data', orec['data'].shape
    #orec['data']=orec['data'].reshape(szo[0],szo[1],szo[2],nspw,szo[3]/nspw,order='F')
    #print 'reshaped rec1 ', orec['data'].shape
    t_rows=orec['data'].shape[3]
    casalog.post('Number of baselines: '+str(b_rows))
    casalog.post('Number of spectral channels: '+str(nchan))
    casalog.post('Number of time pixels: '+str(t_rows))

    casalog.post('Subtracting background defined in timerange: '+timerange)

    for i in range(t_rows):
    	orec['data'][:,:,:,i]-=rec['data'][:,:,:,i]
        if reverse:
    		orec['data'][:,:,:,i]=-orec['data'][:,:,:,i]
    del orec['time']
    del orec['axis_info']
    ms.putdata(orec)
    ms.close()
