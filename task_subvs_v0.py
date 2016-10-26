import os
import shutil
import time
import stat
import numpy as np
from taskinit import ms, cb, tb, qa, casalog, write_history
import pdb

def subvs(vis=None,outputvis=None,timerange=None,spw=None,
          subtime1=None,subtime2=None,splitsel=True,overwrite=False):
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
           default = '' all the spectral channels. Example: spw='0:1~20'
    subtime1 -- Time range 1 of the background to be subtracted from the data 
                 default='' means all times.  format:
                 timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
                 timerange = 'hh:mm:ss~hh:mm:ss'
    subtime2 -- Time range 2 of the backgroud to be subtracted from the data
                 default='' means all times.  examples:
                 timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
                 timerange = 'hh:mm:ss~hh:mm:ss'
    splitsel -- True of False. default = False. If splitsel = False, then the entire input
            measurement set is copied as the output measurement set (outputvis), with 
            background subtracted at selected timerange and spectral channels. 
            If splitsel = True,then only the selected timerange and spectral channels 
            are copied into the output measurement set (outputvis).
    overwrite -- True or False. default = False. If overwrite = True and
                outputvis already exists, the selected subtime and spw in the 
                output measurment set will be replaced with background subtracted 
                visibilities

    """
    #check the visbility ms
    if not outputvis or outputvis.isspace():
        raise ValueError, 'Please specify outputvis'

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

    #define and check the time ranges
    if subtime1 and (type(subtime1)==str):
        [bsubtime1,esubtime1]=subtime1.split('~')
        bsubtime1sec=qa.getvalue(qa.convert(qa.totime(bsubtime1),'s'))
        esubtime1sec=qa.getvalue(qa.convert(qa.totime(esubtime1),'s'))
        timebin1sec=esubtime1sec-bsubtime1sec
        if timebin1sec < 0:
            raise Exception, 'Negative timebin! Please check the "subtime1" parameter.'
        casalog.post('Selected timerange 1: '+subtime1+' as background for uv subtraction.')
    else:
        raise Exception, 'Please enter at least one timerange as the background'

    if subtime2 and (type(subtime2)==str):
        [bsubtime2,esubtime2]=subtime2.split('~')
        bsubtime2sec=qa.getvalue(qa.convert(qa.totime(bsubtime2),'s'))
        esubtime2sec=qa.getvalue(qa.convert(qa.totime(esubtime2),'s'))
        timebin2sec=esubtime2sec-bsubtime2sec
        if timebin2sec < 0:
            raise Exception, 'Negative timebin! Please check the "subtime2" parameter.'
        timebin2=str(timebin2sec)+'s'
        casalog.post('Selected timerange 2: '+subtime2+' as background for uv subtraction.')
        #plus 1s is to ensure averaging over the entire timerange
    else:
        casalog.post('Timerange 2 not selected, using only timerange 1 as background')

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
    else:
        casalog.post('spw not specified, use all frequency channels')

    #Select the background indicated by subtime1
    ms.open(vis,nomodify=True)
    #Select the spw id
    ms.msselect({'time':subtime1})
    if spw and (type(spw)==str):
        ms.selectinit(datadescid=int(spwid))
        nchan=int(echan)-int(bchan)+1
        ms.selectchannel(nchan,int(bchan),1,1)
    rec1=ms.getdata(['data','time','axis_info'],ifraxis=True)
    #print 'shape of the frequency matrix ',rec1['axis_info']['freq_axis']['chan_freq'].shape
    sz1=rec1['data'].shape
    print 'dimension of selected background 1', rec1['data'].shape
    #the data shape is (n_pol,n_channel,n_baseline,n_time), no need to reshape
    #rec1['data']=rec1['data'].reshape(sz1[0],sz1[1],sz1[2],nspw,sz1[3]/nspw,order='F')
    #print 'reshaped rec1 ', rec1['data'].shape
    rec1avg=np.average(rec1['data'],axis=3)
    casalog.post('Averaging the visibilities in subtime1: '+subtime1)
    ms.close()
    if subtime2 and (type(subtime2)==str):
        ms.open(vis,nomodify=True)
        #Select the spw id
        ms.msselect({'time':subtime2})
        if spw and (type(spw)==str):
            ms.selectinit(datadescid=0)
            nchan=int(echan)-int(bchan)+1
            ms.selectchannel(nchan,int(bchan),1,1)
        rec2=ms.getdata(['data','time','axis_info'],ifraxis=True)
        sz2=rec2['data'].shape
        print 'dimension of selected background 2', rec2['data'].shape
        #rec2['data']=rec2['data'].reshape(sz2[0],sz2[1],sz2[2],nspw,sz2[3]/nspw,order='F')
        #print 'reshaped rec1 ', rec2['data'].shape
        rec2avg=np.average(rec2['data'],axis=3)
        ms.close()
        casalog.post('Averaged the visibilities in subtime2: '+subtime2)

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

    if subtime1 and (not subtime2):
        casalog.post('Only "subtime1" is defined, subtracting background defined in subtime1: '+subtime1)
        for i in range(t_rows):
            orec['data'][:,:,:,i]-=rec1avg
    if subtime1 and subtime2 and (type(subtime2)==str):
        casalog.post('Both subtime1 and subtime2 are specified, doing linear interpolation between "subtime1" and "subtime2"') 
        t1=(np.amax(rec1['time'])+np.amin(rec1['time']))/2.
        t2=(np.amax(rec2['time'])+np.amin(rec2['time']))/2.
        touts=orec['time']
        print 't1: ', qa.time(qa.quantity(t1,'s'), form='ymd')
        print 't2: ', qa.time(qa.quantity(t2,'s'), form='ymd')
        for i in range(t_rows):
            tout=touts[i]
            if tout > np.amax([t1,t2]):
                tout=np.amax([t1,t2])
            elif tout < np.amin([t1,t2]):
                tout=np.amin([t1,t2])
            orec['data'][:,:,:,i]-=(rec2avg-rec1avg)*(tout-t1)/(t2-t1)+rec1avg
    
    #orec['data']=orec['data'].reshape(szo[0],szo[1],szo[2],szo[3],order='F')
    #put the modified data back into the output visibility set
    del orec['time']
    del orec['axis_info']
    ms.putdata(orec)
    ms.close()




