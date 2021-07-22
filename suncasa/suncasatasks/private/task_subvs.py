import os
import shutil
import numpy as np
from suncasa.utils import signal_utils as su
import sys

if sys.version_info.major > 2:
    from casatools import ms, quanta, msmetadata
    from casatasks import casalog

    casalog.showconsole(True)
    datams = ms()
    ms_in = ms()
    datamsmd = msmetadata()
    qa = quanta()
else:
    from taskinit import ms, qa, mstool, msmdtool, casalog

    datams = mstool()
    ms_in = mstool()
    datamsmd = msmdtool()


# from taskinit import *
# from callibrary import *
# import pdb

def subvs(vis=None, outputvis=None, timerange='', spw='',
           mode='linear', subtime1='', subtime2='',
           smoothaxis='time', smoothtype='flat', smoothwidth='5',
           splitsel=True, reverse=False, overwrite=False):
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
    mode -- operation mode
            default 'linear' 
                mode = 'linear': use a linear fit for the background to be subtracted
                mode = 'lowpass': act as a lowpass filter---smooth the data using different
                        smooth types and window sizes. Can be performed along either time
                        or frequency axis
                mode = 'highpass': act as a highpass filter---smooth the data first, and 
                        subtract the smoothed data from the original. Can be performed along
                        either time or frequency axis
            mode = 'linear' expandable parameters:
                subtime1 -- Time range 1 of the background to be subtracted from the data 
                             default='' means all times.  format:
                             timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
                             timerange = 'hh:mm:ss~hh:mm:ss'
                subtime2 -- Time range 2 of the backgroud to be subtracted from the data
                             default='' means all times.  examples:
                             timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
                             timerange = 'hh:mm:ss~hh:mm:ss'
            mode = 'lowpass' or 'highpass' expandable parameters:
                smoothaxis -- axis of smooth
                    Default: 'time'
                    smoothaxis = 'time': smooth is along the time axis
                    smoothaxis = 'freq': smooth is along the frequency axis
                smoothtype -- type of the smooth depending on the convolving kernel
                    Default: 'flat'
                    smoothtype = 'flat': convolving kernel is a flat rectangle,
                            equivalent to a boxcar moving smooth
                    smoothtype = 'hanning': Hanning smooth kernel. See numpy.hanning
                    smoothtype = 'hamming': Hamming smooth kernel. See numpy.hamming
                    smoothtype = 'bartlett': Bartlett smooth kernel. See numpy.bartlett
                    smoothtype = 'blackman': Blackman smooth kernel. See numpy.blackman
                smoothwidth -- width of the smooth kernel
                    Default: 5
                    Examples: smoothwidth=5, meaning the width is 5 pixels
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
    # check the visbility ms
    casalog.post('input parameters:')
    casalog.post('vis: ' + vis)
    casalog.post('outputvis: ' + outputvis)
    if mode not in ['linear']:
        casalog.post('smoothaxis: ' + smoothaxis)
        casalog.post('smoothtype: ' + smoothtype)
        casalog.post('smoothwidth: ' + str(smoothwidth))
    if not outputvis or outputvis.isspace():
        raise ValueError('Please specify outputvis')

    if os.path.exists(outputvis):
        if overwrite:
            print("The already existing output measurement set will be updated.")
        else:
            raise ValueError("Output MS {} already exists - will not overwrite.".format(outputvis))
    else:
        if not splitsel:
            shutil.copytree(vis, outputvis)
        else:
            ms_in.open(vis, nomodify=True)
            ms_in.split(outputvis, spw=spw, time=timerange, whichcol='DATA')
            ms_in.close()

    if timerange and (type(timerange) == str):
        [btimeo, etimeo] = timerange.split('~')
        btimeosec = qa.getvalue(qa.convert(qa.totime(btimeo), 's'))
        etimeosec = qa.getvalue(qa.convert(qa.totime(etimeo), 's'))
        timebinosec = etimeosec - btimeosec
        if timebinosec < 0:
            raise Exception('Negative timebin! Please check the "timerange" parameter.')
        casalog.post('Selected timerange: ' + timerange + ' as the time for UV subtraction.')
    else:
        casalog.post('Output timerange not specified, using the entire timerange')

    if spw and (type(spw) == str):
        spwlist = spw.split(';')
    else:
        casalog.post('spw not specified, use all frequency channels')

    # read the output data

    datams.open(outputvis, nomodify=False)

    datamsmd.open(outputvis)
    spwinfod = datams.getspectralwindowinfo()
    spwinfok = spwinfod.keys()
    spwinfok = sorted(spwinfok, key=int)
    spwinfol = [spwinfod[k] for k in spwinfok]
    for s, spi in enumerate(spwinfol):
        print('processing spectral window {}'.format(spi['SpectralWindowId']))
        datams.selectinit(reset=True)
        staql = {'time': '', 'spw': ''}
        if not splitsel:
            # outputvis is identical to input visibility, do the selection
            if timerange and (type(timerange == str)):
                staql['time'] = timerange
            if spw and (type(spw) == str):
                staql['spw'] = spwlist[s]
            if not spw and not timerange:
                # data selection is not made
                print('selecting all spws and times')
                staql['spw'] = str(spi['SpectralWindowId'])
        else:
            # outputvis is splitted, selections have already applied, select all the data
            print('split the selected spws and times')
            staql['spw'] = str(spi['SpectralWindowId'])
        datams.msselect(staql)
        orec = datams.getdata(['data', 'time', 'axis_info'], ifraxis=True)
        npol, nchan, nbl, ntim = orec['data'].shape
        print('dimension of output data', orec['data'].shape)
        casalog.post('Number of baselines: ' + str(nbl))
        casalog.post('Number of spectral channels: ' + str(nchan))
        casalog.post('Number of time pixels: ' + str(ntim))

        try:
            if mode == 'linear':
                # define and check the background time ranges
                if subtime1 and (type(subtime1) == str):
                    [bsubtime1, esubtime1] = subtime1.split('~')
                    bsubtime1sec = qa.getvalue(qa.convert(qa.totime(bsubtime1), 's'))
                    esubtime1sec = qa.getvalue(qa.convert(qa.totime(esubtime1), 's'))
                    timebin1sec = esubtime1sec - bsubtime1sec
                    if timebin1sec < 0:
                        raise Exception('Negative timebin! Please check the "subtime1" parameter.')
                    casalog.post('Selected timerange 1: ' + subtime1 + ' as background for uv subtraction.')
                else:
                    raise Exception('Please enter at least one timerange as the background')
                if subtime2 and (type(subtime2) == str):
                    [bsubtime2, esubtime2] = subtime2.split('~')
                    bsubtime2sec = qa.getvalue(qa.convert(qa.totime(bsubtime2), 's'))
                    esubtime2sec = qa.getvalue(qa.convert(qa.totime(esubtime2), 's'))
                    timebin2sec = esubtime2sec - bsubtime2sec
                    if timebin2sec < 0:
                        raise Exception('Negative timebin! Please check the "subtime2" parameter.')
                    timebin2 = str(timebin2sec) + 's'
                    casalog.post('Selected timerange 2: ' + subtime2 + ' as background for uv subtraction.')
                    # plus 1s is to ensure averaging over the entire timerange
                else:
                    casalog.post('Timerange 2 not selected, using only timerange 1 as background')

                # Select the background indicated by subtime1
                ms_in.open(vis, nomodify=True)
                # Select the spw id
                # ms_in.msselect({'time': subtime1})
                staql0 = {'time': subtime1, 'spw': ''}
                if spw and (type(spw) == str):
                    staql0['spw'] = spwlist[s]
                else:
                    staql0['spw'] = staql['spw']
                ms_in.msselect(staql0)
                rec1 = ms_in.getdata(['data', 'time', 'axis_info'], ifraxis=True)
                # print('shape of the frequency matrix ',rec1['axis_info']['freq_axis']['chan_freq'].shape)
                sz1 = rec1['data'].shape
                print('dimension of selected background 1', rec1['data'].shape)
                # the data shape is (n_pol,n_channel,n_baseline,n_time), no need to reshape
                # rec1['data']=rec1['data'].reshape(sz1[0],sz1[1],sz1[2],nspw,sz1[3]/nspw,order='F')
                # print('reshaped rec1 ', rec1['data'].shape)
                rec1avg = np.average(rec1['data'], axis=3)
                casalog.post('Averaging the visibilities in subtime1: ' + subtime1)
                ms_in.close()
                if subtime2 and (type(subtime2) == str):
                    ms_in.open(vis, nomodify=True)
                    # Select the spw id
                    staql0 = {'time': subtime2, 'spw': ''}
                    if spw and (type(spw) == str):
                        staql0['spw'] = spwlist[s]
                    else:
                        staql0['spw'] = staql['spw']
                    ms_in.msselect(staql0)
                    rec2 = ms_in.getdata(['data', 'time', 'axis_info'], ifraxis=True)
                    sz2 = rec2['data'].shape
                    print('dimension of selected background 2', rec2['data'].shape)
                    # rec2['data']=rec2['data'].reshape(sz2[0],sz2[1],sz2[2],nspw,sz2[3]/nspw,order='F')
                    # print('reshaped rec1 ', rec2['data'].shape)
                    rec2avg = np.average(rec2['data'], axis=3)
                    ms_in.close()
                    casalog.post('Averaged the visibilities in subtime2: ' + subtime2)
                if subtime1 and (not subtime2):
                    casalog.post('Only "subtime1" is defined, subtracting background defined in subtime1: ' + subtime1)
                    t1 = (np.amax(rec1['time']) + np.amin(rec1['time'])) / 2.
                    print('t1: ', qa.time(qa.quantity(t1, 's'), form='ymd', prec=10))
                    for i in range(ntim):
                        orec['data'][:, :, :, i] -= rec1avg
                        if reverse:
                            orec['data'][:, :, :, i] = -orec['data'][:, :, :, i]
                if subtime1 and subtime2 and (type(subtime2) == str):
                    casalog.post(
                        'Both subtime1 and subtime2 are specified, doing linear interpolation between "subtime1" and "subtime2"')
                    t1 = (np.amax(rec1['time']) + np.amin(rec1['time'])) / 2.
                    t2 = (np.amax(rec2['time']) + np.amin(rec2['time'])) / 2.
                    touts = orec['time']
                    print('t1: ', qa.time(qa.quantity(t1, 's'), form='ymd', prec=10))
                    print('t2: ', qa.time(qa.quantity(t2, 's'), form='ymd', prec=10))
                    for i in range(ntim):
                        tout = touts[i]
                        if tout > np.amax([t1, t2]):
                            tout = np.amax([t1, t2])
                        elif tout < np.amin([t1, t2]):
                            tout = np.amin([t1, t2])
                        orec['data'][:, :, :, i] -= (rec2avg - rec1avg) * (tout - t1) / (t2 - t1) + rec1avg
                        if reverse:
                            orec['data'][:, :, :, i] = -orec['data'][:, :, :, i]
            elif mode == 'highpass':
                if smoothtype != 'flat' and smoothtype != 'hanning' and smoothtype != 'hamming' and smoothtype != 'bartlett' and smoothtype != 'blackman':
                    raise Exception('Unknown smoothtype ' + str(smoothtype))
                if smoothaxis == 'time':
                    if smoothwidth <= 0 or smoothwidth >= ntim:
                        raise Exception('Specified smooth width is <=0 or >= the total number of ' + smoothaxis)
                    else:
                        for i in range(orec['data'].shape[0]):
                            for j in range(nchan):
                                for k in range(nbl):
                                    orec['data'][i, j, k, :] -= su.smooth(orec['data'][i, j, k, :],
                                                                          smoothwidth,
                                                                          smoothtype)
                if smoothaxis == 'freq':
                    if smoothwidth <= 0 or smoothwidth >= nchan:
                        raise Exception('Specified smooth width is <=0 or >= the total number of ' + smoothaxis)
                    else:
                        for i in range(orec['data'].shape[0]):
                            for j in range(nbl):
                                for k in range(ntim):
                                    orec['data'][i, :, j, k] -= su.smooth(orec['data'][i, :, j, k],
                                                                          smoothwidth,
                                                                          smoothtype)
            elif mode == 'lowpass':
                if smoothtype != 'flat' and smoothtype != 'hanning' and smoothtype != 'hamming' and smoothtype != 'bartlett' and smoothtype != 'blackman':
                    raise Exception('Unknown smoothtype ' + str(smoothtype))
                if smoothaxis == 'time':
                    if smoothwidth <= 0 or smoothwidth >= ntim:
                        raise Exception('Specified smooth width is <=0 or >= the total number of ' + smoothaxis)
                    else:
                        for i in range(orec['data'].shape[0]):
                            for j in range(nchan):
                                for k in range(nbl):
                                    orec['data'][i, j, k, :] = su.smooth(orec['data'][i, j, k, :],
                                                                         smoothwidth,
                                                                         smoothtype)
                if smoothaxis == 'freq':
                    if smoothwidth <= 0 or smoothwidth >= nchan:
                        raise Exception('Specified smooth width is <=0 or >= the total number of ' + smoothaxis)
                    else:
                        for i in range(orec['data'].shape[0]):
                            for j in range(nbl):
                                for k in range(ntim):
                                    orec['data'][i, :, j, k] = su.smooth(orec['data'][i, :, j, k],
                                                                         smoothwidth,
                                                                         smoothtype)
            else:
                raise Exception('Unknown mode' + str(mode))
        except Exception as instance:
            print('*** Error ***', instance)

        # orec['data']=orec['data'].reshape(szo[0],szo[1],szo[2],szo[3],order='F')
        # put the modified data back into the output visibility set
        del orec['time']
        del orec['axis_info']
        # ms_in.open(outputvis,nomodify=False)
        # if not splitsel:
        # outputvis is identical to input visibility, do the selection
        #    if timerange and (type(timerange==str)):
        #        datams.msselect({'time':timerange})
        #    if spw and (type(spw)==str):
        #        datams.selectinit(datadescid=int(spwid))
        #        nchan=int(echan)-int(bchan)+1
        #        datams.selectchannel(nchan,int(bchan),1,1)
        #    if not spw and not timerange:
        # data selection is not made
        #        datams.selectinit(datadescid=0)
        # else:
        # outputvis is splitted, selections have already applied, select all the data
        #    datams.selectinit(datadescid=0)
        datams.putdata(orec)
    datams.close()
    datamsmd.done()
