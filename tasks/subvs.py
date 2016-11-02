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
import task_subvs
def subvs(vis='', outputvis='', timerange='', spw='', mode='linear', subtime1='', subtime2='', smoothaxis='time', smoothtype='flat', smoothwidth=5, splitsel=True, reverse=False, overwrite=False):

        """Vector-subtraction in UV using selected time ranges and spectral channels as background

    Subvs is a task to do UV vector-subtraction, by selecting time ranges 
    in the data as background. Subvs can be used to subtract the background
    continuum emission to separate the time-dependent emission, e.g. solar 
    coherent radio bursts. 
    
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
                mode = 'lowpass': act as a lowpass filter---smooth the data using different smooth
                        types and smooth window size. Can be performed along either time 
                        or frequency axis
                mode = 'highpass': act as a highpass filter---smooth the data first, and 
                        subtract the smoothed data from the original. Can be performed along either time
                        or frequency axis
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

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['outputvis'] = outputvis
        mytmp['timerange'] = timerange
        mytmp['spw'] = spw
        mytmp['mode'] = mode
        mytmp['subtime1'] = subtime1
        mytmp['subtime2'] = subtime2
        mytmp['smoothaxis'] = smoothaxis
        mytmp['smoothtype'] = smoothtype
        mytmp['smoothwidth'] = smoothwidth
        mytmp['splitsel'] = splitsel
        mytmp['reverse'] = reverse
        mytmp['overwrite'] = overwrite
	pathname="file:///afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks/"
	trec = casac.utils().torecord(pathname+'subvs.xml')

        casalog.origin('subvs')
        if trec.has_key('subvs') and casac.utils().verify(mytmp, trec['subvs']) :
	    result = task_subvs.subvs(vis, outputvis, timerange, spw, mode, subtime1, subtime2, smoothaxis, smoothtype, smoothwidth, splitsel, reverse, overwrite)

	else :
	  result = False
        return result
