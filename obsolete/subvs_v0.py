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
def subvs(vis='', outputvis='', timerange='', spw='', subtime1='', subtime2='', splitsel=True, overwrite=False):

        """Vector-subtraction in UV using selected time ranges and spectral channels as background

    Subvs is a task to do UV vector-subtraction, by selecting time ranges 
    in the data as background. Subvs can be used to subtract the background
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
	subtime1 -- Select the time range 1 of data as the background.
               Visibilities will be vector-averaged in time1 before subtraction. 
               default = '' (will raise an exception); examples,
               timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
               Note: if YYYY/MM/DD is missing date, timerange defaults to the
               first day in the dataset
               timerange='09:14:0~09:54:0' picks 40 min on first day
               timerange='25:00:00~27:30:00' picks 1 hr to 3 hr 30min
               on next day
	subtime2 -- Select the time range 2 of data as the background.
               Visibilities will be vector-averaged in time2 before subtraction.
	       if specified, then linear-interpolated values based on time1 and
	       time2 will be applied to the times specified by "subtime" for 
	       subtraction. 
	       if not specified, then only "time1" is used as the background	
               default = '' (none); examples,
               timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
               Note: if YYYY/MM/DD is missing date, timerange defaults to the
               first day in the dataset
               timerange='09:14:0~09:54:0' picks 40 min on first day
               timerange='25:00:00~27:30:00' picks 1 hr to 3 hr 30min
               on next day
        splitsel -- True of False. default = False. If splitsel = False, then the entire input 
	       measurement set is copied as the output measurement set (outputvis), with 
	       background subtracted at selected timerange and spectral channels. 
	       If splitsel = True,then only the selected timerange and spectral channels 
	       are copied into the output measurement set (outputvis).
	overwrite -- True or False. default = False. If overwrite = True and outputvis
	       already exists, the selected subtime and spw in the already existing 
	       output measurement set will be replaced with subtracted visibilities

        """

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['outputvis'] = outputvis
        mytmp['timerange'] = timerange
        mytmp['spw'] = spw
        mytmp['subtime1'] = subtime1
        mytmp['subtime2'] = subtime2
        mytmp['splitsel'] = splitsel
        mytmp['overwrite'] = overwrite
	pathname="file:///Users/binchen/Dropbox/bc_python/casa_task/"
	trec = casac.utils().torecord(pathname+'subvs.xml')

        casalog.origin('subvs')
        if trec.has_key('subvs') and casac.utils().verify(mytmp, trec['subvs']) :
	    result = task_subvs.subvs(vis, outputvis, timerange, spw, subtime1, subtime2, splitsel, overwrite)

	else :
	  result = False
        return result
