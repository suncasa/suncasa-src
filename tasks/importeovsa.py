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
import task_importeovsa
def importeovsa(idbfiles='', timebin='0s', width=1, visprefix='', nocreatms=True, doconcat=False):

        """Import EOVSA idb file(s) to a measurement set or multiple measurement set

            Imports an arbitrary number of EOVSA idb-format data sets into
            a casa measurement set.  If more than one band is present, they
            will be put in the same measurement set but in a separate spectral
            window.

            Detailed Keyword arguments:

            idbfiles -- Name of input EOVSA idb file(s) or observation time range
                    default: none.  Must be supplied
                    example: idbfiles = 'IDB20160524000518'
                    example: idbfiles=['IDB20160524000518','IDB20160524000528']
                    example: idbfiles=['2016-08-09 12:10:00','2016-08-09 12:50:00']

            visprefix -- Prefix of vis names (may include the path)
            default: none;
            example: visprefix='sun/']


            --- Data Selection ---

            nocreatms -- If copying a new MS file instead of create one from MS simulator.
            default: False

            modelms -- Name of the standard Measurement Set. IF modelms is not provided, use
            '/home/user/sjyu/20160531/ms/sun/SUN/SUN_20160531T142234-10m.1s.ms' as a standard MS.

            doconcat -- If outputing one single MS file


            --- Channel averaging parameter ---

            width -- Number of input channels to average to create an output
            channel. If a list is given, each bin will apply to one spw in
            the selection.
            default: 1 => no channel averaging.
            options: (int) or [int]

            example: chanbin=[2,3] => average 2 channels of 1st selected
            spectral window and 3 in the second one.


            --- Time averaging parameters ---

            timebin -- Bin width for time averaging. When timebin is greater than 0s,
            the task will average data in time. Flagged data will be included
            in the average calculation, unless the parameter keepflags is set to False.
            In this case only partially flagged rows will be used in the average.
            default: '0s'

            combine -- Let the timebin span across scan, state or both.
            State is equivalent to sub-scans. One scan may have several
            state ids. For ALMA MSs, the sub-scans are limited to about
            30s duration each. In these cases, the task will automatically
            add state to the combine parameter. To see the number of states
            in an MS, use the msmd tool. See help msmd.


        
        """

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['idbfiles'] = idbfiles
        mytmp['timebin'] = timebin
        mytmp['width'] = width
        mytmp['visprefix'] = visprefix
        mytmp['nocreatms'] = nocreatms
        mytmp['doconcat'] = doconcat
	pathname="file:///local/software/suncasa/tasks/"
	trec = casac.utils().torecord(pathname+'importeovsa.xml')

        casalog.origin('importeovsa')
        if trec.has_key('importeovsa') and casac.utils().verify(mytmp, trec['importeovsa']) :
	    result = task_importeovsa.importeovsa(idbfiles, timebin, width, visprefix, nocreatms, doconcat)

	else :
	  result = False
        return result
