#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
#from casac import *
import casac
import string
import time
import inspect
import gc
import numpy
from odict import odict
from types import * 
from task_importeovsa import importeovsa
class importeovsa_cli_:
    __name__ = "importeovsa"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (importeovsa_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'idbfiles':None, 'timebin':None, 'width':None, 'visprefix':None, 'nocreatms':None, 'doconcat':None, 'modelms':None, }


    def result(self, key=None):
	    #### and add any that have completed...
	    return None


    def __call__(self, idbfiles=None, timebin=None, width=None, visprefix=None, nocreatms=None, doconcat=None, modelms=None, ):

        """Import EOVSA idb file(s) to a measurement set or multiple measurement set

	Detailed Description: 

            Imports an arbitrary number of EOVSA idb-format data sets into
            a casa measurement set.  If more than one band is present, they
            will be put in the same measurement set but in a separate spectral
            window.
        
	Arguments :
		idbfiles:	Name of input EOVSA idb file(s) or observation time range.
		   Default Value: 

		timebin:	Bin width for time averaging
		   Default Value: 0s

		width:	Width of output channel relative to MS channel (# to average)
		   Default Value: 1

		visprefix:	Prefix of vis names (may include the path).
		   Default Value: 

		nocreatms:	If setting nocreatms True, will simulate a model measurement set for the first idb file and copy the model for the rest of idl files in list. If False, will simulate a new measurement set for every idbfile in list.
		   Default Value: True

		doconcat:	If concatenate multi casa measurement sets to one file.
		   Default Value: False

		modelms:	Name of input model measurement set file. If modelms is assigned, no simulation will start.
		   Default Value: 


	Example :


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
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=sys._getframe(len(inspect.stack())-1).f_globals
	#casac = self.__globals__['casac']
	casalog = self.__globals__['casalog']
	casa = self.__globals__['casa']
	#casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'importeovsa'
        self.__globals__['taskname'] = 'importeovsa'
        ###
        self.__globals__['update_params'](func=self.__globals__['taskname'],printtext=False,ipython_globals=self.__globals__)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults={}
	else:
	    function_signature_defaults=dict(zip(self.__call__.func_code.co_varnames[1:],self.__call__.func_defaults))
	useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
			if (key != 'self') :
			   useLocalDefaults = True

	myparams = {}
	if useLocalDefaults :
	   for item in function_signature_defaults.iteritems():
	       key,val = item
	       keyVal = eval(key)
	       exec('myparams[key] = keyVal')
	       self.parameters[key] = keyVal
	       if (keyVal == None):
	           exec('myparams[key] = '+ key + ' = self.itsdefault(key)')
		   keyVal = eval(key)
		   if(type(keyVal) == dict) :
                      if len(keyVal) > 0 :
		         exec('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
		      else :
		         exec('myparams[key] = ' + key + ' = {}')
	 
        else :
            print ''

            myparams['idbfiles'] = idbfiles = self.parameters['idbfiles']
            myparams['timebin'] = timebin = self.parameters['timebin']
            myparams['width'] = width = self.parameters['width']
            myparams['visprefix'] = visprefix = self.parameters['visprefix']
            myparams['nocreatms'] = nocreatms = self.parameters['nocreatms']
            myparams['doconcat'] = doconcat = self.parameters['doconcat']
            myparams['modelms'] = modelms = self.parameters['modelms']


	result = None

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
        mytmp['modelms'] = modelms
	pathname="file:///local/software/suncasa/tasks/"
	trec = casac.casac.utils().torecord(pathname+'importeovsa.xml')

        casalog.origin('importeovsa')
	try :
          #if not trec.has_key('importeovsa') or not casac.casac.utils().verify(mytmp, trec['importeovsa']) :
	    #return False

          casac.casac.utils().verify(mytmp, trec['importeovsa'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('importeovsa', 'importeovsa.last', myparams, self.__globals__,scriptstr=scriptstr)
          tname = 'importeovsa'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else :
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')
          result = importeovsa(idbfiles, timebin, width, visprefix, nocreatms, doconcat, modelms)
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

	except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
	     tname = 'importeovsa'
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass

        gc.collect()
        return result
#
#
#
    def paramgui(self, useGlobals=True, ipython_globals=None):
        """
        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
        """
        import paramgui
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=sys._getframe(len(inspect.stack())-1).f_globals

        if useGlobals:
	    if ipython_globals == None:
                myf=self.__globals__
            else:
                myf=ipython_globals

            paramgui.setGlobals(myf)
        else:
            paramgui.setGlobals({})

        paramgui.runTask('importeovsa', myf['_ip'])
        paramgui.setGlobals({})

#
#
#
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=sys._getframe(len(inspect.stack())-1).f_globals
        if ipython_globals == None:
            myf=self.__globals__
        else:
            myf=ipython_globals

        a = odict()
        a['idbfiles']  = ''
        a['timebin']  = '0s'
        a['width']  = 1
        a['visprefix']  = ''
        a['nocreatms']  = True
        a['doconcat']  = False
        a['modelms']  = ''


### This function sets the default values but also will return the list of
### parameters or the default value of a given parameter
        if(param == None):
                myf['__set_default_parameters'](a)
        elif(param == 'paramkeys'):
                return a.keys()
        else:
            if(paramvalue==None and subparam==None):
               if(a.has_key(param)):
                  return a[param]
               else:
                  return self.itsdefault(param)
            else:
               retval=a[param]
               if(type(a[param])==dict):
                  for k in range(len(a[param])):
                     valornotval='value'
                     if(a[param][k].has_key('notvalue')):
                        valornotval='notvalue'
                     if((a[param][k][valornotval])==paramvalue):
                        retval=a[param][k].copy()
                        retval.pop(valornotval)
                        if(subparam != None):
                           if(retval.has_key(subparam)):
                              retval=retval[subparam]
                           else:
                              retval=self.itsdefault(subparam)
		     else:
                        retval=self.itsdefault(subparam)
               return retval


#
#
    def check_params(self, param=None, value=None, ipython_globals=None):
      if ipython_globals == None:
          myf=self.__globals__
      else:
          myf=ipython_globals
#      print 'param:', param, 'value:', value
      try :
         if str(type(value)) != "<type 'instance'>" :
            value0 = value
            value = myf['cu'].expandparam(param, value)
            matchtype = False
            if(type(value) == numpy.ndarray):
               if(type(value) == type(value0)):
                  myf[param] = value.tolist()
               else:
                  #print 'value:', value, 'value0:', value0
                  #print 'type(value):', type(value), 'type(value0):', type(value0)
                  myf[param] = value0
                  if type(value0) != list :
                     matchtype = True
            else :
               myf[param] = value
            value = myf['cu'].verifyparam({param:value})
            if matchtype:
               value = False
      except Exception, instance:
         #ignore the exception and just return it unchecked
         myf[param] = value
      return value
#
#
    def description(self, key='importeovsa', subkey=None):
        desc={'importeovsa': 'Import EOVSA idb file(s) to a measurement set or multiple measurement set',
               'idbfiles': 'Name of input EOVSA idb file(s) or observation time range.',
               'timebin': 'Bin width for time averaging',
               'width': 'Width of output channel relative to MS channel (# to average)',
               'visprefix': 'Prefix of vis names (may include the path).',
               'nocreatms': 'If setting nocreatms True, will simulate a model measurement set for the first idb file and copy the model for the rest of idl files in list. If False, will simulate a new measurement set for every idbfile in list.',
               'doconcat': 'If concatenate multi casa measurement sets to one file.',
               'modelms': 'Name of input model measurement set file. If modelms is assigned, no simulation will start.',

              }

#
# Set subfields defaults if needed
#

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['idbfiles']  = ''
        a['timebin']  = '0s'
        a['width']  = 1
        a['visprefix']  = ''
        a['nocreatms']  = True
        a['doconcat']  = False
        a['modelms']  = ''

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if a.has_key(paramname) :
	      return a[paramname]
importeovsa_cli = importeovsa_cli_()
