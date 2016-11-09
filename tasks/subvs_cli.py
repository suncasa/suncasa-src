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
from task_subvs import subvs
class subvs_cli_:
    __name__ = "subvs"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (subvs_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'vis':None, 'outputvis':None, 'timerange':None, 'spw':None, 'mode':None, 'subtime1':None, 'subtime2':None, 'smoothaxis':None, 'smoothtype':None, 'smoothwidth':None, 'splitsel':None, 'reverse':None, 'overwrite':None, }


    def result(self, key=None):
	    #### and add any that have completed...
	    return None


    def __call__(self, vis=None, outputvis=None, timerange=None, spw=None, mode=None, subtime1=None, subtime2=None, smoothaxis=None, smoothtype=None, smoothwidth=None, splitsel=None, reverse=None, overwrite=None, ):

        """Vector-subtraction in UV using selected time ranges and spectral channels as background

	Detailed Description: 



Split is the general purpose program to make a new data set that is a
subset or averaged form of an existing data set.  General selection
parameters are included, and one or all of the various data columns
(DATA, LAG_DATA and/or FLOAT_DATA, and possibly MODEL_DATA and/or
CORRECTED_DATA) can be selected.

Split is often used after the initial calibration of the data to make a
smaller measurement set with only the data that will be used in
further flagging, imaging and/or self-calibration.  split can
average over frequency (channels) and time (integrations).

	Arguments :
		vis:	Name of input measurement set
		   Default Value: 

		outputvis:	Name of output measurement set
		   Default Value: 

		timerange:	Select the time range of the input visbility to be subtracted from
		   Default Value: 

		spw:	Select the spectral channels of the input visibility to be subtracted from
		   Default Value: 

		mode:	Operation: linear, highpass
		   Default Value: linear
		   Allowed Values:
				linear
				lowpass
				highpass

		subtime1:	Select the first time range as the background for uv subtraction 
		   Default Value: 

		subtime2:	Select the second time range as the background for uv subtraction 
		   Default Value: 

		smoothaxis:	Select the axis along which smooth is performed
		   Default Value: time

		smoothtype:	Select the smooth type
		   Default Value: flat

		smoothwidth:	Select the width of the smoothing window
		   Default Value: 5

		splitsel:	Split the selected timerange and spectral channels as outputvis
		   Default Value: True

		reverse:	Reverse the sign of the background-subtracted data (for absorptive structure)
		   Default Value: False

		overwrite:	Overwrite the already existing output measurement set
		   Default Value: False


	Example :


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
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=sys._getframe(len(inspect.stack())-1).f_globals
	#casac = self.__globals__['casac']
	casalog = self.__globals__['casalog']
	casa = self.__globals__['casa']
	#casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'subvs'
        self.__globals__['taskname'] = 'subvs'
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

            myparams['vis'] = vis = self.parameters['vis']
            myparams['outputvis'] = outputvis = self.parameters['outputvis']
            myparams['timerange'] = timerange = self.parameters['timerange']
            myparams['spw'] = spw = self.parameters['spw']
            myparams['mode'] = mode = self.parameters['mode']
            myparams['subtime1'] = subtime1 = self.parameters['subtime1']
            myparams['subtime2'] = subtime2 = self.parameters['subtime2']
            myparams['smoothaxis'] = smoothaxis = self.parameters['smoothaxis']
            myparams['smoothtype'] = smoothtype = self.parameters['smoothtype']
            myparams['smoothwidth'] = smoothwidth = self.parameters['smoothwidth']
            myparams['splitsel'] = splitsel = self.parameters['splitsel']
            myparams['reverse'] = reverse = self.parameters['reverse']
            myparams['overwrite'] = overwrite = self.parameters['overwrite']


	result = None

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
	pathname="file:///local/software/suncasa/tasks/"
	trec = casac.casac.utils().torecord(pathname+'subvs.xml')

        casalog.origin('subvs')
	try :
          #if not trec.has_key('subvs') or not casac.casac.utils().verify(mytmp, trec['subvs']) :
	    #return False

          casac.casac.utils().verify(mytmp, trec['subvs'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('subvs', 'subvs.last', myparams, self.__globals__,scriptstr=scriptstr)
          tname = 'subvs'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else :
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')
          result = subvs(vis, outputvis, timerange, spw, mode, subtime1, subtime2, smoothaxis, smoothtype, smoothwidth, splitsel, reverse, overwrite)
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

	except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
	     tname = 'subvs'
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

        paramgui.runTask('subvs', myf['_ip'])
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
        a['vis']  = ''
        a['outputvis']  = ''
        a['timerange']  = ''
        a['spw']  = ''
        a['mode']  = 'linear'
        a['splitsel']  = True
        a['reverse']  = False
        a['overwrite']  = False

        a['mode'] = {
                    0:odict([{'value':'linear'}, {'subtime1':''}, {'subtime2':''}]), 
                    1:odict([{'value':'lowpass'}, {'smoothaxis':'time'}, {'smoothtype':'hanning'}, {'smoothwidth':5}]), 
                    2:odict([{'value':'highpass'}, {'smoothaxis':'time'}, {'smoothtype':'hanning'}, {'smoothwidth':5}])}

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
    def description(self, key='subvs', subkey=None):
        desc={'subvs': 'Vector-subtraction in UV using selected time ranges and spectral channels as background',
               'vis': 'Name of input measurement set',
               'outputvis': 'Name of output measurement set',
               'timerange': 'Select the time range of the input visbility to be subtracted from',
               'spw': 'Select the spectral channels of the input visibility to be subtracted from',
               'mode': 'Operation: linear, highpass',
               'subtime1': 'Select the first time range as the background for uv subtraction ',
               'subtime2': 'Select the second time range as the background for uv subtraction ',
               'smoothaxis': 'Select the axis along which smooth is performed',
               'smoothtype': 'Select the smooth type',
               'smoothwidth': 'Select the width of the smoothing window',
               'splitsel': 'Split the selected timerange and spectral channels as outputvis',
               'reverse': 'Reverse the sign of the background-subtracted data (for absorptive structure)',
               'overwrite': 'Overwrite the already existing output measurement set',

              }

#
# Set subfields defaults if needed
#

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['vis']  = ''
        a['outputvis']  = ''
        a['timerange']  = ''
        a['spw']  = ''
        a['mode']  = 'linear'
        a['subtime1']  = ''
        a['subtime2']  = ''
        a['smoothaxis']  = 'time'
        a['smoothtype']  = 'flat'
        a['smoothwidth']  = 5
        a['splitsel']  = True
        a['reverse']  = False
        a['overwrite']  = False

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if self.parameters['mode']  == 'linear':
            a['subtime1'] = ''
            a['subtime2'] = ''

        if self.parameters['mode']  == 'lowpass':
            a['smoothaxis'] = 'time'
            a['smoothtype'] = 'hanning'
            a['smoothwidth'] = 5

        if self.parameters['mode']  == 'highpass':
            a['smoothaxis'] = 'time'
            a['smoothtype'] = 'hanning'
            a['smoothwidth'] = 5

        if a.has_key(paramname) :
	      return a[paramname]
subvs_cli = subvs_cli_()
