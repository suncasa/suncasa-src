'''
   Module for converting CASA image analysis outputs to data recognizable by DataBrowser   
'''
# History:
#   2016-Nov-2 BC
#       New module

import numpy as np

def pimfitres2dict(pimfitres):
    ''' parse pimfit results to a dictionary '''
    # check whether the input is a single dictionary (from imfit) or a list (from pimfit)
    if isinstance(imfitres,dict):
        imfitres=[imfitres]
    if (not isinstance(imfitres,list)):
        print 'input "imfitres" is not a list. Abort...'
    out = []
    ntim = len(imfitres)
    for i, res in enumerate(imfitres):
        if 'converged' in res:
            nfreq=len(res['converged'])


         

