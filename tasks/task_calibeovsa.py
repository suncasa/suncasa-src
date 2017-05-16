import os
import gc
import numpy as np
import scipy.constants as constants
import aipy
import eovsapy.chan_util_bc as chan_util_bc
import eovsapy.read_idb as ri
import eovsapy.refcal_anal as ra
from eovsapy.util import Time
from taskinit import *
from gencal_cli import gencal_cli as gencal
from applycal_cli import applycal_cli as applycal

def callibeovsa(vis, caltype=None, docalib=False):
    casalog.origin('eovsacalib')
    if not caltype:
        casalog.post("Caltype not provided. Try to generate all that is applicable.")
        caltype=['refcal']
    if not os.path.exists(vis):
        casalog.post("Input visibility does not exist. Aborting...")
    #if not caltable:
    #    caltable=[os.path.basename(vis).replace('.ms','.'+c) for c in caltype]

    # get band information
    tb.open(vis+'/SPECTRAL_WINDOW')
    nspw=tb.nrows()
    bdname=tb.getcol('NAME')
    bd=[int(b[4:])-1 for b in bdname] #band index from 0 to 33
    nchans=tb.getcol('NUM_CHAN')
    reffreqs=tb.getcol('REF_FREQUENCY')
    cenfreqs=np.zeros((nspw)) 
    tb.close()
    
    # get time stamp, use the beginning of the file
    ms.open(vis)
    summary = ms.summary() 
    ms.close()
    btime = Time(summary['BeginTime'],format='mjd')
    gaintables=[]
    if ('refpha' in caltype) or ('refamp' in caltype) or ('refcal' in caltype):
        refcal = ra.sql2refcal(btime)  
        pha = refcal['pha'] #shape is 15 (nant) x 2 (npol) x 34 (nband)
        amp = refcal['amp']
        ref_t = refcal['timestamp']
        #casalog.post("Reference calibration is derived from observation at "+ref_t.iso)
        print "Reference calibration is derived from observation at "+ref_t.iso

        para_pha=[]
        para_amp=[]
        calpha=np.zeros((nspw,13,2))
        calamp=np.zeros((nspw,13,2))
        for s in range(nspw):
            for n in range(13):
                for p in range(2):
                    calpha[s,n,p]=pha[n,p,bd[s]]
                    calamp[s,n,p]=amp[n,p,bd[s]]
                    para_pha.append(np.degrees(pha[n,p,bd[s]]))
                    para_amp.append(amp[n,p,bd[s]])
        
        if ('refpha' in caltype) or ('refcal' in caltype):
            caltb_pha=os.path.basename(vis).replace('.ms','.refpha')
            gaintables.append(caltb_pha)
            gencal(vis=vis,caltable=caltb_pha,caltype='ph',antenna='0~12',\
               pol='X,Y',spw='0~'+str(nspw-1),parameter=para_pha)
        if ('refamp' in caltype) or ('refcal' in caltype):
            caltb_amp=os.path.basename(vis).replace('.ms','.refamp')
            gaintables.append(caltb_amp)
            gencal(vis=vis,caltable=caltb_amp,caltype='amp',antenna='0~12',\
               pol='X,Y',spw='0~'+str(nspw-1),parameter=para_amp)

    if docalib:
        applycal(vis=vis,gaintable=gaintables)

        
        





