import os
import gc
import numpy as np
import scipy.constants as constants
import aipy
import eovsapy.chan_util_bc as chan_util_bc
import eovsapy.read_idb as ri
import eovsapy.refcal_anal as ra
from eovsapy.util import Time
from taskinit import casalog, tb, ms
from gencal_cli import gencal_cli as gencal
from applycal_cli import applycal_cli as applycal
from eovsa import cal_header as ch
from eovsapy import stateframe as stf


def callibeovsa(vis, caltype=None, docalib=False):
    casalog.origin('eovsacalib')
    if not caltype:
        casalog.post("Caltype not provided. Try to generate all that is applicable.")
        caltype = ['refcal']
    if not os.path.exists(vis):
        casalog.post("Input visibility does not exist. Aborting...")
    if vis.endswith('/'):
        vis= vis[:-1]
    if not vis[-3:] in ['.ms','.MS']:
        casalog.post("Invalid visibility. Aborting...")
    # if not caltable:
    #    caltable=[os.path.basename(vis).replace('.ms','.'+c) for c in caltype]

    # get band information
    tb.open(vis + '/SPECTRAL_WINDOW')
    nspw = tb.nrows()
    bdname = tb.getcol('NAME')
    bd = [int(b[4:]) - 1 for b in bdname]  # band index from 0 to 33
    # nchans = tb.getcol('NUM_CHAN')
    # reffreqs = tb.getcol('REF_FREQUENCY')
    # cenfreqs = np.zeros((nspw))
    tb.close()
    tb.open(vis + '/ANTENNA')
    nant = tb.nrows()
    antname = tb.getcol('NAME')
    antlist = [str(ll) for ll in range(len(antname))]

    # get time stamp, use the beginning of the file
    ms.open(vis)
    summary = ms.summary()
    ms.close()
    btime = Time(summary['BeginTime'], format='mjd')
    gaintables = []
    if ('refpha' in caltype) or ('refamp' in caltype) or ('refcal' in caltype):
        refcal = ra.sql2refcal(btime)
        pha = refcal['pha']  # shape is 15 (nant) x 2 (npol) x 34 (nband)
        amp = refcal['amp']
        ref_t = refcal['timestamp']
        # casalog.post("Reference calibration is derived from observation at "+ref_t.iso)
        print "Reference calibration is derived from observation at " + ref_t.iso

        para_pha = []
        para_amp = []
        calpha = np.zeros((nspw, 13, 2))
        calamp = np.zeros((nspw, 13, 2))
        for s in range(nspw):
            for n in range(13):
                for p in range(2):
                    calpha[s, n, p] = pha[n, p, bd[s]]
                    calamp[s, n, p] = amp[n, p, bd[s]]
                    para_pha.append(np.degrees(pha[n, p, bd[s]]))
                    para_amp.append(amp[n, p, bd[s]])

        if ('refpha' in caltype) or ('refcal' in caltype):
            caltb_pha = os.path.basename(vis).replace('.ms', '.refpha')
            gaintables.append(caltb_pha)
            gencal(vis=vis, caltable=caltb_pha, caltype='ph', antenna='0~12', \
                   pol='X,Y', spw='0~' + str(nspw - 1), parameter=para_pha)
        if ('refamp' in caltype) or ('refcal' in caltype):
            caltb_amp = os.path.basename(vis).replace('.ms', '.refamp')
            gaintables.append(caltb_amp)
            gencal(vis=vis, caltable=caltb_amp, caltype='amp', antenna='0~12', \
                   pol='X,Y', spw='0~' + str(nspw - 1), parameter=para_amp)

        xml, buf = ch.read_calX(4, t=[ref_t, btime], verbose=False)
        if buf:
            dla_t2 = Time(stf.extract(buf[0], xml['Timestamp']), format='lv')
            dlacen_ns2 = stf.extract(buf[0], xml['Delaycen_ns'])
            xml, buf = ch.read_calX(4, t=ref_t)
            dla_t1 = Time(stf.extract(buf, xml['Timestamp']), format='lv')
            dlacen_ns1 = stf.extract(buf, xml['Delaycen_ns'])
            dlacen_ns_diff = dlacen_ns2 - dlacen_ns1
            dlacen_ns_diff[:, 0] -= dlacen_ns_diff[0, 0]
            dlacen_ns_diff[:, 1] -= dlacen_ns_diff[0, 1]
            antennas = ','.join(antlist)
            print 'Multi-band delay is derived from delay center at {} & {}'.format(dla_t1.iso, dla_t2.iso)
            # print '=====Delays relative to Ant 14====='
            # for i, dl in enumerate(dlacen_ns_diff[:, 0] - dlacen_ns_diff[13, 0]):
            #     ant = antlist[i]
            #     print 'Ant eo{0:02d}: x {1:.2f} ns & y {2:.2f} ns'.format(int(ant) + 1, dl,
            #                                                               dlacen_ns_diff[i, 1] - dlacen_ns_diff[13, 1])
            caltb_mbd = os.path.basename(vis).replace('.ms', '.mdb')
            gaintables.append(caltb_mbd)
            gencal(vis=vis, caltable=caltb_mbd, caltype='mbd', pol='X,Y', antenna=antennas,
                   parameter=dlacen_ns_diff.flatten().tolist())
    if docalib:
        applycal(vis=vis, gaintable=gaintables)
