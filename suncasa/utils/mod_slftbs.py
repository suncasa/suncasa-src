import os
import numpy as np
import sys

py3 = sys.version_info.major >= 3
try:
    # from taskinit import gentools, qa, casalog
    from taskinit import tbtool
    # from taskinit import me
    # me = gentools(['me'])[0]
except:
    from casatools import table as tbtool
    # from casatools import ms as mstool
    # from casatools import quanta as qatool
    # from casatools import measures as metool
    # from casatasks import casalog

tb = tbtool()

def cpxx2yy(tb_in=[]):
    if not tb_in:
        print('tb_in not provided. Abort...')
    if type(tb_in) is str:
        tb_in = [tb_in]
    tb.open(tb_in[0] + '/SPECTRAL_WINDOW', nomodify=False)
    nspw = tb.nrows()
    tb.close()
    for ctb in tb_in:
        tb.open(ctb, nomodify=False)
        for s in range(nspw):
            subt = tb.query("DATA_DESC_ID==" + str(s))
            model_d = subt.getcol('MODEL_DATA')
            # cp xx to yy
            model_d[1] = model_d[0]
            subt.putcol('MODEL_DATA', model_d)
            subt.close()
        tb.close()


def concat(tb_in=[], tb_out=None):
    if not tb_in:
        print('tb_in not provided. Abort...')
    if os.path.exists(tb_out):
        os.system('rm -r ' + tb_out)
    # os.system('cp -r '+tb_in[0]+' '+tb_out)
    os.system('cp -r ' + tb_in[0] + ' ' + tb_out)
    tb.open(tb_out + '/SPECTRAL_WINDOW', nomodify=True)
    nspw = tb.nrows()
    tb.close()
    tim = []
    fld = []
    spw = []
    ant1 = []
    ant2 = []
    intv = []
    scan = []
    obid = []
    cpar = []
    para = []
    flag = []
    snr = []
    # wght=[]
    for ctb in tb_in:
        tb.open(ctb, nomodify=True)
        cols = tb.colnames()
        tim0 = tb.getcol(cols[0])
        if len(tim0) == 0:
            continue
        else:
            tim.append(tb.getcol(cols[0]))
            fld.append(tb.getcol(cols[1]))
            spw.append(tb.getcol(cols[2]))
            ant1.append(tb.getcol(cols[3]))
            ant2.append(tb.getcol(cols[4]))
            intv.append(tb.getcol(cols[5]))
            scan.append(tb.getcol(cols[6]))
            obid.append(tb.getcol(cols[7]))
            cpar.append(tb.getcol(cols[8]))
            para.append(tb.getcol(cols[9]))
            flag.append(tb.getcol(cols[10]))
            snr.append(tb.getcol(cols[11]))
            # wght.append(tb.getcol(cols[12]))
        tb.close()

    if len(tim) == 0:
        print('tables have no data. Return')
        return -1
    else:
        tim = np.concatenate(tim)
        fld = np.concatenate(fld)
        spw = np.concatenate(spw)
        ant1 = np.concatenate(ant1)
        ant2 = np.concatenate(ant2)
        intv = np.concatenate(intv)
        scan = np.concatenate(scan)
        obid = np.concatenate(obid)
        cpar = np.concatenate(cpar, axis=2)
        para = np.concatenate(para, axis=2)
        flag = np.concatenate(flag, axis=2)
        snr = np.concatenate(snr, axis=2)
        # wght=np.concatenate(wght)
        tb.open(tb_out, nomodify=False)
        nrows = tb.nrows()
        nrows_new = len(tim)
        tb.addrows(nrows_new - nrows)
        tb.putcol(cols[0], tim)
        tb.putcol(cols[1], fld)
        tb.putcol(cols[2], spw)
        tb.putcol(cols[3], ant1)
        tb.putcol(cols[4], ant2)
        tb.putcol(cols[5], intv)
        tb.putcol(cols[6], scan)
        tb.putcol(cols[7], obid)
        tb.putcol(cols[8], cpar)
        tb.putcol(cols[9], para)
        tb.putcol(cols[10], flag)
        tb.putcol(cols[11], snr)
        tb.close()
        return tb_out
