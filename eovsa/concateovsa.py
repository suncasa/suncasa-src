import os
import numpy as np
from taskinit import tb, casalog
from concat_cli import concat_cli as concat
from clearcal_cli import clearcal_cli as clearcal


def concateovsa(msname, msfiles, visprefix='./', keep_orig_ms=False, cols2rm = ["MODEL_DATA", "CORRECTED_DATA"]):
    concatvis = visprefix + msname
    for ll in msfiles:
        clearcal(ll, addmodel=True)
    concat(vis=msfiles, concatvis=concatvis, timesort=True)
    # Change all observation ids to be the same (zero)
    tb.open(concatvis + '/OBSERVATION', nomodify=False)
    nobs = tb.nrows()
    tim0 = tb.getcell('TIME_RANGE', 0)[0]
    tim1 = tb.getcell('TIME_RANGE', nobs - 1)[1]
    tb.removerows([i + 1 for i in range(nobs - 1)])
    tb.putcell('TIME_RANGE', 0, [tim0, tim1])
    tb.close()
    tb.open(concatvis, nomodify=False)
    obsid = tb.getcol('OBSERVATION_ID')
    newobsid = np.zeros(len(obsid), dtype='int')
    tb.putcol('OBSERVATION_ID', newobsid)
    colnames = tb.colnames()
    for l in range(len(cols2rm)):
        if cols2rm[l] in colnames:
            try:
                tb.removecols(cols2rm[l])
            except:
                pass
    tb.close()
    if not keep_orig_ms:
        for ll in msfiles:
            os.system('rm -rf {}'.format(ll))
