import os
import numpy as np
from taskinit import tb, casalog
from concat_cli import concat_cli as concat
from clearcal_cli import clearcal_cli as clearcal
from split_cli import split_cli as split


def concateovsa(msname, msfiles, visprefix='./', doclearcal=True, keep_orig_ms=False,
                cols2rm=["MODEL_DATA", "CORRECTED_DATA"]):
    concatvis = visprefix + msname
    msfiles_ = []
    for idx, ll in enumerate(msfiles):
        if str(ll).endswith('/'):
            msfiles[idx] = str(ll)[:-1]
    if doclearcal:
        print 'Warning: Corrected column in the input ms file will be cleared!!!'
        for ll in msfiles:
            clearcal(vis=str(ll), addmodel=True)
    else:
        try:
            tmpdir = visprefix + '/tmp_ms/'
            if not os.path.exists(tmpdir):
                os.makedirs(tmpdir)
            for ll in msfiles:
                msfile_ = tmpdir + os.path.basename(str(ll))
                msfiles_.append(msfile_)
                split(vis=str(ll), outputvis=msfile_, datacolumn='corrected')
                clearcal(vis=msfile_, addmodel=True)
        except:
            print 'Warning: Corrected column not found in the input ms file.'
            msfiles_ = msfiles
    if msfiles_:
        concat(vis=msfiles_, concatvis=concatvis, timesort=True)
    else:
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
    if msfiles_ != [] and msfiles_ != msfiles:
        for ll in msfiles_:
            os.system('rm -rf {}'.format(ll))
    if not keep_orig_ms:
        for ll in msfiles:
            os.system('rm -rf {}'.format(ll))
