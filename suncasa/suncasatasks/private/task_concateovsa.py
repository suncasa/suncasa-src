import os
import numpy as np
from taskinit import tb, casalog
from concat_cli import concat_cli as concat
from clearcal_cli import clearcal_cli as clearcal
from split_cli import split_cli as split


def concateovsa(vis, concatvis, datacolumn='corrected', keep_orig_ms=True, cols2rm="model,corrected", freqtol="", dirtol="", respectname=False,
                timesort=True, copypointing=True, visweightscale=[], forcesingleephemfield=""):
    if concatvis[-1] == os.path.sep:
        concatvis = concatvis[:-1]
    if os.path.sep not in concatvis:
        visprefix = './'
    else:
        visprefix = os.path.dirname(concatvis) + os.path.sep
    msfiles = vis
    msfiles_ = []
    for idx, ll in enumerate(msfiles):
        if str(ll).endswith('/'):
            msfiles[idx] = str(ll)[:-1]
    datacolumn = datacolumn.lower()
    if datacolumn == 'data':
        print('DATA columns will be concatenated.')
        for ll in msfiles:
            clearcal(vis=str(ll), addmodel=True)
    elif datacolumn == 'corrected':
        # try:
        print('CORRECTED columns will be concatenated.')
        tmpdir = os.path.join(visprefix, 'tmp_ms') + os.path.sep
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        for ll in msfiles:
            msfile_ = os.path.join(tmpdir, os.path.basename(str(ll)))
            msfiles_.append(msfile_)
            split(vis=str(ll), outputvis=msfile_, datacolumn='corrected')
            clearcal(vis=msfile_, addmodel=True)
    else:
        raise ValueError('Please set datacolumn to be "data" or "corrected"!')

    if msfiles_:
        concat(vis=msfiles_, concatvis=concatvis, freqtol=freqtol, dirtol=dirtol, respectname=respectname, timesort=timesort,
               copypointing=copypointing, visweightscale=visweightscale, forcesingleephemfield=forcesingleephemfield)
        os.system('rm -rf {}'.format(tmpdir))
    else:
        concat(vis=msfiles, concatvis=concatvis, freqtol=freqtol, dirtol=dirtol, respectname=respectname, timesort=timesort,
               copypointing=copypointing, visweightscale=visweightscale, forcesingleephemfield=forcesingleephemfield)
    # Change all observation ids to be the same (zero)
    tb.open(concatvis + '/OBSERVATION', nomodify=False)
    nobs = tb.nrows()
    tim0 = tb.getcell('TIME_RANGE', 0)[0]
    tim1 = tb.getcell('TIME_RANGE', nobs - 1)[1]
    tb.removerows([i + 1 for i in range(nobs - 1)])
    tb.putcell('TIME_RANGE', 0, [tim0, tim1])
    tb.close()

    tb.open(concatvis + '/DATA_DESCRIPTION', nomodify=False)
    nrows = tb.nrows()
    pol_id = tb.getcol('POLARIZATION_ID')
    tb.removerows(np.where(pol_id != 0)[0])
    tb.close()

    tb.open(concatvis, nomodify=False)
    dd_id = tb.getcol('DATA_DESC_ID')
    idx_dd_id, = np.where(dd_id >= nrows / 2)
    dd_id[idx_dd_id] = dd_id[idx_dd_id] - nrows / 2
    tb.putcol('DATA_DESC_ID', dd_id)
    tb.close()

    tb.open(concatvis + '/FIELD', nomodify=False)
    nobs = tb.nrows()
    tb.removerows([i + 1 for i in range(nobs - 1)])
    tb.close()

    tb.open(concatvis + '/SOURCE', nomodify=False)
    nobs = tb.nrows()
    tb.removerows([i + 1 for i in range(nobs - 1)])
    tb.close()

    tb.open(concatvis, nomodify=False)
    obsid = tb.getcol('OBSERVATION_ID')
    newobsid = np.zeros(len(obsid), dtype='int')
    tb.putcol('OBSERVATION_ID', newobsid)
    fldid = tb.getcol('FIELD_ID')
    newfldid = np.zeros(len(fldid), dtype='int')
    tb.putcol('FIELD_ID', newfldid)
    colnames = tb.colnames()

    cols2rm = cols2rm.upper()
    cols2rm = cols2rm.split(',')
    for l in range(len(cols2rm)):
        col = cols2rm[l] + '_DATA'
        if col in colnames:
            try:
                tb.removecols(col)
                print('Column {} removed.'.format(col))
            except:
                pass
    tb.close()

    if msfiles_ != [] and msfiles_ != msfiles:
        for ll in msfiles_:
            os.system('rm -rf {}'.format(ll))
    if not keep_orig_ms:
        for ll in msfiles:
            os.system('rm -rf {}'.format(ll))
            os.system('rm -rf {}.flagversions'.format(ll))
