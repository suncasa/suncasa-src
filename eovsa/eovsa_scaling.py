import glob
import numpy as np
from datetime import datetime
import os
import pytz
from util import Time
import refcal_anal as ra
import multiprocessing as mp
from taskinit import casalog
from importeovsa_cli import importeovsa_cli as importeovsa


def mk_udbms(trange=None, outpath=None, projid='NormalObserving', srcid='Sun', doscaling=True):
    '''
        usage: outfiles = mk_udbms(Time('2017-07-02 15:00'))
    :param trange:
    :param outpath:
    :param projid:
    :param srcid:
    :param doscaling:
    :return:
    '''
    if trange is None:
        trange = Time.now()

    try:
        if len(trange) >= 2:
            trange = Time([trange[0], trange[-1]])
            tdatetime = trange[0].to_datetime()
    except:
        tdatetime = trange.to_datetime()
        local_tz = pytz.timezone('America/Los_Angeles')
        btime = Time(local_tz.localize(tdatetime, is_dst=None).astimezone(pytz.utc))
        etime = Time(btime.mjd + 1.0, format='mjd')
        trange = Time([btime, etime])
    if outpath is None:
        outpath = '/data1/eovsa/fits/UDBms_scl/{}/'.format(tdatetime.strftime("%Y%m"))
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    sclist = ra.findfiles(trange, projid=projid, srcid=srcid)
    ncpu = mp.cpu_count()
    if ncpu > 10:
        ncpu = 10
    if ncpu > len(sclist['scanlist']):
        ncpu = len(sclist['scanlist'])
    if sclist['scanlist']:
        # importeovsa(idbfiles=sclist['scanlist'], ncpu=ncpu, timebin="0s", width=1,
        #             visprefix=outpath, nocreatms=False, doconcat=False, modelms="", doscaling=doscaling)
        importeovsa(idbfiles=sclist['scanlist'], ncpu=ncpu, timebin="0s", width=1,
                    visprefix=outpath, nocreatms=False, doconcat=False, modelms="", doscaling=doscaling,
                    keep_nsclms=False)
        if doscaling:
            # for ll in sclist['scanlist']:
            #     os.system('rm -rf {}.ms'.format(outpath + os.path.basename(ll)))
            outfiles = ['{}_scl.ms'.format(outpath + os.path.basename(ll)) for ll in sclist['scanlist']]
        else:
            outfiles = ['{}.ms'.format(outpath + os.path.basename(ll)) for ll in sclist['scanlist']]
        casalog.post("{} UDB files converted to CASA measurement sets.".format(len(sclist['scanlist'])))
        print "{} UDB files converted to CASA measurement sets.".format(len(sclist['scanlist']))
        return outfiles
    else:
        casalog.post("No UDB files found. Quit.")
        print 'No UDB files found. Quit.'
        return None

        # msfiles = glob.glob('UDB*_scl.ms')
        # ms2concat = []
        # for ll in msfiles:
        #     ms.open(ll)
        #     if ms.summary()['field_0']['name'] == 'Sun':
        #         ms2concat.append(ll)
        #     ms.close()
        #
        # for ll in ms2concat:
        #     clearcal(ll, addmodel=True)
        #
        # concat(vis=ms2concat, concatvis='UDB20170624_SUN.ms', timesort=True)
        #
        # tb.open(concatvis + '/OBSERVATION', nomodify=False)
        # nobs = tb.nrows()
        # tb.removerows([i + 1 for i in range(nobs - 1)])
        # tb.close()
        # tb.open(concatvis, nomodify=False)
        # obsid = tb.getcol('OBSERVATION_ID')
        # newobsid = np.zeros(len(obsid), dtype='int')
        # tb.putcol('OBSERVATION_ID', newobsid)
        # tb.close()
