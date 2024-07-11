from suncasa.suncasatasks import importeovsa
from suncasa.suncasatasks import calibeovsa
import sys
# import importeovsa
# import calibeovsa
from casatasks import split
from eovsapy.read_idb import get_trange_files
from eovsapy.util import Time
from eovsapy import util
import numpy as np
import os
from eovsapy import dump_tsys as dt



def import_calib_idb(trange, workdir=None, ncpu=1, timebin='0s', width=1):

    """
    Script to import and calibrate IDB data based on an input time range
    Parameters
    ----------
    trange: [begin_time, end_time], in eovsa.util Time format.
        Example: Time(['2022-11-12 17:55:00', '2022-11-12 18:10:00'])
    workdir: specify where the working directory is. Default to current path
    Returns
    -------
    vis_out: concatenated CASA measurement set with initial gain, amplitude, and phase calibrations applied
    """
    if not workdir:
        workdir = os.getcwd()
        if workdir[-1] != '/':
            workdir += '/'
    if os.path.exists(workdir) == False:
        os.makedirs(workdir)


    os.chdir(workdir)

    namesuffix = '_' + trange[0].to_datetime().strftime('%H%M') + '-' + trange[1].to_datetime().strftime('%H%M')

    idbdir = util.get_idbdir(trange[0])

    info = dt.rd_fdb(Time(trange[0].to_datetime().strftime('%Y-%m-%d')))
    idxs, = np.where(info['SOURCEID'] == '')
    for k, v in info.items():
        info[k] = info[k][~(info[k] == '')]
    # sidx = np.where(
    #     np.logical_and(info['SOURCEID'] == 'Sun', info['PROJECTID'] == 'NormalObserving') & np.logical_and(
    #         info['ST_TS'].astype(np.float_) >= trange[0].lv,
    #         info['ST_TS'].astype(np.float_) <= trange[1].lv))
    
    term1 = np.logical_and(info['SOURCEID'] == 'Sun', info['PROJECTID'] == 'NormalObserving')
    st_ts_float = info['ST_TS'].astype(np.float_)

    term2_ind = []
    term3_ind = []
    for i in range(len(st_ts_float) - 1):
        if st_ts_float[i] <= trange[0].lv <= st_ts_float[i + 1]:
            term2_ind.append(i)
        if st_ts_float[i] <= trange[1].lv <= st_ts_float[i + 1]:
            term3_ind.append(i)
    if term2_ind == [] or term3_ind == []:
        print("Error: No available IDB list during the indicated time range.")

    sidx = list(range(term2_ind[0], term3_ind[0] + 1))
    # Filter out indices where term1 is not true
    sidx = [i for i in sidx if term1[i]]
    if sidx == []:
        print("Error: No available IDB list for the Sun during the indicated time range.")


    filelist = info['FILE'][sidx]
    # filelist = filelist[np.array([0, 2, 3, 4, 5, 6, 7, 8])]
    # filelist = filelist[np.array([0, 2])]

    outpath = 'msdata/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    inpath = idbdir + '{}/'.format(trange[0].datetime.strftime("%Y%m%d"))
    ncpu = 1

    msfiles = importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin=timebin, width=width,
                                       visprefix=outpath,
                                       nocreatms=False, doconcat=False,
                                       modelms="", doscaling=False, keep_nsclms=False, udb_corr=True,
                                       use_exist_udbcorr=True)
    # msfiles = [outpath + ll + '.ms' for ll in filelist]
    concatvis = os.path.basename(msfiles[0])[:11] + namesuffix + '.ms'
    vis = calibeovsa(msfiles, caltype=['refpha', 'phacal'], interp='nearest', doflag=True, flagant='13~15',
                                doimage=False, doconcat=True,
                                concatvis=concatvis, keep_orig_ms=False)
    
                   

    return concatvis


