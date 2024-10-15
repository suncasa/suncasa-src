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
from suncasa.utils import mstools as mst



def import_calib_idb(trange, workdir=None, ncpu=1, timebin='0s', width=1, udb_corr=True):
    """
    Script to import and calibrate IDB data based on an input time range or list of IDB files.

    :param trange: List containing the begin and end time in eovsa.util Time format or a list of IDB files.
                   Example for time range: Time(['2022-11-12 17:55:00', '2022-11-12 18:10:00'])
                   Example for IDB files: ['IDB202211121804', 'IDB202211121814']
    :type trange: list
    :param workdir: Specify where the working directory is. Defaults to the current path.
    :type workdir: str, optional
    :param ncpu: Number of CPUs to use. Defaults to 1.
    :type ncpu: int, optional
    :param timebin: Time binning interval. Defaults to '0s'.
    :type timebin: str, optional
    :param width: Width of the frequency channels. Defaults to 1.
    :type width: int, optional
    :param udb_corr: Apply correction to input UDB files before import to MS. Defaults to True.
    :type udb_corr: bool, optional
    :returns: Concatenated CASA measurement set with initial gain, amplitude, and phase calibrations applied.
    :rtype: str
    """
    if not workdir:
        workdir = os.getcwd()
        if workdir[-1] != '/':
            workdir += '/'
    if os.path.exists(workdir) == False:
        os.makedirs(workdir)


    os.chdir(workdir)

    try:
        trange = Time(trange)

        # namesuffix = '_' + trange[0].to_datetime().strftime('%H%M') + '-' + trange[1].to_datetime().strftime('%H%M')

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
        print('The timerange corresponds to these files (will take about', len(filelist) * 4, 'minutes to process)')
        for file in filelist:
            print(file)

        ans = input('Do you want to continue? (say no if you want to adjust timerange) [y/n]?').strip().lower()
        if ans == 'n':
            print('Please adjust the timerange and try again.')
            return None

        inpath = idbdir + '{}/'.format(trange[0].datetime.strftime("%Y%m%d"))
        idbfiles = [inpath + ll for ll in filelist]
    except:
        idbfiles = trange

    outpath = 'msdata/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    ncpu = 1

    msfiles = importeovsa(idbfiles=idbfiles, ncpu=ncpu, timebin=timebin, width=width,
                                       visprefix=outpath,
                                       nocreatms=False, doconcat=False,
                                       modelms="", doscaling=False, keep_nsclms=False, udb_corr=udb_corr,
                                       use_exist_udbcorr=True)
    # msfiles = [outpath + ll + '.ms' for ll in filelist]
    if len(msfiles)==0:
        print(f"No measurement set files are imported. Check the input IDB files: {idbfiles}.")
        raise ValueError
    if len(msfiles)==1:
        ms_trange = mst.get_trange(msfiles[0])
    else:
        ms_trange = Time([mst.get_trange(msfiles[0])[0], mst.get_trange(msfiles[-1])[1]])
    namesuffix = '_' + ms_trange[0].to_datetime().strftime('%H%M') + '-' + ms_trange[1].to_datetime().strftime('%H%M')
    concatvis = os.path.basename(msfiles[0])[:11] + namesuffix + '.ms'
    vis = calibeovsa(msfiles, caltype=['refpha', 'phacal'], interp='nearest', doflag=True, flagant='13~15',
                                doimage=False, doconcat=True,
                                concatvis=concatvis, keep_orig_ms=False)
    outputvis = concatvis[:-3] + 'XXYY.ms'
    split(vis=concatvis, outputvis=outputvis, correlation='XX,YY', datacolumn='data')
                   

    return outputvis


