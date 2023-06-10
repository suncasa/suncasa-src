from suncasa.suncasatasks import importeovsa
from suncasa.suncasatasks import calibeovsa
import sys
# import importeovsa
# import calibeovsa
from casatasks import split
from eovsapy.read_idb import get_trange_files
from eovsapy.util import Time
import numpy as np
import os


def import_calib_idb(trange, workdir=None, ncpu=1, timebin='0s', width=1):
    """
    Script to import and calibrate IDB data based on an input time range
    Parameters
    ----------
    trange: [begin_time, end_time], in astropy.time.Time format.
        Example: Time(['2022-11-12T17:55', '2022-11-12T18:10'])
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

    files = get_trange_files(trange)
    msfile = importeovsa(idbfiles=files, ncpu=ncpu, timebin=timebin, width=width, \
                         visprefix=workdir, nocreatms=False, doconcat=True, \
                         modelms='', doscaling=False, keep_nsclms=False, \
                         udb_corr=True)
    vis_out = calibeovsa(msfile, caltype=['refpha', 'phacal'], interp='linear', \
                         doimage=False, doconcat=True, dosplit=True, keep_orig_ms=False)
    return vis_out
