from suncasa.suncasatasks.private import task_importeovsa
from suncasa.suncasatasks.private import task_calibeovsa


def import_calib_idb(trange, workdir=None):
    """
    Script to import and calibrate IDB data based on an input time range
    Parameters
    ----------
    trange: [begin_time, end_time], in astropy.time.Time format.
        Example: Time(['2022-11-12T17:55', '2022-11-12T18:10'])
    workdir: specify where the working directory is. Default to /data1/bchen/flare_pipeline/
    Returns
    -------
    vis_out: concatenated CASA measurement set with initial gain, amplitude, and phase calibrations applied
    """
    if not workdir:
        workdir = '/data1/bchen/flare_pipeline/'
    vis_in = task_importeovsa.importeovsa(trange, visprefix=workdir, ncpu=1, doscaling=False, doconcat=True)
    vis_out = task_calibeovsa.calibeovsa(vis_in, caltype=['refpha', 'phacal'], interp='nearest', doflag=True,
                         flagant='13~15', doimage=False, doconcat=False, 
                         dosplit=True, keep_orig_ms=True)
    return vis_out
