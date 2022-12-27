from suncasa.tasks import task_importeovsa as timporteovsa
from suncasa.tasks import task_calibeovsa as tcalibeovsa


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
    vis_in = timporteovsa.importeovsa(trange, visprefix=workdir, ncpu=1, doscaling=False, doconcat=True)
    vis_out = tcalibeovsa.calibeovsa(vis_in, caltype=['refpha', 'phacal'], interp='nearest', doflag=True,
                                     flagant='13~15', doimage=False, doconcat=False, keep_orig_ms=True)
    return vis_out
