import copy
import os
import platform
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import drms
import hvpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import sunpy
import sunpy.map as smap
from astropy import units as u
from astropy.time import Time
from hvpy.datasource import DataSource
from sunpy.net import Fido, attrs as a

from ..dspec import dspec as ds
# from config import get_and_create_download_dir
from ..utils import helioimage2fits as hf
from ..utils import mstools

data_sources_aia = {
    131: DataSource.AIA_131,
    193: DataSource.AIA_193,
    4500: DataSource.AIA_4500,
    1600: DataSource.AIA_1600,
    211: DataSource.AIA_211,
    94: DataSource.AIA_94,
    1700: DataSource.AIA_1700,
    304: DataSource.AIA_304,
    171: DataSource.AIA_171,
    335: DataSource.AIA_335
}

systemname = platform.system()

sunpy1 = sunpy.version.major >= 1
sunpy3 = sunpy.version.major >= 3
py3 = sys.version_info.major >= 3
if py3:
    # For Python 3.0 and later
    from urllib.request import urlopen
else:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

from ..suncasatasks import ptclean6 as ptclean
from ..casa_compat import import_casatools, import_casatasks

tasks = import_casatasks('split', 'tclean', 'casalog')
split = tasks.get('split')
tclean = tasks.get('tclean')
casalog = tasks.get('casalog')

tools = import_casatools(['tbtool', 'mstool', 'qatool'])

tbtool = tools['tbtool']
mstool = tools['mstool']
qatool = tools['qatool']
tb = tbtool()
ms = mstool()
qa = qatool()

c_external = False
from matplotlib.dates import DateFormatter
from astropy.coordinates import SkyCoord
import matplotlib as mpl

import matplotlib.colors as colors
import matplotlib.patches as patches
from ..utils import DButil
from mpl_toolkits.axes_grid1 import make_axes_locatable
from suncasa.utils import plot_mapX as pmX
from suncasa.io import ndfits
from tqdm import tqdm

sunpy1 = sunpy.version.major >= 1
if sunpy1:
    from sunpy.coordinates import sun
else:
    from sunpy import sun
    import sunpy.cm.cm as cm_sunpy

polmap = {'RR': 0, 'LL': 1, 'I': 0, 'V': 1, 'XX': 0, 'YY': 1}


def validate_and_reset_restoringbeam(restoringbm):
    """
    Validates the format of the restoringbeam string. If the format is incorrect,
    it prints an error message and resets the value to an empty string.

    Parameters:
    - restoringbeam (str): The restoring beam size to validate.

    Returns:
    - str: The original restoringbeam if valid, or an empty string if invalid.
    """
    # Pattern to match: optionally a plus sign, followed by one or more digits,
    # optionally followed by a decimal point with more digits, and ending with 'arcsec'
    pattern = r'^\+?\d+(\.\d+)?arcsec$'

    # Check if the string matches the pattern
    if not re.match(pattern, restoringbm):
        print(
            'Error in the format of the provided restoringbeam, which should be in the format of "100arcsec". Resetting to "".')
        return ''  # Reset to an empty string if the format is incorrect
    else:
        return restoringbm  # Return the original string if the format is correct


def get_normalization(vmin, vmax, scale):
    """
    Returns a normalization object based on the given scaling type.

    :param vmin: Minimum value for normalization
    :type vmin: float
    :param vmax: Maximum value for normalization
    :type vmax: float
    :param scale: Type of scaling, either 'linear' or 'log'
    :type scale: str
    :return: The normalization object based on the given scaling
    :rtype: colors.Normalize or colors.LogNorm
    :raises ValueError: If `scale` is not 'linear' or 'log'
    """
    # print(f"vmin: {vmin}, vmax: {vmax}, scale: {scale}")
    if scale.lower() == 'linear':
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    elif scale.lower() == 'log':
        if vmin <= 0:
            vmin = None
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        raise ValueError('Only "linear" and "log" are accepted for scaling.')

    return norm


def read_imres(imresfile):
    from astropy.time import Time
    py3 = sys.version_info.major >= 3

    def uniq(lst):
        last = object()
        nlst = []
        for item in lst:
            if item == last:
                continue
            nlst.append(item)
            last = item
        return nlst

    imres = np.load(imresfile, encoding="latin1", allow_pickle=True)
    imres = imres['imres'].item()

    if not py3:
        iterop_ = imres.iteritems()
    else:
        iterop_ = imres.items()
    for k, v in iterop_:
        imres[k] = list(np.array(v))
    Spw = sorted(list(set(imres['Spw'])))
    # Spw = [str(int(sp)) for sp in Spw]
    Spw = [str(int(sp)) if '~' not in sp else sp for sp in Spw]
    nspw = len(Spw)
    imres['Freq'] = [list(ll) for ll in imres['Freq']]
    Freq = sorted(uniq(imres['Freq']))

    plttimes = list(set(imres['BeginTime']))
    plttimes = sorted(plttimes)
    ntime = len(plttimes)
    # sort the imres according to time
    images = np.array(imres['ImageName'])
    btimes = Time(imres['BeginTime'])
    etimes = Time(imres['EndTime'])
    spws = np.array(imres['Spw'])
    obs = np.array(imres['Obs'])[0]
    vis = np.array(imres['Vis'])[0]
    inds = btimes.argsort()
    images_sort = images[inds].reshape(ntime, nspw)
    btimes_sort = btimes[inds].reshape(ntime, nspw)
    etimes_sort = etimes[inds].reshape(ntime, nspw)
    spws_sort = spws[inds].reshape(ntime, nspw)
    return {'images': images_sort,
            'btimes': btimes_sort,
            'etimes': etimes_sort,
            'spws_sort': spws_sort,
            'spws': Spw,
            'freq': Freq,
            'plttimes': Time(plttimes),
            'obs': obs,
            'vis': vis}


def checkspecnan(spec):
    import matplotlib
    from pkg_resources import parse_version
    if parse_version(matplotlib.__version__) < parse_version('1.5.0'):
        spec[np.isnan(spec)] = 0.0
    return spec


def get_goes_data(t=None, sat_num=None):
    ''' Reads GOES data from https://umbra.nascom.nasa.gov/ repository, for date
        and satellite number provided.  If sat_num is None, data for all available
        satellites are downloaded, with some sanity check used to decide the best.
        If the Time() object t is None, data for the day before the current date
        are read (since there is a delay of 1 day in availability of the data).
        Returns:
           goes_t    GOES time array in plot_date format
           goes_data GOES 1-8 A lightcurve
        '''
    from sunpy.util.config import get_and_create_download_dir
    import shutil
    from astropy.io import fits
    import ssl
    if t is None:
        t = Time(Time.now().mjd - 1, format='mjd')
    yr = t.iso[:4]
    datstr = t.iso[:10].replace('-', '')
    context = ssl._create_unverified_context()
    if sat_num is None:
        f = urlopen.urlopen('https://umbra.nascom.nasa.gov/goes/fits/' + yr, context=context)
        lines = f.readlines()
        sat_num = []
        for line in lines:
            idx = line.find(datstr)
            if idx != -1:
                sat_num.append(line[idx - 2:idx])
    if type(sat_num) is int:
        sat_num = [str(sat_num)]
    filenames = []
    for sat in sat_num:
        filename = 'go' + sat + datstr + '.fits'
        url = 'https://umbra.nascom.nasa.gov/goes/fits/' + yr + '/' + filename
        f = urlopen.urlopen(url, context=context)
        with open(get_and_create_download_dir() + '/' + filename, 'wb') as g:
            shutil.copyfileobj(f, g)
        filenames.append(get_and_create_download_dir() + '/' + filename)
    pmerit = 0
    for file in filenames:
        gfits = fits.open(file)
        data = gfits[2].data['FLUX'][0][:, 0]
        good, = np.where(data > 1.e-8)
        tsecs = gfits[2].data['TIME'][0]
        merit = len(good)
        date_elements = gfits[0].header['DATE-OBS'].split('/')
        if merit > pmerit:
            print('File:', file, 'is best')
            pmerit = merit
            goes_data = data
            goes_t = Time(date_elements[2] + '-' + date_elements[1] + '-' + date_elements[0]).plot_date + tsecs / 86400.
    try:
        return goes_t, goes_data
    except:
        print('No good GOES data for', datstr)
        return None, None


def ms_clearhistory(msfile):
    from taskinit import tb
    tb_history = msfile + '/HISTORY'
    tb_history_bk = msfile + '/HISTORY_bk'
    if os.path.exists(tb_history_bk):
        os.system('rm -rf {0}'.format(tb_history_bk))
    os.system('cp -r {0} {1}'.format(tb_history, tb_history_bk))
    tb.open(tb_history, nomodify=False)
    nrows = tb.nrows()
    if nrows > 0:
        tb.removerows(range(nrows))
    tb.close()


def ms_restorehistory(msfile):
    tb_history = msfile + '/HISTORY'
    os.system('rm -rf {0}'.format(tb_history))
    os.system('mv {0}_bk {0}'.format(tb_history))


aiadir_default = '/srg/data/sdo/aia/level1/'


def get_mapcube_time(mapcube):
    from astropy.time import Time
    t = []
    for idx, mp in enumerate(mapcube):
        if mp.meta.has_key('t_obs'):
            tstr = mp.meta['t_obs']
        else:
            tstr = mp.meta['date-obs']
        t.append(tstr)
    return Time(t)


def uniq(lst):
    last = object()
    nlst = []
    for item in lst:
        if item == last:
            continue
        nlst.append(item)
        last = item
    return nlst


def get_colorbar_params(fbounds, stepfactor=1):
    cfq = fbounds['cfreqs']
    nspws = len(cfq)
    freqmask = []
    if 'bounds_lo' in fbounds.keys():
        # bounds = np.hstack((fbounds['bounds_lo'], fbounds['bounds_hi'][-1]))
        bounds = np.hstack((fbounds['bounds_lo'], fbounds['bounds_hi']))
        bounds = np.sort(uniq(list(set(bounds))))
        # print(bounds)
        # bounds2 = np.hstack((fbounds['bounds_lo'], fbounds['bounds_hi']))
        # bounds2 = np.array(uniq(list(set(bounds2))))
        # print(bounds2)
        if bounds[0] > fbounds['bounds_all'][0]:
            bounds = np.hstack((fbounds['bounds_all'][0], bounds))
        if bounds[-1] < fbounds['bounds_all'][-1]:
            bounds = np.hstack((bounds, fbounds['bounds_all'][-1]))
        if fbounds['bounds_lo'][0] > fbounds['bounds_all'][0]:
            fbd_lo, fbd_hi = fbounds['bounds_all'][0], fbounds['bounds_lo'][0]
        else:
            fbd_lo, fbd_hi = None, None
        freqmask.append((fbd_lo, fbd_hi))
        for sidx, s in enumerate(fbounds['bounds_hi']):
            if sidx == nspws - 1:
                if fbounds['bounds_hi'][sidx] < fbounds['bounds_all'][-1]:
                    fbd_lo, fbd_hi = fbounds['bounds_hi'][sidx], fbounds['bounds_all'][-1]
                else:
                    fbd_lo, fbd_hi = None, None
            else:
                fbd_lo, fbd_hi = s, fbounds['bounds_lo'][sidx + 1]
            freqmask.append((fbd_lo, fbd_hi))
    else:
        bounds = fbounds['bounds_all']
    step = int(nspws // int(50 / stepfactor))
    if step > 1:
        ticks = cfq[::step]
    else:
        ticks = cfq
    fmin = fbounds['bounds_all'][0]
    fmax = fbounds['bounds_all'][-1]
    # if vmin not in ticks:
    #     ticks = np.hstack((vmin, ticks))
    # if vmax not in ticks:
    #     ticks = np.hstack((ticks, vmax))
    return ticks, bounds, fmax, fmin, freqmask


def parse_trange(trange):
    """
    Parses the time range input and returns start and end times.

    :param trange: Time range for the query.
    :type trange: list or astropy.time.Time
    :raises ValueError: If `trange` is empty or the start time is later than the end time.
    :return: Start and end times as astropy Time objects.
    :rtype: tuple of astropy.time.Time
    """
    if not trange:
        raise ValueError("trange cannot be empty. Provide at least one timestamp.")

    trange = Time(trange) if not isinstance(trange, Time) else trange

    if isinstance(trange, Time) and trange.isscalar:
        delta = np.array([-24, 24]) / 86400.0  # seconds to days
        tst, ted = Time(trange.jd + delta, format='jd')
    elif len(trange) == 2:
        tst, ted = trange[0], trange[1]
        if ted < tst:
            raise ValueError("Start time must occur earlier than end time!")
    else:
        tst, ted = trange[0], trange[-1]
        if ted < tst:
            raise ValueError("Start time must occur earlier than end time!")

    return tst, ted


def download_aia_data(trange, wavelengths=[171], cadence=None, outdir='./'):
    """
    Downloads AIA data for specified time range and wavelengths.

    :param trange: Time range for the query.
    :type trange: list or astropy.time.Time
    :param wavelengths: List of wavelengths to download, defaults to [171].
    :type wavelengths: list, optional
    :param cadence: Desired data cadence, defaults to None.
    :type cadence: Quantity, optional
    :param outdir: Output directory for the downloaded data, defaults to './'.
    :type outdir: str, optional
    :return: List of downloaded files.
    :rtype: list
    """
    trange = Time(trange) if not isinstance(trange, Time) else trange
    outdir = Path(outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    if isinstance(wavelengths, (int, float)):
        wavelengths = [wavelengths]
    elif isinstance(wavelengths, str) and wavelengths.lower() == 'all':
        wavelengths = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]

    print(f"{len(wavelengths)} passbands to download")
    downloaded_files = []

    if trange.isscalar:
        # Attempt downloading using JP2 single downloader first
        try:
            for wave in wavelengths:
                downloaded_files.append(
                    download_single_jp2(trange.to_datetime(), wave, outdir, {wave: DataSource.AIA_171}))
            return downloaded_files
        except Exception as e:
            print(f"Single JP2 download failed: {e}")

    # If single download fails or we have a range, use the full range downloader
    tst, ted = parse_trange(trange)

    # Attempt downloading using JP2 downloader
    try:
        downloaded_files += download_jp2(tst, ted, wavelengths, outdir, cadence)
    except Exception as e:
        print(f"JP2 download failed: {e}")
        # Attempt downloading using JSOC
        try:
            downloaded_files += download_using_jsoc(tst, ted, wavelengths, cadence, outdir)
        except Exception as e:
            print(f"JSOC download failed: {e}")
            # Attempt downloading using SunPy Fido
            try:
                downloaded_files += download_using_fido(tst, ted, wavelengths, cadence, outdir)
            except Exception as e:
                print(f"Fido download failed: {e}")

    return downloaded_files


def download_jp2(tstart, tend, wavelengths, outdir, cadence=None):
    """
    Download AIA data in JP2 format using hvpy.

    :param tstart: Start time for the data.
    :type tstart: astropy.time.Time
    :param tend: End time for the data.
    :type tend: astropy.time.Time
    :param wavelengths: List of wavelengths to download.
    :type wavelengths: list
    :param outdir: Output directory for the downloaded data.
    :type outdir: str
    :param cadence: Desired data cadence, defaults to None.
    :type cadence: Quantity, optional
    :return: List of downloaded files.
    :rtype: list
    """

    downloaded_files = []
    for wave in wavelengths:
        if wave not in data_sources_aia:
            print(f"Wavelength {wave} not supported for JP2 download.")
            continue

        tdt = timedelta(seconds=12 if wave not in [1600, 1700] else 24) if cadence is None else timedelta(
            seconds=cadence.value)

        st = tstart.to_datetime()
        et = tend.to_datetime()
        timestamps = [st + i * tdt for i in range(int((et - st) / tdt) + 1)]
        for timestamp in tqdm(timestamps, desc=f"Downloading {wave}"):
            downloaded_files.append(download_single_jp2(timestamp, wave, outdir, data_sources_aia))
    return downloaded_files


def download_single_jp2(timestamp, wave, outdir, data_sources):
    """
    Download a single JP2 file for the given timestamp and wavelength.

    :param timestamp: Timestamp for the data.
    :type timestamp: datetime
    :param wave: Wavelength of the data.
    :type wave: int
    :param outdir: Output directory for the downloaded data.
    :type outdir: str
    :param data_sources: Dictionary of data sources for different wavelengths.
    :type data_sources: dict
    :return: Path to the downloaded JP2 file.
    :rtype: str
    """
    from datetime import datetime
    try:
        product = 'aia.lev1_euv_12s' if wave not in [1600, 1700] else 'aia.lev1_uv_24s'
        res = hvpy.getClosestImage(date=timestamp, sourceId=data_sources[wave])
        timestamp_real = datetime.strptime(res['date'], '%Y-%m-%d %H:%M:%S')
        timestr = timestamp_real.strftime("%Y-%m-%dT%H%M%S")
        filename = f"{product}.{timestr}Z.{wave:.0f}.image_lev1.jp2"
        outfile = os.path.join(outdir, filename)
        if os.path.exists(outfile):
            print(f"File {outfile} already exists.")
        else:
            print(f"Downloading {filename}")
            hvpy.save_file(hvpy.getJP2Image(timestamp_real, sourceId=data_sources[wave]), filename=outfile,
                           overwrite=False)
        return outfile
    except Exception as e:
        print(f"Failed to download for time {timestamp} and wavelength {wave}: {e}")
        return None


def download_using_jsoc(tst, ted, wavelengths, cadence, outdir):
    """
    Download AIA data using the JSOC API.

    :param tst: Start time for the data.
    :type tst: astropy.time.Time
    :param ted: End time for the data.
    :type ted: astropy.time.Time
    :param wavelengths: List of wavelengths to download.
    :type wavelengths: list
    :param cadence: Desired data cadence, defaults to None.
    :type cadence: Quantity, optional
    :param outdir: Output directory for the downloaded data.
    :type outdir: pathlib.Path
    :raises drms.ExportError: If JSOC export fails.
    :return: List of downloaded files.
    :rtype: list
    """
    client = drms.Client(email=os.environ.get("JSOC_EMAIL"))
    downloaded_files = []
    for wave in wavelengths:
        query_str = build_jsoc_query(wave, tst, ted, cadence)
        print(f"JSOC query: {query_str}")
        try:
            result = client.export(query_str, method="url", protocol="fits", email='suncasa-group@njit.edu')
            files = result.download(outdir.as_posix())
            downloaded_files.extend(files)
            print(f"JSOC: {files}")
        except drms.ExportError as e:
            print(f"JSOC export failed for wavelength {wave}: {e}")
    return downloaded_files


def build_jsoc_query(wave, tst, ted, cadence):
    """
    Build JSOC query string for data export.

    :param wave: Wavelength of the data.
    :type wave: int
    :param tst: Start time for the data.
    :type tst: astropy.time.Time
    :param ted: End time for the data.
    :type ted: astropy.time.Time
    :param cadence: Desired data cadence, defaults to None.
    :type cadence: Quantity, optional
    :return: JSOC query string.
    :rtype: str
    """
    product = 'aia.lev1_euv_12s' if wave not in [1600, 1700] else 'aia.lev1_uv_24s'
    cadence_str = f"@{cadence.value:.0f}s" if cadence else ""
    return f"{product}[{tst.isot}Z-{ted.isot}Z{cadence_str}][?WAVELNTH={int(wave)}?]{{image}}"


def download_using_fido(tst, ted, wavelengths, cadence, outdir):
    """
    Download AIA data using SunPy Fido.

    :param tst: Start time for the data.
    :type tst: astropy.time.Time
    :param ted: End time for the data.
    :type ted: astropy.time.Time
    :param wavelengths: List of wavelengths to download.
    :type wavelengths: list
    :param cadence: Desired data cadence, defaults to None.
    :type cadence: Quantity, optional
    :param outdir: Output directory for the downloaded data.
    :type outdir: pathlib.Path
    :return: List of downloaded files.
    :rtype: list
    """
    downloaded_files = []
    for wave in wavelengths:
        wave_range = a.Wavelength(wave * u.angstrom, wave * u.angstrom)
        try:
            result = Fido.search(a.Time(tst.iso, ted.iso), wave_range, a.Instrument('AIA'))
            files = Fido.fetch(result, path=str(outdir / '{file}.fits'))
            downloaded_files.extend(files)
            print(f"Fido: Downloaded {len(files)} files for wavelength {wave}")
        except Exception as e:
            print(f"Fido download failed for wavelength {wave}: {e}")
    return downloaded_files


def trange2aiafits(trange, aiawave, aiadir):
    """
    Retrieve or download AIA FITS files for the specified time range and wavelength.

    :param trange: Time range for the query.
    :type trange: list or astropy.time.Time
    :param aiawave: Wavelength of the AIA data.
    :type aiawave: int
    :param aiadir: Directory to search for the data.
    :type aiadir: str
    :return: Path to the AIA FITS files.
    :rtype: str or None
    """
    trange = parse_trange(trange)
    if (trange[1] - trange[0]).jd < 12.0 / 24.0 / 3600.0:
        trange = Time(np.mean(trange.jd) + np.array([-1.0, 1.0]) * 6.0 / 24.0 / 3600.0, format='jd')
    aiafits = DButil.readsdofile(datadir=aiadir_default, wavelength=aiawave, trange=trange, isexists=True)
    if not aiafits:
        aiafits = DButil.readsdofileX(datadir='./', wavelength=aiawave, trange=trange, isexists=True)
    if not aiafits:
        aiafits = DButil.readsdofileX(datadir=aiadir, wavelength=aiawave, trange=trange, isexists=True)
    if not aiafits:
        download_aia_data(trange, wavelengths=[aiawave])
        aiafits = DButil.readsdofileX(datadir='./', wavelength=aiawave, trange=trange, isexists=True)
    return aiafits


def parse_rdata(rdata, meta, icmap=None, stokes='I,V', sp=None, show_warnings=False):
    '''
    rdata, meta: (required) data and header of a fits file readed using ndfits.read
    icmap: (optional) colormap for plotting radio images
    stokes: (optional) polarizations to visualizing
    sp: (optional) the spectral window to plot, if there are multiple spectral windows in the fits file
    Returns
    -------
    cmaps: A dictionary contains colormaps for the selected polarizations
    datas: A dictionary contains image data for the selected polarizations
    '''
    if not show_warnings:
        import warnings
        warnings.filterwarnings("ignore")
    ndim = meta['naxis']
    pol_axis = meta['pol_axis']
    npol_fits = meta['npol']
    pols = stokes.split(',')
    npol_in = len(pols)
    if npol_fits != npol_in:
        warnings.warn("The input stokes setting is not matching those in the provided fitsfile.")
    freq_axis = meta['freq_axis']
    nfreq = meta['nfreq']
    cmap_I = plt.get_cmap(icmap)
    cmap_V = plt.get_cmap('RdBu')
    datas = {}
    cmaps = {}
    slc1 = [slice(None)] * ndim
    slc2 = [slice(None)] * ndim
    if freq_axis is not None:
        if sp is None:
            pass
        else:
            slc1[freq_axis] = slice(sp, sp + 1)
    if npol_in > 1:
        if npol_fits > 1:
            if stokes == 'I,V':
                slc1[pol_axis] = slice(0, 1)
                slc2[pol_axis] = slice(1, 2)
                datas['I'] = (rdata[tuple(slc1)] + rdata[tuple(slc2)]) / 2.0
                datas['V'] = (rdata[tuple(slc1)] - rdata[tuple(slc2)]) / (rdata[tuple(slc1)] + rdata[tuple(slc2)])
                cmaps['I'] = cmap_I
                cmaps['V'] = cmap_V
            else:
                slc1[pol_axis] = slice(polmap[pols[0]], polmap[pols[0]] + 1)
                slc2[pol_axis] = slice(polmap[pols[1]], polmap[pols[1]] + 1)
                datas[pols[0]] = rdata[tuple(slc1)]
                datas[pols[1]] = rdata[tuple(slc2)]
                cmaps[pols[0]] = cmap_I
                cmaps[pols[1]] = cmap_I
        else:
            if pol_axis is None:
                datas[pols[0]] = rdata[tuple(slc1)]
                cmaps[pols[0]] = cmap_I
            else:
                slc1[pol_axis] = slice(0, 1)
                datas[pols[0]] = rdata[tuple(slc1)]
                cmaps[pols[0]] = cmap_I
    else:
        if npol_fits > 1:
            if pols[0] in ['I', 'V']:
                slc1[pol_axis] = slice(0, 1)
                slc2[pol_axis] = slice(1, 2)
                if pols[0] == 'I':
                    datas['I'] = (rdata[tuple(slc1)] + rdata[tuple(slc2)]) / 2.0
                    cmaps['I'] = cmap_I
                else:
                    datas['V'] = (rdata[tuple(slc1)] - rdata[tuple(slc2)]) / (rdata[tuple(slc1)] + rdata[tuple(slc2)])
                    cmaps['V'] = cmap_V
            else:
                slc1[pol_axis] = slice(polmap[pols[0]], polmap[pols[0]] + 1)
                datas[pols[0]] = rdata[tuple(slc1)]
                cmaps[pols[0]] = cmap_I
        else:
            if pol_axis is None:
                datas[pols[0]] = rdata[tuple(slc1)]
                cmaps[pols[0]] = cmap_I
            else:
                slc1[pol_axis] = slice(0, 1)
                datas[pols[0]] = rdata[tuple(slc1)]
                cmaps[pols[0]] = cmap_I
    if pol_axis is not None:
        if not py3:
            iterop_ = datas.iteritems()
        else:
            iterop_ = datas.items()
        for k, v in iterop_:
            datas[k] = np.squeeze(v, axis=pol_axis)
    return cmaps, datas


def mk_qlook_image(vis, ncpu=1, timerange='', twidth=12, stokes='I,V', antenna='', imagedir=None, spws=[], toTb=True,
                   sclfactor=1.0, overwrite=True, doslfcal=False, datacolumn='data',
                   phasecenter='', robust=0.0, niter=500, gain=0.1, imsize=[512], cell=['5.0arcsec'], pbcor=True,
                   reftime='', restoringbeam=[''],
                   refbmsize=70., reffreq=1.0, minbmsize=4.0,
                   mask='', docompress=False, wrapfits=True,
                   uvrange='', subregion='', c_external=True, show_warnings=False):
    vis = [vis]
    subdir = ['/']
    outfits_list = []

    for idx, f in enumerate(vis):
        if f[-1] == '/':
            vis[idx] = f[:-1]

    if not imagedir:
        imagedir = './'
    msfile = vis[0]
    ms.open(msfile)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    imres = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': [], 'Freq': [],
             'Obs': []}
    # axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    # freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    # freqInfo_ravel = freqInfo.ravel()
    ms.close()
    nspw = len(spwInfo)

    if not spws:
        if observatory == 'EVLA':
            spws = list(np.arange(nspw).astype(str))
        if observatory == 'EOVSA':
            spws = ['1~5', '6~10', '11~15', '16~25']
    if observatory == 'EOVSA':
        sto = stokes.split(',')
        try:
            next(s for s in ['XX', 'YY'] if s in stokes)
        except:
            nsto = len(sto)
            if nsto == 0:
                stokes = 'XX'
            elif nsto == 1:
                stokes = 'XX'
            elif nsto == 2:
                stokes = 'XX,YY'
            elif nsto == 4:
                stokes = 'XX,YY,XY,YX'
            else:
                stokes = 'XX,YY'
            print('Provide stokes: {}. However EOVSA has linear feeds. Force stokes to be {}'.format(','.join(sto),
                                                                                                     stokes))

    # pdb.set_trace()
    msfilebs = os.path.basename(msfile)
    imdir = imagedir + subdir[0]
    if not os.path.exists(imdir):
        os.makedirs(imdir)
    if doslfcal:
        pass
        # slfcalms = './' + msfilebs + '.rr'
        # split(msfile, outputvis=slfcalms, datacolumn='corrected', correlation='RR')

    cfreqs = mstools.get_bandinfo(msfile, spws)
    if restoringbeam == ['']:
        restoringbms = [''] * nspw
    else:
        if observatory == 'EOVSA':
            restoringbms = mstools.get_bmsize(cfreqs, refbmsize=refbmsize, reffreq=reffreq, minbmsize=minbmsize)
        else:
            restoringbms = [''] * nspw
    for sp, spw in enumerate(spws):
        spwran = [s.zfill(2) for s in spw.split('~')]

        spw_ = spw.split('~')
        if len(spw_) == 2:
            freqran = [(spwInfo['{}'.format(s)]['RefFreq'] + spwInfo['{}'.format(s)]['TotalWidth'] / 2.0) / 1.0e9 for s
                       in spw.split('~')]
        elif len(spw_) == 1:
            s = spw_[0]
            freqran = np.array([0, spwInfo['{}'.format(s)]['TotalWidth']]) + spwInfo['{}'.format(s)]['RefFreq']
            freqran = freqran / 1.0e9
            freqran = list(freqran)
        else:
            raise ValueError("Keyword 'spw' in wrong format")

        if restoringbms[sp] == '':
            restoringbm = ['']
        else:
            restoringbm = ['{:.1f}arcsec'.format(restoringbms[sp])]
        cfreq = cfreqs[sp]
        if cell == ['5.0arcsec'] and imsize == [512]:
            if cfreq < 10.:
                imsize = 512
                cell = ['5arcsec']
            else:
                imsize = 1024
                cell = ['2.5arcsec']
        if len(spwran) == 2:
            spwstr = spwran[0] + '~' + spwran[1]
        else:
            spwstr = spwran[0]

        imagesuffix = '.spw' + spwstr.replace('~', '-')
        # if cfreq > 10.:
        #     antenna = antenna + ';!0&1;!0&2'  # deselect the shortest baselines
        sto = stokes.replace(',', '')
        if c_external:
            cleanscript = os.path.join(imdir, 'ptclean_external.py')
            resfile = os.path.join(imdir, os.path.basename(msfile) + '.res.npz')
            os.system('rm -rf {}'.format(cleanscript))
            inpdict = {'vis': msfile,
                       'imageprefix': imdir,
                       'imagesuffix': imagesuffix,
                       'timerange': timerange,
                       'twidth': twidth,
                       'spw': spw,
                       'uvrange': uvrange,
                       'restoringbeam': restoringbm,
                       'mask': mask,
                       'ncpu': ncpu,
                       'niter': niter,
                       'gain': gain,
                       'antenna': antenna,
                       'imsize': imsize,
                       'cell': cell,
                       'stokes': sto,
                       'doreg': True,
                       'usephacenter': True,
                       'phasecenter': phasecenter,
                       'docompress': docompress,
                       'reftime': reftime,
                       'overwrite': overwrite,
                       'toTb': toTb,
                       'sclfactor': sclfactor,
                       'datacolumn': datacolumn,
                       'pbcor': pbcor,
                       'subregion': subregion,
                       'weighting': 'briggs',
                       'robust': robust}
            for key, val in inpdict.items():
                if type(val) is str:
                    inpdict[key] = '"{}"'.format(val)
            fi = open(cleanscript, 'wb')
            fi.write('from ptclean3_cli import ptclean3_cli as ptclean3 \n')
            fi.write('import numpy as np \n')
            if not show_warnings:
                fi.write('import warnings \n')
                fi.write('warnings.filterwarnings("ignore") \n')
            ostrs = []
            if not py3:
                iterop_ = inpdict.iteritems()
            else:
                iterop_ = inpdict.items()
            for k, v in iterop_:
                ostrs.append('{}={}'.format(k, v))
            ostr = ','.join(ostrs)
            fi.write('res = ptclean3({}) \n'.format(ostr))
            fi.write('np.savez("{}",res=res) \n'.format(resfile))
            fi.close()

            os.system('casa --nologger -c {}'.format(cleanscript))
            res = np.load(resfile)
            res = res['res'].item()
        else:
            res = ptclean(vis=msfile,
                          imageprefix=imdir,
                          imagesuffix=imagesuffix,
                          timerange=timerange,
                          twidth=twidth,
                          spw=spw,
                          uvrange=uvrange,
                          restoringbeam=restoringbm,
                          mask=mask,
                          ncpu=ncpu,
                          niter=niter,
                          gain=gain,
                          antenna=antenna,
                          imsize=imsize,
                          cell=cell,
                          stokes=sto,
                          doreg=True,
                          usephacenter=True,
                          phasecenter=phasecenter,
                          docompress=docompress,
                          reftime=reftime,
                          overwrite=overwrite,
                          toTb=toTb,
                          sclfactor=sclfactor,
                          datacolumn=datacolumn,
                          pbcor=pbcor,
                          subregion=subregion,
                          weighting='briggs',
                          robust=robust)
            ## perform a sanity check if all intended images are created.
            ntrial = 0
            while np.count_nonzero(np.array(res['Succeeded']) == False) > 0 and ntrial < 2:
                print('some intended images are not created. Trying to create them again.... trial #{}'.format(ntrial))
                res = ptclean(vis=msfile,
                              imageprefix=imdir,
                              imagesuffix=imagesuffix,
                              timerange=timerange,
                              twidth=twidth,
                              spw=spw,
                              uvrange=uvrange,
                              restoringbeam=restoringbm,
                              mask=mask,
                              ncpu=ncpu,
                              niter=niter,
                              gain=gain,
                              antenna=antenna,
                              imsize=imsize,
                              cell=cell,
                              stokes=sto,
                              doreg=True,
                              usephacenter=True,
                              phasecenter=phasecenter,
                              docompress=docompress,
                              reftime=reftime,
                              overwrite=overwrite,
                              toTb=toTb,
                              sclfactor=sclfactor,
                              datacolumn=datacolumn,
                              pbcor=pbcor,
                              subregion=subregion,
                              weighting='briggs',
                              robust=robust)
                ntrial += 1
            if np.count_nonzero(np.array(res['Succeeded']) == False) > 0:
                print('Some intended images can not be created after {} trails. Preceed with missing images.'.format(
                    ntrial))

        if res:
            imres['Succeeded'] += res['Succeeded']
            imres['BeginTime'] += res['BeginTime']
            imres['EndTime'] += res['EndTime']
            imres['ImageName'] += res['ImageName']
            imres['Spw'] += [spwstr] * len(res['ImageName'])
            imres['Vis'] += [msfile] * len(res['ImageName'])
            imres['Freq'] += [freqran] * len(res['ImageName'])
            imres['Obs'] += [observatory] * len(res['ImageName'])
        else:
            return None

    # save it for debugging purposes
    imresfile = os.path.join(imagedir, '{}.imres.npz'.format(os.path.basename(msfile)))
    np.savez(imresfile, imres=imres)

    if wrapfits:
        imres = read_imres(imresfile)
        nt = len(imres['plttimes'])
        imresnew = {'images': [], 'btimes': [], 'etimes': [], 'spw': imres['spws'], 'vis': imres['vis'],
                    'freq': imres['freq'], 'obs': imres['obs']}
        imresnew['btimes'] = imres['btimes'][:, 0].iso
        imresnew['etimes'] = imres['etimes'][:, 0].iso
        qlookallbdfitsdir = imagedir.replace('qlookfits', 'qlookallbdfits')
        if not os.path.exists(qlookallbdfitsdir):
            os.makedirs(qlookallbdfitsdir)
        for tidx in tqdm(range(nt)):
            fitsfiles = [os.path.join(imagedir, os.path.basename(fitsfile)) for fitsfile in imres['images'][tidx, :]]
            outname = imres['btimes'][tidx, 0].strftime('{}.%Y%m%dT%H%M%S.%f.allbd.fits'.format(observatory))
            imresnew['images'].append(os.path.join(imagedir, outname))
            outfits = os.path.join(qlookallbdfitsdir, outname)
            ndfits.wrap(fitsfiles, outfitsfile=outfits, docompress=True, imres=imres)
            outfits_list.append(outfits)
        imresnewfile = os.path.join(qlookallbdfitsdir, os.path.basename(imresfile))
        np.savez(imresnewfile, imres=imresnew)
        os.system('rm -rf {}'.format(imagedir))
        os.system('mv {} {}'.format(qlookallbdfitsdir, imagedir))
        outfits_list = [item.replace(qlookallbdfitsdir, imagedir) for item in outfits_list]
        return imresnew, outfits_list
    else:
        return imres, outfits_list


def radio_image_clevels(imres, clevels=[0.2, 0.4, 0.6, 0.8], snr=2., timerange_bkg=None, overbright=1e12):
    '''
    Function to calculate contour levels for plotting from a timeseries of the input radio images.
    The contour levels are fixed to those relative to the flare peak.
    :param imres: a dictionary contains the information of the FITS files (produced by qlookplot)
    :param clevels: relative contour levels with regard to the flare peak
    :snr: minimum enhancement during the peak time against the background level. Default 3
    :timerange_bkg: timerange of the background. Example: ['2024-01-03T17:46', '2024-01-03T17:47']. 
                    Default to None. If None, use the median of all FITS files
    :overbright: Brightness temperature threshold in K, above which the peak values are deemed corrupted and not considered
                    Default to 1e12 K.
    '''
    max_fs = []
    for i, f in enumerate(imres['images']):
        meta, rdata = ndfits.read(f)
        max_f = np.nanmax(rdata, axis=(1, 2))
        max_fs.append(max_f)
        if i == 0:
            peaks = np.zeros_like(max_f)
            idx0, = np.where((max_f < overbright))
            peaks[idx0] = max_f[idx0]
        idx, = np.where((max_f > peaks) & (max_f < overbright))
        if len(idx) > 0:
            peaks[idx] = max_f[idx]

    if timerange_bkg is None:
        peaks_bkg = max_fs[0]
    else:
        idx0 = max(np.argmin(np.abs(Time(imres['btimes']).mjd - Time(timerange_bkg[0]).mjd)), 0)
        idx1 = min(np.argmin(np.abs(Time(imres['btimes']).mjd - Time(timerange_bkg[1]).mjd)), len(imres['btimes']) - 1)
        peaks_bkg = np.nanmedian(max_fs[idx0:idx1], axis=0)

    clevelsfix = []
    for p, p_bkg in zip(peaks, peaks_bkg):
        if p / p_bkg > snr:
            clevelsfix.append(np.array(clevels) * p)
        else:
            clevelsfix.append(np.array([2.]) * p)

    return clevelsfix


def get_data_limits(data, dmin, dmax, dnorm):
    if dmax is None:
        dmax = np.nanpercentile(data, 99)
    if dmin is None:
        dmin = np.nanpercentile(data, 1)
        if dnorm == 'log' and dmin <= 0:
            dmin = np.nanpercentile(data[data > 0], 1)
    return dmin, dmax


def setup_axes(fig, npols, nspw, plotaia):
    if npols == 1:
        hnspw = max(nspw // 2, 1)
        ncols = hnspw
        nrows = 4
        gs = gridspec.GridSpec(nrows, ncols, height_ratios=[4, 4, 1, 1])

        if nspw <= 1 or plotaia:
            axs = [fig.add_subplot(gs[:2, :hnspw])]
        else:
            axs = [fig.add_subplot(gs[0, 0])]
            for ll in range(1, nspw):
                axs.append(fig.add_subplot(gs[ll // hnspw, ll % hnspw], sharex=axs[0], sharey=axs[0]))

        axs_dspec = [fig.add_subplot(gs[2:, :])]

    elif npols == 2:
        hnspw = max(nspw // 2, 1)
        ncols = hnspw + 2
        nrows = 4
        gs = gridspec.GridSpec(nrows, ncols, height_ratios=[1, 1, 1, 1])

        if nspw <= 1 or plotaia:
            axs = [fig.add_subplot(gs[:2, 2:]), fig.add_subplot(gs[2:, 2:])]
        else:
            axs = [fig.add_subplot(gs[0, 2])]
            for ll in range(1, nspw):
                axs.append(fig.add_subplot(gs[ll // hnspw, ll % hnspw + 2], sharex=axs[0], sharey=axs[0]))
            for ll in range(nspw):
                axs.append(fig.add_subplot(gs[ll // hnspw + 2, ll % hnspw + 2], sharex=axs[0], sharey=axs[0]))

        axs_dspec = [fig.add_subplot(gs[:2, :2]), fig.add_subplot(gs[2:, :2])]

    return axs, axs_dspec, gs


def plt_qlook_image(imres, timerange='', spwplt=None, figdir='./qlookimgs/', specdata=None,
                    verbose=False, stokes='I,V', fov=None,
                    imax=None, imin=None, icmap='RdYlBu', inorm='linear',
                    amax=None, amin=None, acmap='gray_r', anorm='log',
                    dmax=None, dmin=None, dcmap='viridis', dnorm='linear',
                    sclfactor=1.0,
                    nclevels=3,
                    clevels=None, clevelsfix=None, aiafits='', aiadir=None, aiawave=171, plotaia=True,
                    freqbounds=None, moviename='',
                    alpha_cont=1.0, custom_mapcubes=[], opencontour=False, movieformat='html', ds_normalised=False,
                    minsnr=5, timtol=10. / 60. / 24., overwrite=True):
    """
    Generate quick-look images of solar radio data with optional AIA overlays.

    This function plots dynamic spectra, AIA overlays, and radio images for a specified time range and field of view.
    It can handle different plotting configurations, including stokes parameters and plotting custom contour levels.

    :param imres: Input image results dictionary from suncasa.utils.qlookplot.
    :type imres: dict or str
    :param timerange: Time range for plotting, defaults to ''. If '', use the entire timerange of the input imres.
                      Can be a string or `astropy.time.Time` object, e.g., '2024-01-03T17:46~2024-01-03T17:47'.
    :type timerange: str or `astropy.time.Time`, optional
    :param spwplt: List of spectral windows to plot, defaults to None.
    :type spwplt: list, optional
    :param figdir: Directory to save the generated plots, defaults to './qlookimgs/'.
    :type figdir: str, optional
    :param specdata: Spectral data for plotting dynamic spectra.
    :type specdata: object, optional
    :param verbose: If True, print verbose output, defaults to True.
    :type verbose: bool, optional
    :param stokes: Comma-separated stokes parameters to plot, defaults to 'I,V'.
    :type stokes: str, optional
    :param fov: Field of view for the plot, defaults to None.
    :type fov: list, optional
    :param imax: Maximum value for radio image plot, defaults to None.
    :type imax: float, optional
    :param imin: Minimum value for radio image plot, defaults to None.
    :type imin: float, optional
    :param icmap: Colormap for radio image, defaults to 'RdYlBu'.
    :type icmap: str, optional
    :param inorm: Normalization for radio images. Must be 'linear' or 'log'.
    :type inorm: str, optional
    :param amax, amin: Max and min values for AIA data normalization.
    :type amax: float, optional
    :type amin: float, optional
    :param acmap: Colormap for AIA images, defaults to None.
    :type acmap: str, optional
    :param anorm: Normalization for AIA data. Must be 'linear' or 'log'.
    :type anorm: str, optional
    :param dmax, dmin: Max and min values for dynamic spectra normalization.
    :type dmax: float, optional
    :type dmin: float, optional
    :param dcmap: Colormap for dynamic spectra, defaults to 'viridis'.
    :type dcmap: str, optional
    :param dnorm: Normalization for dynamic spectra. Must be 'linear' or 'log'.
    :type dnorm: str, optional
    :param sclfactor: Scaling factor for spectral data, defaults to 1.0.
    :type sclfactor: float, optional
    :param nclevels: Number of contour levels for radio images, calculated between `imin` and `imax`.
                     Used if `clevels` and `clevelsfix` are not provided.
    :type nclevels: int, optional
    :param clevels: Contour levels for radio images, defined as percentages of the maximum value in the image.
                    Overrides `nclevels` if provided. For example, `[0.3, 0.5, 0.8]` plots contours at 30%, 50%,
                    and 80% of the maximum value for each spectral window.
    :type clevels: list, optional
    :param clevelsfix: Fixed contour levels for each spectral window, calculated using `radio_image_clevels`.
                       These levels are relative to the flare peak and defined as percentages of the peak value.
                       Overrides both `nclevels` and `clevels` if provided.
    :type clevelsfix: list, optional
    :param aiafits: Path to AIA FITS files, defaults to ''.
    :type aiafits: str, optional
    :param aiadir: Directory to search for AIA FITS files, defaults to None.
    :type aiadir: str, optional
    :param aiawave: AIA wavelength to use, defaults to 171.
    :type aiawave: int, optional
    :param plotaia: If True, plot AIA data, defaults to True.
    :type plotaia: bool, optional
    :param freqbounds: Frequency bounds for plotting, defaults to None.
    :type freqbounds: dict, optional
    :param moviename: Name for the output movie, defaults to ''.
    :type moviename: str, optional
    :param alpha_cont: Alpha value for contours, defaults to 1.0.
    :type alpha_cont: float, optional
    :param custom_mapcubes: List of custom mapcubes to plot, defaults to [].
    :type custom_mapcubes: list, optional
    :param opencontour: If True, use open contours, defaults to False.
    :type opencontour: bool, optional
    :param movieformat: Format for the output movie, defaults to 'html'.
    :type movieformat: str, optional
    :param ds_normalised: If True, normalize the dynamic spectra, defaults to False.
    :type ds_normalised: bool, optional
    :param minsnr: Minimum signal-to-noise ratio for plotting, defaults to 5.
    :type minsnr: int, optional
    :param timtol: Time tolerance for matching AIA data, defaults to 10./60./24. (in days).
    :type timtol: float, optional
    :param overwrite: If True, overwrite existing plots, defaults to True.
    :type overwrite: bool, optional

    :raises ValueError: If the input parameters are not valid.

    :return: Path to the generated movie file.
    :rtype: str

    Example usage:
    --------------
    ```python
    from suncasa.utils import qlookplot as ql
    from astropy.time import Time
    import numpy as np
    import matplotlib.pyplot as plt

    # Load imres from the npz file from a prior qlookplot run
    imres = np.load('./qlookfits/IDB20240514_1645-1738.XXYY.cal.1s.ms.slfcaled.imres.npz', allow_pickle=True)['imres'].item()

    # Load dynamic spectrum
    from suncasa.dspec import dspec
    d = dspec.Dspec()
    d.read('eovsa.spec.flare_id_20240514164700.fits', source='eovsa')

    # Calculate contour levels
    # Define the background time to determine bands with sufficient flare enhancement
    timerange_bkg = Time(['2024-05-14T16:39', '2024-05-14T16:42'])
    clevelsfix = ql.radio_image_clevels(imres, timerange_bkg=timerange_bkg, snr=2)

    # Define timerange, fov, dmin, dmax (for spectrogram), and AIA normalization for plotting images
    fov = [[760., 1060.], [-430., -130.]]
    fov = None
    ql.plt_qlook_image(imres, stokes='I', figdir='./qlookimgs/',
                       specdata=d, icmap=plt.get_cmap('RdYlBu'), aiadir='./', fov=fov,
                       opencontour=True, dcmap='viridis', dnorm='log', clevelsfix=clevelsfix, movieformat='mp4')
    ```
    """
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    if sunpy1:
        from sunpy.coordinates import sun
    else:
        from sunpy import sun
    import astropy.units as u

    if isinstance(icmap, str):
        icmap = plt.get_cmap(icmap)
    if plotaia:
        nclevels = 2

    pols = stokes.split(',')
    npols = len(pols)
    outmovie = ''

    if 'images' in imres.keys():
        wrapfits = True
    else:
        wrapfits = False
    if freqbounds is None:
        freqbounds = mstools.get_bandinfo(imres['vis'], spw=imres['spw'], returnbdinfo=True)
    cfreqs = freqbounds['cfreqs']
    cfreqs_all = freqbounds['cfreqs_all']
    freq_dist = (cfreqs - cfreqs_all[0]) / (cfreqs_all[-1] - cfreqs_all[0])
    if wrapfits:
        # imresnew = {'images': [], 'btimes': [], 'etimes': [], 'spw': imres['spws'], 'vis': imres['vis'],
        # 'freq': imres['freq'], 'obs': imres['obs']}
        btimes = Time(imres['btimes'])
        etimes = Time(imres['etimes'])
        # tpltidxs, = np.where(np.logical_and(btimes.jd >= t_ran[0].jd, etimes.jd <= t_ran[1].jd))
        # btimes = btimes[tpltidxs]
        # etimes = etimes[tpltidxs]
        plttimes = btimes
        ntime = len(plttimes)
        if 'obs' in imres.keys():
            observatory = imres['obs']
        else:
            observatory = ''
        Spw = imres['spw']
        Freq = imres['freq']
        nspw = len(Spw)
        images_sort = imres['images']
    else:
        # btimes = Time(imres['BeginTime'])
        # etimes = Time(imres['EndTime'])
        # tpltidxs, = np.where(np.logical_and(btimes.jd >= t_ran[0].jd, etimes.jd <= t_ran[1].jd))
        # if not py3:
        #     iterop_ = imres.iteritems()
        # else:
        #     iterop_ = imres.items()
        # for k, v in iterop_:
        #     imres[k] = list(np.array(v)[tpltidxs])
        if 'Obs' in imres.keys():
            observatory = imres['Obs'][0]
        else:
            observatory = ''

        Spw = sorted(list(set(imres['Spw'])))
        nspw = len(Spw)
        imres['Freq'] = [list(ll) for ll in imres['Freq']]
        Freq = sorted(uniq(imres['Freq']))

        plttimes = list(set(imres['BeginTime']))
        plttimes = Time(sorted(plttimes))
        ntime = len(plttimes)
        # sort the imres according to time
        images = np.array(imres['ImageName'])
        btimes = Time(imres['BeginTime'])
        etimes = Time(imres['EndTime'])
        spws = np.array(imres['Spw'])
        suc = np.array(imres['Succeeded'])
        inds = btimes.argsort()
        images_sort = images[inds].reshape(ntime, nspw)
        btimes_sort = btimes[inds].reshape(ntime, nspw)
        suc_sort = suc[inds].reshape(ntime, nspw)
        spws_sort = spws[inds].reshape(ntime, nspw)
        btimes = btimes_sort[:, 0]

    if isinstance(timerange, str):
        if timerange == '':
            t_ran = [btimes[0], etimes[-1]]
        else:
            tstart, tend = timerange.split('~')
            t_ran = Time([qa.quantity(tstart, 'd')['value'], qa.quantity(tend, 'd')['value']], format='mjd')
    else:
        t_ran = Time(timerange)

    dt = np.nanmean(np.diff(plttimes.mjd)) * 24. * 60 * 60  ## time cadence in seconds
    print(f'time cadence: {dt:.1f} s')
    if custom_mapcubes:
        cmpc_plttimes_mjd = []
        for cmpc in custom_mapcubes['mapcube']:
            cmpc_plttimes_mjd.append(get_mapcube_time(cmpc).mjd)

    # if verbose:
    #     print('{0:d} figures to plot'.format(len(tpltidxs)))
    plt.ioff()
    if np.iscomplexobj(specdata.data):
        spec = np.abs(specdata.data)
    else:
        spec = specdata.data
    spec = spec * sclfactor
    spec = checkspecnan(spec)
    if spec.ndim == 2:
        (nfreq, ntim) = spec.shape
        nbl = 1
        npol = 2
        spec_ = np.zeros((npol, nbl, nfreq, ntim))
        spec_[0, 0, :, :] = spec
        spec_[0, 0, :, :] = spec
        spec = spec_
    else:
        (npol, nbl, nfreq, ntim) = spec.shape
    # tidx = range(ntim)
    fidx = range(nfreq)
    freq = specdata.freq_axis
    freqghz = freq / 1e9
    pol = ''.join(pols)
    spec_tim = specdata.time_axis
    tidx, = np.where(np.logical_and(spec_tim > t_ran[0], spec_tim < t_ran[1]))
    spec_tim_plt = spec_tim.plot_date

    spec_plt = None
    polstr = None
    cmaps = None
    dranges = None
    iranges = None

    if npols == 1:
        if pol in ['RR', 'XX']:
            spec_plt = spec[0, 0, :, :]
        elif pol in ['LL', 'YY']:
            spec_plt = spec[1, 0, :, :]
        elif pol == 'I':
            spec_plt = (spec[0, 0, :, :] + spec[1, 0, :, :]) / 2
        elif pol == 'V':
            spec_plt = (spec[0, 0, :, :] - spec[1, 0, :, :]) / (spec[0, 0, :, :] + spec[1, 0, :, :])
        spec_plt = [spec_plt]

        dmin, dmax = get_data_limits(np.array(spec_plt).flatten(), dmin, dmax, dnorm)
        dranges = [[dmin, dmax]]
        iranges = [[imin, imax]]
        cmaps = [plt.get_cmap(dcmap)]
        print(f'Plotting the dynamic spectrum in polarization {pol}')
    elif npols == 2:
        R_plot = np.abs(spec[0, 0, :, :])
        L_plot = np.abs(spec[1, 0, :, :])

        if pol == 'RRLL':
            spec_plt = [R_plot, L_plot]
            polstr = ['RR', 'LL']
        elif pol == 'XXYY':
            spec_plt = [R_plot, L_plot]
            polstr = ['XX', 'YY']
        elif pol == 'IV':
            I_plot = (R_plot + L_plot) / 2
            V_plot = (R_plot - L_plot) / (2 * I_plot)
            spec_plt = [I_plot, V_plot]
            polstr = ['I', 'V']
            cmaps = [plt.get_cmap(dcmap), plt.get_cmap('RdBu')]
            dranges = [[dmin, dmax], [-1, 1]]
            iranges = [[imin, imax], [-1, 1]]
        else:
            raise ValueError(f'Invalid polarization: {pol}')

        dmin, dmax = get_data_limits(np.array(spec_plt).flatten(), dmin, dmax, dnorm)
        if dranges is None:
            dranges = [[dmin, dmax]] * len(spec_plt)
        if iranges is None:
            iranges = [[imin, imax]] * len(spec_plt)
        if cmaps is None:
            cmaps = [plt.get_cmap(dcmap)] * len(spec_plt)

        print(f'Plotting the dynamic spectrum in polarization {pol}')

    else:
        raise ValueError(f'Unsupported value for npols: {npols}')

    fig_size = (10, 12)
    fig = plt.figure(figsize=fig_size)
    axs, axs_dspec, gs = setup_axes(fig, npols, nspw, plotaia)

    for ax in axs + axs_dspec:
        ax.tick_params(direction='out', axis='both')
    # fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    # pdb.set_trace()
    if plotaia:
        '''check if aiafits files exist'''
        if aiafits == '' or aiafits is None:
            if aiadir is not None:
                aiafiles = []
                for i in tqdm(range(ntime)):
                    plttime = btimes[i]
                    aiafile = DButil.readsdofileX(datadir=aiadir, wavelength=aiawave, trange=plttime, isexists=True,
                                                  timtol=timtol)
                    if not aiafile:
                        aiafile = DButil.readsdofile(datadir=aiadir_default, wavelength=aiawave, trange=plttime,
                                                     isexists=True,
                                                     timtol=timtol)
                    if not aiafile:
                        aiafile = DButil.readsdofileX(datadir='./', wavelength=aiawave, trange=plttime, isexists=True,
                                                      timtol=timtol)
                    if aiafile == []:
                        if verbose:
                            print('No AIA fits files found. Downloading AIA data...')
                        aiafile = download_aia_data(trange=plttime, wavelengths=aiawave, cadence=dt * u.second)[0]
                        # print(f'download jp2: {aiafile}')
                    aiafiles.append(aiafile)
                    # print(i, aiafile)

                # if np.count_nonzero(aiafiles) <ntime:
                #     aiafiles = []
                #     print('No AIA fits files found. Downloading AIA data...')
                #     download_aia_data(trange=t_ran, wavelength=aiawave, cadence=dt * u.second)
                #     aiadir = './'
                #     for i in tqdm(range(ntime)):
                #         plttime = btimes[i]
                #         aiafile = DButil.readsdofileX(datadir=aiadir, wavelength=aiawave, trange=plttime,
                #                                       isexists=True,
                #                                       timtol=timtol)
                #         if aiafile == []:
                #             aiafiles.append(None)
                #         else:
                #             aiafiles.append(aiafile)
        else:
            # from suncasa.utils import stackplotX as stp
            # st = stp.Stackplot(aiafits)
            # tjd_aia = st.tplt.jd
            aiafiles = aiafits

    for i, plttime in enumerate(tqdm(plttimes)):
        plt.ioff()
        # plt.clf()
        for ax in axs:
            ax.cla()
        plttime = btimes[i]
        figname = f'{observatory}_qlimg_{plttime.isot.replace(":", "").replace("-", "")[:19]}.png'
        fignameful = os.path.join(figdir, figname)
        if i > 0 and os.path.exists(fignameful) and not overwrite:
            continue
        # tofd = plttime.mjd - np.fix(plttime.mjd)

        # if tofd < 16. / 24. or sum(
        #         suci) < nspw - 2:  # if time of the day is before 16 UT (and 24 UT), skip plotting (because the old antennas are not tracking)
        #     continue
        # fig=plt.figure(figsize=(9,6))
        # fig.suptitle('EOVSA @ '+plttime.iso[:19])
        # if verbose:
        #     print('Plotting image at: ', plttime.iso)

        if plttime == plttimes[0]:
            dspecvspans = []
            for pol in range(npols):
                ax = axs_dspec[pol]
                _dnorm = get_normalization(dranges[pol][0], dranges[pol][1], dnorm)
                cmap_plt = copy.copy(cmaps[pol])
                cmap_plt.set_bad(cmap_plt(0.0))
                im_spec = ax.pcolormesh(spec_tim_plt[tidx], freqghz, spec_plt[pol][:, tidx], cmap=cmap_plt,
                                        norm=_dnorm, rasterized=True)
                ax.set_xlim(spec_tim_plt[tidx[0]], spec_tim_plt[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
                ax.set_ylabel('Frequency [GHz]')
                for idx, freq in enumerate(Freq):  ##dotted lines indicating spws
                    if nspw <= 10:
                        # ax.axhspan(freq[0], freq[1], linestyle='dotted', edgecolor='w', alpha=0.7, facecolor='none')
                        xtext, ytext = ax.transAxes.inverted().transform(
                            ax.transData.transform([spec_tim_plt[tidx[0]], np.mean(freq)]))
                        # ax.text(xtext + 0.01, ytext, 'spw ' + Spw[idx], color='w', transform=ax.transAxes,
                        #         fontweight='bold', ha='left', va='center',
                        #         fontsize=8, alpha=0.5)
                        colors_spws = icmap(freq_dist)
                        ax.annotate('', xy=(xtext - 0.02, ytext), xycoords='axes fraction', xytext=(xtext, ytext),
                                    textcoords='axes fraction',
                                    arrowprops=dict(arrowstyle='-', color=colors_spws[idx], linewidth=2),
                                    horizontalalignment='right')

                ax.text(0.01, 0.98, 'Stokes ' + pols[pol], color='w', transform=ax.transAxes, fontweight='bold',
                        ha='left', va='top')
                dspecvspans.append(ax.axvspan(btimes[i].plot_date, etimes[i].plot_date, color='w', alpha=0.4))
                ax_pos = ax.get_position().extents
                x0, y0, x1, y1 = ax_pos
                h, v = x1 - x0, y1 - y0
                x0_new = x0 + 0.10 * h
                y0_new = y0 + 0.20 * v
                x1_new = x1 - 0.03 * h
                y1_new = y1 - 0.00 * v
                # ax.set_position(mpl.transforms.Bbox([[x0_new, y0_new], [x1_new, y1_new]]))
                if pol == npols - 1:
                    ax.xaxis_date()
                    ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
                    ax.set_xlabel('Time [UT]', fontsize=9)
                    for xlabel in ax.get_xmajorticklabels():
                        xlabel.set_rotation(30)
                        xlabel.set_horizontalalignment("right")
                else:
                    ax_pos = ax.get_position().extents
                    ax_pos2 = axs_dspec[-1].get_position().extents
                    x0, y0, x1, y1 = ax_pos
                    h, v = x1 - x0, y1 - y0
                    x0_new = x0
                    y0_new = ax_pos2[-1]
                    x1_new = x1
                    y1_new = y0_new + v
                    # ax.set_position(mpl.transforms.Bbox([[x0_new, y0_new], [x1_new, y1_new]]))
                    ax.xaxis.set_visible(False)
                if ds_normalised == False:
                    divider = make_axes_locatable(ax)
                    cax_spec = divider.append_axes('right', size='3.0%', pad=0.05)
                    cax_spec.tick_params(direction='out')
                    clb_spec = plt.colorbar(im_spec, ax=ax, cax=cax_spec)
                    clb_spec.set_label('Flux [sfu]')
        else:
            for pol in range(npols):
                xy = dspecvspans[pol].get_xy()
                xy[:, 0][np.array([0, 1, 4])] = btimes[i].plot_date
                xy[:, 0][np.array([2, 3])] = etimes[i].plot_date
                dspecvspans[pol].set_xy(xy)
        if plotaia:
            # pdb.set_trace()
            aia_jp2 = False
            if np.count_nonzero(aiafiles) > 0:
                try:
                    aiafits = aiafiles[i]
                    if verbose:
                        print(f'plotting AIA image at {plttime.iso}: {aiafits}')
                    aiamap = smap.Map(aiafits)
                    if aiafits.endswith('.jp2'):
                        aia_jp2 = True
                    if not aia_jp2:
                        aiamap = DButil.normalize_aiamap(aiamap)
                        data = aiamap.data
                        data[data < 1.0] = 1.0
                        aiamap = smap.Map(data, aiamap.meta)
                except Exception as e:
                    aiamap = None
                    if verbose:
                        if hasattr(e, 'message'):
                            print(e.message)
                        else:
                            print(e)
                    print('error in reading aiafits. Proceed without AIA')
            else:
                aiamap = None
                pass
                # aiamap = st.mapcube[np.nanargmin(np.abs(tjd_aia - plttime.jd))]
        else:
            aiamap = None

        if wrapfits:
            suci_ = True
            meta, rdata = ndfits.read(images_sort[i])
        else:
            suci_ = suc_sort[i]
            meta, rdata = None, None

        colors_spws = icmap(freq_dist)
        spwpltCounts = 0
        # print(f'starting plotting radio images. spwpltCounts: {spwpltCounts} initiated')
        for s, sp in enumerate(Spw):
            if spwplt is not None:
                if sp not in spwplt:
                    continue
            # print(f'spwpltCounts: {spwpltCounts}, s: {s}, sp: {sp}')
            # print(image)
            for pidx, pol in enumerate(pols):
                if wrapfits:
                    suci = suci_
                else:
                    suci = suci_[s]
                if suci:
                    # try:
                    if wrapfits:
                        pass
                    else:
                        image = images_sort[i, s]
                        meta, rdata = ndfits.read(image)
                    cmaps, datas = parse_rdata(rdata, meta, icmap=icmap,
                                               stokes=stokes, sp=s)
                    # except:
                    #     continue

                    rmap = smap.Map(np.squeeze(datas[pol][:, :]), meta['header'])
                else:
                    # make an empty map
                    data = np.zeros((512, 512))
                    hgln_obs = 0.
                    rsun_ref = sun.constants.radius.value
                    if sunpy1:
                        dsun_obs = sun.earth_distance(Time(plttime)).to(u.meter).value
                        rsun_obs = sun.angular_radius(Time(plttime)).value
                        hglt_obs = sun.B0(Time(plttime)).value
                    else:
                        dsun_obs = sun.sunearth_distance(Time(plttime)).to(u.meter).value
                        rsun_obs = sun.solar_semidiameter_angular_size(Time(plttime)).value
                        hglt_obs = sun.heliographic_solar_center(Time(plttime))[1].value
                    header = {"DATE-OBS": plttime.isot, "EXPTIME": 0., "CDELT1": 5., "NAXIS1": 512, "CRVAL1": 0.,
                              "CRPIX1": 257, "CUNIT1": "arcsec",
                              "CTYPE1": "HPLN-TAN", "CDELT2": 5., "NAXIS2": 512, "CRVAL2": 0., "CRPIX2": 257,
                              "CUNIT2": "arcsec",
                              "CTYPE2": "HPLT-TAN", "HGLT_OBS": hglt_obs,
                              "HGLN_OBS": hgln_obs,
                              "RSUN_OBS": rsun_obs,
                              "RSUN_REF": rsun_ref,
                              "DSUN_OBS": dsun_obs}
                    rmap = smap.Map(data, header)
                # resample the image for plotting
                if fov is not None:
                    fov = [np.array(ll) for ll in fov]
                    try:
                        pad = max(np.diff(fov[0])[0], np.diff(fov[1])[0])
                        rmap = rmap.submap((fov[0] + np.array([-1.0, 1.0]) * pad) * u.arcsec,
                                           (fov[1] + np.array([-1.0, 1.0]) * pad) * u.arcsec)
                    except:
                        pad = max(fov[0][1] - fov[0][0], fov[1][1] - fov[1][0])
                        bl = SkyCoord((fov[0][0] - pad) * u.arcsec, (fov[1][0] - pad) * u.arcsec,
                                      frame=rmap.coordinate_frame)
                        tr = SkyCoord((fov[0][1] + pad) * u.arcsec, (fov[1][1] + pad) * u.arcsec,
                                      frame=rmap.coordinate_frame)
                        if sunpy3:
                            rmap = rmap.submap(bl, top_right=tr)
                        else:
                            rmap = rmap.submap(bl, tr)

                else:
                    dim = u.Quantity([256, 256], u.pixel)
                    rmap = rmap.resample(dim)
                if plotaia:
                    if npols > 1:
                        ax = axs[pidx]
                    else:
                        ax = axs[0]
                else:
                    if npols > 1:
                        ax = axs[pidx + s * 2]
                    else:
                        ax = axs[s]

                # Check the quality of the data
                # (ny, nx) = data.shape
                max_pix = np.nanmax(rmap.data)
                min_pix = np.nanmin(rmap.data)
                # data[int(ny * 0.25):int(ny * 0.75), int(nx * 0.5):int(nx * 0.75)] = np.nan
                rms = np.nanstd(rmap.data)
                dyn_range = max_pix / rms

                rmap_ = pmX.Sunmap(rmap)
                if dyn_range > minsnr and max_pix > 2 * np.abs(min_pix):
                    rmap_flag = False
                else:
                    rmap_flag = True
                    # print(f'Bad data quality at {plttime.iso} for spw {sp} pol {pol}. dyn_range: {dyn_range:.1f}, max_pix: {max_pix:.1f}, rms: {rms:.1f}. Skip plotting.')
                    # continue

                if plotaia:
                    if aiamap:
                        if amax is None:
                            amax = np.nanmax(aiamap.data)
                        if amin is None:
                            amin = 1.0
                        if aia_jp2:
                            _anorm = get_normalization(0, 255, 'linear')
                        else:
                            _anorm = get_normalization(amin, amax, anorm)
                        if acmap is None:
                            acmap = 'gray_r'

                        if nspw > 1:
                            # print(f'adding radio images at s: {s}, sp: {sp}: spwpltCounts: {spwpltCounts}')
                            if spwpltCounts == 0:
                                aiamap_ = pmX.Sunmap(aiamap)
                                aiamap_.imshow(axes=ax, cmap=acmap,
                                               norm=_anorm,
                                               interpolation='nearest')
                                # print(
                                #     f'radio image at sp:{sp} pol:{pol} at {plttime.iso} aiamap.data.max: {np.nanmax(aiamap.data)}')
                        else:
                            aiamap_ = pmX.Sunmap(aiamap)
                            aiamap_.imshow(axes=ax, cmap=acmap,
                                           norm=_anorm,
                                           interpolation='nearest')
                    else:
                        rmap_blank = sunpy.map.Map(np.full_like(rmap.data, np.nan), rmap.meta)
                        rmap_blank_ = pmX.Sunmap(rmap_blank)
                        if nspw > 1:
                            if spwpltCounts == 0:
                                rmap_blank_.imshow(axes=ax)
                        else:
                            rmap_blank_.imshow(axes=ax)

                    if not rmap_flag:
                        try:
                            clevels1 = np.array(clevelsfix[s])  # firstly try contour with clevelsfix if given
                        except:
                            try:
                                clevels1 = np.linspace(iranges[pidx][0], iranges[pidx][1], nclevels)
                            except:
                                try:
                                    clevels1 = np.array(clevels) * np.nanmax(rmap.data)
                                except:
                                    clevels1 = np.linspace(0.5, 1.0, 2) * np.nanmax(rmap.data)
                        if np.any(clevels1):
                            if nspw > 1:
                                if opencontour:
                                    rmap_.contour(axes=ax, levels=clevels1,
                                                  colors=[colors_spws[s]] * len(clevels1),
                                                  alpha=alpha_cont, linewidths=2)
                                else:
                                    if clevels1[0]<np.nanmax(rmap.data):
                                        rmap_.contourf(axes=ax, levels=[clevels1[0], np.nanmax(rmap.data)],
                                                       colors=[colors_spws[s]] * 2,
                                                       alpha=alpha_cont)
                            else:
                                rmap_.contour(axes=ax, levels=clevels1, cmap=icmap)
                else:
                    if rmap_flag:
                        rmap_blank = sunpy.map.Map(np.full_like(rmap.data, np.nan), rmap.meta)
                        rmap_blank_ = pmX.Sunmap(rmap_blank)
                        if nspw > 1:
                            if spwpltCounts == 0:
                                rmap_blank_.imshow(axes=ax)
                        else:
                            rmap_blank_.imshow(axes=ax)
                    else:
                        _inorm = get_normalization(iranges[pidx][0], iranges[pidx][1], inorm)
                        rmap_.imshow(axes=ax, norm=_inorm, cmap=cmaps[pol],
                                     interpolation='nearest')
                        rmap_.draw_limb(axes=ax)
                        rmap_.draw_grid(axes=ax)
                if custom_mapcubes:
                    for cmpcidx, cmpc in enumerate(custom_mapcubes['mapcube']):
                        dtcmpc = np.mean(np.diff(cmpc_plttimes_mjd[cmpcidx]))
                        timeline = cmpc_plttimes_mjd[cmpcidx] - Time(plttime).mjd
                        if np.min(np.abs(timeline)) <= dtcmpc:
                            if 'levels' in custom_mapcubes.keys():
                                levels = np.array(custom_mapcubes['levels'][cmpcidx])
                            else:
                                levels = np.linspace(0.2, 0.9, 3)
                            if 'color' in custom_mapcubes.keys():
                                color = custom_mapcubes['color'][cmpcidx]
                            else:
                                color = None
                            cmpidx = np.argmin(np.abs(timeline))
                            cmp = cmpc[cmpidx]
                            if 'label' in custom_mapcubes.keys():
                                label = custom_mapcubes['label'][cmpcidx]
                            else:
                                label = '-'.join(['{:.0f}'.format(ll) for ll in cmp.measurement.value]) + ' {}'.format(
                                    cmp.measurement.unit)
                            cmp_ = pmX.Sunmap(cmp)
                            cmp_.contour(axes=ax, levels=np.array(levels) * np.nanmax(cmp.data), colors=color)
                            ax.text(0.97, (len(custom_mapcubes['mapcube']) - cmpcidx - 1) * 0.06 + 0.03, label,
                                    horizontalalignment='right',
                                    verticalalignment='bottom', transform=ax.transAxes, color=color)
                # ax.set_autoscale_on(True)
                if fov:
                    ax.set_xlim(fov[0])
                    ax.set_ylim(fov[1])
                else:
                    # ax.set_xlim([-1220, 1220])
                    # ax.set_ylim([-1220, 1220])
                    ax.set_xlim(rmap_.xrange.value)
                    ax.set_ylim(rmap_.yrange.value)
                if spwpltCounts == 0 and pidx == 0:
                    timetext = ax.text(0.99, 0.98, '', color='w', fontweight='bold', fontsize=9, ha='right', va='top',
                                       transform=ax.transAxes)
                    timetext.set_text(plttime.iso[:19])
                if nspw <= 1 or plotaia == False:
                    try:
                        ax.text(0.98, 0.01, '{1} @ {0:.1f} GHz'.format(cfreqs[s], pol),
                                color='w', transform=ax.transAxes, fontweight='bold', ha='right')
                    except:
                        ax.text(0.98, 0.01, '{1} @ {0:.1f} GHz'.format(0., pol), color='w',
                                transform=ax.transAxes, fontweight='bold',
                                ha='right')

                ax.set_title(' ')
                # ax.xaxis.set_visible(False)
                # ax.yaxis.set_visible(False)
            spwpltCounts += 1
        if i == 0 and plotaia:
            if nspw > 1:
                import matplotlib.colorbar as colorbar
                ticks, bounds, fmax, fmin, freqmask = get_colorbar_params(freqbounds)

                for pidx in range(npols):
                    ax = axs[pidx]
                    divider = make_axes_locatable(ax)
                    cax_freq = divider.append_axes('right', size='6.0%', pad=0.1)
                    # cax_freq.tick_params(direction='out')
                    cb = colorbar.ColorbarBase(cax_freq, norm=colors.Normalize(vmin=fmin, vmax=fmax), cmap=icmap,
                                               orientation='vertical', boundaries=bounds, spacing='proportional',
                                               ticks=ticks, format='%4.1f', alpha=alpha_cont)
                    # Freqs = [np.mean(fq) for fq in Freq]
                    # mpl.colorbar.ColorbarBase(cax_freq, cmap=icmap, norm=colors.Normalize(vmax=Freqs[-1], vmin=Freqs[0]))
                    for fbd_lo, fbd_hi in freqmask:
                        if fbd_hi is not None:
                            cax_freq.axhspan(fbd_lo, fbd_hi, hatch='//', edgecolor='k', facecolor='#BBBBBB')
                    cax_freq.set_ylabel('Frequency [GHz]')
                    cax_freq.tick_params(axis="y", pad=-20., length=0, colors='k', labelsize=8)
                    cax_freq.axhline(fmin, xmin=1.0, xmax=1.2, color='k', clip_on=False)
                    cax_freq.axhline(fmax, xmin=1.0, xmax=1.2, color='k', clip_on=False)
                    cax_freq.text(1.25, 0.0, '{:.1f}'.format(fmin), fontsize=9, transform=cax_freq.transAxes,
                                  va='center',
                                  ha='left')
                    cax_freq.text(1.25, 1.0, '{:.1f}'.format(fmax), fontsize=9, transform=cax_freq.transAxes,
                                  va='center',
                                  ha='left')
        # figname = observatory + '_qlimg_' + plttime.isot.replace(':', '').replace('-', '')[:19] + '.png'
        # fig_tdt = plttime.to_datetime())
        # fig_subdir = fig_tdt.strftime("%Y/%m/%d/")
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        if verbose:
            print('Saving plot to: ' + fignameful)
        if i == 0:
            gs.tight_layout(fig, rect=[0.08, 0, 0.98, 1.0])
        fig.savefig(fignameful)
    plt.close(fig)
    if not moviename:
        moviename = 'movie'
    if movieformat.lower() == 'html':
        DButil.img2html_movie(figdir, outname=moviename)
    else:
        DButil.img2movie(figdir, outname=moviename)
    outmovie = figdir + moviename + '.' + movieformat
    return outmovie


def dspec_external(vis, workdir='./', specfile=None, ds_normalised=False):
    dspecscript = os.path.join(workdir, 'dspec.py')
    if not specfile:
        specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
    os.system('rm -rf {}'.format(dspecscript))
    fi = open(dspecscript, 'wb')
    fi.write('from suncasa.dspec import dspec as ds \n')
    if ds_normalised == False:
        fi.write(
            'specdata = ds.Dspec("{0}", specfile="{1}", domedian=True, verbose=True, savespec=True, usetbtool=True) \n'.format(
                vis,
                specfile))
    else:
        fi.write(
            'specdata = ds.Dspec("{0}", specfile="{1}", domedian=True, verbose=True, savespec=True, usetbtool=True, ds_normalised=True) \n'.format(
                vis,
                specfile))
    fi.close()
    os.system('casa --nologger -c {}'.format(dspecscript))


def qlookplot(vis, timerange=None, spw='', spwplt=None,
              workdir='./', specfile=None,
              xycen=None, fov=[500., 500.], xyrange=None,
              restoringbeam=[''],
              refbmsize=70., reffreq=1.0, minbmsize=4.0,
              antenna='', uvrange='', stokes='RR,LL',
              robust=0.0,
              weighting='briggs', niter=500,
              imsize=[512], cell=['5.0arcsec'], mask='', gain=0.1, pbcor=True,
              interactive=False, datacolumn='data',
              reftime='',
              toTb=True, sclfactor=1.0, subregion='',
              usemsphacenter=True,
              imagefile=None, outfits='',
              docompress=True,
              wrapfits=True,
              nclevels=3,
              clevels=None, clevelsfix=None, calpha=0.5, opencontour=False,
              imax=None, imin=None, icmap=None, inorm='linear',
              dmin=None, dmax=None, dcmap=None, dnorm='linear',
              plotaia=True, aiawave=171, aiafits=None, aiadir=None,
              amax=None, amin=None, acmap=None, anorm='log',
              goestime=None,
              mkmovie=False, ncpu=1, twidth=1, movieformat='html',
              cleartmpfits=True, overwrite=True,
              clearmshistory=False, show_warnings=False, verbose=False, quiet=False, ds_normalised=False):
    '''
    Generate quick-look plots and dynamic spectra for solar radio observations.
    Required inputs:
               :param vis: Path to the calibrated CASA measurement set.
    Important optional inputs:
            timerange: Timerange for analysis in standard CASA format. Defaults to entire range, which can be slow.
            spw: spectral window (SPW) selection following the CASA syntax.
                 Examples: spw='1:2~60' (spw id 1, channel range 2-60); spw='*:1.2~1.3GHz' (selects all channels within 1.2-1.3 GHz; note the *)
                 spw can be a list of spectral windows, i.e, ['0', '1', '2', '3', '4', '5', '6', '7']
            spwplt: Subset of SPW to display, defaults to all specified in `spw`.
            workdir: Working directory for temporary files, defaults to current directory.
            specfile: Path to a saved dynamic spectrum file (from suncasa.dspec.dspec.Dspec())
                        6or generate a median dynamic spectrum on the fly if not provided.
    Optional inputs:
            goestime: goes plot time, example ['2016/02/18 18:00:00','2016/02/18 23:00:00']
            xycen: center of the image in helioprojective coordinates (HPLN/HPLT), in arcseconds. Example: [900, -150.]
            fov: field of view in arcsecs. Example: [500., 500.]
            xyrange: field of view in solar XY coordinates. Format: [[x1,x2],[y1,y2]]. Example: [[900., 1200.],[0,300]]
                     ***NOTE: THIS PARAMETER OVERWRITES XYCEN AND FOV***

            Restoring Beam Parameters:
                restoringbeam: A list specifying the sizes of the restoring beam explicitly in arc seconds (e.g., ['100arcsec', '80arcsec', '50arcsec', ...]).
                                Must match the number of SPWs if not ['']. If specified, these values override automatic beam size calculations for each SPW.

                refbmsize:      The reference beam size in arc seconds. This parameter is used in conjunction with `reffreq` to calculate the beam size for all SPWs,
                                assuming the beam size is inversely proportional to the frequency.
                                The parameters `refbmsize`,`reffreq` and `minbmsize` are only used if `restoringbeam` is set to [''].

                reffreq:        The reference frequency in GHz, used together with `refbmsize` to calculate the beam sizes for SPWs.

                minbmsize:      Minimum beam size in arcseconds, overrides smaller calculated sizes.

            CASA tclean parameters: refer to CASA tclean documentation for more details.
                antenna: baseline to generate dynamic spectrum
                uvrange: uvrange to select baselines for generating dynamic spectrum
                stokes: polarization of the clean image, can be 'RR,LL' or 'I,V'
                robust:
                weighting:
                niter:
                imsize:
                cell:
                mask: only accept CASA region format (https://casaguides.nrao.edu/index.php/CASA_Region_Format)
                gain:
                pbcor:
                interactive:
                datacolumn:
            image registration parameters:
                reftime: Reference time for image alignment.
                toTb: Bool. Convert the default Jy/beam to brightness temperature?
                sclfactor: scale the image values by its value (e.g., sclfactor = 100 to compensate VLA 20 dB attenuator)
                subregion: only write the data within the sub-region selection. See 'help par.region' for details.
                usephacenter: Bool -- if True, correct for the RA and DEC in the ms file based on solar empheris.
                     Otherwise assume the phasecenter is correctly pointed to the solar disk center
                     (EOVSA case)
                imagefile: Use specified CASA radio image file for registration; otherwise, generate anew.
                outfits: Use specified FITS file of a radio image for output; otherwise, generate anew.
                docompress: if compress the outfits
                wrapfits: if wrap the fits files of multiple spectral windows at one given time interval into a combined fits file.
                overwrite: if overwrite the existed outfits file (default: True).

            radio image plotting parameters:
                nclevels: Number of contour levels for radio image plots.
                clevels: Specific contour levels for radio image plots.
                opencontour: Boolean. Plots open contours if True; filled contours otherwise.
                icmap: Color map (string or Colormap object) for radio images/contours.
                imax, imin: Color scale range, defining normalization before color mapping.
                inorm: Normalization method (string or Normalize object), overriding imax and imin.

            radio dynamic spectrum plotting parameters:
                dcmap: Color map (string or Colormap object) for the dynamic spectrum.
                dmin, dmax: Color scale range for dynamic spectrum normalization before color mapping.
                dnorm: Normalization method (string or Normalize object), overriding dmax and dmin.

            SDO/AIA image plotting parameters:
                plotaia: Boolean. Downloads and plots AIA image at specified aiawave if True.
                aiawave: AIA image passband to download and display.
                aiafits: Directly plots AIA image from provided FITS file, skipping download. (note: users can provide any solar image FITS file for plotting).
                aiadir: Searches this directory for AIA image files to skip download.
                acmap: Color map (string or Colormap object) for AIA images.
                amin, amax: Color scale range for AIA image normalization before color mapping.
                anorm: Normalization method (string or Normalize object), overriding amax and amin.
            movie parameters:
                mkmovie: Boolean. Generates a movie from radio images over multiple time intervals if True.
                ncpu: Number of CPUs for parallel clean operations with ptclean.
                twidth: Time pixel averaging width (default: 1).
                movieformat: Output movie format, either 'html' or 'mp4'.
    '''

    outfits_list = None
    outmovie = None
    from importlib import reload
    reload(mstools)
    if not show_warnings:
        import warnings
        warnings.filterwarnings("ignore")

    if aiadir == None:
        aiadir = './'
    if xycen:
        xc, yc = xycen
        if len(fov) == 1:
            fov = fov * 2
        xlen, ylen = fov
        # if parse_version(sunpy.__version__) > parse_version('0.8.0'):
        #     xyrange = [[xc - xlen / 2.0, yc - ylen / 2.0], [xc + xlen / 2.0, yc + ylen / 2.0]]
        # else:
        xyrange = [[xc - xlen / 2.0, xc + xlen / 2.0], [yc - ylen / 2.0, yc + ylen / 2.0]]
    stokes_allowed = ['RR,LL', 'I,V', 'RRLL', 'IV', 'XXYY', 'XX,YY', 'RR', 'LL', 'I', 'V', 'XX', 'YY']
    if not stokes in stokes_allowed:
        print('Error: wrong stokes parameter ' + str(stokes) + '. Allowed values are ' + ';  '.join(stokes_allowed))
        return -1
    if stokes == 'RRLL':
        stokes = 'RR,LL'
    elif stokes == 'XXYY':
        stokes = 'XX,YY'
    elif stokes == 'IV':
        stokes = 'I,V'

    if dcmap is None:
        dcmap = plt.get_cmap('afmhot')

    polmap = {'RR': 0, 'LL': 1, 'I': 0, 'V': 1, 'XX': 0, 'YY': 1}
    pols = stokes.split(',')
    npol_in = len(pols)

    if vis[-1] == '/':
        vis = vis[:-1]
    if not os.path.exists(vis):
        print('Error: input measurement not exist')
        return -1
    if clearmshistory:
        ms_clearhistory(vis)
    if aiafits is None:
        aiafits = ''
    # split the data
    # generating dynamic spectrum
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    os.chdir(workdir)
    workdir = './'
    if specfile:
        try:
            spec_inst = ds.Dspec()
            try:
                spec_inst.read(specfile)
            except:
                sources = ["eovsa", "suncasa", "lwa", "general"]
                for source in sources:
                    try:
                        spec_inst.read(specfile, source=source)
                        break
                    except:
                        pass
        except:
            print('Provided dynamic spectrum file not numpy npz. Generating one from the visibility data')
            specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
            if c_external:
                dspec_external(vis, workdir=workdir, specfile=specfile)
                spec_inst = ds.Dspec(specfile)
            else:
                spec_inst = ds.Dspec(vis, specfile=specfile, domedian=True, verbose=True, usetbtool=True)

    else:
        print('Dynamic spectrum file not provided; Generating one from the visibility data')
        # specdata = ds.get_dspec(vis, domedian=True, verbose=True)
        specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
        if c_external:
            dspec_external(vis, workdir=workdir, specfile=specfile)
            spec_inst = ds.Dspec(specfile)
        else:
            spec_inst = ds.Dspec(vis, specfile=specfile, domedian=True, verbose=True, usetbtool=True)

    try:
        tb.open(vis + '/POINTING')
        starttim = Time(tb.getcell('TIME_ORIGIN', 0) / 24. / 3600., format='mjd')
        endtim = Time(tb.getcell('TIME_ORIGIN', tb.nrows() - 1) / 24. / 3600., format='mjd')
    except:
        tb.open(vis)
        starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
        endtim = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
    tb.close()
    datstr = starttim.iso[:10]

    if timerange is None or timerange == '':
        starttim1 = starttim
        endtim1 = endtim
        timerange = '{0}~{1}'.format(starttim.iso.replace('-', '/').replace(' ', '/'),
                                     endtim.iso.replace('-', '/').replace(' ', '/'))
    else:
        try:
            (tstart, tend) = timerange.split('~')
            if tstart[2] == ':':
                starttim1 = Time(datstr + 'T' + tstart)
                endtim1 = Time(datstr + 'T' + tend)
                timerange = '{0}/{1}~{0}/{2}'.format(datstr.replace('-', '/'), tstart, tend)
            else:
                starttim1 = Time(qa.quantity(tstart, 'd')['value'], format='mjd')
                endtim1 = Time(qa.quantity(tend, 'd')['value'], format='mjd')
        except ValueError:
            print("keyword 'timerange' in wrong format")

    midtime_mjd = (starttim1.mjd + endtim1.mjd) / 2.

    if vis.endswith('/'):
        vis = vis[:-1]
    visname = os.path.basename(vis)
    bt = starttim1.plot_date
    et = endtim1.plot_date

    # find out min and max frequency for plotting in dynamic spectrum
    ms.open(vis)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    spwInfo = ms.getspectralwindowinfo()
    nspwall = len(spwInfo)

    if not spw:
        if observatory == 'EOVSA':
            # if nspwall == 31:
            #     spw = list((np.arange(30) + 1).astype(str))
            #     spwselec = '1~' + str(30)
            #     spw = [str(sp) for sp in spw]
            # else:
            spw = list(np.arange(nspwall).astype(str))
            spwselec = '0~' + str(nspwall - 1)
            spw = [str(sp) for sp in spw]
        else:
            spwselec = '0~' + str(nspwall - 1)
            spw = [spwselec]
    else:
        if type(spw) is list:
            spwselec = ';'.join(spw)
        else:
            spwselec = spw
            if ';' in spw:
                spw = spw.split(';')
            else:
                spw = [spw]  # spw=spw.split(';')

    nspws = len(spw)

    if icmap is None:
        if plotaia:
            if nspws > 1:
                icmap = plt.get_cmap('RdYlBu')
            else:
                icmap = plt.get_cmap('gist_heat')
        else:
            icmap = plt.get_cmap('gist_heat')
    else:
        icmap = plt.get_cmap(icmap)

    print('spw is', spw)
    bdinfo = mstools.get_bandinfo(vis, spw=spw, returnbdinfo=True)
    # print(freqbounds)
    cfreqs = bdinfo['cfreqs']
    cfreqs_all = bdinfo['cfreqs_all']
    freq_dist = lambda fq: (fq - cfreqs_all[0]) / (cfreqs_all[-1] - cfreqs_all[0])
    staql = {'timerange': timerange, 'spw': spwselec}
    if ms.msselect(staql, onlyparse=True):
        ndx = ms.msselectedindices()
        chan_sel = ndx['channel']
        bspw = chan_sel[0, 0]
        bchan = chan_sel[0, 1]
        espw = chan_sel[-1, 0]
        echan = chan_sel[-1, 2]
        bfreq = spwInfo[str(bspw)]['Chan1Freq'] + spwInfo[str(bspw)]['ChanWidth'] * bchan
        efreq = spwInfo[str(espw)]['Chan1Freq'] + spwInfo[str(espw)]['ChanWidth'] * echan
        bfreqghz = bfreq / 1e9
        efreqghz = efreq / 1e9
        if verbose:
            print('selected timerange {}'.format(timerange))
            print('selected frequency range {0:6.3f} to {1:6.3f} GHz'.format(bfreqghz, efreqghz))
    else:
        print("spw or timerange selection failed. Aborting...")
        ms.close()
        return -1
    ms.close()

    if observatory == 'EOVSA':
        if npol_in == 2:
            if stokes == 'RRLL' or stokes == 'RR,LL':
                print('Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be XXYY')
                stokes = 'XX,YY'
        else:
            if stokes == 'RR':
                print('Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be XX')
                stokes = 'XX'
            elif stokes == 'LL':
                print('Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be YY')
                stokes = 'YY'
            else:
                pass

    if mkmovie:
        plt.ioff()
        # fig = plt.figure(figsize=(12, 7.5), dpi=100)
        if outfits:
            pass
        else:
            eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
            if observatory == 'EOVSA' or (not usemsphacenter):
                print('This is EOVSA data')
                # use RA and DEC from FIELD ID 0
                tb.open(vis + '/FIELD')
                phadir = tb.getcol('PHASE_DIR').flatten()
                tb.close()
                ra0 = phadir[0]
                dec0 = phadir[1]
            else:
                ra0 = eph['ra'][0]
                dec0 = eph['dec'][0]

            if not xycen:
                # use solar disk center as default
                phasecenter = 'J2000 ' + str(ra0) + 'rad ' + str(dec0) + 'rad'
            else:
                x0 = np.radians(xycen[0] / 3600.)
                y0 = np.radians(xycen[1] / 3600.)
                p0 = np.radians(eph['p0'][0])  # p angle in radians
                raoff = -((x0) * np.cos(p0) - y0 * np.sin(p0)) / np.cos(eph['dec'][0])
                decoff = (x0) * np.sin(p0) + y0 * np.cos(p0)
                newra = ra0 + raoff
                newdec = dec0 + decoff
                phasecenter = 'J2000 ' + str(newra) + 'rad ' + str(newdec) + 'rad'
            print('use phasecenter: ' + phasecenter)
            qlookfitsdir = os.path.join(workdir, 'qlookfits/')
            # qlookallbdfitsdir = os.path.join(workdir, 'qlookallbdfits/')
            qlookfigdir = os.path.join(workdir, 'qlookimgs/')
            imresfile = os.path.join(qlookfitsdir, '{}.imres.npz'.format(os.path.basename(vis)))
            if overwrite:
                imres, outfits_list = mk_qlook_image(vis, timerange=timerange, spws=spw, twidth=twidth, ncpu=ncpu,
                                                     imagedir=qlookfitsdir, phasecenter=phasecenter, stokes=stokes,
                                                     mask=mask,
                                                     uvrange=uvrange, robust=robust, niter=niter, gain=gain,
                                                     imsize=imsize, cell=cell,
                                                     pbcor=pbcor,
                                                     reftime=reftime, restoringbeam=restoringbeam, sclfactor=sclfactor,
                                                     docompress=docompress,
                                                     wrapfits=wrapfits,
                                                     overwrite=overwrite,
                                                     c_external=c_external,
                                                     subregion=subregion,
                                                     show_warnings=show_warnings)
            else:
                if os.path.exists(imresfile):
                    imres = np.load(imresfile, allow_pickle=True)
                    imres = imres['imres'].item()
                else:
                    print('Image results file not found; Creating new images.')
                    imres, outfits_list = mk_qlook_image(vis, timerange=timerange, spws=spw, twidth=twidth, ncpu=ncpu,
                                                         imagedir=qlookfitsdir, phasecenter=phasecenter, stokes=stokes,
                                                         mask=mask,
                                                         uvrange=uvrange, robust=robust, niter=niter, gain=gain,
                                                         imsize=imsize,
                                                         cell=cell, pbcor=pbcor,
                                                         reftime=reftime, restoringbeam=restoringbeam,
                                                         sclfactor=sclfactor,
                                                         docompress=docompress,
                                                         wrapfits=wrapfits,
                                                         overwrite=overwrite,
                                                         c_external=c_external,
                                                         subregion=subregion,
                                                         show_warnings=show_warnings)
            if not os.path.exists(qlookfigdir):
                os.makedirs(qlookfigdir)

            # clevelsfix = radio_image_clevels(imres, snr=2)
            outmovie = plt_qlook_image(imres, timerange=timerange, spwplt=spwplt, figdir=qlookfigdir,
                                       specdata=spec_inst,
                                       verbose=verbose,
                                       stokes=stokes, fov=xyrange,
                                       amax=amax, amin=amin, acmap=acmap, anorm=anorm,
                                       imax=imax, imin=imin, icmap=icmap, inorm=inorm,
                                       dmax=dmax, dmin=dmin, dcmap=dcmap, dnorm=dnorm,
                                       sclfactor=sclfactor,
                                       aiafits=aiafits, aiawave=aiawave, aiadir=aiadir, plotaia=plotaia,
                                       freqbounds=bdinfo, alpha_cont=calpha,
                                       opencontour=opencontour,
                                       nclevels=nclevels,
                                       # clevelsfix = clevelsfix,
                                       movieformat=movieformat, ds_normalised=ds_normalised)

    else:
        if np.iscomplexobj(spec_inst.data):
            spec = np.abs(spec_inst.data)
        else:
            spec = spec_inst.data
        spec = spec * sclfactor
        spec = checkspecnan(spec)
        if spec.ndim == 2:
            (nfreq, ntim) = spec.shape
            nbl = 1
            npol_fits = 2
            spec_ = np.zeros((npol_fits, nbl, nfreq, ntim))
            spec_[0, 0, :, :] = spec
            spec_[0, 0, :, :] = spec
            spec = spec_
        else:
            (npol_fits, nbl, nfreq, ntim) = spec.shape
        fidx = range(nfreq)
        freq = spec_inst.freq_axis
        freqghz = freq / 1e9
        spec_tim = spec_inst.time_axis
        spec_tim_plt = spec_tim.plot_date
        plt.ion()
        if not quiet:
            # fig = plt.figure(figsize=(11.65, 8.74), dpi=100)
            fig = plt.figure(figsize=(11.80, 8.80), dpi=80)
            ax1 = plt.subplot2grid((6, 8), (0, 0), rowspan=2, colspan=2)
            ax2 = plt.subplot2grid((6, 8), (2, 0), rowspan=2, colspan=2, sharex=ax1, sharey=ax1)
            ax3 = plt.subplot2grid((6, 8), (4, 0), rowspan=2, colspan=2)
            ax4 = plt.subplot2grid((6, 8), (0, 2), rowspan=3, colspan=3)
            ax5 = plt.subplot2grid((6, 8), (3, 2), rowspan=3, colspan=3)
            ax6 = plt.subplot2grid((6, 8), (0, 5), rowspan=3, colspan=3, sharex=ax4, sharey=ax4)
            ax7 = plt.subplot2grid((6, 8), (3, 5), rowspan=3, colspan=3, sharex=ax5, sharey=ax5)

            specs = {}
            if npol_in > 1:
                if npol_fits > 1:
                    if stokes == 'I,V':
                        specs['I'] = (np.absolute(spec[0, 0, :, :]) + np.absolute(spec[1, 0, :, :])) / 2.0
                        specs['V'] = (np.absolute(spec[0, 0, :, :]) - np.absolute(spec[1, 0, :, :])) / 2.0
                    else:
                        specs[pols[0]] = np.absolute(spec[0, 0, :, :])
                        specs[pols[1]] = np.absolute(spec[1, 0, :, :])
                else:
                    warnings.warn(
                        "The provided specfile only provides one polarization. The polarization of the dynamic spectrum could be wrong.")
                    specs[pols[0]] = np.absolute(spec[0, 0, :, :])
                    specs[pols[1]] = np.zeros_like(spec[0, 0, :, :])
            else:
                if npol_fits > 1:
                    if stokes == 'I':
                        specs['I'] = (np.absolute(spec[0, 0, :, :]) + np.absolute(spec[1, 0, :, :])) / 2.0
                    elif stokes == 'V':
                        specs['V'] = (np.absolute(spec[0, 0, :, :]) - np.absolute(spec[1, 0, :, :])) / 2.0
                    else:
                        specs[pols[0]] = np.absolute(spec[polmap[pols[0]], 0, :, :])
                else:
                    specs[pols[0]] = np.absolute(spec[0, 0, :, :])

            print('plot the dynamic spectrum in pol ' + ' & '.join(pols))

            if dnorm is 'linear':
                dnorm = colors.Normalize(vmax=dmax, vmin=dmin)
            elif dnorm is 'log':
                dnorm = colors.LogNorm(vmax=dmax, vmin=dmin)



            axs = [ax1, ax2]
            for axidx, ax in enumerate(axs):
                if axidx < npol_in:
                    ax.pcolormesh(spec_tim_plt, freqghz, specs[pols[axidx]], cmap=dcmap, norm=dnorm,
                                  rasterized=True)
                    ax.set_title(observatory + ' ' + datstr + ' ' + pols[axidx], fontsize=9)
                ax.set_autoscale_on(True)
                ax.add_patch(patches.Rectangle((bt, bfreqghz), et - bt, efreqghz - bfreqghz, ec='w', fill=False))
                ax.plot([(bt + et) / 2.], [(bfreqghz + efreqghz) / 2.], '*w', ms=12)
                for tick in ax.get_xticklabels():
                    tick.set_rotation(30)
                    tick.set_fontsize(8)
                ax.set_ylabel('Frequency (GHz)', fontsize=9)
                if axidx == 1:
                    ax.xaxis_date()
                    ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
                    locator = mpl.dates.AutoDateLocator()
                    ax.xaxis.set_major_locator(locator)
                    ax.set_xlim(spec_tim_plt[0], spec_tim_plt[-1])
                    ax.set_ylim(freqghz[0], freqghz[-1])

            # import pdb
            # pdb.set_trace()
            # Second part: GOES plot
            if goestime:
                btgoes = goestime[0]
                etgoes = goestime[1]
            else:
                # datstrg = datstr.replace('-', '/')
                tdur = spec_tim[-1].jd - spec_tim[0].jd
                # btgoes = datstr + ' ' + qa.time(qa.quantity(tim[0] - tdur, 's'), form='clean', prec=9)[0]
                # etgoes = datstr + ' ' + qa.time(qa.quantity(tim[-1] + tdur, 's'), form='clean', prec=9)[0]
                btgoes, etgoes = Time([spec_tim[0].jd - tdur, spec_tim[-1].jd + tdur], format='jd').iso
            if verbose:
                print('Acquire GOES soft X-ray data in from ' + btgoes + ' to ' + etgoes)

            # ax3 = plt.subplot(gs1[2])

            try:
                import socket
                socket.setdefaulttimeout(60)
                from sunpy.timeseries import TimeSeries
                from sunpy.time import TimeRange, parse_time
                from sunpy.net import Fido, attrs as a
                results = Fido.search(a.Time(TimeRange(btgoes, etgoes)), a.Instrument('XRS'))
                files = Fido.fetch(results)
                goest = TimeSeries(files)
                if isinstance(goest, list):
                    import pandas as pd
                    gdata = [g.data for g in goest]
                    gdata = pd.concat(gdata, join="inner")
                else:
                    gdata = goest.data
                try:
                    goes_dates = mpl.dates.date2num(parse_time(gdata.index))
                except:
                    goes_dates = Time(gdata.index).plot_date
                if np.abs(gdata['xrsb'].mean()) > 1e-9:
                    goesdata = gdata['xrsb']
                    goesdif = np.diff(gdata['xrsb'])
                else:
                    goes_dates, goesdata = get_goes_data(Time((spec_tim[-1].mjd + spec_tim[0].mjd) / 2.0, format='mjd'))
                    goesdif = np.diff(goesdata)

                gmax = np.nanmax(goesdif)
                gmin = np.nanmin(goesdif)
                ran = gmax - gmin
                db = 2.8 / ran
                goesdifp = goesdif * db + gmin + (-6)
                ax3.step(goes_dates, np.log10(goesdata), '-', label='1.0--8.0 $\AA$', color='red', lw=1.0)
                ax3.step(goes_dates[0:-1], goesdifp, '-', label='Derivative', color='blue', lw=0.5)

                ax3.set_ylim([-8, -3])
                ax3.set_yticks([-8, -7, -6, -5, -4, -3])
                ax3.set_yticklabels(
                    [r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$'])
                ax3.set_title('Goes Soft X-ray', fontsize=9)
                ax3.set_ylabel('Watts m$^{-2}$')
                ax3.set_xlabel(Time(spec_tim_plt[0], format='plot_date').iso[0:10])
                ax3.axvspan(spec_tim_plt[0], spec_tim_plt[-1], alpha=0.2)
                ax3.set_xlim(Time([btgoes, etgoes]).plot_date)

                for tick in ax3.get_xticklabels():
                    tick.set_fontsize(8)
                    tick.set_rotation(30)

                ax3_2 = ax3.twinx()
                # ax3_2.set_yscale("log")
                ax3_2.set_ylim([-8, -3])
                ax3_2.set_yticks([-8, -7, -6, -5, -4, -3])
                ax3_2.set_yticklabels(['A', 'B', 'C', 'M', 'X', ''])

                ax3.yaxis.grid(True, 'major')
                ax3.xaxis.grid(False, 'major')
                ax3.legend(prop={'size': 6})

                formatter = mpl.dates.DateFormatter('%H:%M')
                ax3.xaxis.set_major_formatter(formatter)
                locator = mpl.dates.AutoDateLocator()
                ax3.xaxis.set_major_locator(locator)

                ax3.fmt_xdata = mpl.dates.DateFormatter('%H:%M')
            except:
                print('Error in downloading GOES soft X-ray data. Proceeding with out soft X-ray plot.')
                ax3.set_title('Goes Soft X-ray', fontsize=9)

            # third part
            # start to download the fits files
            if plotaia:
                if acmap is None:
                    if sunpy1:
                        cmap_aia = plt.get_cmap('sdoaia{}'.format(aiawave))
                    else:
                        cmap_aia = cm_sunpy.get_cmap('sdoaia{}'.format(aiawave))
                else:
                    cmap_aia = plt.get_cmap(acmap)
                cmap_aia.set_bad(cmap_aia(0.0))
                if not aiafits:
                    try:
                        if int(aiawave) in [171, 131, 94, 335, 304, 211, 193]:
                            tdf = 6. / 24 / 3600
                        else:
                            tdf = 12. / 24 / 3600
                        newlist = trange2aiafits(Time([midtime_mjd - tdf, midtime_mjd + tdf], format='mjd'), aiawave,
                                                 aiadir)
                    except:
                        newlist = [-1]
                else:
                    newlist = [aiafits]

                try:
                    aiafits = newlist[0]
                    aiamap = smap.Map(aiafits)
                    aia_jp2 = False
                    if aiafits.endswith('.jp2'):
                        aia_jp2 = True
                    if not aia_jp2:
                        aiamap = DButil.normalize_aiamap(aiamap)
                        data = aiamap.data
                        data[data < 1.0] = 1.0
                        aiamap = smap.Map(data, aiamap.meta)
                except:
                    print('error in reading aiafits. Proceed without AIA')

        if (os.path.exists(outfits)) and (not overwrite):
            pass
        else:
            if not imagefile:
                eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
                if observatory == 'EOVSA' or (not usemsphacenter):
                    print('This is EOVSA data')
                    # use RA and DEC from FIELD ID 0
                    tb.open(vis + '/FIELD')
                    phadir = tb.getcol('PHASE_DIR').flatten()
                    tb.close()
                    ra0 = phadir[0]
                    dec0 = phadir[1]
                    if stokes == 'RRLL' or stokes == 'RR,LL':
                        print('Provide stokes: ' + str(
                            stokes) + '. However EOVSA has linear feeds. Force stokes to be IV')
                        stokes = 'I,V'
                else:
                    ra0 = eph['ra'][0]
                    dec0 = eph['dec'][0]

                if not xycen:
                    # use solar disk center as default
                    phasecenter = 'J2000 ' + str(ra0) + 'rad ' + str(dec0) + 'rad'
                else:
                    x0 = np.radians(xycen[0] / 3600.)
                    y0 = np.radians(xycen[1] / 3600.)
                    p0 = np.radians(eph['p0'][0])  # p angle in radians
                    raoff = -((x0) * np.cos(p0) - y0 * np.sin(p0)) / np.cos(eph['dec'][0])
                    decoff = (x0) * np.sin(p0) + y0 * np.cos(p0)
                    newra = ra0 + raoff
                    newdec = dec0 + decoff
                    phasecenter = 'J2000 ' + str(newra) + 'rad ' + str(newdec) + 'rad'

                if nspws > 1:
                    imagefiles, fitsfiles = [], []

                    if restoringbeam == ['']:
                        if observatory == 'EOVSA':
                            restoringbms = mstools.get_bmsize(cfreqs, refbmsize=refbmsize, reffreq=reffreq,
                                                              minbmsize=minbmsize)
                        else:
                            restoringbms = [''] * nspws
                    else:
                        try:
                            restoringbms = [float(b.replace('arcsec', '')) for b in restoringbeam]
                        except:
                            print('Error encountered while processing the provided restoring beam sizes. '
                                  'They should be specified in the format "numberarcsec" (e.g., "100arcsec").'
                                  'Falling back to using circular beams calculated as proportional to 1/freq,'
                                  'based on provided reference beam size, reference frequency,'
                                  'and minimum beam size settings.')
                            restoringbms = mstools.get_bmsize(cfreqs, refbmsize=refbmsize, reffreq=reffreq,
                                                              minbmsize=minbmsize)
                    print(f'restoringbms: {restoringbms}')
                    sto = stokes.replace(',', '')
                    print('Original phasecenter: ' + str(ra0) + str(dec0))
                    print('use phasecenter: ' + phasecenter)
                    print('do clean for ' + timerange + ' stokes ' + sto)

                    for s, sp in enumerate(tqdm(spw, desc="Processing spectral window")):
                        if restoringbms[s] == '':
                            restoringbm = ['']
                        else:
                            restoringbm = ['{:.1f}arcsec'.format(restoringbms[s])]
                        spwran = [s_.zfill(2) for s_ in sp.split('~')]
                        if len(spwran) == 2:
                            spstr = spwran[0] + '~' + spwran[1]
                        else:
                            spstr = spwran[0]
                        imagename = os.path.join(workdir, visname + '_s' + spstr + '.outim')
                        junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.pb', '.sumwt', '.image',
                                 '.image.pbcor']
                        for junk in junks:
                            if os.path.exists(imagename + junk):
                                os.system('rm -rf ' + imagename + junk + '*')
                        if verbose:
                            print('use beamsize {}'.format(restoringbm))
                        tclean(vis=vis,
                               imagename=imagename,
                               selectdata=True,
                               spw=sp,
                               timerange=timerange,
                               stokes=sto,
                               niter=niter, gain=gain,
                               antenna=antenna,
                               interactive=interactive,
                               mask=mask,
                               uvrange=uvrange,
                               pbcor=True,
                               imsize=imsize,
                               cell=cell,
                               datacolumn=datacolumn,
                               restoringbeam=restoringbm,
                               weighting=weighting,
                               robust=robust,
                               phasecenter=phasecenter)

                        if pbcor:
                            junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.image', '.pb', '.sumwt']
                            imagefile = imagename + '.image.pbcor'
                        else:
                            junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.image.pbcor', '.pb', '.sumwt']
                            imagefile = imagename + '.image'
                        for junk in junks:
                            if os.path.exists(imagename + junk):
                                os.system('rm -rf ' + imagename + junk)

                        ofits = imagefile + '.fits'
                        imagefiles.append(imagefile)
                        fitsfiles.append(ofits)
                    hf.imreg(vis=vis, imagefile=imagefiles, timerange=[timerange] * len(imagefiles),
                             fitsfile=fitsfiles, verbose=verbose, overwrite=True, sclfactor=sclfactor, toTb=toTb,
                             docompress=False)
                    # print('fits file ' + ','.join(fitsfiles) + ' selected')
                    if not outfits:
                        outfits = mstools.time2filename(vis, timerange=timerange) + '.image.fits'

                    ndfits.wrap(fitsfiles, outfitsfile=outfits, docompress=docompress)
                    if cleartmpfits:
                        for junk in imagefiles + fitsfiles:
                            os.system('rm -rf {}'.format(junk))
                    warnings.warn(
                        "If the provided spw is not equally spaced, the frequency information of the fits file {} that combining {} could be a wrong. Use it with caution!".format(
                            outfits, ','.join(fitsfiles)))
                else:
                    # single image
                    if restoringbeam == ['']:
                        if observatory == 'EOVSA':
                            # Don't use CASA's default. Trying my best to get a good estimate of the restoring beam
                            restoringbms = mstools.get_bmsize(cfreqs, refbmsize=refbmsize, reffreq=reffreq,
                                                              minbmsize=minbmsize)
                            restoringbeam = str(np.mean(restoringbms[0])) + 'arcsec'
                    else:
                        restoringbeam = restoringbeam[0]
                        restoringbeam = validate_and_reset_restoringbeam(restoringbeam)

                    imagename = os.path.join(workdir, visname + '.outim')
                    junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.pb', '.sumwt', '.image', '.image.pbcor']
                    for junk in junks:
                        if os.path.exists(imagename + junk):
                            os.system('rm -rf ' + imagename + junk + '*')
                    sto = stokes.replace(',', '')
                    print('do clean for ' + timerange + ' in spw ' + ';'.join(spw) + ' stokes ' + sto)
                    print('Original phasecenter: ' + str(ra0) + str(dec0))
                    print('use phasecenter: ' + phasecenter)
                    print('use restoringbeam {}'.format(restoringbeam))

                    tclean(vis=vis,
                           imagename=imagename,
                           selectdata=True,
                           spw=';'.join(spw),
                           timerange=timerange,
                           stokes=sto,
                           antenna=antenna,
                           niter=niter, gain=gain,
                           interactive=interactive,
                           mask=mask,
                           uvrange=uvrange,
                           pbcor=True,
                           imsize=imsize,
                           cell=cell,
                           datacolumn=datacolumn,
                           restoringbeam=restoringbeam,
                           weighting=weighting,
                           robust=robust,
                           phasecenter=phasecenter)

                    if pbcor:
                        junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.image', '.pb', '.sumwt']
                        imagefile = imagename + '.image.pbcor'
                    else:
                        junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.image.pbcor', '.pb', '.sumwt']
                        imagefile = imagename + '.image'
                    for junk in junks:
                        if os.path.exists(imagename + junk):
                            os.system('rm -rf ' + imagename + junk)
                    if not outfits:
                        outfits = mstools.time2filename(vis, timerange=timerange) + '.image.fits'
                    hf.imreg(vis=vis, imagefile=imagefile, timerange=timerange, reftime=reftime,
                             fitsfile=outfits, verbose=verbose, overwrite=True, sclfactor=sclfactor, toTb=toTb,
                             docompress=docompress)
                    print('fits file ' + outfits + ' selected')
            else:
                if not outfits:
                    outfits = mstools.time2filename(vis, timerange=timerange) + '.image.fits'
                hf.imreg(vis=vis, imagefile=imagefile, timerange=timerange, reftime=reftime,
                         fitsfile=outfits, verbose=verbose, overwrite=True, sclfactor=sclfactor, toTb=toTb,
                         docompress=docompress)
                print('fits file ' + outfits + ' selected')
        if verbose:
            print('vis', vis, 'imagefile', imagefile, 'timerange', timerange, 'reftime', reftime, 'fitsfile', outfits,
                  'verbose', verbose, 'overwrite', True, 'sclfactor', sclfactor, 'toTb', toTb, 'docompress', docompress)

        if not quiet:
            ax4.cla()
            ax5.cla()
            ax6.cla()
            ax7.cla()

            rfits = outfits
            # if nspws>1:
            #     pass
            # else:
            if isinstance(rfits, list):
                rfits = rfits[0]

            meta, rdata = ndfits.read(rfits)
            rmap = smap.Map(np.squeeze(rdata), meta['header'])
            if rmap is None:
                print('radio fits file not recognized by sunpy.map. Aborting...')
                return -1

            cmaps, datas = parse_rdata(rdata, meta, icmap=icmap,
                                       stokes=stokes)

            if not xyrange:
                if xycen:
                    x0 = xycen[0] * u.arcsec
                    y0 = xycen[1] * u.arcsec
                if not xycen:
                    row, col = rmap.data.shape
                    positon = np.nanargmax(rmap.data)
                    m, n = divmod(positon, col)
                    if sunpy1:
                        x0 = rmap.bottom_left_coord.Tx + rmap.scale[1] * (n + 0.5) * u.pix
                        y0 = rmap.bottom_left_coord.Ty + rmap.scale[0] * (m + 0.5) * u.pix
                    else:
                        x0 = rmap.xrange[0] + rmap.scale[1] * (n + 0.5) * u.pix
                        y0 = rmap.yrange[0] + rmap.scale[0] * (m + 0.5) * u.pix
                if len(fov) == 1:
                    fov = [fov] * 2
                sz_x = fov[0] * u.arcsec
                sz_y = fov[1] * u.arcsec
                x1 = x0 - sz_x / 2.
                x2 = x0 + sz_x / 2.
                y1 = y0 - sz_y / 2.
                y2 = y0 + sz_y / 2.
                xyrange = [[x1.to(u.arcsec).value, x2.to(u.arcsec).value],
                           [y1.to(u.arcsec).value, y2.to(u.arcsec).value]]
            else:
                sz_x = (xyrange[0][1] - xyrange[0][0]) * u.arcsec
                sz_y = (xyrange[1][1] - xyrange[1][0]) * u.arcsec

            clvls = {}
            if nspws < 2:
                for pol in pols:
                    if pol == 'V':
                        clvls[pol] = np.array([0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8])
                    else:
                        if clevels is None:
                            clvls[pol] = np.linspace(0.2, 0.9, 5)
                        else:
                            clvls[pol] = np.array(clevels)
            else:
                for pol in pols:
                    if pol == 'V':
                        clvls[pol] = np.array([0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8])
                    else:
                        if clevels is None:
                            clvls[pol] = np.linspace(0.3, 1, 2)
                        else:
                            clvls[pol] = np.array(clevels)

            if 'aiamap' in vars():
                if amax is None:
                    amax = np.nanmax(aiamap.data)
                if amin is None:
                    amin = 1.0

                if aia_jp2:
                    _anorm = get_normalization(0, 255, 'linear')
                else:
                    _anorm = get_normalization(amin, amax, anorm)

                title0 = 'AIA {0:.0f} Å'.format(aiamap.wavelength.value)
                aiamap_ = pmX.Sunmap(aiamap)

                axs = [ax4, ax6]
                aiamap_.draw_limb(axes=axs)
                aiamap_.draw_grid(axes=axs)
                aiamap_.imshow(axes=axs, cmap=cmap_aia, norm=_anorm, interpolation='nearest')
                for axidx, ax in enumerate(axs):
                    ax.set_title(title0, fontsize=9)
                    rect = mpl.patches.Rectangle((xyrange[0][0], xyrange[1][0]), sz_x.value, sz_y.value, edgecolor='w',
                                                 facecolor='none')
                    ax.add_patch(rect)

                axs = [ax5, ax7]
                aiamap_.draw_limb(axes=axs)
                aiamap_.draw_grid(axes=axs)
                aiamap_.imshow(axes=axs, cmap=cmap_aia, norm=_anorm, interpolation='nearest')

                axs = [[ax4, ax5], [ax6, ax7]]
                draw_limb_grid_flag = True
                for s, sp in enumerate(spw):
                    if spwplt is not None:
                        if sp not in spwplt:
                            continue

                    for pidx, pol in enumerate(pols):
                        rcmap = [icmap(freq_dist(cfreqs[s]))] * len(clvls[pol])
                        if meta['naxis'] > 2:
                            rmap_plt = smap.Map(np.squeeze(datas[pol][s, :, :]), meta['header'])
                        else:
                            rmap_plt = smap.Map(np.squeeze(datas[pol]), meta['header'])
                        rmap_plt_ = pmX.Sunmap(rmap_plt)
                        rmap_data_max = np.nanmax(rmap_plt.data)
                        if rmap_data_max <= 0.0 or rmap_data_max is np.nan:
                            print(f'Warning: the max value of the map is {rmap_data_max}. Skip plotting this map.')
                            continue
                        if nspws > 1:
                            if opencontour:
                                rmap_plt_.contour(axes=[axs[pidx][0], axs[pidx][1]], colors=rcmap,
                                                  levels=clvls[pol][:1] * np.nanmax(rmap_plt.data), alpha=calpha)
                            else:
                                rmap_plt_.contourf(axes=[axs[pidx][0], axs[pidx][1]], colors=rcmap,
                                                   levels=clvls[pol] * np.nanmax(rmap_plt.data), alpha=calpha)
                        else:
                            rmap_plt_.contour(axes=[axs[pidx][0], axs[pidx][1]], cmap=cmaps[pol],
                                              levels=clvls[pol] * np.nanmax(rmap_plt.data), alpha=calpha)
                        if draw_limb_grid_flag:
                            rmap_plt_.draw_limb(axes=[axs[pidx][0], axs[pidx][1]])
                            rmap_plt_.draw_grid(axes=[axs[pidx][0], axs[pidx][1]])
                            draw_limb_grid_flag = False
                            if nspws < 2:
                                title = title0 + ' + {0} {1:6.3f} GHz'.format(observatory, (bfreqghz + efreqghz) / 2.0)
                            else:
                                title = title0 + ' + {0} multi spws'.format(observatory)
                            axs[pidx][0].set_title(title + ' ' + pols[pidx], fontsize=9)
                            rect = mpl.patches.Rectangle((xyrange[0][0], xyrange[1][0]), sz_x.value, sz_y.value,
                                                         edgecolor='w',
                                                         facecolor='none')
                            axs[pidx][0].add_patch(rect)

                ax4.text(0.02, 0.02, 'AIA {0:.0f} '.format(aiamap.wavelength.value) + aiamap.date.strftime('%H:%M:%S'),
                         verticalalignment='bottom',
                         horizontalalignment='left', transform=ax4.transAxes, color='w', fontsize=9)
                ax6.text(0.02, 0.02, 'AIA {0:.0f} '.format(aiamap.wavelength.value) + aiamap.date.strftime('%H:%M:%S'),
                         verticalalignment='bottom',
                         horizontalalignment='left', transform=ax6.transAxes, color='w', fontsize=9)
            else:
                axs = [[ax4, ax5], [ax6, ax7]]
                if nspws < 2:
                    title = '{0} {1:6.3f} GHz'.format(observatory, (bfreqghz + efreqghz) / 2.0)
                    for pidx, pol in enumerate(pols):
                        if meta['naxis'] > 2:
                            rmap_plt = smap.Map(datas[pol][0, :, :], meta['header'])
                        else:
                            rmap_plt = smap.Map(datas[pol], meta['header'])
                        rmap_plt_ = pmX.Sunmap(rmap_plt)
                        rmap_plt_.imshow(axes=[axs[pidx][0], axs[pidx][1]], cmap=cmaps[pol], interpolation='nearest')
                        axs[pidx][0].set_title(title + ' ' + pols[pidx], fontsize=9)
                        rmap_plt_.draw_limb(axes=[axs[pidx][0], axs[pidx][1]])
                        rmap_plt_.draw_grid(axes=[axs[pidx][0], axs[pidx][1]])
                        rect = mpl.patches.Rectangle((xyrange[0][0], xyrange[1][0]), sz_x.value, sz_y.value,
                                                     edgecolor='w',
                                                     facecolor='none')
                        axs[pidx][0].add_patch(rect)
                        # rmap_plt_.imshow(axes=axs[pidx][1], cmap=cmaps[pol],interpolation = 'nearest')
                        # rmap_plt_.draw_limb(axes=axs[pidx][1])
                        # rmap_plt_.draw_grid(axes=axs[pidx][1])
                else:
                    title = '{0} multi spw'.format(observatory, (bfreqghz + efreqghz) / 2.0)
                    for s, sp in enumerate(spw):
                        if spwplt is not None:
                            if sp not in spwplt:
                                continue
                        for pidx, pol in enumerate(pols):
                            rcmap = [cmaps[pol](freq_dist(cfreqs[s]))] * len(clvls[pol])
                            rmap_plt = smap.Map(np.squeeze(datas[pol][s, :, :]), meta['header'])
                            rmap_plt_ = pmX.Sunmap(rmap_plt)
                            if opencontour:
                                rmap_plt_.contour(axes=[axs[pidx][0], axs[pidx][1]], colors=rcmap,
                                                  levels=clvls[pol][:1] * np.nanmax(rmap_plt.data), alpha=calpha)
                            else:
                                rmap_plt_.contourf(axes=[axs[pidx][0], axs[pidx][1]], colors=rcmap,
                                                   levels=clvls[pol] * np.nanmax(rmap_plt.data), alpha=calpha)
                            axs[pidx][0].set_title(title + ' ' + pols[pidx], fontsize=9)
                            rmap_plt_.draw_limb(axes=[axs[pidx][0], axs[pidx][1]])
                            rmap_plt_.draw_grid(axes=[axs[pidx][0], axs[pidx][1]])
                            if s == 0:
                                rect = mpl.patches.Rectangle((xyrange[0][0], xyrange[1][0]), sz_x.value, sz_y.value,
                                                             edgecolor='w',
                                                             facecolor='none')
                                axs[pidx][0].add_patch(rect)
                            # rmap_plt_.contourf(axes=axs[pidx][1], colors=rcmap,
                            #                    levels=clvls[pol] * np.nanmax(rmap_plt.data), alpha=calpha)
                            # rmap_plt_.draw_limb(axes=axs[pidx][1])
                            # rmap_plt_.draw_grid(axes=axs[pidx][1])

            ax6.set_xlim(-1220, 1220)
            ax6.set_ylim(-1220, 1220)
            ax7.set_xlim(xyrange[0])
            ax7.set_ylim(xyrange[1])
            ax4.set_ylabel('')
            # ax6.set_yticklabels([])
            ax5.set_ylabel('')
            # ax7.set_yticklabels([])
            ax5.text(0.02, 0.02, observatory + ' ' + rmap.date.strftime('%H:%M:%S.%f'), verticalalignment='bottom',
                     horizontalalignment='left',
                     transform=ax5.transAxes, color='k', fontsize=9)
            ax7.text(0.02, 0.02, observatory + ' ' + rmap.date.strftime('%H:%M:%S.%f'), verticalalignment='bottom',
                     horizontalalignment='left',
                     transform=ax7.transAxes, color='k', fontsize=9)

            axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
            try:
                axs = axs + [ax3_2]
            except:
                pass
            for ax in axs:
                for tick in ax.get_xticklabels():
                    tick.set_fontsize(8)
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(8)
                ax.set_xlabel(ax.get_xlabel(), fontsize=9)
                ax.set_ylabel(ax.get_ylabel(), fontsize=9)

            fig.subplots_adjust(top=0.94, bottom=0.07, left=0.06, right=0.93, hspace=0.80, wspace=0.88)

            if nspws >= 2:
                # try:
                import matplotlib.colorbar as colorbar
                axs = [ax4, ax7]
                ax1_pos = axs[0].get_position().extents
                ax2_pos = axs[1].get_position().extents
                caxcenter = (ax1_pos[2] + ax2_pos[0]) / 2.0 - ax1_pos[2] + ax2_pos[2]
                caxwidth = (ax2_pos[0] - ax1_pos[2]) / 2.0
                cayheight = ax1_pos[3] - 0.05 - ax2_pos[1]
                cax = plt.axes((caxcenter - caxwidth / 2.0, ax2_pos[1], caxwidth, cayheight))

                ticks, bounds, vmax, vmin, freqmask = get_colorbar_params(bdinfo)

                cb = colorbar.ColorbarBase(cax, norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=icmap,
                                           orientation='vertical', boundaries=bounds, spacing='proportional',
                                           ticks=ticks, format='%4.1f', alpha=calpha)

                for fbd_lo, fbd_hi in freqmask:
                    if fbd_hi is not None:
                        cax.axhspan(fbd_lo, fbd_hi, hatch='//', edgecolor='k', facecolor='#BBBBBB')

                ax.text(0.5, 1.04, 'MW', ha='center', va='bottom', transform=cax.transAxes, color='k',
                        fontweight='normal')
                ax.text(0.5, 1.01, '[GHz]', ha='center', va='bottom', transform=cax.transAxes, color='k',
                        fontweight='normal')
                cax.xaxis.set_visible(False)
                cax.tick_params(axis="y", pad=-20., length=0, colors='k', labelsize=8)
                cax.axhline(vmin, xmin=1.0, xmax=1.2, color='k', clip_on=False)
                cax.axhline(vmax, xmin=1.0, xmax=1.2, color='k', clip_on=False)
                cax.text(1.25, 0.0, '{:.1f}'.format(vmin), fontsize=9, transform=cax.transAxes, va='center', ha='left')
                cax.text(1.25, 1.0, '{:.1f}'.format(vmax), fontsize=9, transform=cax.transAxes, va='center', ha='left')
                # cax2 = cax.twiny()
                # cax2.set_visible(False)
                # cax2.tick_params(axis="y", pad=0., length=10, colors='k', labelsize=8)
                # cax2.set_yticks([vmin,vmax])
                # except:
                #     print('Failed to plot SPW colorbar')

            fig.canvas.draw_idle()
            fig.show()
    if clearmshistory:
        ms_restorehistory(vis)
    return outfits, outfits_list, outmovie
