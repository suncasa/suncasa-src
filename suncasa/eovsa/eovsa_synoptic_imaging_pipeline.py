import numpy as np
import numpy.ma as ma
from datetime import datetime, timedelta, time
import os
from glob import glob
import matplotlib.pyplot as plt
from suncasa.utils import helioimage2fits as hf
from sunpy import map as smap
import astropy.units as u
from scipy import ndimage
from astropy.time import Time
import pandas as pd
import shutil
from suncasa.utils import mstools as mstl
from suncasa.casa_compat import import_casatasks
from tqdm import tqdm
from suncasa.io import ndfits
import logging
import inspect
import argparse
import socket

hostname = socket.gethostname()
if hostname in ['pipeline', 'inti.hpcnet.campus.njit.edu']:
    is_on_server = True
else:
    is_on_server = False

logging.basicConfig(level=logging.INFO, format='EOVSA pipeline: [%(levelname)s] - %(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

tasks = import_casatasks('gaincal', 'applycal', 'clearcal', 'delmod', 'ft', 'uvsub', 'split', 'concat', 'flagmanager',
                         'flagdata', 'tclean', 'hanningsmooth', 'imhead')
gaincal = tasks.get('gaincal')
applycal = tasks.get('applycal')
clearcal = tasks.get('clearcal')
delmod = tasks.get('delmod')
ft = tasks.get('ft')
uvsub = tasks.get('uvsub')
split = tasks.get('split')
concat = tasks.get('concat')
flagmanager = tasks.get('flagmanager')
flagdata = tasks.get('flagdata')
tclean = tasks.get('tclean')
hanningsmooth = tasks.get('hanningsmooth')
imhead = tasks.get('imhead')

from suncasa.casa_compat import import_casatools

tools = import_casatools(['qatool', 'iatool', 'cltool', 'mstool', 'tbtool'])
qatool = tools['qatool']
iatool = tools['iatool']
cltool = tools['cltool']
mstool = tools['mstool']
tbtool = tools['tbtool']
qa = qatool()
ia = iatool()
ms = mstool()
tb = tbtool()


def log_print(level, message):
    # Get the caller's stack frame and extract information
    frame = inspect.stack()[1]
    module_name = inspect.getmodulename(frame[1])
    function_name = frame[3]

    # Format the custom context information
    context = f"{module_name}::{function_name}"

    # Prepare the full log message
    full_message = f"{context} - {message}"

    # Fetch the logger
    logger = logging.getLogger()

    # Log the message with the specified level
    if level.upper() == 'WARNING':
        logger.warning(full_message)
    elif level.upper() == 'INFO':
        logger.info(full_message)
    elif level.upper() == 'DEBUG':
        logger.debug(full_message)
    elif level.upper() == 'ERROR':
        logger.error(full_message)
    elif level.upper() == 'CRITICAL':
        logger.critical(full_message)
    else:
        logger.info(full_message)  # Default to INFO if an unsupported level is given


def is_factor_of_60_minutes(tdt):
    """
    Check if tdt is a factor of 60 minutes or a harmonic of 60 minutes.

    :param tdt: Time duration to check.
    :type tdt: timedelta
    :return: True if tdt is a factor or harmonic of 60 minutes, False otherwise.
    :rtype: bool
    """
    minutes = tdt.total_seconds() / 60
    return 60 % minutes == 0 or minutes % 60 == 0


def generate_trange_series(tbg, ted, tdt, snap_to_full_hour=False):
    """
    Generate a list of time ranges using pandas.date_range, with options to snap the ranges to full hours and
    adjust the first and last time range based on specific conditions.

    :param tbg: The start time.
    :type tbg: str or datetime-like
    :param ted: The end time.
    :type ted: str or datetime-like
    :param tdt: The duration of each time range.
    :type tdt: timedelta
    :param snap_to_full_hour: Whether to snap the time series to full hours, defaults to False.
    :type snap_to_full_hour: bool
    :return: A list of tuples, each representing the start and end time of a range.
    :rtype: list of tuple
    """

    if not is_factor_of_60_minutes(tdt) and snap_to_full_hour == True:
        snap_to_full_hour = False
        log_print('WARNING',
                  'snap_to_full_hour is set to False because tdt is not a factor of 60 minutes or a harmonic of 60 minutes.')

    tbg_datetime = pd.to_datetime(tbg)
    ted_datetime = pd.to_datetime(ted)

    time_ranges = []

    # Handle the first range differently to potentially snap its end to the next full hour
    first_end = tbg_datetime + tdt
    if snap_to_full_hour:
        if first_end.hour > tbg_datetime.hour or first_end > ted_datetime:
            first_end = first_end.replace(minute=0, second=0, microsecond=0)
            if first_end <= tbg_datetime:  # In case tdt < 1 hour
                first_end += pd.Timedelta(hours=1)
        first_end = min(first_end, ted_datetime)
    time_ranges.append((tbg_datetime, first_end))

    # Generate subsequent ranges, ensuring they start on the full hour when snap_to_full_hour is True
    current_time = first_end
    if snap_to_full_hour and current_time != ted_datetime:
        current_time = current_time.replace(minute=0, second=0, microsecond=0)

    while current_time < ted_datetime:
        end_time = min(current_time + tdt, ted_datetime)
        if snap_to_full_hour and end_time.hour > current_time.hour:
            end_time = end_time.replace(minute=0, second=0, microsecond=0)
        time_ranges.append((current_time, end_time))
        current_time = end_time

    # Check for and handle a short first or last range
    if len(time_ranges) > 1:
        # If the first range is too short, merge it with the next one
        if (time_ranges[1][0] - time_ranges[0][0]) < (0.5 * tdt):
            time_ranges[1] = (time_ranges[0][0], time_ranges[1][1])
            time_ranges.pop(0)
        # If the last range is too short, merge it with the previous one
        if (time_ranges[-1][1] - time_ranges[-1][0]) < (0.5 * tdt):
            time_ranges[-2] = (time_ranges[-2][0], time_ranges[-1][1])
            time_ranges.pop()

    return time_ranges


def trange2timerange(trange):
    """
    Convert a time range tuple in datetime format to a string representation.

    :param trange: A tuple containing start and end times as datetime objects.
    :type trange: tuple
    :return: A string representation of the time range.
    :rtype: str
    """

    sttime, edtime = trange
    timerange = f"{sttime.strftime('%Y/%m/%d/%H:%M:%S')}~{edtime.strftime('%Y/%m/%d/%H:%M:%S')}"
    return timerange


def rotateimage(data, xc_centre, yc_centre, p_angle):
    """
    Rotate an image around a specified point (xc_centre, yc_centre) by a given angle.

    :param data: The image data.
    :type data: numpy.ndarray
    :param xc_centre: The x-coordinate of the rotation center.
    :type xc_centre: int
    :param yc_centre: The y-coordinate of the rotation center.
    :type yc_centre: int
    :param p_angle: The rotation angle in degrees.
    :type p_angle: float
    :return: The rotated image.
    :rtype: numpy.ndarray
    """

    padX = [data.shape[1] - xc_centre, xc_centre]
    padY = [data.shape[0] - yc_centre, yc_centre]
    imgP = np.pad(data, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, p_angle, reshape=False, order=0, prefilter=False)
    return imgR[padY[0]:-padY[1], padX[0]:-padX[1]]


def sunpymap2helioimage(sunmap, out_image):
    """
    Rotate a SunPy map from helioprojective to RA-DEC coordinates and write it to a CASA image format.

    :param sunmap: The input SunPy Map object to be rotated.
    :type sunmap: sunpy.map.Map
    :param out_image: The filepath for the output CASA image.
    :type out_image: str
    :return: The filepath to the output CASA image format.
    :rtype: str
    """
    p_ang = sunmap.meta['p_angle'] * u.deg
    ia.open(out_image)
    data_rot = rotateimage(sunmap.data, int(sunmap.reference_pixel.x.value), int(sunmap.reference_pixel.y.value),
                           -p_ang.to('deg').value)
    # data = ia.getchunk()
    data_rot = data_rot[np.newaxis, np.newaxis, :, :].T
    # indics = data_rot == 0.00
    # data_rot[indics] = data[indics]
    ia.putchunk(data_rot)
    ia.close()
    return out_image


def solar_diff_rot_image(in_map, newtime, out_image, showplt=False):
    """
    Reproject a SunPy map to account for solar differential rotation to a new observation time.

    :param in_map: The input SunPy Map object to be reprojected.
    :type in_map: sunpy.map.Map
    :param newtime: The new time to which the map is reprojected.
    :type newtime: astropy.time.Time
    :param out_image: The path for the output image file in CASA format.
    :type out_image: str
    :param showplt: Boolean flag to show plots of the original and reprojected maps, defaults to False.
    :type showplt: bool
    :return: The path to the output CASA image format.
    :rtype: str
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    from sunpy.coordinates import Helioprojective, propagate_with_solar_surface
    ## sunpy reproject the map.date as the reference time. so we have to use the diff between the time of the each map and the ref map
    reftime = in_map.date + in_map.exposure_time / 2
    out_time = in_map.date - (reftime - Time(newtime))
    # out_frame = Helioprojective(observer='earth', obstime=out_time,
    #                             rsun=eomap_ref.coordinate_frame.rsun)
    out_frame = Helioprojective(observer=in_map.observer_coordinate, obstime=out_time,
                                rsun=in_map.coordinate_frame.rsun)
    out_center = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=out_frame)
    out_ref_pixel = [in_map.reference_pixel.x.value, in_map.reference_pixel.y.value] * in_map.reference_pixel.x.unit
    out_header = smap.make_fitswcs_header(in_map.data.shape,
                                          out_center,
                                          reference_pixel=out_ref_pixel,
                                          scale=u.Quantity(in_map.scale))
    out_wcs = WCS(out_header)
    with propagate_with_solar_surface():
        out_map, footprint = in_map.reproject_to(out_wcs, return_footprint=True)

    out_data = out_map.data

    out_data[footprint == 0] = in_map.data[footprint == 0]
    out_map = smap.Map(out_data, out_map.meta)
    out_map.meta['p_angle'] = in_map.meta['p_angle']
    if showplt:
        fig = plt.figure(figsize=(12, 4))

        ax1 = fig.add_subplot(121, projection=in_map)
        in_map.plot(axes=ax1, title='Original map')
        plt.colorbar()

        ax2 = fig.add_subplot(122, projection=out_map)
        out_map.plot(axes=ax2,
                     title=f"Reprojected to an Earth observer {(Time(newtime) - reftime).to('day')} later")
        plt.colorbar()

        plt.show()
    sunpymap2helioimage(out_map, out_image)

    return out_image


def get_bmsize(cfreq, refbmsize=70.0, reffreq=1.0, minbmsize=4.0):
    """
    Calculate the beam size at given frequencies based on a reference beam size at a reference frequency.
    This function supports both single frequency values and lists of frequencies.

    :param cfreq: Input frequencies in GHz, can be a float or a list of floats.
    :type cfreq: float or list
    :param refbmsize: Reference beam size in arcsec, defaults to 70.0.
    :type refbmsize: float, optional
    :param reffreq: Reference frequency in GHz, defaults to 1.0.
    :type reffreq: float, optional
    :param minbmsize: Minimum beam size in arcsec, defaults to 4.0.
    :type minbmsize: float, optional
    :return: Beam size at the given frequencies, same type as input (float or numpy array).
    :rtype: float or numpy.ndarray
    """
    # Ensure cfreq is an array for uniform processing
    cfreq = np.array(cfreq, dtype=float)

    # Calculate beam size
    bmsize = refbmsize * reffreq / cfreq

    # Enforce minimum beam size
    bmsize = np.maximum(bmsize, minbmsize)

    # If the original input was a single float, return a single float
    if bmsize.size == 1:
        return bmsize.item()  # Convert numpy scalar to Python float
    else:
        return bmsize


def calc_diskmodel(slashdate, nbands, freq, defaultfreq):
    from astropy.time import Time
    # Default disk size measured for 2019/09/03
    # todo add monthly fitting procedure for the disk size and flux density
    defaultsize = np.array([990.6, 989.4, 988.2, 987.1, 986.0, 984.9, 983.8, 982.7, 981.7, 980.7,
                            979.7, 978.8, 977.8, 976.9, 976.0, 975.2, 974.3, 973.5, 972.7, 972.0,
                            971.2, 970.5, 969.8, 969.1, 968.5, 967.8, 967.2, 966.7, 966.1, 965.6,
                            965.1, 964.6, 964.1, 963.7, 963.3, 962.9, 962.5, 962.1, 961.8, 961.5,
                            961.3, 961.0, 960.8, 960.6, 960.4, 960.2, 960.1, 960.0, 959.9, 959.8])

    # Get current solar distance and modify the default size accordingly
    try:
        from sunpy.coordinates.sun import earth_distance
        fac = earth_distance('2019/09/03') / earth_distance(slashdate)
    except:
        import sunpy.coordinates.ephemeris as eph
        fac = eph.get_sunearth_distance('2019/09/03') / eph.get_sunearth_distance(slashdate)

    newsize = defaultsize * fac.to_value()
    if nbands == 34:
        if Time(slashdate.replace('/', '-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectral windows (bands 4-34 -> spw 0-30)
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[[0] + list(range(2, 34))])
    dsize = np.array([str(i)[:5] + 'arcsec' for i in newsize], dtype='U12')

    # These are nominal flux densities * 2, determined on 2019/09/03
    defaultfdens = np.array([891282, 954570, 1173229, 1245433, 1373730, 1506802,
                             1613253, 1702751, 1800721, 1946756, 2096020, 2243951,
                             2367362, 2525968, 2699795, 2861604, 3054829, 3220450,
                             3404182, 3602625, 3794312, 3962926, 4164667, 4360683,
                             4575677, 4767210, 4972824, 5211717, 5444632, 5648266,
                             5926634, 6144249, 6339863, 6598018, 6802707, 7016012,
                             7258929, 7454951, 7742816, 7948976, 8203206, 8411834,
                             8656720, 8908130, 9087766, 9410760, 9571365, 9827078,
                             10023598, 8896671])
    fdens = defaultfdens
    if nbands == 34:
        if Time(slashdate.replace('/', '-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectal windows (bands 4-34 -> spw 0-30)
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[[0] + list(range(2, 34))])
    return dsize, fdens


def writediskxml(dsize, fdens, freq, xmlfile='SOLDISK.xml'):
    import xml.etree.ElementTree as ET
    # create the file structure
    sdk = ET.Element('SOLDISK')
    sdk_dsize = ET.SubElement(sdk, 'item')
    sdk_fdens = ET.SubElement(sdk, 'item')
    sdk_freqs = ET.SubElement(sdk, 'item')
    sdk_dsize.set('disk_size', ','.join(dsize))
    sdk_fdens.set('flux_dens', ','.join(['{:.1f}Jy'.format(s) for s in fdens]))
    sdk_freqs.set('freq', ','.join(freq))

    # create a new XML file with the results
    mydata = ET.tostring(sdk)
    if isinstance(mydata, bytes):
        mydata = mydata.decode()
    if os.path.exists(xmlfile):
        os.system('rm -rf ' + xmlfile)
    with open(xmlfile, 'w') as sf:
        sf.write(mydata)
    return xmlfile


def readdiskxml(xmlfile):
    import astropy.units as u
    import xml.etree.ElementTree as ET
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    diskinfo = {}
    for elem in root:
        d = elem.attrib
        for k, v in d.items():
            v_ = v.split(',')
            v_ = [u.Unit(f).to_string().split(' ') for f in v_]
            diskinfo[k] = []
            for val, uni in v_:
                diskinfo[k].append(float(val))
            diskinfo[k] = np.array(diskinfo[k]) * u.Unit(uni)

    return diskinfo


def gaussian2d(x, y, amplitude, x0, y0, sigma_x, sigma_y, theta):
    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = amplitude * np.exp(- (a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) ** 2)))
    return g


def image_adddisk(eofile, diskinfo, edgeconvmode='frommergeddisk', caltbonly=False, bmfactor=2.0, overwrite=True):
    '''

    :param eofile: input image FITS file
    :param diskinfo: disk information file
    :param edgeconvmode: edge convolve mode, 'frommergeddisk' or 'frombeam'
    :param caltbonly: calculate the Tb of the disk and return the value
    :param bmfactor:  factor to multiply the beam major and minor axes to get the sigma of the Gaussian kernel
    :param overwrite: whether to overwrite the output fits file if it already exists
    :return:
    '''

    from sunpy import map as smap
    from suncasa.utils import plot_mapX as pmX
    from scipy import constants
    import astropy.units as u
    from sunpy import io as sio
    dsize = diskinfo['disk_size']
    fdens = diskinfo['flux_dens']
    freqs = diskinfo['freq']
    eomap = smap.Map(eofile)
    data = eomap.data  # remember the data order is reversed due to the FITS convension
    if np.all(data == 0):
        log_print('WARNING', 'The input image is empty. Skipping...')
        return None, None, None
    eomap_ = pmX.Sunmap(eomap)
    header = eomap.meta
    bmaj = header['bmaj'] * 3600 * u.arcsec
    bmin = header['bmin'] * 3600 * u.arcsec
    cell = (header['cdelt1'] * u.Unit(header['cunit1']) + header['cdelt2'] * u.Unit(header['cunit2'])) / 2.0
    bmsize = (bmaj + bmin) / 2.0
    keys = list(header.keys())
    values = list(header.values())
    mapx, mapy = eomap_.map2wcsgrids(cell=False)
    mapx = mapx[:-1, :-1]
    mapy = mapy[:-1, :-1]
    rdisk = np.sqrt(mapx ** 2 + mapy ** 2)

    k_b = constants.k
    c_l = constants.c
    const = 2. * k_b / c_l ** 2
    pix_area = (cell.to(u.rad).value) ** 2
    jy_to_si = 1e-26
    factor2 = 1.
    try:
        ## this works for uncompressed FITS
        faxis = keys[values.index('FREQ')][-1]
        crval_freq = header['CRVAL' + faxis]
        cdelt_freq = header['CDELT' + faxis]
        nu = crval_freq + cdelt_freq * (1 - header['CRPIX' + faxis])
        nul = crval_freq + cdelt_freq * (1 - header['CRPIX' + faxis])
        nuh = crval_freq + cdelt_freq * (header['NAXIS' + faxis] - header['CRPIX' + faxis])
        cunit_freq = header['CUNIT' + faxis]
    except:
        ## this works for compressed FITS
        crval_freq = header['REFFRQ']
        cdelt_freq = header['DLTFRQ']
        try:
            cunit_freq = header['FQUNIT']
        except:
            cunit_freq = 'Hz'
        nu = crval_freq
        nul = crval_freq - cdelt_freq / 2.0
        nuh = crval_freq + cdelt_freq / 2.0

    if caltbonly:
        edgeconvmode = ''
    if edgeconvmode == 'frommergeddisk':
        ## get the frequency range of the image
        nu_bound = (np.array([nul, nuh]) + 0.5 * np.array([-1, 1]) * cdelt_freq) * u.Unit(
            cunit_freq)
        nu_bound = nu_bound.to(u.GHz)
        ## get the frequencies of the disk models
        fidxs = np.logical_and(freqs > nu_bound[0], freqs < nu_bound[1])
        ny, nx = rdisk.shape
        freqs_ = freqs[fidxs]
        fdens_ = fdens[fidxs] / 2.0  # divide by 2 because fdens is 2x solar flux density
        dsize_ = dsize[fidxs]
        fdisk_ = np.empty((len(freqs_), ny, nx))
        fdisk_[:] = np.nan
        for fidx, freq in enumerate(freqs_):
            fdisk_[fidx, ...][rdisk <= dsize_[fidx].value] = 1.0
            # nu = crval_freq + cdelt_freq * (1 - header['CRPIX' + faxis])
            factor = const * freq.to(u.Hz).value ** 2  # SI unit
            jy2tb = jy_to_si / pix_area / factor * factor2
            fdisk_[fidx, ...] = fdisk_[fidx, ...] / np.nansum(fdisk_[fidx, ...]) * fdens_[fidx].value
            fdisk_[fidx, ...] = fdisk_[fidx, ...] * jy2tb
        #         # fdisk_[np.isnan(fdisk_)] = 0.0
        tbdisk = np.nanmean(fdisk_, axis=0)
        tbdisk[np.isnan(tbdisk)] = 0.0

        sig2fwhm = 2.0 * np.sqrt(2 * np.log(2)) * bmfactor
        x0, y0 = 0, 0
        sigx, sigy = bmaj.value / sig2fwhm, bmin.value / sig2fwhm
        theta = -(90.0 - header['bpa']) * u.deg
        x = (np.arange(31) - 15) * cell.value
        y = (np.arange(31) - 15) * cell.value
        x, y = np.meshgrid(x, y)
        kernel = gaussian2d(x, y, 1.0, x0, y0, sigx, sigy, theta.to(u.radian).value)
        kernel = kernel / np.nansum(kernel)
        from scipy import signal
        tbdisk = signal.fftconvolve(tbdisk, kernel, mode='same')
    else:
        freqghz = nu / 1.0e9
        factor = const * nu ** 2  # SI unit
        jy2tb = jy_to_si / pix_area / factor * factor2
        p_dsize = np.poly1d(np.polyfit(freqs.value, dsize.value, 15))
        p_fdens = np.poly1d(
            np.polyfit(freqs.value, fdens.value, 15)) / 2.  # divide by 2 because fdens is 2x solar flux density
        if edgeconvmode == 'frombeam':
            from scipy.special import erfc
            factor_erfc = 2.0  ## erfc function ranges from 0 to 2
            fdisk = erfc((rdisk - p_dsize(freqghz)) / bmsize.value) / factor_erfc
        else:
            fdisk = np.zeros_like(rdisk)
            fdisk[rdisk <= p_dsize(freqghz)] = 1.0
        fdisk = fdisk / np.nansum(fdisk) * p_fdens(freqghz)
        tbdisk = fdisk * jy2tb

    tb_disk = np.nanmax(tbdisk)
    if caltbonly:
        return tb_disk
    else:
        datanew = data + tbdisk
        # datanew[np.isnan(data)] = 0.0
        header['TBDISK'] = tb_disk
        header['TBUNIT'] = 'K'
        eomap_disk = smap.Map(datanew, header)
        nametmp = eofile.split('.')
        nametmp.insert(-1, 'disk')
        outfits = '.'.join(nametmp)
        datanew = datanew.astype(np.float32)
        if os.path.exists(outfits):
            if overwrite:
                os.system('rm -rf ' + outfits)
            else:
                eomap_disk, tb_disk, outfits
        sio.write_file(outfits, datanew, header)
        return eomap_disk, tb_disk, outfits


def mk_diskmodel(outname='disk', direction='J2000 10h00m00.0s 20d00m00.0s',
                 reffreq='2.8GHz', flux=660000.0, eqradius='16.166arcmin', polradius='16.166arcmin',
                 pangle='21.1deg', overwrite=True):
    ''' Create a blank solar disk model image (or optionally a data cube)
        outname       String to use for part of the image and fits file names (default 'disk')
        direction     String specifying the position of the Sun in RA and Dec.  Default
                        means use the standard string "J2000 10h00m00.0s 20d00m00.0s"
        reffreq       The reference frequency to use for the disk model (the frequency at which
                        the flux level applies). Default is '2.8GHz'.
        flux          The flux density, in Jy, for the entire disk. Default is 66 sfu.
        eqradius      The equatorial radius of the disk.  Default is
                        16 arcmin + 10" (for typical extension of the radio limb)
        polradius     The polar radius of the disk.  Default is
                        16 arcmin + 10" (for typical extension of the radio limb)
        pangle        The solar P-angle (geographic position of the N-pole of the Sun) in
                        degrees E of N.  This only matters if eqradius != polradius
        index         The spectral index to use at other frequencies.  Default None means
                        use a constant flux density for all frequencies.
        cell          The cell size (assumed square) to use for the image.  The image size
                        is determined from a standard radius of 960" for the Sun, divided by
                        cell size, increased to nearest power of 512 pixels. The default is '2.0arcsec',
                        which results in an image size of 1024 x 1024.
        Note that the frequency increment used is '325MHz', which is the width of EOVSA bands
          (not the width of individual science channels)
    '''

    diskcl = outname + reffreq + '.cl'
    if os.path.exists(diskcl):
        if overwrite:
            os.system('rm -rf ' + diskcl)
        else:
            return diskcl

    cl = cltool()

    try:
        aspect = 1.01  # Enlarge the equatorial disk by 1%
        eqradius = qa.quantity(eqradius)
        diamajor = qa.quantity(2 * aspect * eqradius['value'], eqradius['unit'])
        polradius = qa.quantity(polradius)
        diaminor = qa.quantity(2 * polradius['value'], polradius['unit'])
    except:
        print('Radius', eqradius, polradius,
              'does not have the expected format, number + unit where unit is arcmin or arcsec')
        return

    # Add 90 degrees to pangle, due to angle definition in addcomponent() -- it puts the majoraxis vertical
    pangle = qa.add(qa.quantity(pangle), qa.quantity('90deg'))
    # Flux density is split between XX and YY
    cl.addcomponent(dir=direction, flux=flux / 2.0, fluxunit='Jy', freq=reffreq, shape='disk',
                    majoraxis=diamajor, minoraxis=diaminor, positionangle=pangle)
    cl.setrefdirframe(0, 'J2000')
    cl.rename(diskcl)
    cl.done()
    return diskcl


def insertdiskmodel(vis, sizescale=1.0, fdens=None, dsize=None, xmlfile='SOLDISK.xml', writediskinfoonly=False,
                    active=False, overwrite=False):
    # Apply size scale adjustment (default is no adjustment)
    for i in range(len(dsize)):
        dsz = dsize[i]
        if isinstance(dsz, bytes):
            dsz = dsz.decode()
        num, unit = dsz.split('arc')
        dsize[i] = str(float(num) * sizescale)[:6] + 'arc' + unit

    msfile = vis

    ms.open(msfile)
    spwinfo = ms.getspectralwindowinfo()
    nspw = len(spwinfo.keys())
    ms.done()
    diskcldir = 'diskcl/'
    if not os.path.exists(diskcldir):
        os.makedirs(diskcldir)
    frq = []
    spws = range(nspw)
    for sp in spws:
        spw = spwinfo[str(sp)]
        frq.append('{:.4f}GHz'.format((spw['RefFreq'] + spw['TotalWidth'] / 2.0) / 1e9))
    frq = np.array(frq)
    writediskxml(dsize, fdens, frq, xmlfile=xmlfile)

    if not writediskinfoonly:
        tb.open(msfile + '/FIELD')
        phadir = tb.getcol('PHASE_DIR').flatten()
        tb.close()
        ra = phadir[0]
        dec = phadir[1]
        direction = 'J2000 ' + str(ra) + 'rad ' + str(dec) + 'rad'

        diskcl = []
        for sp in tqdm(spws, desc='Generating {} disk models'.format(nspw), ascii=True):
            diskcl.append(
                mk_diskmodel(outname=diskcldir + 'disk{:02d}_'.format(sp),
                             direction=direction, reffreq=frq[sp],
                             flux=fdens[sp], eqradius=dsize[sp], polradius=dsize[sp], overwrite=overwrite))

        if not active:
            delmod(msfile, otf=True, scr=True)
            for sp in tqdm(spws, desc='Inserting disk model', ascii=True):
                ft(vis=msfile, spw=str(sp), field='', model="", nterms=1,
                   reffreq="", complist=str(diskcl[sp]), incremental=False, usescratch=True)
        else:
            for sp in tqdm(spws, desc='Inserting disk model', ascii=True):
                model_ft = mstl.getmodel(msfile, spw=str(sp))
                ft(vis=msfile, spw=str(sp), field='', model="", nterms=1,
                   reffreq="", complist=str(diskcl[sp]), incremental=False, usescratch=True)
                model_disk = mstl.getmodel(msfile, spw=str(sp))
                mstl.putmodel(msfile, spw=str(sp), model=model_ft + model_disk)

        return msfile, diskcl


def uvrange_uplim_from_freq(x, x0, x1, y0, y1):
    """
    Calculate the upper limit of the UV range from frequency using linear interpolation.

    :param x: The frequency for which the UV range upper limit is calculated.
    :type x: float
    :param x0: The start frequency of the interpolation range.
    :type x0: float
    :param x1: The end frequency of the interpolation range.
    :type x1: float
    :param y0: The UV range upper limit at the start frequency.
    :type y0: float
    :param y1: The UV range upper limit at the end frequency.
    :type y1: float
    :return: The calculated UV range upper limit for the given frequency.
    :rtype: float
    """
    return y0 + (x - x0) * ((y1 - y0) / (x1 - x0))


def disk_slfcal(msfile, tbg, ted, disk_params, workdir='./', overwrite=True, iterbands=False):
    """
    Perform disk self-calibration on measurement set data.

    :param msfile: The measurement set file to be calibrated.
    :type msfile: str
    :param tbg: The beginning time of the calibration range.
    :type tbg: datetime
    :param ted: The end time of the calibration range.
    :type ted: datetime
    :param disk_params: A dictionary containing disk model parameters.
    :type disk_params: dict
    :param iterbands: Boolean flag to run gaincal iterating over frequency bands, defaults to False. iterbands = True is only useful when uvrange_uplim is used, which is currently not implemented.
    :type iterbands: bool
    :return: The filepath of the self-calibrated measurement set.
    :rtype: str
    """

    # Add the reference model back to the visibility data
    # ft(vis=msfile, model=img_ref_model, usescratch=True)
    msfile_diskscl = msfile.replace('.ms', '.disk_slfcal.ms')
    msfile_diskscled = msfile.replace('.ms', '.disk_slfcaled.ms')
    if overwrite:
        if os.path.isdir(msfile_diskscled):
            shutil.rmtree(msfile_diskscled, ignore_errors=True)
    else:
        if os.path.isdir(msfile_diskscled):
            return msfile_diskscled

    freq = disk_params['freq']
    dsize = disk_params['dsize']
    fdens = disk_params['fdens']
    diskxmlfile = disk_params['diskxmlfile']
    msfile, diskcl = insertdiskmodel(msfile, dsize=dsize, fdens=fdens, xmlfile=diskxmlfile, active=True)

    log_print('INFO', 'First round of phase disk selfcal on the disk using solution interval inf')
    caltbs_dk = []
    caltb = os.path.join(workdir, f"caltb_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}_dk1.pha")
    if os.path.isdir(caltb):
        shutil.rmtree(caltb, ignore_errors=True)
    if iterbands:
        appendmode = False
        for sp, fghz in enumerate(freq):
            if os.path.isdir(caltb):
                appendmode = True
            # uvrange_uplim = uvrange_uplim_from_freq(fghz, x0=2.8875, x1=17.1875, y0=0.2, y1=0.8) * 3
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    spw=f'{sp}',
                    # uvrange=f"<{uvrange_uplim:.1f}Klambda",
                    uvrange="",
                    antenna="0~12&0~12",
                    # solint=tdtstr,
                    solint='inf',
                    combine="scan", interp="linear",
                    refant="0", refantmode="flex", minsnr=1.0,
                    gaintype="G", calmode="p", append=appendmode)
    else:
        gaincal(vis=msfile, caltable=caltb, selectdata=True,
                uvrange="",
                antenna="0~12&0~12",
                # solint=tdtstr,
                solint='inf',
                combine="scan", interp="linear",
                refant="0", refantmode="flex", minsnr=1.0,
                gaintype="G", calmode="p")
    if os.path.isdir(caltb):
        caltbs_dk.append(caltb)

    # plotms(vis='temp_20211124/caltb_1545-2334_dk1.pha', xaxis='freq', yaxis='phase',
    #        iteraxis='antenna', coloraxis='',
    #        plotrange=[-1, -1, -180, 180])

    if len(caltbs_dk) > 0:
        log_print('INFO', 'Second round of phase disk selfcal on the disk using solution interval 1min')
        caltb = os.path.join(workdir, f"caltb_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}_dk2.pha")
        if os.path.isdir(caltb):
            shutil.rmtree(caltb, ignore_errors=True)
        if iterbands:
            appendmode = False
            for sp, fghz in enumerate(freq):
                if os.path.isdir(caltb):
                    appendmode = True
                # uvrange_uplim = uvrange_uplim_from_freq(fghz, x0=2.8875, x1=17.1875, y0=0.2, y1=0.8) * 3
                gaincal(vis=msfile, caltable=caltb, selectdata=True,
                        spw=f'{sp}',
                        # uvrange=f"<{uvrange_uplim:.1f}Klambda",
                        uvrange="",
                        antenna="0~12&0~12",
                        solint='int',
                        gaintable=caltbs_dk,
                        combine="scan", interp="linear",
                        refant="0", refantmode="flex", minsnr=1.0,
                        gaintype="G", calmode="p", append=appendmode)
        else:
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    uvrange="",
                    antenna="0~12&0~12",
                    solint='int',
                    gaintable=caltbs_dk,
                    combine="scan", interp="linear",
                    refant="0", refantmode="flex", minsnr=1.0,
                    gaintype="G", calmode="p")
        if os.path.isdir(caltb):
            caltbs_dk.append(caltb)

    # plotms(vis='temp_20211124/caltb_1545-2334_dk2.pha', xaxis='time', yaxis='phase',
    #        iteraxis='spw', coloraxis='antenna1',
    #        plotrange=[-1, -1, -180, 180])

    if len(caltbs_dk) > 0:
        log_print('INFO', 'Final round of amplitude disk selfcal using solution interval inf')
        caltb = os.path.join(workdir, f"caltb_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}_dk.amp")
        if os.path.isdir(caltb):
            shutil.rmtree(caltb, ignore_errors=True)
        if iterbands:
            appendmode = False
            for sp, fghz in enumerate(freq):
                if os.path.isdir(caltb):
                    appendmode = True
                # uvrange_uplim = uvrange_uplim_from_freq(fghz, x0=2.8875, x1=17.1875, y0=0.2, y1=0.8) * 3
                gaincal(vis=msfile, caltable=caltb, selectdata=True,
                        spw=f'{sp}',
                        # uvrange=f"<{uvrange_uplim:.1f}Klambda",
                        uvrange="",
                        antenna="0~12&0~12",
                        solint='inf',
                        gaintable=caltbs_dk,
                        combine="scan", interp="linear",
                        refant="10", refantmode="flex",
                        minsnr=1.0,
                        gaintype="G", calmode="a", append=appendmode)
        else:
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    uvrange="",
                    antenna="0~12&0~12",
                    solint='inf',
                    gaintable=caltbs_dk,
                    combine="scan", interp="linear",
                    refant="10", refantmode="flex",
                    minsnr=1.0,
                    gaintype="G", calmode="a")
        if os.path.isdir(caltb):
            caltbs_dk.append(caltb)

    # plotms(vis='temp_20211124/caltb_1545-2334_dk.amp', xaxis='time', yaxis='amp',
    #        iteraxis='spw', coloraxis='antenna1',
    #        plotrange=[-1, -1, -1, -1])

    # clearcal(msfile)
    if len(caltbs_dk) > 0:
        log_print('INFO', 'Applying disk self-calibration solutions to the data and split to a new ms')
        applycal(vis=msfile, selectdata=True, antenna="0~12", gaintable=caltbs_dk, interp="linear", calwt=False,
                 applymode="calonly")

        if os.path.isdir(msfile_diskscl):
            shutil.rmtree(msfile_diskscl, ignore_errors=True)
        split(msfile, outputvis=msfile_diskscl, datacolumn="corrected")

        for sp, dkcl in tqdm(enumerate(diskcl), desc='Inserting disk model', ascii=True):
            ft(vis=msfile_diskscl, spw=str(sp), field='', model="", nterms=1,
               reffreq="", complist=str(dkcl), incremental=False, usescratch=True)
        log_print('INFO', 'Removing the disk model from the data. The residual will be ARs only')
        uvsub(vis=msfile_diskscl)
        split(msfile_diskscl, outputvis=msfile_diskscled, datacolumn="corrected")
        shutil.rmtree(msfile_diskscl, ignore_errors=True)
    else:
        log_print('WARNING',
                  'Failed to create calibration table for disk self-calibration. Proceeding without self-calibration...')
        msfile_diskscled = msfile

    return msfile_diskscled


def shift_corr(mmsfiles, trange_series, spws, imagemodel, imagemodel_fits, reftime_master, workdir='./',
               pols='XX', overwrite=False, do_featureslfcal=False, do_diskslfcal=True, disk_params={}, do_sbdcal=True,
               verbose=True):
    """
    Corrects smearing effects in solar observation data by aligning them with a model image.

    :param mmsfiles: Paths to the measurement set files to be corrected.
    :type mmsfiles: list of str
    :param trange_series: Each tuple contains start and end times for the measurement set files.
    :type trange_series: list of tuple
    :param spws: Spectral window indices to be considered.
    :type spws: list of string
    :param imagemodel: Path to the model image in CASA measurement set format.
    :type imagemodel: str
    :param imagemodel_fits: Path to the model image in FITS format.
    :type imagemodel_fits: str
    :param reftime_master: Reference time for the rotation of the model image.
    :type reftime_master: astropy.time.Time
    :param workdir: Working directory for output files, defaults to './'.
    :type workdir: str
    :param do_featureslfcal: Boolean flag to perform feature self-calibration, defaults to False.
    :type do_featureslfcal: bool
    :param pols: Polarization types to be considered, defaults to 'XX'.
    :type pols: str
    :param overwrite: Boolean flag to overwrite existing files, defaults to False.
    :type overwrite: bool
    :param do_diskslfcal: Boolean flag to perform disk self-calibration, defaults to True.
    :type do_diskslfcal: bool
    :param disk_params: Dictionary containing parameters for disk self-calibration.
    :type disk_params: dict
    :param do_sbdcal: Boolean flag to perform single-band delay calibration, defaults to True.
    :type do_sbdcal: bool
    :param verbose: Boolean flag to enable verbose output, defaults to True.
    :type verbose: bool
    :return: Paths to the corrected measurement set files.
    :rtype: list of str
    """

    # Load the model image maps and calculate the reference time
    imagemodel_map = smap.Map(imagemodel_fits)

    # Calculate mid-points of time ranges for rotation
    trange_series_mid = [(start + (end - start) / 2) for start, end in trange_series]
    mmsfiles_rot = []
    total_slots = len(trange_series)
    if verbose:
        print("=" * 64)
        print("-" * 64)
        print(f"Processing {total_slots} sub-MSs for solar rotation correction.")

    trange_series_bg = trange_series[0][0]
    trange_series_ed = trange_series[-1][1]
    trange_series_str = f"{trange_series_bg.strftime('%H%M')}-{trange_series_ed.strftime('%H%M')}UT"
    reftime_master_str = f"{reftime_master.to_datetime().strftime('%H%M')}UT"
    imagemodel_rot_master = []
    for sp, immodel_map, immodel, immodel_fits in zip(spws, imagemodel_map, imagemodel, imagemodel_fits):
        spwstr = format_spw(sp)
        immodel_rot_master = os.path.join(workdir,
                                          f"eovsa_{trange_series_str}.s{spwstr}.rot2_{reftime_master_str}.master.model")
        if overwrite:
            if os.path.isdir(immodel_rot_master):
                shutil.rmtree(immodel_rot_master, ignore_errors=True)
        if not os.path.isdir(immodel_rot_master):
            shutil.copytree(immodel, immodel_rot_master, dirs_exist_ok=True)
            solar_diff_rot_image(immodel_map, reftime_master, immodel_rot_master)
        imagemodel_rot_master.append(immodel_rot_master)

    for tidx, (mmsfile, time_range, time_mid) in enumerate(zip(mmsfiles, trange_series, trange_series_mid)):
        tbg, ted = time_range
        trange = trange2timerange(time_range)
        if mmsfile is None:
            mmsfiles_rot.append(None)
            if verbose:
                log_print('INFO',
                          f"Corrected {tidx + 1} of {total_slots}: No sub-MS file found for this timerange. Skipping...")
            continue
        mmsfile_basename = os.path.splitext(mmsfile.rstrip('/'))[0]
        msfile_rot = f"{mmsfile_basename}.shift_corrected.ms"
        if overwrite:
            if os.path.isdir(msfile_rot):
                shutil.rmtree(msfile_rot, ignore_errors=True)
        else:
            if os.path.exists(msfile_rot):
                mmsfiles_rot.append(msfile_rot)
                if verbose:
                    log_print('INFO',
                              f"Corrected {tidx + 1} of {total_slots}: {mmsfile} → {msfile_rot} (already exists and was not re-created)")
                continue
        for sp, immodel_map, immodel in zip(spws, imagemodel_map, imagemodel):
            spwstr = format_spw(sp)
            imfile_model_rot = os.path.join(workdir,
                                            f"eovsa_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}UT.s{spwstr}.rot.model")
            if not os.path.isdir(msfile_rot):
                if overwrite:
                    if os.path.isdir(imfile_model_rot):
                        shutil.rmtree(imfile_model_rot, ignore_errors=True)
                if not os.path.isdir(imfile_model_rot):
                    shutil.copytree(immodel, imfile_model_rot, dirs_exist_ok=True)
                    if verbose:
                        log_print('INFO',
                                  f"Rotate the model image ({os.path.basename(immodel)}) to {time_mid.strftime('%H:%M')}UT")
                    solar_diff_rot_image(immodel_map, time_mid, imfile_model_rot)
                if verbose:
                    log_print('INFO',
                              f"Add the rotated model ({os.path.basename(imfile_model_rot)}) to the visibility data")
                ft(vis=mmsfile, model=imfile_model_rot, spw=sp, usescratch=True)
        if do_featureslfcal:
            if verbose:
                log_print('INFO', 'Performing feature self-calibration')
            mmsfile_slfcaled = f"{mmsfile_basename}.slfcaled.ms"
            caltb = os.path.join(workdir, f"caltb_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}.pha")
            if overwrite:
                if os.path.isdir(mmsfile_slfcaled):
                    shutil.rmtree(mmsfile_slfcaled, ignore_errors=True)
                if os.path.isdir(caltb):
                    shutil.rmtree(caltb, ignore_errors=True)
            if verbose:
                log_print('INFO', 'Attempt feature self-calibration starting with a uvrange > 1.5 klambda')
            if not os.path.isdir(caltb):
                if pols == 'XXYY':
                    mstl.gaincalXY(vis=mmsfile, caltable=caltb, pols=pols, selectdata=True,
                                   timerange=trange,
                                   uvrange='>1.5klambda',
                                   combine="scan", antenna='0~12&0~12', refant='0', solint='inf',
                                   refantmode="strict",
                                   gaintype='G', minsnr=1.0, calmode='p', append=False)
                else:
                    gaincal(vis=mmsfile, caltable=caltb, selectdata=True,
                            timerange=trange,
                            uvrange='>1.5klambda',
                            combine="scan", antenna='0~12&0~12', refant='0', solint='inf', refantmode="strict",
                            gaintype='G',
                            minsnr=1.0, calmode='p', append=False)
            if not os.path.isdir(caltb):
                if verbose:
                    log_print('INFO',
                              'The first feature self-calibration attempt fails, try again with no uvrange condition')
                if pols == 'XXYY':
                    mstl.gaincalXY(vis=mmsfile, caltable=caltb, pols=pols, selectdata=True,
                                   timerange=trange,
                                   uvrange='',
                                   combine="scan", antenna='0~12&0~12', refant='0', solint='inf',
                                   refantmode="strict",
                                   gaintype='G', minsnr=1.0, calmode='p', append=False)
                else:
                    gaincal(vis=mmsfile, caltable=caltb, selectdata=True,
                            timerange=trange,
                            uvrange='',
                            combine="scan", antenna='0~12&0~12', refant='0', solint='inf', refantmode="strict",
                            gaintype='G',
                            minsnr=1.0, calmode='p', append=False)
            if do_sbdcal:
                caltb_k = os.path.join(workdir, f"caltb_{tbg.strftime('%H%M')}-{ted.strftime('%H%M')}.sbd")
                if os.path.isdir(caltb_k):
                    os.system('rm -rf ' + caltb_k)
                if pols == 'XXYY':
                    mstl.gaincalXY(vis=mmsfile, caltable=caltb_k, pols=pols, selectdata=True,
                                   timerange=trange,
                                   uvrange='',
                                   combine="scan", antenna='0~12&0~12', refant='0', solint='inf',
                                   refantmode="strict",
                                   gaintype='K', calmode='p',
                                   minblperant=4, minsnr=2, append=False)
                else:
                    gaincal(vis=mmsfile, caltable=caltb_k, selectdata=True,
                            timerange=trange,
                            uvrange='',
                            combine="scan", antenna='0~12&0~12', refant='0', solint='inf', refantmode="strict",
                            gaintype='K',
                            minsnr=2, calmode='p',
                            minblperant=4, append=False)
            if os.path.isdir(caltb):
                if do_sbdcal:
                    if os.path.isdir(caltb_k):
                        caltb = [caltb, caltb_k]
                log_print('INFO', 'Applying feature self-calibration solutions to the data and split to a new ms')
                applycal(vis=mmsfile, selectdata=True, antenna="0~12", gaintable=caltb, interp="linear",
                         calwt=False, applymode="calonly")
                # split the corrected data to a new ms
                mstl.splitX(mmsfile, outputvis=mmsfile_slfcaled, datacolumn="corrected",
                            datacolumn2="model_data")
                shutil.rmtree(mmsfile, ignore_errors=True)
                mmsfile = mmsfile_slfcaled
            else:
                if verbose:
                    log_print('INFO',
                              f'Feature self-calibration failed for {mmsfile}. Proceeding without self-calibration...')
        uvsub(vis=mmsfile)  # Leaves only the residual in the corrected data column
        delmod(vis=mmsfile)

        # Add the reference model back to the visibility data
        for sp, immodel_rot_master in zip(spws, imagemodel_rot_master):
            ft(vis=mmsfile, model=immodel_rot_master, spw=sp, usescratch=True)
        uvsub(vis=mmsfile, reverse=True)  # Now the viz data contains residuals + corrected models
        # Generate the corrected MS file
        if os.path.isdir(msfile_rot):
            shutil.rmtree(msfile_rot, ignore_errors=True)
        mstl.splitX(vis=mmsfile, outputvis=msfile_rot, datacolumn='corrected', datacolumn2='model_data')
        shutil.rmtree(mmsfile, ignore_errors=True)

        if do_diskslfcal:
            msfile_rot = disk_slfcal(msfile_rot, tbg, ted, disk_params, workdir)
        mmsfiles_rot.append(msfile_rot)
        if verbose:
            log_print('INFO', f"Corrected {tidx + 1} of {total_slots}: {mmsfile} → {msfile_rot}")
    return mmsfiles_rot


def split_mms(msname, timerange_series, spw='', workdir='./', overwrite=False, verbose=True):
    """
    Splits a measurement set into multiple subsets based on specified time ranges.

    :param msname: Path to the original measurement set.
    :type msname: str
    :param timerange_series: Time ranges for each split, with each tuple containing start and end datetime objects.
    :type timerange_series: list of tuple
    :param spw: Spectral window specification for the splitting process, defaults to ''.
    :type spw: str
    :param workdir: Target directory for saving the split measurement sets, defaults to './'.
    :type workdir: str
    :param overwrite: Boolean flag to overwrite existing files, defaults to False.
    :type overwrite: bool
    :param verbose: Boolean flag to enable verbose output, defaults to True.
    :type verbose: bool
    :return: Paths to the created split measurement sets, corresponding to each specified time range.
    :rtype: list of str
    """
    tseries_beg = timerange_series[0][0]
    tseries_end = timerange_series[-1][1]
    tseries_range = f"{tseries_beg.strftime('%Y/%m/%d/%H:%M:%S')}-{tseries_end.strftime('%Y/%m/%d/%H:%M:%S')}UT"
    mmsfiles = []
    total_slots = len(timerange_series)
    if verbose:
        log_print('INFO', f"Initiating split of {msname} at {tseries_range} into {total_slots} segments")

    for tidx, (start, end) in enumerate(timerange_series):
        timerange = f"{start.strftime('%Y/%m/%d/%H:%M:%S')}~{end.strftime('%Y/%m/%d/%H:%M:%S')}"
        outputvis = os.path.join(workdir, f"eovsa_{start.strftime('%H%M')}-{end.strftime('%H%M')}.ms")
        if overwrite:
            if os.path.isdir(outputvis):
                shutil.rmtree(outputvis, ignore_errors=True)
        else:
            if os.path.isdir(outputvis):
                mmsfiles.append(outputvis)
                if verbose:
                    log_print('INFO',
                              f"Sub-MS {tidx + 1}/{total_slots} created: {outputvis} (Time range: {timerange}) (already exists and was not re-created)")
                continue
        try:
            split(vis=msname, timerange=timerange, spw=spw, datacolumn='data', correlation='XX',
                  outputvis=outputvis)
            mmsfiles.append(outputvis)
            if verbose:
                log_print('INFO', f"Sub-MS {tidx + 1}/{total_slots} created: {outputvis} (Time range: {timerange})")
        except Exception as e:
            if verbose:
                log_print('ERROR', f"Failed to create sub-MS {tidx + 1} for {timerange} due to {e}")
            mmsfiles.append(None)
    if verbose:
        processed_count = sum(1 for f in mmsfiles if f is not None)
        log_print('INFO', f"Completed: {processed_count}/{total_slots} sub-MSs successfully processed")

    return mmsfiles


def all_paths_exist(paths):
    count = sum(1 for path in paths if os.path.exists(path))
    return count == len(paths)


def ant_trange(vis):
    ''' Figure out nominal times for tracking of old EOVSA antennas, and return time
        range in CASA format
    '''
    from eovsapy import eovsa_array as ea
    from astropy.time import Time

    # Get timerange from the visibility file
    # msinfo = dict.fromkeys(['vis', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs', 'observatory'])
    ms.open(vis)
    # metadata = ms.metadata()
    scans = ms.getscansummary()
    ms.done()
    sk = sorted(list(scans.keys()))
    vistrange = np.array([scans[sk[0]]['0']['BeginTime'], scans[sk[-1]]['0']['EndTime']])

    # Get the Sun transit time, based on the date in the vis file name (must have UDByyyymmdd in the name)
    aa = ea.eovsa_array()
    date = vis.split('UDB')[-1][:8]
    slashdate = date[:4] + '/' + date[4:6] + '/' + date[6:8]
    aa.date = slashdate
    sun = aa.cat['Sun']
    mjd_transit = Time(aa.next_transit(sun).datetime(), format='datetime').mjd
    # Construct timerange limits based on +/- 3h55m from transit time (when all dishes are nominally tracking)
    # and clip the visibility range not to exceed those limits
    mjdrange = np.clip(vistrange, mjd_transit - 0.1632, mjd_transit + 0.1632)
    trange = Time(mjdrange[0], format='mjd').iso[:19] + '~' + Time(mjdrange[1], format='mjd').iso[:19]
    trange = trange.replace('-', '/').replace(' ', '/')
    return trange


def format_spw(spw):
    """
    Formats the spectral window (spw) string for file naming, ensuring start and end values are separated by a dash and zero-padded to two digits.

    :param str spw: The spectral window in the format "start~end", where start and end are integers.
    :return: A formatted string with start and end values zero-padded to two digits, separated by a dash.
    :rtype: str
    """
    return '-'.join(['{:02d}'.format(int(sp_)) for sp_ in spw.split('~')])


def rm_imname_extensions(imname, keep_ext=[], verbose=False):
    """
    Remove directories and files that match the base image name (`imname`) with specific extensions.

    :param imname: The base name of images/directories to remove specific extensions for.
    :type imname: str
    :param keep_ext: A list of extensions to keep. Default is an empty list.
    :type keep_ext: list
    :param verbose: If True, print detailed information about the removal process. Default is False.
    :type verbose: bool
    """
    # List of known extensions to match and remove
    extensions = ['.psf', '.pb', '.residual', '.mask', '.prev.mask', '.sumwt', '.model', '.image', '.image.pbcor']
    for ext in keep_ext:
        extensions.remove(ext)

    # Loop through each extension and remove the matching directories or files
    for ext in extensions:
        # Construct the full pattern for each extension, including directories and files
        pattern_dir = f"{imname}{ext}/"  # Pattern for directories
        pattern_file = f"{imname}{ext}"  # Pattern for files

        # Find all matches for the current extension
        matching_dirs = glob(pattern_dir)
        matching_files = glob(pattern_file)

        # Remove each matched directory
        for dir_path in matching_dirs:
            if os.path.isdir(dir_path):  # Double-check if it's a directory
                if verbose: log_print('INFO', f"Removing directory: {dir_path}")
                shutil.rmtree(dir_path)
            else:
                if verbose: log_print('INFO', f"Expected a directory but found none: {dir_path}")

        # Remove each matched file
        for file_path in matching_files:
            if os.path.isfile(file_path):  # Double-check if it's a file
                if verbose: log_print('INFO', f"Removing file: {file_path}")
                os.remove(file_path)
            else:
                if verbose: log_print('INFO', f"Expected a file but found none: {file_path}")


def check_image_zeros(imname):
    """
    Check if the image contains non-zero data.
    """
    ia = iatool()
    ia.open(imname)
    data = ia.getchunk(getmask=True)
    ia.close()
    return np.count_nonzero(data) <= 1


def format_param(param):
    if isinstance(param, str):
        return f'"{param}"'
    else:
        return str(param)


def run_tclean_automasking(vis, sp, trange, uvrange, datacolumn, imname,
                           imsize, cell, stokes, scales, niter, reffreq,
                           pbcor,
                           savemodel,
                           usemask,
                           restoringbeam,
                           sidelobethreshold,
                           noisethreshold,
                           lownoisethreshold,
                           minbeamfrac,
                           negativethreshold,
                           growiterations):
    """
    Wrapper function for the tclean task in CASA.

    :param str vis: Name of the visibility data set to be imaged.
    :param str sp: Spectral window selection.
    :param str trange: Time range selection.
    :param str uvrange: UV range selection.
    :param str datacolumn: Specifies which data column to use.
    :param str imname: Name of the output image.
    :param list imsize: Image size in pixels.
    :param list cell: Cell size specification.
    :param str stokes: Stokes parameters to image.
    :param list scales: Scales for multi-scale clean.
    :param int niter: Maximum number of iterations.
    :param float reffreq: Reference frequency for the image.
    :param bool pbcor: Specifies whether to perform primary beam correction.
    :param str savemodel: Specifies whether to save the model.
    :param bool usemask: Specifies whether to use auto-masking.
    :param list restoringbeam: Restoring beam parameters.
    :param float sidelobethreshold: Side lobe threshold for auto-masking.
    :param float noisethreshold: Noise threshold for auto-masking.
    :param float lownoisethreshold: Low noise threshold for auto-masking.
    :param float negativethreshold: Negative threshold for auto-masking.
    :param int growiterations: Number of grow iterations for auto-masking.
    """
    params_str = ', '.join([
        f'vis={format_param(vis)}',
        f'spw={format_param(sp)}',
        f'timerange={format_param(trange)}',
        f'uvrange={format_param(uvrange)}',
        f'datacolumn={format_param(datacolumn)}',
        f'imagename={format_param(imname)}',
        f'imsize={format_param(imsize)}',
        f'cell={format_param(cell)}',
        f'stokes={format_param(stokes)}',
        f'scales={format_param(scales)}',
        f'niter={format_param(niter)}',
        f'reffreq={format_param(reffreq)}',
        f'pbcor={format_param(pbcor)}',
        f'savemodel={format_param(savemodel)}',
        f'usemask={format_param(usemask)}',
        f'restoringbeam={format_param(restoringbeam)}',
        f'sidelobethreshold={format_param(sidelobethreshold)}',
        f'noisethreshold={format_param(noisethreshold)}',
        f'lownoisethreshold={format_param(lownoisethreshold)}',
        f'minbeamfrac={format_param(minbeamfrac)}',
        f'negativethreshold={format_param(negativethreshold)}',
        f'growiterations={format_param(growiterations)}'
    ])

    log_print('INFO', f'Running tclean with parameters: {params_str}')

    tclean(vis=vis, selectdata=True, spw=sp, timerange=trange,
           uvrange=uvrange, antenna="0~12&0~12",
           datacolumn=datacolumn, imagename=imname,
           imsize=imsize, cell=cell,
           stokes=stokes, projection="SIN", specmode="mfs",
           interpolation="linear", deconvolver="multiscale",
           scales=scales, nterms=2, smallscalebias=0.6,
           restoration=True, weighting="briggs",
           robust=0.0, niter=niter, gain=0.05,
           reffreq=reffreq,
           savemodel=savemodel,
           pbcor=pbcor,
           veltype='radio',
           outframe='TOPO',
           restoringbeam=restoringbeam,
           usemask=usemask, pbmask=0.0,
           sidelobethreshold=sidelobethreshold,
           noisethreshold=noisethreshold,
           lownoisethreshold=lownoisethreshold,
           minbeamfrac=minbeamfrac,
           negativethreshold=negativethreshold,
           smoothfactor=1.0, cutthreshold=0.01,
           growiterations=growiterations, dogrowprune=True,
           minpercentchange=-1.0)


class FrequencySetup:
    """
    Manages frequency setup based on observation date for radio astronomy imaging.

    This class is initialized with an observation time and calculates
    essential frequency parameters such as effective observing frequencies (eofreq)
    and spectral windows (spws) based on the observation date. It provides methods
    to calculate reference frequency and bandwidth for given spectral windows.

    :param tim: Observation time used to determine the frequency setup.
    :type tim: astropy.time.Time

    Attributes:
    - tim (astropy.time.Time): The observation time.
    - spw2band (numpy.ndarray): An array mapping spectral window indices to band numbers.
    - bandwidth (float): The bandwidth in GHz.
    - defaultfreq (numpy.ndarray): The default effective observing frequencies in GHz.
    - nbands (int): Number of bands.
    - eofreq (numpy.ndarray): Effective observing frequencies based on the observation time.
    - spws (list of str): Spectral window selections based on the observation time.

    :Example:

    >>> from astropy.time import Time
    >>> tim = Time('2022-01-01T00:00:00', format='isot')
    >>> freq_setup = FrequencySetup(tim)
    >>> crval, cdelt = freq_setup.get_reffreq_and_cdelt('5~10')
    >>> print(crval, cdelt)
    """

    def __init__(self, tim=None):
        if tim is None:
            tim = Time.now()
        self.tim = tim
        self.spw2band = np.array([0, 1] + list(range(4, 52)))
        self.bandwidth = 0.325  # 325 MHz
        self.defaultfreq = 1.1 + self.bandwidth * (self.spw2band + 0.5)

        if self.tim.mjd > 58536:
            self.nbands = 52
            self.eofreq = self.defaultfreq
            self.spws = ['0~1', '2~4', '5~10', '11~20', '21~30', '31~40', '41~49']
        else:
            self.bandwidth = 0.5  # 500 MHz
            self.nbands = 34
            self.eofreq = 1.419 + np.arange(self.nbands) * self.bandwidth
            self.spws = ['1~3', '4~9', '10~16', '17~24', '25~30']

    def get_reffreq_and_cdelt(self, spw):
        """
        Calculates the reference frequency (CRVAL) and the frequency delta (CDELT)
        for a given spectral window range.

        This method takes a spectral window selection and computes the mean of the effective
        observing frequencies (eofreq) within that range as the reference frequency. It also
        calculates the bandwidth covered by the spectral window range as the frequency delta.

        :param spw: Spectral window selection, specified as a range 'start~end' or a single value.
        :type spw: str

        :return: A tuple containing the reference frequency and frequency delta, both in GHz.
        :rtype: (str, str)

        :Example:

        >>> crval, cdelt = freq_setup.get_reffreq_and_cdelt('5~10')
        >>> print(crval, cdelt)
        """
        sp_st, sp_ed = (int(s) for s in spw.split('~')) if '~' in spw else (int(spw), int(spw))
        crval = f'{np.mean(self.eofreq[sp_st:sp_ed + 1]):.4f}GHz'
        nband = (self.eofreq[sp_ed] - self.eofreq[sp_st]) / self.bandwidth + 1
        cdelt = f'{nband * self.bandwidth:.4f}GHz'
        return crval, cdelt


# def get_eofreq(vis):
#     spw2band = np.array([0, 1] + list(range(4, 52)))
#     bandwidth = 0.325  ## 325 MHz
#     defaultfreq = 1.1 + bandwidth * (spw2band + 0.5)
#     if mstl.get_trange(vis)[0].mjd > 58536:
#         # After 2019 Feb 22, the band numbers changed to 1-52, and spw from 0-49
#         nbands = 52
#         eofreq = defaultfreq
#         spws = ['0~1', '2~4', '5~10', '11~20', '21~30', '31~40', '41~49']
#     else:
#         # Before 2019 Feb 22, the band numbers were 1-34, and spw from 0-30
#         bandwidth = 0.5  ## 500 MHz
#         nbands = 34
#         eofreq = 1.419 + np.arange(nbands) * bandwidth
#         spws = ['1~3', '4~9', '10~16', '17~24', '25~30']
#     return eofreq, spws, bandwidth
#
#
# def get_reffreq(spw, reffreqs, bandwidth):
#     if '~' in spw:
#         sp_st, sp_ed = [int(i) for i in spw.split('~')]
#     else:
#         sp_st = sp_ed = int(spw)
#     nband = (reffreqs[sp_ed] - reffreqs[sp_st]) / bandwidth + 1
#     crval = f'{np.mean(reffreqs[sp_st:sp_ed + 1]):.4f}GHz'
#     cdelt = f'{nband * bandwidth:.4f}GHz'
#     return crval, cdelt


def fd_images(vis,
              cleanup=False,
              image_marker='',
              timerange='',
              niter=None,
              cell=['2.5arcsec'],
              imsize=[1024],
              spws=['0~1', '2~4', '5~10', '11~20', '21~30', '31~43'],
              imgoutdir='./',
              bright=None,
              stokes="XX",
              uvrange='',
              toTb=True,
              pbcor=True,
              datacolumn='data',
              savemodel="none",
              usemask='auto-multithresh',
              overwrite=False,
              compress=False,
              dryrun=False):
    """
    Generates full-disk images, optionally cleans up interim images, and performs image registration.

    This function creates full-disk solar images based on visibility data, performs an optional cleanup of interim images,
    and aligns the resulting images to a standard solar disk model. It supports dynamic image size, cell size, and spectral
    window (SPW) selection. The function also handles the creation of FITS files from the generated images, ensuring the reference
    frequency in the FITS header is calculated as the middle of the selected frequency range.

    :param str vis: Path to the visibility data.
    :param bool cleanup: If True, deletes interim images after processing. Default is False.
    :param str image_marker: Additional identifier for the image name. Default is ''.
    :param str timerange: Range of time to select from data. Default is ''.
    :param int niter: Number of iterations for the cleaning algorithm. If None, defaults to 5000.
    :param list cell: Size of the image cell, e.g., ['2.5arcsec'].
    :param list imsize: Dimensions of the output image, e.g., [1024].
    :param list spws: Spectral windows to process, specified as ranges, e.g., ['0~1', '2~5'].
    :param str imgoutdir: Output directory for the resulting images and FITS files. Defaults to './'.
    :param list bright: Flags to indicate brightness processing per SPW. Defaults to all True.
    :param str stokes: Stokes parameter to use. Default is "XX".
    :param str uvrange: UV range to select from data. Default is ''.
    :param bool toTb: If True, converts image to temperature scale. Default is True.
    :param bool pbcor: If True, applies primary beam correction. Default is True.
    :param str datacolumn: Data column to use from the visibility data. Default is 'data'.
    :param str savemodel: Options for saving the model visibilities. Choices are 'none', 'virtual', and 'modelcolumn'. Default is 'none'.
    :param bool overwrite: If True, overwrites existing FITS files. Default is False.
    :param bool dryrun: If True, only returns the paths to the generated images. No actual image processing is performed. Default is False.

    :return: A tuple containing two lists: paths to the generated FITS files and paths to the CASA image files.
    :rtype: (list of str, list of str)

    .. note::
       1. The reference frequency in the FITS header is calculated as the middle of the selected frequency range. This ensures
       accurate representation of the frequency information in the generated images.
       2. The function internally manages directory creation for images, applies multiscale cleaning based on initial beam size
       determination, and uses default or specified scales for image cleaning. Errors in beam size determination lead to fallback
       on default scales. The function finally aligns and converts the images to FITS format, with options for temperature scale
       conversion and phase center alignment.
    """
    # Check if "images" directory exists (if not, create it and mark it for later deletion)
    if pbcor:
        imext = '.image.pbcor'
    else:
        imext = '.image'
    imgtmpdir = os.path.join(imgoutdir, 'images')
    try:
        if os.stat(imgtmpdir):
            rm_images = False  # Mark as not removeable
    except:
        os.mkdir(imgtmpdir)
        if cleanup:
            rm_images = True  # Mark as removeable
        else:
            rm_images = False  # Mark as not removeable

    if timerange == '':
        trange = ant_trange(vis)
        (tstart, tend) = timerange.split('~')
        tbg = Time(qa.quantity(tstart, 's')['value'], format='mjd').to_datetime()
        fitsname_prefix = f'eovsa_{tbg.strftime("%Y%m%d")}'
    else:
        trange = timerange
        (tstart, tend) = timerange.split('~')
        tbg = Time(qa.quantity(tstart, 's')['value'], format='mjd').to_datetime()
        ted = Time(qa.quantity(tend, 's')['value'], format='mjd').to_datetime()
        tdt = ted - tbg
        if tdt.total_seconds() >= 60:
            tdtstr = f'{int(tdt.total_seconds() / 60)}min'
        else:
            tdtstr = f'{int(tdt.total_seconds())}sec'
        fitsname_prefix = f'eovsa_{tbg.strftime("%Y%m%dT%H%M%S")}_{tdtstr}'

    freq_setup = FrequencySetup(Time(tbg))

    # eofreq, spws, bandwidth = get_eofreq(vis)
    if image_marker != '':
        image_marker = image_marker.lstrip('.')
        image_marker = f'.{image_marker}'
    if toTb:
        fitsname_suffix = '.tb.fits'
    else:
        fitsname_suffix = '.fits'
    if niter is None:
        niter = 5000
    if bright is None:
        bright = [True] * len(spws)
    imagefile = []
    fitsfile = []
    cellvalue = np.float_(cell[0].replace('arcsec', ''))
    for s, sp in enumerate(spws):
        if bright[s]:
            spwstr = format_spw(sp)
            imname = os.path.join(imgtmpdir, f"{fitsname_prefix}.s{spwstr}{image_marker}")
            outfits = os.path.join(imgoutdir, f"{fitsname_prefix}.s{spwstr}{image_marker}{fitsname_suffix}")
            if overwrite:
                if os.path.exists(outfits):
                    os.remove(outfits)
            else:
                if os.path.exists(outfits):
                    imagefile.append(imname + imext)
                    fitsfile.append(outfits)
                    log_print('INFO',
                              f'Processing tclean for SPW {spwstr} -- (timerange: {timerange}) image_marker: {image_marker} -- {outfits} already exists. Skipping.')
                    continue
            if dryrun:
                imagefile.append(imname + imext)
                fitsfile.append(outfits)
                continue
            log_print('INFO',
                      f'Processing tclean for SPW {spwstr} -- (timerange: {timerange}) image_marker: {image_marker}')
            if os.path.exists(imname + imext):
                rm_imname_extensions(imname)
            # try:
            #     if os.path.isdir(imname + '.init.image'):
            #         ## if the init.image already exists, get the beam size from the existing image
            #         bmaj, bmin, bpa, beamunit, bpaunit = hf.getbeam(imagefile=[imname + '.init.image'])
            #     else:
            #         ## use multiscale to make high quality images
            #         ## do an initial hogbom clean with niter=0 to determine the beam size
            #         tclean(vis=vis, selectdata=True, spw=sp, timerange=trange,
            #                uvrange='',
            #                antenna="0~12&0~12",
            #                datacolumn=datacolumn, imagename=imname, imsize=imsize, cell=cell,
            #                stokes=stokes, projection="SIN", specmode="mfs", interpolation="linear",
            #                deconvolver="hogbom",
            #                weighting="briggs", robust=0,
            #                niter=0, gain=0.05, savemodel="none")
            #         bmaj, bmin, bpa, beamunit, bpaunit = hf.getbeam(imagefile=[imname + '.image'])
            #     bmsz = np.nanmean([bmaj, bmin])
            #     if bmsz < 5.0:
            #         bmsz = 5.0
            #     scalesfactor = [0, 1, 3]  ## [ 0, 1xbeam, 3xbeam]
            #     scales = [np.ceil(l * bmsz / cellvalue) for l in scalesfactor]
            #     log_print('INFO',
            #               f'Applying scales {scales} for multiscale clean, determined by initial beam size assessment.')
            # except:
            #     scales = [0, 5, 15]
            #     log_print('WARNING:',
            #               f'Unable to determine beam size from initial clean. Applying default scales {scales} for multiscale clean.')
            reffreq, cdelt4_real = freq_setup.get_reffreq_and_cdelt(sp)
            bmsz = get_bmsize(float(reffreq.rstrip('GHz')), minbmsize=5.0)
            scalesfactor = [0, 1, 3]  ## [ 0, 1xbeam, 3xbeam]
            scales = [np.ceil(l * bmsz / cellvalue) for l in scalesfactor]
            log_print('INFO',
                      f'Applying scales {scales} for multiscale clean, determined by pre-defined beam size.')
            restoringbeam = [f'{bmsz:.3f}arcsec']
            try:
                initial_params = {'sidelobethreshold': 1.5, 'noisethreshold': 2.5, 'lownoisethreshold': 1.5,
                                  'minbeamfrac': 0.3,
                                  'negativethreshold': 5.0,
                                  'growiterations': 75}
                secondary_params = {'sidelobethreshold': 1.0, 'noisethreshold': 2.0, 'lownoisethreshold': 1.0,
                                    'minbeamfrac': 0.2,
                                    'negativethreshold': 3.0,
                                    'growiterations': 75}
                final_params = {'sidelobethreshold': 0.5, 'noisethreshold': 1.0, 'lownoisethreshold': 0.5,
                                'minbeamfrac': 0.2,
                                'negativethreshold': 0.0,
                                'growiterations': 75}

                log_print('INFO', 'Initial tclean run')
                rm_imname_extensions(imname)
                # if len(eofreq) > 0:
                run_tclean_automasking(vis, sp, trange, uvrange, datacolumn, imname, imsize, cell, stokes, scales,
                                       niter, reffreq,
                                       pbcor,
                                       savemodel,
                                       usemask,
                                       restoringbeam,
                                       **initial_params)

                if check_image_zeros(imname + imext):
                    log_print('INFO', 'Second tclean run with adjusted parameters')
                    rm_imname_extensions(imname)
                    run_tclean_automasking(vis, sp, trange, uvrange, datacolumn, imname, imsize, cell, stokes, scales,
                                           niter, reffreq,
                                           pbcor,
                                           savemodel,
                                           usemask,
                                           restoringbeam,
                                           **secondary_params)

                if check_image_zeros(imname + imext):
                    log_print('INFO', 'Final tclean run with further adjusted parameters')
                    rm_imname_extensions(imname)
                    run_tclean_automasking(vis, sp, trange, uvrange, datacolumn, imname, imsize, cell, stokes, scales,
                                           niter, reffreq,
                                           pbcor,
                                           usemask,
                                           restoringbeam,
                                           savemodel, **final_params)
                if check_image_zeros(imname + imext):
                    log_print('INFO', 'tclean run failed with automasking. Trying without automasking...')
                    rm_imname_extensions(imname)
                    tclean(vis=vis, selectdata=True, spw=sp, timerange=trange,
                           uvrange=uvrange, antenna="0~12&0~12",
                           datacolumn=datacolumn, imagename=imname,
                           imsize=imsize, cell=cell,
                           stokes=stokes, projection="SIN", specmode="mfs",
                           interpolation="linear", deconvolver="multiscale",
                           pbcor=True,
                           restoration=True, weighting="briggs",
                           robust=0.0, niter=niter, gain=0.05, savemodel=savemodel,
                           usemask='user', pbmask=0.0)
                if check_image_zeros(imname + imext):
                    log_print('WARNING',
                              f'Final tclean run for SPW {spwstr} produced an image with all zeros. Skipping...')
                    imagefile.append(None)
                    fitsfile.append(None)

                ## the cdelt value of the frequency axis is not correct in the image header, so we need to correct it.
                ctype4 = imhead(imname + imext, mode='get', hdkey='ctype4')
                if ctype4 == 'Frequency':
                    cdelt4 = qa.convert(imhead(imname + imext, mode='get', hdkey='cdelt4'), "GHz")
                    if cdelt4_real != '':
                        log_print('INFO', f"Correcting cdelt4 from {cdelt4['value']:.4f}GHz to {cdelt4_real}")
                        cdelt4_new = qa.quantity(cdelt4_real)
                        imhead(imname + imext, mode='put', hdkey='cdelt4', hdvalue=cdelt4_new)
                if os.path.exists(outfits):
                    os.remove(outfits)

                imagefile.append(imname + imext)
                fitsfile.append(outfits)
            except Exception as e:
                log_print('ERROR', f'Failed to generate image for SPW {spwstr} due to {e}')
                imagefile.append(None)
                fitsfile.append(None)
    if dryrun:
        return fitsfile, imagefile
    # Check if all elements in imagefile are None
    if all(item is None for item in imagefile):
        log_print('WARNING', "All elements in imagefile are None. Skipping imreg.")
    else:
        # print(
        #     f"vis={vis}, "
        #     f"imagefile={imagefile}, "
        #     f"fitsfile={fitsfile}, "
        #     f"timerange={[trange] * len(fitsfile)}, "
        #     f"toTb={toTb}, "
        #     f"usephacenter={False}, "
        #     f"overwrite={overwrite}, "
        #     f"docompress={compress}"
        # )
        hf.imreg(vis=vis, imagefile=imagefile, fitsfile=fitsfile, timerange=[trange] * len(fitsfile), toTb=toTb,
                 usephacenter=False, overwrite=overwrite, docompress=compress)
        if image_marker == '.init':
            rm_imname_extensions(imname, keep_ext=[imext, '.model'])
        else:
            rm_imname_extensions(imname)
    if rm_images:
        shutil.rmtree(imgtmpdir)  # Remove all images and the folder named images
    return fitsfile, imagefile


def merge_FITSfiles(fitsfilesin, outfits, exptime_weight=False, suppress_ondiskres=False, suppress_thrshd=0.3,
                    overwrite=True):
    """
    Merges multiple FITS files into a single output file by calculating the mean of stacked data.

    This function stacks the data from a list of input FITS files along a new axis, optionally applying exposure time
    weighting and suppression of on-disk residuals based on a threshold relative to the disk brightness temperature.
    The mean of the stacked data is then written to a specified output FITS file.

    Parameters
    ----------
    fitsfilesin : list of str
        List of file paths to the input FITS files to be merged.
    outfits : str
        File path for the output FITS file where the merged data will be saved.
    exptime_weight : bool, optional
        If True, the exposure time is used as a weight for calculating the mean of the stacked data. Defaults to False.
    suppress_ondiskres : bool, optional
        If True, suppresses on-disk residuals based on `suppress_thrshd`. Defaults to False.
    suppress_thrshd : float, optional
        Threshold for suppression of on-disk residuals, expressed as a fraction of the disk brightness temperature.
        May suppress weak but real emission. Defaults to 0.3.
    overwrite : bool, optional
        If True, the output file is overwritten if it already exists. Defaults to True.

    Raises
    ------
    ValueError
        If any of the input FITS files cannot be read.

    Returns
    -------
    None

    Notes
    -----
    The suppression of on-disk residuals is achieved by applying a sigmoid function to pixels within a specified range
    of the disk brightness temperature, effectively reducing the contribution of residuals while preserving the overall
    structure.

    Examples
    --------
    Merge three FITS files without exposure time weighting and with on-disk residual suppression:

    >>> merge_FITSfiles(['sun_01.fits', 'sun_02.fits', 'sun_03.fits'], 'sun_merged.fits',
    ...                 exptime_weight=False, suppress_ondiskres=True, suppress_thrshd=0.3, overwrite=True)
    """
    data_stack = []
    exptimes = []
    for file in fitsfilesin:
        meta_, data_ = ndfits.read(file)
        if data_ is not None:
            data = data_
            meta = meta_
            exptimes.append(meta['header']['EXPTIME'])
            data_stack.append(np.squeeze(data))
    data_stack = np.dstack(data_stack)
    mdata_stack = ma.masked_invalid(data_stack)
    if suppress_ondiskres:
        import copy
        from sunpy.map.maputils import all_coordinates_from_map, coordinate_is_on_solar_disk
        def sigmoid(x):
            return 1 / (1 + np.exp(-x))

        mdata_stack_ = copy.deepcopy(mdata_stack)
        tbdisk = meta['header']['tbdisk']
        # Step 1: Identify the pixels in the range tbdisk to 1.1*tbdisk
        threshd = suppress_thrshd  # 30% of the disk brightness temperature
        mask = np.logical_and(mdata_stack_ >= tbdisk, mdata_stack_ <= (1 + threshd) * tbdisk)
        # Step 2: Calculate the suppression factor
        # Define the transition boundaries
        lower_bound = tbdisk
        upper_bound = (1 + threshd) * tbdisk
        # Normalize the data to the range [0, 1]
        normalized_data = (mdata_stack_ - lower_bound) / (upper_bound - lower_bound)
        # Apply the sigmoid function to the normalized data
        suppression_factor = sigmoid(10 * (normalized_data - 0.5))
        # Step 3: Apply the suppression factor to the pixels
        mdata_stack_[mask] *= suppression_factor[mask]
        eomap = meta['refmap']
        hpc_coords = all_coordinates_from_map(eomap)
        mask_disk = (coordinate_is_on_solar_disk(hpc_coords))
        # Expand mask_disk along the third axis to match the shape of mdata_stack
        mask_disk_expanded = np.repeat(mask_disk[..., np.newaxis], mdata_stack.shape[2], axis=2)
        mdata_stack[mask_disk_expanded] = mdata_stack_[mask_disk_expanded]

    if exptime_weight:
        ## use the exposure time as the weight to calculate the mean
        weights = np.array(exptimes)
        # Adjust the weights of the first and last images by halving them to mitigate edge effects.
        weights[[0, -1]] /= 2
        weights = weights / np.sum(weights)
        weights = weights[np.newaxis, np.newaxis, :]
        data_stack = np.nansum(mdata_stack * weights, axis=-1)
    else:
        data_stack = np.nanmean(mdata_stack, axis=-1)
    meta['header']['EXPTIME'] = np.nansum(exptimes)
    date_obs = Time(meta['header']['date-obs']).datetime.replace(hour=20, minute=0, second=0) - timedelta(
        seconds=meta['header']['EXPTIME'] / 2.0)
    meta['header']['date-obs'] = date_obs.strftime('%Y-%m-%dT%H:%M:%S.%f')
    ndfits.write(outfits, data_stack, meta['header'], overwrite=overwrite)
    log_print('INFO', f'Merging {outfits}')


from multiprocessing import Pool
from functools import partial


def process_time_block(tidx_ted_tbg, msfile_in, msname, subdir, total_blocks, tdt, tdtstr, spws, niter_init,
                       reftime_master, do_diskslfcal, disk_params, pols='XX', do_sbdcal=False, overwrite=False):
    tidx, (tbg, ted) = tidx_ted_tbg
    timerange = trange2timerange([tbg, ted])

    if isinstance(msfile_in, list):
        msfile = msfile_in[tidx]
        log_print('INFO',
                  f"msfile_in is a list. Processing time block {tidx + 1} of {total_blocks} (Time range: {timerange}) ... ")
    else:
        msfile = msfile_in
        log_print('INFO',
                  f"msfile_in is a file. Processing time block {tidx + 1} of {total_blocks} (Time range: {timerange}) ... ")
    combined_vis_sub = os.path.join(subdir, f'{msname}_shift_corrected.block{tidx + 1}.seg{tdtstr}.ms')

    tr_series = generate_trange_series(tbg, ted, tdt)
    mmsfiles = split_mms(msfile, tr_series, workdir=subdir, overwrite=overwrite, verbose=True)

    log_print('INFO', f"Dry run to check if the final FITS images exist for {timerange} ...")
    fitsfile, imagefile = fd_images(msfile,
                                    timerange=timerange,
                                    cleanup=False,
                                    niter=5000, spws=spws,
                                    pbcor=False,
                                    imgoutdir=subdir,
                                    dryrun=True)
    if all_paths_exist(fitsfile):
        log_print('INFO', f"Final FITS images exist for {timerange}. Skip this block.")
        return None  # Skip processing if images exist

    fits_ref, img_ref = fd_images(msfile,
                                  timerange=timerange,
                                  cleanup=False,
                                  overwrite=overwrite,
                                  pbcor=False,
                                  image_marker='init',
                                  niter=niter_init, spws=spws,
                                  imgoutdir=subdir)
    img_ref_name = [l.rstrip('.image') for l in img_ref]
    img_ref_model = [f'{l}.model' for l in img_ref_name]
    img_ref_model_fits = [f'{l}.model.fits' for l in img_ref_name]

    hf.imreg(vis=msfile, imagefile=img_ref, fitsfile=fits_ref, timerange=[timerange] * len(spws), toTb=False,
             usephacenter=False, overwrite=overwrite, docompress=True)

    hf.imreg(vis=msfile, imagefile=img_ref_model, fitsfile=img_ref_model_fits,
             timerange=[timerange] * len(spws),
             toTb=False, usephacenter=False, overwrite=overwrite, docompress=True)

    mmsfiles_rot = shift_corr(mmsfiles, tr_series, spws, img_ref_model, img_ref_model_fits, reftime_master,
                              workdir=subdir, do_featureslfcal=True, pols=pols,
                              overwrite=overwrite, do_diskslfcal=do_diskslfcal, disk_params=disk_params,
                              do_sbdcal=do_sbdcal,
                              verbose=True)

    if os.path.isdir(combined_vis_sub):
        shutil.rmtree(combined_vis_sub, ignore_errors=True)
    concat(vis=[l for l in mmsfiles_rot if l is not None], concatvis=combined_vis_sub)

    log_print('INFO', f"Time block {tidx + 1} of {total_blocks} completed.")
    return combined_vis_sub


def process_imaging_timerange(tbg_ted, msfile_in, spws, subdir, overwrite):
    tidx, (tbg, ted) = tbg_ted
    if isinstance(msfile_in, list):
        ## this is for running the code on pipeline server
        msfile = msfile_in[tidx]
    else:
        msfile = msfile_in
    if msfile is None:
        return None, None
    timerange = trange2timerange((tbg, ted))
    fitsfile, imagefile = fd_images(msfile,
                                    timerange=timerange,
                                    pbcor=False,
                                    overwrite=overwrite,
                                    cleanup=False,
                                    niter=5000, spws=spws,
                                    imgoutdir=subdir,
                                    compress=True)
    return fitsfile, imagefile


def pipeline_run(vis, outputvis='', workdir=None, slfcaltbdir=None, imgoutdir=None, figoutdir=None, clearcache=False,
                 pols='XX', mergeFITSonly=False, verbose=True, do_diskslfcal=True, overwrite=False, niter_init=200,
                 ncpu='auto',
                 tr_series_imaging=None,
                 spws_imaging=None, hanning=False, do_sbdcal=False):
    """
    Executes the EOVSA data processing pipeline for solar observation data.

    :param vis: Path to the input measurement set (MS) data.
    :type vis: str
    :param outputvis: Output path for the processed MS, defaults to an empty string.
    :type outputvis: str, optional
    :param workdir: Working directory for intermediate data, defaults to None which sets it to '/data1/workdir'.
    :type workdir: str, optional
    :param slfcaltbdir: Directory for storing calibration tables, defaults to None.
    :type slfcaltbdir: str, optional
    :param imgoutdir: Output directory for image files, defaults to None.
    :type imgoutdir: str, optional
    :param figoutdir: Output directory for figures, defaults to None.
    :type figoutdir: str, optional
    :param clearcache: If True, clears the cache after processing, defaults to False.
    :type clearcache: bool, optional
    :param pols: Polarization types to process, defaults to 'XX'.
    :type pols: str, optional
    :param mergeFITSonly: If True, skips processing and only merges FITS files, defaults to False.
    :type mergeFITSonly: bool, optional
    :param verbose: Enables verbose output during processing, defaults to True.
    :type verbose: bool, optional
    :param do_diskslfcal: If True, performs disk self-calibration, defaults to True.
    :type do_diskslfcal: bool, optional
    :param overwrite: If True, overwrites existing files, defaults to False.
    :type overwrite: bool, optional
    :param niter_init: Initial number of iterations for imaging, defaults to 200.
    :type niter_init: int, optional
    :param ncpu: Specifies the number of CPUs for parallel processing, defaults to 'auto'.
    :type ncpu: str or int, optional
    :param tr_series_imaging: Time ranges for imaging, defaults to None.
    :type tr_series_imaging: list of tuple, optional
    :param spws_imaging: Spectral windows selected for imaging, defaults to None.
    :type spws_imaging: list of str, optional
    :param hanning: If True, applies Hanning smoothing to the data, defaults to False.
    :type hanning: bool, optional
    :param do_sbdcal: Boolean flag to perform single-band delay calibration, defaults to False.
    :type do_sbdcal: bool
    :return: Path to the processed visibility data.
    :rtype: str

    :example:
    ## if you want to specify the spectral windows for imaging
    >>> spws_imaging = ['5~10', '11~20', '21~30']
    ## if you want to specify the time range for imaging
    >>> from datetime import datetime, timedelta
    >>> from suncasa.eovsa import eovsa_synoptic_imaging_pipeline as esip
    >>> tbg_imaging = datetime(2024, 4, 8, 16, 58, 0)
    >>> ted_imaging = datetime(2024, 4, 8, 19, 00, 0)
    >>> tdt_imaging = timedelta(minutes=2)
    >>> tr_series_imaging = esip.generate_trange_series(tbg_imaging, ted_imaging, tdt_imaging)
    """

    msfile = vis.rstrip('/')
    if workdir is None:
        workdir = './'
    os.chdir(workdir)
    if slfcaltbdir is None:
        slfcaltbdir = workdir + '/'
    if imgoutdir is None:
        imgoutdir = workdir + '/'
    ## todo figoutdir is obsolete. will be removed once the pipeline is stable.
    if figoutdir is None:
        figoutdir = workdir + '/'

    ## the msfile we use 60 min model image to correct the data in 10 min interval. the model image is shifted to the reftime_master (20:00 UT of each day).
    tdt_master = timedelta(hours=1)
    # tdt_master = timedelta(minutes=30)
    tdt = timedelta(minutes=10)
    viz_timerange = ant_trange(msfile)
    (tstart, tend) = viz_timerange.split('~')
    tbg_master = Time(qa.quantity(tstart, 's')['value'], format='mjd').to_datetime()
    ted_master = Time(qa.quantity(tend, 's')['value'], format='mjd').to_datetime()
    reftime_master = Time(datetime.combine(tbg_master.date(), time(20, 0)))
    tr_series_master = generate_trange_series(tbg_master, ted_master, tdt_master, snap_to_full_hour=True)
    tdtmststr = f'{int(tdt_master.total_seconds() / 60)}min'
    tdtstr = f'{int(tdt.total_seconds() / 60)}min'
    date_str = tbg_master.strftime('%Y%m%d')
    tdtmst_str = f'{int(tdt_master.total_seconds() / 60)}m'
    subdir = tbg_master.strftime('temp_%Y%m%d')
    freq_setup = FrequencySetup(Time(tbg_master))
    spws = freq_setup.spws
    defaultfreq = freq_setup.defaultfreq
    freq = defaultfreq
    nbands = freq_setup.nbands
    slashdate = tbg_master.strftime('%Y/%m/%d')
    dsize, fdens = calc_diskmodel(slashdate, nbands, freq, defaultfreq)

    if os.path.isdir(subdir) == False:
        os.makedirs(subdir)

    msname, _ = os.path.splitext(msfile)
    msname = os.path.basename(msname)

    msfile_copy = os.path.join(subdir, f'{msname}.ms')
    if not os.path.exists(msfile_copy):
        shutil.copytree(msfile, msfile_copy)
    msfile = msfile_copy
    diskxmlfile = msfile + '.SOLDISK.xml'
    disk_params = {'dsize': dsize, 'fdens': fdens, 'freq': freq, 'diskxmlfile': diskxmlfile}

    combined_vis = os.path.join(subdir, f'{msname}_shift_corrected.b{tdtmststr}.s{tdtstr}.ms')
    if outputvis == '':
        outputvis = os.path.join(workdir, f'{msname}.b{tdtmststr}.s{tdtstr}.shift_corr.ms')
    outputvis = outputvis.rstrip('/')
    # combined_vis_disk = os.path.join(subdir, f'{msname}_shift_corrected.b{tdtmststr}.s{tdtstr}.disk.ms')

    ### --------------------------------------------------------------###
    ## pre-processing block. flagging, hanning smoothing...
    ### --------------------------------------------------------------###
    run_start_time_pre_proc = datetime.now()
    if os.path.isdir(msfile + '.flagversions') == True:
        if verbose:
            log_print('INFO', f'Flagversion of {msfile} already exists. Skipped...')
        # flagmanager(msfile, mode='restore', versionname='pipeline_init')
        # shutil.rmtree(msfile + '.flagversions',ignore_errors=True)
    else:
        if verbose:
            log_print('INFO', f'Flagging any high amplitudes viz from flares or RFIs in {msfile}')
        flagmanager(msfile, mode='save', versionname='pipeline_init')
        flagdata(vis=msfile, mode="tfcrop", spw='', action='apply', display='',
                 timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
        flagmanager(msfile, mode='save', versionname='pipeline_remove_RFI-and-BURSTS')

    if hanning:
        msfile_hanning = msfile + '.hanning'
        if not os.path.exists(msfile_hanning):
            if verbose:
                log_print('INFO', f'Perform hanningsmooth for {msfile}...')
            hanningsmooth(vis=msfile, datacolumn='data', outputvis=msfile_hanning)
        else:
            if verbose:
                log_print('INFO', f'The hanningsmoothed {msfile} already exists. Skipped...')
        msfile = msfile_hanning

    run_end_time_pre_proc = datetime.now()
    elapsed_time = run_end_time_pre_proc - run_start_time_pre_proc
    elapsed_time_pre_proc = elapsed_time.total_seconds() / 60
    log_print('INFO', f"Pre-processing completed. Elapsed time: {elapsed_time_pre_proc:.1f} minutes")

    ### --------------------------------------------------------------###
    ## shift correction block.
    ### --------------------------------------------------------------###

    total_blocks = len(tr_series_master)
    if verbose:
        log_print('INFO', f"Processing {msfile} into {total_blocks} time blocks based on the viz time ranges...")
    run_start_time_shift_corr = datetime.now()

    if ncpu == 'auto':
        ncpu = total_blocks
    else:
        try:
            ncpu = int(ncpu)
        except ValueError:
            raise ValueError("ncpu must be an integer or 'auto'.")

    mmsfiles_rot_all = []
    if not mergeFITSonly:
        if ncpu == 1:
            log_print('INFO', f"Using 1 CPU for serial processing ...")
            for tidx, (tbg, ted) in enumerate(tr_series_master):
                combined_vis_sub = process_time_block((tidx, (tbg, ted)),
                                                      msfile_in=msfile,
                                                      msname=msname,
                                                      subdir=subdir,
                                                      total_blocks=total_blocks,
                                                      tdt=tdt,
                                                      tdtstr=tdtstr,
                                                      spws=spws,
                                                      niter_init=niter_init,
                                                      reftime_master=reftime_master,
                                                      do_diskslfcal=do_diskslfcal,
                                                      disk_params=disk_params, pols=pols,
                                                      do_sbdcal=do_sbdcal, overwrite=overwrite)
                mmsfiles_rot_all.append(combined_vis_sub)
        else:
            log_print('INFO', f"Using {ncpu} CPUs for parallel processing ...")
            if is_on_server:
                log_print('INFO', f"Running on pipeline server. Using 1 CPU for splitting ...")
                mmsfiles = split_mms(msfile, tr_series_master, workdir=subdir, overwrite=overwrite, verbose=True)
                msfile2proc = mmsfiles
            else:
                msfile2proc = msfile
            worker = partial(process_time_block,
                             msfile_in=msfile2proc,
                             msname=msname,
                             subdir=subdir,
                             total_blocks=total_blocks,
                             tdt=tdt,
                             tdtstr=tdtstr,
                             spws=spws,
                             niter_init=niter_init,
                             reftime_master=reftime_master,
                             do_diskslfcal=do_diskslfcal,
                             disk_params=disk_params, pols=pols,
                             do_sbdcal=do_sbdcal, overwrite=overwrite)

            with Pool(ncpu) as pool:
                results = pool.map(worker, enumerate(tr_series_master))
            # with Pool(ncpu) as pool:
            #     # Using map_async instead of map
            #     result_object = pool.map_async(worker, enumerate(tr_series_master))
            #     results = result_object.get()
            mmsfiles_rot_all = [res for res in results if res is not None]

        if os.path.isdir(combined_vis) == True:
            shutil.rmtree(combined_vis, ignore_errors=True)

        concat(vis=[l for l in mmsfiles_rot_all if l is not None], concatvis=combined_vis)

        add_disk_before_imaging = False
        if add_disk_before_imaging:
            ## insert disk back
            combined_vis, diskcl = insertdiskmodel(combined_vis, dsize=dsize, fdens=fdens, xmlfile=diskxmlfile,
                                                   active=False)
            ## add disk model back to the viz
            uvsub(vis=combined_vis, reverse=True)
            combined_vis_disk = combined_vis + '.disk'
            split(vis=combined_vis, outputvis=combined_vis_disk, datacolumn='corrected')
            combined_vis = combined_vis_disk

    run_end_time_shift_corr = datetime.now()
    elapsed_time = run_end_time_shift_corr - run_start_time_shift_corr
    elapsed_time_shift_corr = elapsed_time.total_seconds() / 60
    log_print('INFO', f"Self-calibration processes completed. Elapsed time: {elapsed_time_shift_corr:.1f} minutes")

    ### --------------------------------------------------------------###
    ## imaging block.
    ### --------------------------------------------------------------###
    run_start_time_imaging = datetime.now()
    ## imaging the final combined ms file
    if spws_imaging is None:
        spws_imaging = spws
    if tr_series_imaging is None:
        tr_series_imaging = tr_series_master
    total_blocks_imaging = len(tr_series_imaging)
    if not mergeFITSonly:
        if ncpu == 1:
            for tidx, (tbg, ted) in enumerate(tr_series_imaging):
                timerange = trange2timerange((tbg, ted))
                log_print('INFO',
                          f"Imaging time block {tidx + 1} of {total_blocks_imaging} (Time range: {timerange}) ... ")
                fitsfile, imagefile = fd_images(combined_vis,
                                                timerange=timerange,
                                                pbcor=False,
                                                overwrite=overwrite,
                                                cleanup=False,
                                                niter=5000, spws=spws_imaging,
                                                # usemask='user', ## toggle this for single band imaging
                                                imgoutdir=subdir)
        else:
            if is_on_server:
                msfiles_in = mmsfiles_rot_all
            else:
                msfiles_in = combined_vis
            # Prepare partial function with pre-filled arguments
            process_with_params = partial(process_imaging_timerange, msfile_in=msfiles_in, spws=spws_imaging,
                                          subdir=subdir, overwrite=overwrite)

            with Pool(ncpu) as p:
                p.map(process_with_params, enumerate(tr_series_imaging))

    run_end_time_imaging = datetime.now()
    elapsed_time = run_end_time_imaging - run_start_time_imaging
    elapsed_time_imaging = elapsed_time.total_seconds() / 60
    log_print('INFO', f"Imaging processes completed. Elapsed time: {elapsed_time_imaging:.1f} minutes")

    ### --------------------------------------------------------------###
    ## post-processing block. adding disk back...
    ### --------------------------------------------------------------###
    ## merge the daily synoptic images to a single daily synoptic image
    run_start_time_post_proc = datetime.now()
    diskinfo = readdiskxml(diskxmlfile)
    syndaily_fitsfiles = []
    ## filename convention
    ## eovsa.dailysynoptic.20211124T200000_UTC.s02-05.tb.fits
    ## eovsa.synoptic_180m.20211124T180000_UTC.s02-05.tb.fits
    for spw in spws_imaging:
        spwstr = format_spw(spw)
        eofiles_rot = sorted(glob(f"{subdir}/eovsa_{date_str}T*_??min.s{spwstr}.tb.fits"))

        eofiles_rot_disk = []
        for eofile in eofiles_rot:
            datetimestr = os.path.basename(eofile).split('_')[1]
            synfitsfile = os.path.join(imgoutdir, f"eovsa.synoptic_{tdtmst_str}.{datetimestr}_UTC.s{spwstr}.tb.fits")
            eomap_disk, tb_disk, eofile_disk = image_adddisk(eofile, diskinfo)
            if eofile_disk is not None:
                shutil.move(eofile_disk, synfitsfile)
                log_print('INFO', f'Adding solar disk back to {synfitsfile}')
                eofiles_rot_disk.append(synfitsfile)

        syndaily_fitsfile = os.path.join(imgoutdir, f"eovsa.synoptic_daily.{date_str}T200000_UTC.s{spwstr}.tb.fits")
        merge_FITSfiles(eofiles_rot_disk, syndaily_fitsfile, exptime_weight=True)
        syndaily_fitsfiles.append(syndaily_fitsfiles)
        ## todo synoptic_60m FITS images currently align to 20 UT. Need to rotate them to better reflect actual observation times. Also, the FITS files need to be compressed.

    ## remove the intermediate ms files
    for i in mmsfiles_rot_all:
        if os.path.isdir(i):
            shutil.rmtree(i, ignore_errors=True)

    if outputvis:
        if os.path.exists(outputvis):
            shutil.rmtree(outputvis)
        shutil.move(combined_vis, outputvis)
        combined_vis = outputvis

        newdiskxmlfile = '{}.SOLDISK.xml'.format(outputvis)
        if os.path.exists(newdiskxmlfile):
            shutil.rmtree(newdiskxmlfile)
        shutil.move(diskxmlfile, newdiskxmlfile)

    if clearcache:
        shutil.rmtree(subdir, ignore_errors=True)

    run_end_time_post_proc = datetime.now()
    elapsed_time = run_end_time_post_proc - run_start_time_post_proc
    elapsed_time_post_proc = elapsed_time.total_seconds() / 60
    log_print('INFO', f"Merging completed. Elapsed time: {elapsed_time_post_proc:.1f} minutes")

    # Create a dictionary to store the time taken for each step
    time_dict = {
        "Pre-processing": elapsed_time_pre_proc,
        "Self-calibration": elapsed_time_shift_corr,
        "Imaging": elapsed_time_imaging,
        "Post-processing": elapsed_time_post_proc
    }

    # Calculate the total processing time
    total_time = sum(time_dict.values())

    # Log the total processing time
    log_print('INFO', f"Total processing time: {total_time:.1f} minutes.")

    # Log the time taken for each step
    for step, tim in time_dict.items():
        log_print('INFO', f"    --> {step}: {tim:.1f} minutes.")
    return combined_vis


if __name__ == '__main__':
    '''
    this code is trying to address the issue of smearing effect in all-day synthesis images.
    description of the issue: https://www.ovsa.njit.edu/wiki/index.php/All-Day_Synthesis_Issues
    This code is built on python 3.8.
    '''
    description = 'this code is trying to address the issue of smearing effect in all-day synthesis images. Description of the issue: https://www.ovsa.njit.edu/wiki/index.php/All-Day_Synthesis_Issues. This code is built on python 3.8.'
    parser = argparse.ArgumentParser(
        description="Executes the EOVSA data processing pipeline for solar observation data.")
    parser.add_argument('vis', type=str, help='Path to the input measurement set (MS) data.')
    parser.add_argument('--outputvis', type=str, default='', help='Output path for the processed MS.')
    parser.add_argument('--workdir', type=str, help='Working directory for intermediate data.')
    parser.add_argument('--slfcaltbdir', type=str, help='Directory for storing calibration tables.')
    parser.add_argument('--imgoutdir', type=str, help='Output directory for image files.')
    parser.add_argument('--figoutdir', type=str, help='Output directory for figures.')
    parser.add_argument('--clearcache', action='store_true', help='Clears the cache after processing.')
    parser.add_argument('--pols', type=str, default='XX', help='Polarization types to process.')
    parser.add_argument('--mergeFITSonly', action='store_true', help='Skips processing and only merges FITS files.')
    # parser.add_argument('--verbose', action='store_true', help='Enables verbose output during processing.')
    # parser.add_argument('--do_diskslfcal', action='store_true', help='Performs disk self-calibration.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrites existing files.')
    parser.add_argument('--niter_init', type=int, default=200, help='Initial number of iterations for imaging.')
    parser.add_argument('--ncpu', type=str, default='auto',
                        help="Specifies the number of CPUs for parallel processing.")
    parser.add_argument('--tr_series_imaging', nargs='*', help='Time ranges for imaging, expects a list of tuples.')
    parser.add_argument('--spws_imaging', nargs='*', help='Spectral windows selected for imaging.')
    parser.add_argument('--hanning', action='store_true', help='Applies Hanning smoothing to the data.')
    parser.add_argument('--do_sbdcal', action='store_true', help='Perform single-band delay calibration.')

    args = parser.parse_args()

    pipeline_run(
        vis=args.vis,
        outputvis=args.outputvis,
        workdir=args.workdir,
        slfcaltbdir=args.slfcaltbdir,
        imgoutdir=args.imgoutdir,
        figoutdir=args.figoutdir,
        clearcache=args.clearcache,
        pols=args.pols,
        mergeFITSonly=args.mergeFITSonly,
        # verbose=args.verbose,
        # do_diskslfcal=args.do_diskslfcal,
        overwrite=args.overwrite,
        niter_init=args.niter_init,
        ncpu=args.ncpu,
        tr_series_imaging=args.tr_series_imaging,
        spws_imaging=args.spws_imaging,
        hanning=args.hanning,
        do_sbdcal=args.do_sbdcal
    )
