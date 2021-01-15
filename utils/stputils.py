import numpy as np
import glob
import os

# import json
# import pickle
# from functools import wraps
# import astropy.units as u

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"


def insertchar(source_str, insert_str, pos):
    return source_str[:pos] + insert_str + source_str[pos:]


def get_curve_grad(x, y):
    '''
    get the grad of at data point
    :param x:
    :param y:
    :return: grad,posang
    '''
    deltay = np.roll(y, -1) - np.roll(y, 1)
    deltay[0] = y[1] - y[0]
    deltay[-1] = y[-1] - y[-2]
    deltax = np.roll(x, -1) - np.roll(x, 1)
    deltax[0] = x[1] - x[0]
    deltax[-1] = x[-1] - x[-2]
    grad = deltay / deltax
    posang = np.arctan2(deltay, deltax)
    return {'grad': grad, 'posang': posang}


def findDist(x, y):
    dx = np.diff(x)
    dy = np.diff(y)
    dist = np.hypot(dx, dy)
    return np.insert(dist, 0, 0.0)


def paramspline(x, y, length, s=0):
    from scipy.interpolate import splev, splprep
    tck, u = splprep([x, y], s=s)
    unew = np.linspace(0, u[-1], length)
    out = splev(unew, tck)
    xs, ys = out[0], out[1]
    grads = get_curve_grad(xs, ys)
    return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def polyfit(x, y, length, deg, keepxorder=False):
    T = False
    if len(x) == 2:
        if np.abs(np.diff(y) / np.diff(x)) < 1:
            xproc = x
            yproc = y
        else:
            xproc = y
            yproc = x
            T = True
    else:
        xproc = x
        yproc = y
    if keepxorder:
        xs = np.linspace(xproc[0], xproc[-1], length)
    else:
        xs = np.linspace(np.nanmin(xproc), np.nanmax(xproc), length)
    z = np.polyfit(x=xproc, y=yproc, deg=deg)
    p = np.poly1d(z)
    ys = p(xs)
    grads = get_curve_grad(xs, ys)
    if T:
        return {'xs': ys, 'ys': xs, 'grads': 1/grads['grad'], 'posangs': np.pi/2-grads['posang']}
    else:
        return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def improfile(z, xi, yi, interp='cubic'):
    '''
    Pixel-value cross-section along line segment in an image
    :param z: an image array
    :param xi and yi: equal-length vectors specifying the pixel coordinates of the endpoints of the line segment
    :param interp: interpolation type to sampling, 'nearest' or 'cubic'
    :return: the intensity values of pixels along the line
    '''
    import scipy.ndimage
    imgshape = z.shape
    if len(xi) != len(yi):
        raise ValueError('xi and yi must be equal-length!')
    if len(xi) < 2:
        raise ValueError('xi or yi must contain at least two elements!')
    for idx, ll in enumerate(xi):
        if not 0 < ll < imgshape[1]:
            return np.nan
        if not 0 < yi[idx] < imgshape[0]:
            return np.nan
    if len(xi) == 2:
        length = np.floor(np.hypot(np.diff(xi), np.diff(yi))[0]).astype(np.int)
        x, y = np.linspace(xi[0], xi[1], length), np.linspace(yi[0], yi[1], length)
    else:
        x, y = xi, yi
    if interp == 'cubic':
        zi = scipy.ndimage.map_coordinates(z, np.vstack((x, y)))
    else:
        zi = z[np.floor(y).astype(np.int), np.floor(x).astype(np.int)]

    return zi


def map2wcsgrids(snpmap, cell=True, antialiased=True):
    '''

    :param snpmap:
    :param cell: if True, return the coordinates of the pixel centers. if False, return the coordinates of the pixel boundaries
    :return:
    '''
    # embed()
    import astropy.units as u
    ny, nx = snpmap.data.shape
    x0, x1 = snpmap.xrange.to(u.arcsec).value
    y0, y1 = snpmap.yrange.to(u.arcsec).value
    dx = snpmap.scale.axis1.to(u.arcsec / u.pix).value
    dy = snpmap.scale.axis2.to(u.arcsec / u.pix).value

    if cell:
        mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)
        mapx = np.tile(mapx, ny).reshape(ny, nx)
        mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
    else:
        nx += 1
        ny += 1
        mapx, mapy = np.linspace(x0 - dx / 2.0, x1 + dx / 2.0, nx), np.linspace(y0 - dy / 2.0, y1 + dy / 2.0, ny)
        mapx = np.tile(mapx, ny).reshape(ny, nx)
        mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
    return mapx, mapy


def readsdofile(datadir=None, wavelength=None, trange=None, isexists=False, timtol=1):
    '''
    read sdo file from local database
    :param datadir:
    :param wavelength:
    :param trange: the timestamp or timerange in Julian days. if is timerange, return a list of files in the timerange
    :param isexists: check if file exist. if files exist, return file name
    :param timtol: time difference tolerance in days for considering data as the same timestamp
    :return:
    '''
    from astropy.time import Time
    import sunpy.map
    from datetime import date
    from datetime import timedelta as td

    trange = Time(trange)
    wavelength = str(wavelength)
    wavelength = wavelength.lower()
    if timtol < 12. / 3600 / 24:
        timtol = 12. / 3600 / 24

    if len(trange.iso) == 2:
        if trange[1] < trange[0]:
            raise ValueError('start time must be occur earlier than end time!')
        else:
            sdofitspath = []
            jdtimestr = [Time(ll, format='jd').iso for ll in trange]
            ymd = [ll.split(' ')[0].split('-') for ll in jdtimestr]
            d1 = date(int(ymd[0][0]), int(ymd[0][1]), int(ymd[0][2]))
            d2 = date(int(ymd[1][0]), int(ymd[1][1]), int(ymd[1][2]))
            delta = d2 - d1
            for i in range(delta.days + 1):
                ymd = d1 + td(days=i)
                sdofitspathtmp = glob.glob(
                    datadir + '/{:04d}/{:02d}/{:02d}/aia.lev1_*Z.{}.image_lev1.fits'.format(ymd.year, ymd.month,
                                                                                            ymd.day, wavelength))
                if len(sdofitspathtmp) > 0:
                    sdofitspath = sdofitspath + sdofitspathtmp
            if len(sdofitspath) == 0:
                if isexists:
                    return sdofitspath
                else:
                    raise ValueError(
                        'No SDO file found under {} at the time range of {} to {}. Download the data with EvtBrowser first.'.format(
                            datadir,
                            jdtimestr[0],
                            jdtimestr[1]))
            sdofits = [os.path.basename(ll) for ll in sdofitspath]
            sdotimeline = Time(
                [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
                 sdofits],
                format='iso', scale='utc')
            sdofitspathnew = [x for (y, x) in sorted(zip(sdotimeline.jd, sdofitspath))]
            sdofitsnew = [os.path.basename(ll) for ll in sdofitspathnew]
            sdotimelinenew = Time(
                [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
                 sdofitsnew], format='iso',
                scale='utc')
            sdofile = list(
                np.array(sdofitspathnew)[
                    np.where(np.logical_and(trange[0].jd <= sdotimelinenew.jd, sdotimelinenew.jd <= trange[1].jd))[0]])
            return sdofile
    else:
        jdtimstr = trange.iso
        ymd = jdtimstr.split(' ')[0].split('-')
        sdofitspath = glob.glob(
            datadir + '/{}/{}/{}/aia.lev1_*Z.{}.image_lev1.fits'.format(ymd[0], ymd[1], ymd[2], wavelength))
        if len(sdofitspath) == 0:
            return []  # raise ValueError('No SDO file found under {}.'.format(datadir))
        sdofits = [os.path.basename(ll) for ll in sdofitspath]
        sdotimeline = Time(
            [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
             sdofits],
            format='iso', scale='utc')
        if timtol < np.min(np.abs(sdotimeline.jd - trange.jd)):
            return []  # raise ValueError('No SDO file found at the select timestamp. Download the data with EvtBrowser first.')
        idxaia = np.argmin(np.abs(sdotimeline.jd - trange.jd))
        sdofile = sdofitspath[idxaia]
        if isexists:
            return sdofile
        else:
            try:
                sdomap = sunpy.map.Map(sdofile)
                return sdomap
            except:
                raise ValueError('File not found or invalid input')


def get_map_corner_coord(sunpymap):
    bottom_left_coord = sunpymap.bottom_left_coord
    top_right_coord = sunpymap.top_right_coord
    x0, y0 = bottom_left_coord.Tx, bottom_left_coord.Ty
    x1, y1 = top_right_coord.Tx, top_right_coord.Ty
    unit = x0.unit
    return np.array([x0.value, x1.value, y0.value, y1.value]) * unit
