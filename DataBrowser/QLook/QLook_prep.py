import numpy as np
from bokeh.models import (
    ColumnDataSource)
import pandas as pd
import jdutil
import os
from astropy.io import fits
import sunpy.map
import scipy.optimize as opt
import pylab as plt
import astropy.units as u
import time
from functools import partial
import pickle
import matplotlib as mplt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import sunpy.cm as spcm

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

# define the colormaps
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]


def get_contour_data(X, Y, Z):
    cs = plt.contour(X, Y, Z, levels=(np.arange(5, 10, 2) / 10.0 * Z.max()).tolist())
    xs = []
    ys = []
    xt = []
    yt = []
    col = []
    text = []
    isolevelid = 0
    for isolevel in cs.collections:
        isocol = isolevel.get_color()[0]
        thecol = 3 * [None]
        theiso = '{:.0f}%'.format(cs.get_array()[isolevelid] / Z.max() * 100)
        isolevelid += 1
        for i in range(3):
            thecol[i] = int(255 * isocol[i])
        thecol = '#%02x%02x%02x' % (thecol[0], thecol[1], thecol[2])

        for path in isolevel.get_paths():
            v = path.vertices
            x = v[:, 0]
            y = v[:, 1]
            xs.append(x.tolist())
            ys.append(y.tolist())
            xt.append(x[len(x) / 2])
            yt.append(y[len(y) / 2])
            text.append(theiso)
            col.append(thecol)

    source = ColumnDataSource(data={'xs': xs, 'ys': ys, 'line_color': col, 'xt': xt, 'yt': yt, 'text': text})
    return source


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                       + c * ((y - yo) ** 2)))
    return g.ravel()


def maxfit(image, plot=False):
    # NAME:
    #    maxfit
    # PURPOSE:
    #    find maximum position of an image using 2D Gaussian fit
    # INPUTS:
    #    image: sunpy

    data = image.data
    nx, ny = data.shape
    dx, dy = image.scale.x.value, image.scale.y.value
    xc, yc = image.center.x.value, image.center.y.value
    mapx, mapy = (np.linspace(0, nx - 1, nx) - image.reference_pixel.x.value + 1 + 0.5) * dx + xc, (
    np.linspace(0, ny - 1, ny) - image.reference_pixel.y.value + 1 + 0.5) * dy + yc
    mapx, mapy = np.meshgrid(mapx, mapy)
    idxmax = np.where(data == np.amax(data))
    idxhm = np.where(data >= np.amax(data) / 2)
    hmfw = np.sqrt(len(idxhm[0]) / np.pi) * dx
    xmax, ymax = mapx[idxmax[0][0], idxmax[1][0]], mapy[idxmax[0][0], idxmax[1][0]]
    # add some noise to the data and try to fit the data generated beforehand
    theta = np.arctan(
        abs(idxhm[1] - image.reference_pixel.y.value).sum() / abs(idxhm[0] - image.reference_pixel.x.value).sum())
    initial_guess = (np.amax(data), xmax, ymax, hmfw, hmfw, theta, data.std())

    try:
        popt, pcov = opt.curve_fit(twoD_Gaussian, (mapx, mapy), data.ravel(), p0=initial_guess)
    except:
        popt = np.empty((7))
        popt[:] = np.nan

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = plt.subplot(121, projection=image)
        im = image.plot()
        # plt.colorbar()
        ax1.set_autoscale_on(False)
        data_fitted = twoD_Gaussian((mapx, mapy), *popt).reshape(nx, ny)
        # image_fitted = copy.deepcopy(image)
        # image_fitted.data = data_fitted
        # image_fitted.draw_contours([10,20,30,40,50,60,70,80,90] * u.percent)
        plt.plot((popt[1] * u.arcsec).to(u.deg), (popt[2] * u.arcsec).to(u.deg), 'o',
            transform=ax1.get_transform('world'))
        plt.contour((mapx * u.arcsec).to(u.deg), (mapy * u.arcsec).to(u.deg), data_fitted,
            levels=(np.arange(5, 10, 2) / 10.0 * np.amax(data_fitted)).tolist(),
            transform=ax1.get_transform('world'))
    # from sunpy.net.helioviewer import HelioviewerClient
    # import matplotlib.colors as colors
    # hv = HelioviewerClient()
    # filepath = hv.download_jp2(image.date, observatory='SDO', instrument='AIA', detector='AIA', measurement='171',directory='/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/output/U01/database/J2000/', overwrite=True)
    # image_aia = sunpy.map.Map(filepath)
    # image_aia.plot_settings['norm'] = colors.LogNorm(30, image_aia.max())
    # # fig = plt.figure()
    # lengthx = nx*dx/2.0*u.arcsec
    # lengthy = ny*dy/2.0*u.arcsec
    # x0 = xc*u.arcsec
    # y0 = yc*u.arcsec
    # image_aia_submap = image_aia.submap(u.Quantity([x0 - lengthx, x0 + lengthx]),
    # 	                     u.Quantity([y0 - lengthy, y0 + lengthy]))
    # ax2 = plt.subplot(122,projection=image_aia_submap)
    # im = image_aia_submap.plot()
    # # plt.colorbar()
    # ax2.set_autoscale_on(False)
    # plt.plot((popt[1]*u.arcsec).to(u.deg), (popt[2]*u.arcsec).to(u.deg), 'o',
    # 	transform=ax2.get_transform('world'))
    # plt.contour((mapx*u.arcsec).to(u.deg),(mapy*u.arcsec).to(u.deg),data_fitted,levels=(np.arange(5,10,2)/10.0*data_fitted.max()).tolist(),transform=ax2.get_transform('world'))

    return popt


prep_posmaxfit = True
prep_thumbnail = False

if prep_posmaxfit:
    # Gauss_params=[]

    # start_timestamp = time.time()
    # for f in fits_local:
    # 	if os.path.exists(vla_local_fitspath+f):
    # 		hdulist = fits.open(vla_local_fitspath+f)
    # 		hdu=hdulist[0]
    # 		vlamap = sunpy.map.Map((hdu.data[0,0,:,:], hdu.header))
    # 		gauss_params = maxfit(vlamap,plot=False)
    # 		fexist = True
    # 	else:
    # 		gauss_params = np.empty((7))
    # 		gauss_params[:] = np.nan
    # 		fexist = False
    # 	Gauss_params.append([f,gauss_params,fexist])

    # print("--- %s seconds ---" % (time.time() - start_timestamp))
    structure_id = 'U01'
    specdata = np.load('U01_2014-Nov-01T164613.125~164629.475.spec.npz')
    spec = specdata['spec']
    npol = specdata['npol']
    nbl = specdata['nbl']
    ntim = specdata['ntim']
    nfreq = specdata['nfreq']
    tim = specdata['tim']
    freq = specdata['freq'] / 1e9
    bl = specdata['bl'].item()
    b = 0
    pol = 'LL'
    sz_spec = spec.shape
    spec_plt = np.zeros((4, sz_spec[2], sz_spec[3]))
    spec_plt[0, :, :] = spec[1, b, :, :]
    dtim = tim - tim[0]
    tim_map = ((np.tile(tim, nfreq).reshape(nfreq, ntim) / 3600. / 24. + 2400000.5)) * 86400.
    freq_map = np.tile(freq, ntim).reshape(ntim, nfreq).swapaxes(0, 1)
    xx = tim_map.flatten()
    yy = freq_map.flatten()
    vla_local_fitspath = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/output/U01/database/vla/local/'
    vla_local_thumbnailpath = './QLook/static/thumbnail/U01/local/'
    vla_global_fitspath = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/output/U01/database/vla/global/'
    fits_local = []
    fits_global = []

    StructureIdList = dict(
        time=[[tim[0], tim[-1], tim[-1], tim[0]]],
        freq=[[freq[0], freq[0], freq[-1], freq[-1]]],
        str_id=[['U01']]
    )

    with open('structure_id_list', 'wb') as f:
        pickle.dump(StructureIdList, f)

    # pb = progressbar(idx_tim.size, "*")
    # for ii in range(8405,len(idx_selec[0])):
    # for ii in range(0,idx_tim.size):

    timestrs=[]
    for ii in range(len(xx)):
        t_int = 0.05
        f_int = 2.0  # MHz
        t0 = xx[ii]  # -0.5*t_int
        timestr0 = jdutil.jd_to_datetime(t0 / 3600. / 24.)
        timestr = timestr0.strftime('%Y-%m-%dT%H%M%S') + '.{:03d}'.format(int(round(timestr0.microsecond / 1e3)))
        timestrs.append(timestr0.strftime('%H:%M:%S') + '.{:03d}'.format(int(round(timestr0.microsecond / 1e3))))
        f0 = yy[ii] * 1e3
        freqstr = '{:d}MHz'.format(int(round(f0)))
        fits_local.append(structure_id + '_' + timestr + '_' + freqstr + '.slfcal.image.cutout.fits')
        fits_global.append(structure_id + '_' + timestr + '_' + freqstr + '.slfcal.image.fits')

    Gauss_params = []
    start_timestamp = time.time()


    def func(vla_local_fitspath=None, iterable=None):
        if os.path.exists(vla_local_fitspath + iterable):
            hdulist = fits.open(vla_local_fitspath + iterable)
            hdu = hdulist[0]
            vlamap = sunpy.map.Map((hdu.data[0, 0, :, :], hdu.header))
            gauss_params = maxfit(vlamap, plot=False)
            fexist = True
        else:
            gauss_params = np.empty((7))
            gauss_params[:] = np.nan
            fexist = False
        return [gauss_params, fexist]


    func_partial = partial(func, vla_local_fitspath)
    iterable = fits_local
    import multiprocessing as mp

    ncpu = mp.cpu_count() / 2
    pool = mp.Pool(processes=ncpu)
    result = []
    # dummpy = pool.map_async(func_partial,iterable)
    for iterable in fits_local:
        result.append(pool.apply_async(func_partial, (iterable,)))
    pool.close()
    pool.join()

    print("--- %s seconds ---" % (time.time() - start_timestamp))
    for res in result:
        Gauss_params.append((res.get()))

    dspecDF = pd.DataFrame({'time': xx - xx[0],
                            'freq': yy,
                            'timestr': timestrs,
                            'dspec': spec_plt[0, :, :].flatten(),
                            'fits_local': fits_local,
                            'fits_global': fits_global,
                            'fits_exist': [ll[1] for ll in Gauss_params],
                            'thumbnail': [vla_local_thumbnailpath + ll[0:-4] + 'small.jpg' for ll in
                                          fits_local],
                            # 'thumbnail':['http://bokeh.pydata.org/static/snake3D.png' for ll in fits_local],
                            'x_pos': [ll[0][1] for ll in Gauss_params],
                            'y_pos': [ll[0][2] for ll in Gauss_params],
                            'x_width': [ll[0][3] for ll in Gauss_params],
                            'y_width': [ll[0][4] for ll in Gauss_params],
                            'amp_gaus': [ll[0][0] for ll in Gauss_params],
                            'theta': [ll[0][5] for ll in Gauss_params],
                            'amp_offset': [ll[0][6] for ll in Gauss_params]})

    with open('dspecDF-save', 'wb') as f:
        pickle.dump(dspecDF, f)
else:
    with open('dspecDF-save', 'rb') as f:
        dspecDF = pickle.load(f)

# prep the thumbnails of the vla images
if prep_thumbnail:
    start_timestamp = time.time()
    plt.ioff()
    for ll in range(len(dspecDF)):
        # for ll in range(408,409):
        # for ll in range(7812,7813):
        if not os.path.exists(vla_local_thumbnailpath + dspecDF.loc[ll]['fits_local'][0:-4] + 'small.jpg'):
            if not np.isnan(dspecDF.loc[ll]['x_pos']) and dspecDF.loc[ll]['amp_gaus'] > 0:
                hdulist = fits.open(vla_local_fitspath + dspecDF.loc[ll]['fits_local'])
                hdu = hdulist[0]
                vlamap = sunpy.map.Map((hdu.data[0, 0, :, :], hdu.header))
                nx, ny = vlamap.data.shape
                dx, dy = vlamap.scale.x.value, vlamap.scale.y.value
                xc, yc = vlamap.center.x.value, vlamap.center.y.value
                mapx, mapy = (np.linspace(0, nx - 1, nx) - vlamap.reference_pixel.x.value + 1 + 0.5) * dx + xc, (
                np.linspace(0, ny - 1, ny) - vlamap.reference_pixel.y.value + 1 + 0.5) * dy + yc
                mapx, mapy = np.meshgrid(mapx, mapy)
                fig = plt.figure(figsize=(4, 4))
                ax1 = plt.subplot(111, projection=vlamap)
                im = vlamap.plot()
                # plt.colorbar()
                ax1.set_autoscale_on(False)
                popt = [dspecDF.loc[ll, :]['amp_gaus'],
                        dspecDF.loc[ll, :]['x_pos'],
                        dspecDF.loc[ll, :]['y_pos'],
                        dspecDF.loc[ll, :]['x_width'],
                        dspecDF.loc[ll, :]['y_width'],
                        dspecDF.loc[ll, :]['theta'],
                        dspecDF.loc[ll, :]['amp_offset']]
                data_fitted = twoD_Gaussian((mapx, mapy), *popt).reshape(nx, ny)
                plt.plot((popt[1] * u.arcsec).to(u.deg), (popt[2] * u.arcsec).to(u.deg), 'o',
                    transform=ax1.get_transform('world'), markersize=3)
                print np.amax(data_fitted), ll
                if popt[0] + popt[6] > 0:
                    plt.contour((mapx * u.arcsec).to(u.deg), (mapy * u.arcsec).to(u.deg), data_fitted,
                        levels=(np.arange(5, 10, 2) / 10.0 * (popt[0] + popt[6])).tolist(),
                        transform=ax1.get_transform('world'))
                fig.tight_layout()
                fig.savefig(vla_local_thumbnailpath + dspecDF.loc[ll]['fits_local'][0:-4] + 'small.jpg')
                plt.close(fig)
    print("--- %s seconds ---" % (time.time() - start_timestamp))
