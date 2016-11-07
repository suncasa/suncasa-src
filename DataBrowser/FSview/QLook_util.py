import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pylab as plt
import scipy.optimize as opt
from bokeh.models import (ColumnDataSource)
import matplotlib.cm as cm

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"


def get_contour_data(X, Y, Z):
    cs = plt.contour(X, Y, Z, levels=(np.arange(5, 10, 2) / 10.0 * Z.max()).tolist(), cmap=cm.Greys_r)
    xs = []
    ys = []
    xt = []
    yt = []
    col = []
    text = []
    isolevelid = 0
    for isolevel in cs.collections:
        # isocol = isolevel.get_color()[0]
        # thecol = 3 * [None]
        theiso = '{:.0f}%'.format(cs.get_array()[isolevelid] / Z.max() * 100)
        isolevelid += 1
        # for i in xrange(3):
        # thecol[i] = int(255 * isocol[i])
        thecol = '#%02x%02x%02x' % (220, 220, 220)
        # thecol = '#03FFF9'

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
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
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
