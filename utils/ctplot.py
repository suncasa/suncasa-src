import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import astropy.units as u
from astropy.io import fits
import numpy as np
import sunpy.map
from suncasa.utils import DButil
from suncasa.utils.puffin import PuffinMap


def maxfit(smap, width=[5, 5], mapxy=None):
    '''
    find maximum position of an image using 2D Gaussian fit
    :param smap: sunpy map
    :param width:
    :param mapxy: [mapx,mapy]
    :return:
    '''

    data = smap.data
    if mapxy:
        mapx, mapy = mapxy
    else:
        XX, YY = np.meshgrid(np.arange(smap.data.shape[0]), np.arange(smap.data.shape[1]))
        mapx, mapy = smap.pixel_to_data(XX * u.pix, YY * u.pix)
        mapx, mapy = mapx.value, mapy.value

    idxmax = np.where(data == np.nanmax(data))

    try:
        z = np.polyfit(x=mapx[idxmax[0][0], idxmax[1][0] - width[0]:idxmax[1][0] + width[0] + 1],
                       y=data[idxmax[0][0], idxmax[1][0] - width[0]:idxmax[1][0] + width[0] + 1], deg=2)
        xmax = -z[1] / z[0] / 2.0
        z = np.polyfit(x=mapy[idxmax[0][0] - width[1]:idxmax[0][0] + width[1] + 1, idxmax[1][0]],
                       y=data[idxmax[0][0] - width[1]:idxmax[0][0] + width[1] + 1, idxmax[1][0]], deg=2)
        ymax = -z[1] / z[0] / 2.0
        return [xmax, ymax]
    except:
        print '2D polynomial fitting failed.'
        return None


def contour1chn(vlafile, aiafile, chn=0, pol=0, x_range=[], y_range=[], levels=[0.2, 0.4, 0.6, 0.8]):
    aiamap = sunpy.map.Map(aiafile)
    if x_range and y_range:
        aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec),
                               u.Quantity(y_range * u.arcsec))
    hdulist = fits.open(vlafile)
    hdu = hdulist[0]
    vladata = hdu.data[pol, chn, :, :]
    vlamap = sunpy.map.Map((vladata, hdu.header))
    cmap = sunpy.map.CompositeMap(aiamap)
    cmap.add_map(vlamap, levels=np.array(levels) * np.nanmax(vlamap.data))
    cmap.peek()


def contourchans(vlafile, aiafile, pol=0, chans=[], x_range=[], y_range=[], levels=[0.9], **kwargs):
    aiamap = sunpy.map.Map(aiafile)
    if x_range and y_range:
        aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec),
                               u.Quantity(y_range * u.arcsec))
    hdulist = fits.open(vlafile)
    hdu = hdulist[0]
    vladata = hdu.data[pol, 0, :, :]
    vlamap = sunpy.map.Map((vladata, hdu.header))
    vlamap_pfmap = PuffinMap(smap=vlamap)
    mapx, mapy = vlamap_pfmap.meshgrid()

    fig = plt.figure()
    ax = plt.subplot(projection=aiamap)
    im = aiamap.plot()
    # # Prevent the image from being re-scaled while overplotting.
    ax.set_autoscale_on(False)
    nfreq = hdu.data[pol, :, :, :].shape[1]
    if len(chans) == 0:
        chans = range(0, nfreq)
    nchan = len(chans)
    for idx, chan in enumerate(chans):
        vladata = hdu.data[pol, chan, :, :]
        vlamap = sunpy.map.Map((vladata, hdu.header))
        SRC_vlamap_contour = DButil.get_contour_data(mapx.value, mapy.value, vlamap.data, levels=levels)
        if SRC_vlamap_contour.data['xs']:
            xs = SRC_vlamap_contour.data['xs'][0] * u.arcsec
            ys = SRC_vlamap_contour.data['ys'][0] * u.arcsec
            plt.plot(xs.to(u.deg), ys.to(u.deg), transform=ax.get_transform('world'),
                     color=cm.jet(int(float(chan) / nfreq * 255)), zorder=nchan - idx, **kwargs)
    freqs = (hdu.header['CRVAL3'] + np.arange(nfreq) * hdu.header['CDELT3']) / 1e9
    fig.subplots_adjust(right=0.8)
    ax1 = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cmap = cm.jet
    norm = colors.Normalize(vmin=freqs[0], vmax=freqs[-1])
    cb = colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
    cb.set_label('Frequency [GHz]')


def pltcentriod(vlafile, aiafile, pol=0, chans=[], x_range=[], y_range=[], **kwargs):
    aiamap = sunpy.map.Map(aiafile)
    if x_range and y_range:
        aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec),
                               u.Quantity(y_range * u.arcsec))
    hdulist = fits.open(vlafile)
    hdu = hdulist[0]
    vladata = hdu.data[pol, 0, :, :]
    vlamap = sunpy.map.Map((vladata, hdu.header))
    XX, YY = np.meshgrid(np.arange(vlamap.data.shape[0]), np.arange(vlamap.data.shape[1]))
    mapx, mapy = vlamap.pixel_to_data(XX * u.pix, YY * u.pix)
    mapx, mapy = mapx.value, mapy.value
    fig = plt.figure()
    ax = plt.subplot(projection=aiamap)
    im = aiamap.plot()
    # # Prevent the image from being re-scaled while overplotting.
    ax.set_autoscale_on(False)
    nfreq = hdu.data[pol, :, :, :].shape[1]
    if len(chans) == 0:
        chans = range(0, nfreq)
    nchan = len(chans)
    for idx, chan in enumerate(chans):
        vladata = hdu.data[pol, chan, :, :]
        if np.nanmax(vladata):
            vlamap.data = vladata
            maxxy = maxfit(vlamap, mapxy=[mapx, mapy])
            if maxxy:
                x, y = maxxy * u.arcsec
                plt.plot(x.to(u.deg), y.to(u.deg), '+', transform=ax.get_transform('world'),
                         color=cm.jet(int(float(chan) / nfreq * 255)), zorder=nchan - idx, **kwargs)
    freqs = (hdu.header['CRVAL3'] + np.arange(nfreq) * hdu.header['CDELT3']) / 1e9
    fig.subplots_adjust(right=0.8)
    ax1 = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cmap = cm.jet
    norm = colors.Normalize(vmin=freqs[0], vmax=freqs[-1])
    cb = colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
    cb.set_label('Frequency [GHz]')