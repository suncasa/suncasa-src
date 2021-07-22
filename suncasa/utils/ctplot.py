import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as colorbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
import numpy as np
import sunpy.map
from suncasa.utils import DButil


def add_scalebar(length, x, y, ax, align='right', label='', label_position='below', yoff=-0.005, xoff=None, fontsize=None, **kwargs):
    '''
    :param length : length of scalebar in data coordinate
    :param x,y : in the coordinate system of the Axes; (0,0) is bottom left of the axes, and (1,1) is top right of the axes.
    :param ax:
    :param text:
    :param yoff : Yoffset of the label in the coordinate system of the Axes
    :return:
    '''
    x, y = ax.transData.inverted().transform(ax.transAxes.transform([x, y]))
    if align == 'right':
        x0, y0 = x - length, y
        x1, y1 = x, y
    elif align == 'left':
        x0, y0 = x, y
        x1, y1 = x + length, y
    elif align == 'middle':
        x0, y0 = x - length / 2.0, y
        x1, y1 = x + length / 2.0, y
    kwargs_text = {}
    if 'c' or 'color' in kwargs.keys():
        try:
            kwargs_text['color'] = kwargs['c']
        except:
            kwargs_text['color'] = kwargs['color']
    else:
        kwargs['c'] = 'white'
        kwargs_text['color'] = 'white'
    if fontsize is not None:
        kwargs_text['fontsize'] = fontsize
    ax.plot([x0, x1], [y0, y1], **kwargs)
    xtext, ytext = ax.transAxes.inverted().transform(ax.transData.transform([(x0 + x1) / 2.0, y0]))
    if label_position == 'above':
        verticalalign = 'bottom'
        ytext = ytext + yoff
    if label_position == 'below':
        verticalalign = 'top'
        ytext = ytext + yoff
    if xoff:
        xtext = xtext + xoff
    ax.text(xtext, ytext, label, horizontalalignment='center', verticalalignment=verticalalign, transform=ax.transAxes, **kwargs_text)


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
        XX, YY = np.meshgrid(np.arange(smap.data.shape[1]), np.arange(smap.data.shape[0]))
        # mapx, mapy = smap.pixel_to_data(XX * u.pix, YY * u.pix)
        mesh = smap.pixel_to_world(XX * u.pix, YY * u.pix)
        mapx, mapy = mesh.Tx.value, mesh.Tx.value

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
        print('2D polynomial fitting failed.')
        return None


def contour1chn(vlafile, aiafile, chn=0, pol=0, x_range=[], y_range=[], levels=[0.2, 0.4, 0.6, 0.8], cmap='jet', *args, **kargs):
    aiamap = sunpy.map.Map(aiafile)
    if x_range and y_range:
        aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec), u.Quantity(y_range * u.arcsec))
    ax = plt.subplot()
    aiamap.plot(axes=ax)
    if type(vlafile) is not list:
        vlafile = [vlafile]
    for idx, vf in enumerate(vlafile):
        hdulist = fits.open(vf)
        hdu = hdulist[0]
        dims = hdu.data.shape
        if len(dims) == 2:
            vladata = hdu.data
        elif len(dims) == 4:
            vladata = hdu.data[pol, chn, :, :]
        else:
            raise ValueError('check you import vla fits file')
        vlamap = sunpy.map.Map((vladata, hdu.header))
        sy, sx = vlamap.data.shape
        XX, YY = np.meshgrid(np.arange(sy), np.arange(sx))
        vlamapx = vlamap.pixel_to_data(XX * u.pix, YY * u.pix).Tx
        vlamapy = vlamap.pixel_to_data(XX * u.pix, YY * u.pix).Ty
        ax.contour(vlamapx.value, vlamapy.value, vlamap.data, levels=np.array(levels) * np.nanmax(vlamap.data),
                   colors=[cm.get_cmap(cmap)(float(idx) / (len(vlafile) - 1))] * len(levels), *args, **kargs)


def plot_compsite_map(vlafile, aiafile, outfile='', label='', pol=0, chans=[], chan_mask=None, x_range=[], y_range=[], levels=[0.9],
                      plotstyle='centroid', figsize=(10, 8), figdpi=100, width=[5, 5], zorder=1, maponly=False, vlaonly=False, dspecdata={},
                      thrshd=None, cmap='jet', plt_clbar=True, aiaplt_args={'axes': None, 'cmap': None, 'vmax': None, 'vmin': None}, **kwargs):
    cmap = cm.get_cmap(cmap)
    if aiaplt_args['cmap'] is None:
        aiaplt_args.pop('cmap', None)
    if outfile:
        plt.ioff()
    if not vlaonly:
        if type(aiafile) != sunpy.map.sources.sdo.AIAMap and type(aiafile) != sunpy.map.compositemap.CompositeMap:
            aiamap = sunpy.map.Map(aiafile)
        else:
            aiamap = aiafile
        if x_range and y_range:
            try:
                aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec), u.Quantity(y_range * u.arcsec))
            except:
                from astropy.coordinates import SkyCoord
                bl = SkyCoord(x_range[0] * u.arcsec, y_range[0] * u.arcsec, frame=aiamap.coordinate_frame)
                tr = SkyCoord(x_range[1] * u.arcsec, y_range[1] * u.arcsec, frame=aiamap.coordinate_frame)
                aiamap = aiamap.submap(bl, tr)

    if type(vlafile) == dict:
        if dspecdata:
            if dspecdata['stack'] == 'Vstack':
                fig = plt.figure(figsize=(9, 12))
                plt.subplot(2, 1, 1)
            elif dspecdata['stack'] == 'Hstack':
                fig = plt.figure(figsize=(20, 8))
                plt.subplot(1, 2, 1)
        else:
            if aiaplt_args['axes'] is None:
                fig = plt.figure(figsize=figsize)
                plt.subplot()
            else:
                fig = plt.gcf()
        if not vlaonly:
            aiamap.plot(**aiaplt_args)
            ax1 = plt.gca()
            ax1.set_autoscale_on(False)
        else:
            ax1 = aiaplt_args['axes']
            ax1.set_autoscale_on(False)

        if label:
            ax1.text(0.98, 0.98, label, horizontalalignment='right', verticalalignment='top', color='white', transform=ax1.transAxes, fontsize=14)
        clrange = vlafile['ColorMapper']['crange']
        if not maponly:
            cargsort = np.argsort(np.array(vlafile['ColorMapper']['c']))
            cargsort = cargsort[::zorder]
            vlafile['x'] = np.array(vlafile['x'])[cargsort]
            vlafile['y'] = np.array(vlafile['y'])[cargsort]
            vlafile['ColorMapper']['c'] = np.array(vlafile['ColorMapper']['c'])[cargsort]
            vlafile['s'] = np.array(vlafile['s'])[cargsort]
            im1 = ax1.scatter(vlafile['x'], vlafile['y'], c=vlafile['ColorMapper']['c'], s=vlafile['s'], vmin=clrange[0], vmax=clrange[-1], cmap=cmap,
                              **kwargs)
        else:
            im1 = ax1.scatter([], [], c=[], s=[], vmin=clrange[0], vmax=clrange[-1], cmap=cmap, **kwargs)
        if not dspecdata:
            if plt_clbar:
                cb1 = plt.colorbar(im1, orientation='vertical', ax=ax1)
                if type(vlafile) == dict:
                    if 'title' in vlafile['ColorMapper'].keys():
                        cb1.set_label(vlafile['ColorMapper']['title'])
                    else:
                        cb1.set_label('Frequency [GHz]')
                else:
                    cb1.set_label('Frequency [GHz]')
        else:
            tim = dspecdata['time']
            dt = np.mean(np.diff(tim))
            cmapspec = cmap
            cmapspec.set_bad('white', 1.0)
            normspec = colors.Normalize(vmin=dspecdata['drange'][0], vmax=dspecdata['drange'][1])
            if dspecdata['stack'] == 'Vstack':
                plt.subplot(2, 1, 2)
            elif dspecdata['stack'] == 'Hstack':
                plt.subplot(1, 2, 2)
            ax2 = plt.gca()
            im2 = plt.pcolormesh(tim, dspecdata['freq'], dspecdata['peak'], cmap=cmapspec, norm=normspec)
            ax2.add_patch(patches.Rectangle((vlafile['t'], dspecdata['freq'][0]), dt, dspecdata['freq'][-1] - dspecdata['freq'][0], facecolor='black',
                                            edgecolor='white', alpha=0.3))
            ax2.set_xlim(tim[0], tim[-1])
            ax2.set_ylim(dspecdata['freq'][0], dspecdata['freq'][-1])
            ax2.set_title('Vector Dynamic spectrum')
            labels = ax2.get_xticks().tolist()
            newlabels = [Time(lb / 24. / 3600., format='jd').iso.split(' ')[1] for lb in labels]
            ax2.set_xticklabels(newlabels, rotation=45)
            if plt_clbar:
                cb1 = plt.colorbar(im1, orientation='vertical', ax=ax1)
                if type(vlafile) == dict:
                    if 'title' in vlafile['ColorMapper'].keys():
                        cb1.set_label(vlafile['ColorMapper']['title'])
                    else:
                        cb1.set_label('Frequency [GHz]')
                else:
                    cb1.set_label('Frequency [GHz]')
                cb2 = plt.colorbar(im2, orientation='vertical', ax=ax2)
                cb2.set_label('Max intensity [Jy/beam]')
    else:
        if aiaplt_args['axes'] is None:
            fig = plt.figure(figsize=figsize)
            plt.subplot()
        else:
            fig = plt.gcf()
        if not vlaonly:
            aiamap.plot(**aiaplt_args)
            ax1 = plt.gca()
            ax1.set_autoscale_on(False)
        else:
            ax1 = aiaplt_args['axes']
            ax1.set_autoscale_on(False)
        # ax1.xaxis.set_ticks_position('top')
        # ax1.yaxis.set_ticks_position('right')
        if label:
            ax1.text(0.98, 0.98, label, horizontalalignment='right', verticalalignment='top', color='white', transform=ax1.transAxes, fontsize=14)
        # # Prevent the image from being re-scaled while overplotting.
        hdulist = fits.open(vlafile)
        hdu = hdulist[0]
        nfreq = hdu.data.shape[1]
        if not maponly:
            vladata = hdu.data[pol, 0, :, :]
            vlamap = sunpy.map.Map((vladata, hdu.header))
            XX, YY = np.meshgrid(np.arange(vlamap.data.shape[1]), np.arange(vlamap.data.shape[0]))
            mapxy = vlamap.pixel_to_data(XX * u.pix, YY * u.pix)
            mapx, mapy = mapxy.Tx.value, mapxy.Ty.value
            clrange = (hdu.header['CRVAL3'] + np.arange(nfreq) * hdu.header['CDELT3']) / 1e9
            if len(chans) == 0:
                chans = range(0, nfreq)
            if plotstyle == 'centroid':
                pltdata = {'x': [], 'y': [], 's': [], 'c': []}
                for idx, chan in enumerate(chans):
                    vladata = hdu.data[pol, chan, :, :]
                    if np.nanmax(vladata):
                        vlamap.data = vladata
                        maxxy = maxfit(vlamap, mapxy=[mapx, mapy], width=width)
                        if maxxy:
                            x, y = maxxy
                            pltdata['x'].append(x)
                            pltdata['y'].append(y)
                            pltdata['c'].append((hdu.header['CRVAL3'] + chan * hdu.header['CDELT3']) / 1e9)
                            pltdata['s'].append(np.nanmax(vladata))
                cargsort = np.argsort(np.array(pltdata['c']))
                cargsort = cargsort[::zorder]
                pltdata['x'] = np.array(pltdata['x'])[cargsort]
                pltdata['y'] = np.array(pltdata['y'])[cargsort]
                pltdata['c'] = np.array(pltdata['c'])[cargsort]
                pltdata['s'] = np.array(pltdata['s'])[cargsort]
                im1 = ax1.scatter(pltdata['x'], pltdata['y'], c=pltdata['c'], s=pltdata['s'], vmin=clrange[0], vmax=clrange[-1], cmap=cmap, **kwargs)
                if plt_clbar:
                    cb1 = plt.colorbar(im1, orientation='vertical', ax=ax1)
                    cb1.set_label('Frequency [GHz]')
            elif plotstyle == 'contour':
                nchan = len(chans)
                for idx, chan in enumerate(chans):
                    if chan_mask is not None:
                        freq = (hdu.header['CRVAL3'] + chan * hdu.header['CDELT3']) / 1e9
                        if freq not in chan_mask:
                            continue
                    vladata = hdu.data[pol, chan, :, :]
                    vlamap = sunpy.map.Map((vladata, hdu.header))
                    SRC_vlamap_contour = DButil.get_contour_data(mapx, mapy, vlamap.data, levels=levels)
                    if SRC_vlamap_contour.data['xs']:
                        for ii, xs in enumerate(SRC_vlamap_contour.data['xs']):
                            x, y = xs, SRC_vlamap_contour.data['ys'][ii]
                            if not thrshd or np.nanmax(vladata) >= thrshd:
                                plt.plot(x, y, color=cmap(int(float(chan) / nfreq * 255)), zorder=nchan + zorder * idx, **kwargs)
                if plt_clbar:
                    fig.subplots_adjust(right=0.8)
                    cax1 = fig.add_axes([0.85, 0.1, 0.01, 0.8])
                    norm = colors.Normalize(vmin=clrange[0], vmax=clrange[-1])
                    cb1 = colorbar.ColorbarBase(cax1, cmap=cmap, norm=norm, orientation='vertical')
                    cb1.set_label('Frequency [GHz]')
    try:
        cb1.ax.set_aspect(40)
        cb1.ax.tick_params(direction='in')
    except:
        pass
    try:
        cb2.ax.set_aspect(40)
        cb2.ax.tick_params(direction='in')
    except:
        pass

    if aiaplt_args['axes'] is None:
        fig.tight_layout(pad=3)
    if outfile:
        if outfile.endswith('.eps'):
            fig.savefig(outfile, format='eps', dpi=int(figdpi))
        elif outfile.endswith('.pdf'):
            fig.savefig(outfile, format='pdf', dpi=int(figdpi))
        elif outfile.endswith('.png'):
            fig.savefig(outfile, format='png', dpi=int(figdpi))
        elif outfile.endswith('.jpeg'):
            fig.savefig(outfile, format='jpeg', dpi=int(figdpi))
        elif outfile.endswith('.jpg'):
            fig.savefig(outfile, format='jpg', dpi=int(figdpi))
        else:
            raise ValueError('Can not save {}. Provide a filename with extension of eps, pdf, png or jpeg.'.format(outfile))
