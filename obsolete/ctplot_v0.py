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
        XX, YY = np.meshgrid(np.arange(smap.data.shape[1]), np.arange(smap.data.shape[0]))
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
    dims = hdu.data.shape
    if len(dims) == 2:
        vladata = hdu.data
    elif len(dims) == 4:
        vladata = hdu.data[pol, chn, :, :]
    else:
        raise ValueError('check you import vla fits file')
    vlamap = sunpy.map.Map((vladata, hdu.header))
    cmap = sunpy.map.CompositeMap(aiamap)
    cmap.add_map(vlamap, levels=np.array(levels) * np.nanmax(vlamap.data))
    cmap.peek()


def plotmap(vlafile, aiafile, outfile='', label='', pol=0, chans=[], x_range=[], y_range=[], levels=[0.9],
            plotstyle='centroid', figsize=(10, 8), figdpi=100,
            zorder=-1, maponly=False, dspecdata={}, **kwargs):
    import matplotlib.gridspec as gridspec
    if outfile:
        plt.ioff()
    if type(aiafile) != sunpy.map.sources.sdo.AIAMap:
        aiamap = sunpy.map.Map(aiafile)
    else:
        aiamap = aiafile
    if x_range and y_range:
        aiamap = aiamap.submap(u.Quantity(x_range * u.arcsec),
                               u.Quantity(y_range * u.arcsec))

    if type(vlafile) == dict:
        if dspecdata:
            if dspecdata['stack'] == 'Vstack':
                gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.02], height_ratios=[1, 0.5])
                fig = plt.figure(figsize=(9, 12))
            elif dspecdata['stack'] == 'Hstack':
                gs = gridspec.GridSpec(1, 4, width_ratios=[1, 0.02, 1.5, 0.02], height_ratios=[1])
                fig = plt.figure(figsize=(18, 8))
            ax1 = plt.subplot(gs[0], projection=aiamap)
            ax2 = plt.subplot(gs[1])
            ax3 = plt.subplot(gs[2])
            ax4 = plt.subplot(gs[3])
        else:
            gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.05], height_ratios=[1])
            fig = plt.figure(figsize=figsize)
            ax1 = plt.subplot(gs[0], projection=aiamap)
            ax2 = plt.subplot(gs[1])
        aiamap.plot(axes=ax1)
        ax1.set_autoscale_on(False)

        if label:
            ax1.text(0.98, 0.98, label, horizontalalignment='right', verticalalignment='top', color='white',
                     transform=ax1.transAxes, fontsize=14)
        clrange = vlafile['ColorMapper']['crange']
        if not maponly:
            nchan = len(vlafile['ColorMapper']['c'])
            for idx, chan in enumerate(vlafile['ColorMapper']['c']):
                x, y = [vlafile['x'][idx], vlafile['y'][idx]] * u.arcsec
                s = vlafile['s'][idx]
                ax1.scatter(x.to(u.deg), y.to(u.deg), transform=ax1.get_transform('world'),
                            facecolor=cm.jet(int((chan - clrange[0]) / (clrange[-1] - clrange[0]) * 255)), edgecolor='black',
                            s=s, zorder=nchan + zorder * idx, **kwargs)
        if not dspecdata:
            cmap = cm.jet
            norm = colors.Normalize(vmin=clrange[0], vmax=clrange[-1])
            cb1 = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')
            if type(vlafile) == dict:
                cb1.set_label(vlafile['ColorMapper']['title'])
            else:
                cb1.set_label('Frequency [GHz]')
        else:
            tim = dspecdata['time']
            dt = np.mean(np.diff(tim))
            cmapspec = cm.jet
            cmapspec.set_bad('white', 1.0)
            normspec = colors.Normalize(vmin=dspecdata['drange'][0], vmax=dspecdata['drange'][1])
            ax3.pcolormesh(tim, dspecdata['freq'], dspecdata['peak'], cmap=cmapspec, norm=normspec)
            ax3.add_patch(patches.Rectangle((vlafile['t'], dspecdata['freq'][0]), dt,
                                            dspecdata['freq'][-1] - dspecdata['freq'][0], facecolor='black',
                                            edgecolor='white', alpha=0.3))
            ax3.set_xlim(tim[0], tim[-1] + dt)
            ax3.set_ylim(dspecdata['freq'][0], dspecdata['freq'][-1])
            ax3.set_title('Vector Dynamic spectrum')
            labels = ax3.get_xticks().tolist()
            newlabels = [Time(lb / 24. / 3600., format='jd').iso.split(' ')[1] for lb in labels]
            ax3.set_xticklabels(newlabels, rotation=45)
            if dspecdata['stack'] == 'Vstack':
                fig.tight_layout(h_pad=3, pad=3)
                ax3.set_ylabel('Frequency [GHz]')
                pos1 = ax1.get_position()
                pos2 = ax2.get_position()
                pos4 = ax4.get_position()
                pos2new = [pos2.x0 - (1.0 - pos1.x1) / 2, pos2.y0, pos2.width, pos2.height]
                ax2.set_position(pos2new)
                pos4new = [pos4.x0 - (1.0 - pos1.x1) / 2, pos4.y0, pos4.width, pos4.height]
                ax4.set_position(pos4new)
            elif dspecdata['stack'] == 'Hstack':
                fig.tight_layout(w_pad=0, pad=5)
                ax3.set_yticklabels([])
                pos2 = ax2.get_position()
                pos4 = ax4.get_position()
                pos2new = [pos2.x0 - 0.03, pos2.y0, pos2.width, pos2.height]
                ax2.set_position(pos2new)
                pos4new = [pos4.x0 - 0.01, pos4.y0, pos4.width, pos4.height]
                ax4.set_position(pos4new)
            cmap = cm.jet
            norm = colors.Normalize(vmin=clrange[0], vmax=clrange[-1])
            cb1 = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')
            cb2 = colorbar.ColorbarBase(ax4, cmap=cmapspec, norm=normspec, orientation='vertical')
            if type(vlafile) == dict:
                cb1.set_label(vlafile['ColorMapper']['title'])
            else:
                cb1.set_label('Frequency [GHz]')
            cb2.set_label('Max intensity [Jy/beam]')
    else:
        fig = plt.figure(figsize=figsize)
        ax1 = plt.subplot(projection=aiamap)
        aiamap.plot()
        if label:
            ax1.text(0.98, 0.98, label, horizontalalignment='right', verticalalignment='top', color='white',
                     transform=ax1.transAxes, fontsize=14)
        # # Prevent the image from being re-scaled while overplotting.
        ax1.set_autoscale_on(False)
        hdulist = fits.open(vlafile)
        hdu = hdulist[0]
        if not maponly:
            vladata = hdu.data[pol, 0, :, :]
            vlamap = sunpy.map.Map((vladata, hdu.header))
            XX, YY = np.meshgrid(np.arange(vlamap.data.shape[1]), np.arange(vlamap.data.shape[0]))
            mapx, mapy = vlamap.pixel_to_data(XX * u.pix, YY * u.pix)
            mapx, mapy = mapx.value, mapy.value
            nfreq = hdu.data[pol, :, :, :].shape[1]
            if len(chans) == 0:
                chans = range(0, nfreq)
            nchan = len(chans)
            if plotstyle == 'centroid':
                for idx, chan in enumerate(chans):
                    vladata = hdu.data[pol, chan, :, :]
                    if np.nanmax(vladata):
                        vlamap.data = vladata
                        maxxy = maxfit(vlamap, mapxy=[mapx, mapy])
                        if maxxy:
                            x, y = maxxy * u.arcsec
                            plt.scatter(x.to(u.deg), y.to(u.deg), transform=ax1.get_transform('world'),
                                        facecolor=cm.jet(int(float(chan) / nfreq * 255)), edgecolor='black',
                                        s=np.nanmax(vladata),
                                        zorder=nchan + zorder * idx, **kwargs)
            elif plotstyle == 'contour':
                for idx, chan in enumerate(chans):
                    vladata = hdu.data[pol, chan, :, :]
                    vlamap = sunpy.map.Map((vladata, hdu.header))
                    SRC_vlamap_contour = DButil.get_contour_data(mapx, mapy, vlamap.data, levels=levels)
                    if SRC_vlamap_contour.data['xs']:
                        xs = SRC_vlamap_contour.data['xs'][0] * u.arcsec
                        ys = SRC_vlamap_contour.data['ys'][0] * u.arcsec
                        plt.plot(xs.to(u.deg), ys.to(u.deg), transform=ax1.get_transform('world'),
                                 color=cm.jet(int(float(chan) / nfreq * 255)), zorder=nchan + zorder * idx, **kwargs)

        clrange = (hdu.header['CRVAL3'] + np.arange(nfreq) * hdu.header['CDELT3']) / 1e9
        fig.subplots_adjust(right=0.8)
        cax1 = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        cmap = cm.jet
        norm = colors.Normalize(vmin=clrange[0], vmax=clrange[-1])
        cb1 = colorbar.ColorbarBase(cax1, cmap=cmap, norm=norm, orientation='vertical')
        cb1.set_label('Frequency [GHz]')

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
            raise ValueError(
                'Can not save {}. Provide a filename with extension of eps, pdf, png or jpeg.'.format(outfile))
