from astropy.io import fits
import astropy.table
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_dspec(filename, doplot=False, vmax=None, vmin=None, norm=None, cmap=None):
    """
    Read EOVSA Dynamic Spectrum FITS file <filename> and return a spectrogram dictionary.
    Optionally show an overview plot if doplot switch is set.

    Example:
    --------
    >>> from suncasa.eovsa import eovsa_dspec as ds
    >>> from astropy.time import Time
    >>> from matplotlib.colors import LogNorm
    ## Read EOVSA Dynamic Spectrum FITS file <filename>
    >>> filename = 'EOVSA_TPall_20170713.fts'
    >>> s = ds.get_dspec(filename, doplot=True, cmap='gist_heat', norm=LogNorm(vmax=2.1e3, vmin=40))
    ## To access the data in the spectrogram object, use
    >>> spec = s['spectrogram']                    ## (Array of amplitudes in SFU, of size nfreq,ntimes)
    >>> fghz = s['spectrum_axis']                  ## (Array of frequencies in GHz, of size nfreq)
    >>> tim = Time(s['time_axis'], format='mjd')   ## (Array of UT times in astropy.time object, of size ntimes)

    Parameters
    ----------
    filename :   filename of the spectrogram fits file

    doplot : Boolean, optional

    vmin, vmax : scalar, optional
    When using scalar data and no explicit norm, vmin and vmax
    define the data range that the colormap covers. By default,
    the colormap covers the complete value range of the supplied data.
    vmin, vmax are ignored if the norm parameter is used.

    norm : `matplotib.colors Normalization object` or str, optional
    The Normalize instance used to scale scalar data to the [0, 1]
    range before mapping to colors using cmap. By default, a linear
    scaling mapping the lowest value to 0 and the highest to 1 is used.
    This parameter is ignored for RGB(A) data.

    cmap : `matplotib.colors.Colormap` or str
    A colormap instance or the name of a registered colormap.

    Returns
    -------
    spectrogram : dictionary
    """

    hdulist = fits.open(filename)
    spec = hdulist[0].data
    fghz = np.array(astropy.table.Table(hdulist[1].data)['sfreq'])
    tim = astropy.table.Table(hdulist[2].data)
    tmjd= np.array(tim['mjd']) + np.array(tim['time']) / 24. / 3600 / 1000
    tim = Time(tmjd, format='mjd')
    timplt = tim.plot_date
    ntim = len(timplt)
    nfreq = len(fghz)
    if doplot:
        fig, ax = plt.subplots(figsize=(7, 4))
        # if vmax is None:
        #     vmax = np.nanmax(spec)
        # if vmin is None:
        #     vmin = 0.1
        # if norm is None:
        #     norm = colors.Normalize(vmax=vmax, vmin=vmin)
        pcm = ax.pcolormesh(timplt, fghz, spec, norm=norm, vmax=vmax, vmin=vmin, cmap=cmap,rasterized='True')
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        ax.set_ylim(fghz[0], fghz[-1])
        ax.set_xlim(tim[0].plot_date, tim[-1].plot_date)
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('EOVSA Dynamic Spectrum for ' + tim[0].datetime.strftime('%Y-%m-%d'))
        divider = make_axes_locatable(ax)
        cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
        clb_spec = plt.colorbar(pcm, ax=ax, cax=cax_spec, label='Flux density [sfu]')
        cmap = pcm.get_cmap()
        bkg_facecolor=cmap(0.0)
        ax.set_facecolor(bkg_facecolor)
        dtim = np.diff(tmjd)
        idxs_tgap = np.where(dtim*24*60 >= 0.5)[0]
        if len(idxs_tgap)>0:
            for idx in idxs_tgap:
                if idx+1 < len(tmjd):
                    ax.axvspan(tim[idx].plot_date,tim[idx+1].plot_date,facecolor=bkg_facecolor,edgecolor='none')



        def format_coord(x, y):
            col = np.argmin(np.absolute(timplt - x))
            row = np.argmin(np.absolute(fghz - y))
            if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                timstr = tim[col].isot
                flux = spec[row, col]
                return 'time {0}, freq = {1:.3f} GHz, flux = {2:.2f} sfu'.format(timstr, y, flux)
            else:
                return 'x = {0}, y = {1:.3f}'.format(x, y)

        ax.format_coord = format_coord
        fig.tight_layout()
        plt.show()
    return {'spectrogram': spec, 'spectrum_axis': fghz, 'time_axis': tmjd}
