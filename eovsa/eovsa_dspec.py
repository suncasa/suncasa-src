from astropy.io import fits
import astropy.table
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def dspec(tpfits, vmax=None, vmin=None, norm=None, cmap=None):
    hdulist = fits.open(tpfits)
    spec = hdulist[0].data
    fghz = np.array(astropy.table.Table(hdulist[1].data)['sfreq'])
    tim = astropy.table.Table(hdulist[2].data)
    tim = Time(np.array(tim['mjd']) + np.array(tim['time']) / 24. / 3600 / 1000, format='mjd')
    timplt = tim.plot_date
    ntim = len(timplt)
    nfreq = len(fghz)
    fig, ax = plt.subplots(figsize=(9, 4))
    if vmax is None:
        vmax = np.nanmax(spec)
    if vmin is None:
        vmin = 0.1
    if norm is None:
        norm = colors.Normalize(vmax=vmax, vmin=vmin)
    pcm = ax.pcolormesh(timplt, fghz, spec, norm=norm, cmap=cmap)
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
    return fig, ax
