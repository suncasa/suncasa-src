from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_dspec(filenames, doplot=False, vmax=None, vmin=None, norm=None, cmap='viridis'):
    """
    Read one or more e-Callisto spectrogram FITS files and return a combined spectrogram dictionary.
    Optionally display a quicklook plot.

    Parameters
    ----------
    filenames : str or list of str
        Path(s) to e-Callisto FITS file(s).
    doplot : bool
        Whether to plot the dynamic spectrum.
    vmax, vmin : float
        Color scale limits.
    norm : matplotlib.colors.Normalize
        Custom normalization for the color scale.
    cmap : str or colormap
        Colormap name or object.

    Returns
    -------
    dict
        {
            'spectrogram': 2D array [nfreq, total_ntimes],
            'spectrum_axis': 1D array of frequencies [MHz],
            'time_axis': 1D array of MJD values
        }
    """
    if isinstance(filenames, str):
        filenames = [filenames]

    all_specs = []
    all_times = []
    freqs = None

    for i, filename in enumerate(filenames):
        hdulist = fits.open(filename)
        spec = hdulist[0].data  # shape: (nfreq, ntimes)
        header = hdulist[0].header
        info = hdulist[1].data[0]

        time_sec = np.array(info['TIME'])  # seconds of the day
        freqs_file = np.array(info['FREQUENCY'])  # in MHz

        # Convert to MJD using DATE-OBS + TIME-OBS in header
        date_obs_str = header.get('DATE-OBS', '').replace('/', '-') + 'T' + header.get('TIME-OBS', '00:00:00.000')
        t0 = Time(date_obs_str, format='isot', scale='utc')
        tmjd = t0.mjd + time_sec / 86400.0

        if i == 0:
            freqs = freqs_file
        else:
            if not np.allclose(freqs_file, freqs, atol=1e-3):
                print(f"Warning: frequency axis mismatch in file {filename}")

        all_specs.append(spec)
        all_times.append(tmjd)
        hdulist.close()

    # Concatenate along time axis
    spec_combined = np.concatenate(all_specs, axis=1)  # concatenate time axis
    time_combined = np.concatenate(all_times)

    # Optional plot
    if doplot:
        tim = Time(time_combined, format='mjd')
        timplt = tim.plot_date
        ntim = len(timplt)
        nfreq = len(freqs)

        fig, ax = plt.subplots(figsize=(9, 4))
        pcm = ax.pcolormesh(timplt, freqs, spec_combined, shading='auto', cmap=cmap,
                            norm=norm, vmax=vmax, vmin=vmin, rasterized=True)
        ax.set_title("e-Callisto Dynamic Spectrum")
        ax.set_xlabel("Time [UT]")
        ax.set_ylabel("Frequency [MHz]")
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        ax.set_ylim(freqs[0], freqs[-1])
        ax.set_xlim(tim[0].plot_date, tim[-1].plot_date)
        ax.set_facecolor(pcm.get_cmap()(0.0))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.05)
        plt.colorbar(pcm, cax=cax, label="Intensity [arb. units]")

        def format_coord(x, y):
            col = np.argmin(np.abs(timplt - x))
            row = np.argmin(np.abs(freqs - y))
            if 0 <= col < ntim and 0 <= row < nfreq:
                timstr = tim[col].isot
                flux = spec_combined[row, col]
                return f'time: {timstr}, freq: {y:.1f} MHz, intensity: {flux:.2f}'
            else:
                return f'x = {x:.3f}, y = {y:.1f}'

        ax.format_coord = format_coord
        fig.tight_layout()
        plt.show()

    return {
        'spectrogram': spec_combined,
        'spectrum_axis': freqs,
        'time_axis': time_combined
    }
