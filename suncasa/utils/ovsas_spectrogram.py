import os
from datetime import datetime, timedelta

def plot(timestamp=None, timerange=None, figdir='/common/lwa/spec_v2/daily/', combine=True,
         clip=[10, 99.995], add_logo=False, fast_plot=True, interactive=False):
    """
    Plot the OVRO-LWA and EOVSA spectrograms along with STIX and GOES light curves for a given timestamp or time range.

    :param timestamp: The timestamp for the data to be plotted. If not provided, the current UTC time is used.
    :type timestamp: datetime.datetime, optional
    :param timerange: The time range to be plotted. If provided, it overrides the timestamp.
    :type timerange: list of datetime.datetime, optional
    :param figdir: Directory where the figures will be saved, defaults to '/common/lwa/spec_v2/daily/'.
    :type figdir: str, optional
    :param combine: If True, combine all plots into a single figure, defaults to True. Otherwise, save each plot (OVRO-LWA, EOVSA, STIX, GOES) separately.
    :type combine: bool, optional
    :param clip: The percentile values for clipping the color scale of EOVSA and OVRO-LWA spectrograms, defaults to [10, 99.995].
    :type clip: list of float, optional
    :param add_logo: If True, add logos to the plots, defaults to False.
    :type add_logo: bool, optional
    :param fast_plot: If True, use fast plotting methods, defaults to True. If you want to plot the full resolution, set it to False. But it will take longer to plot. If the time range is less than 30 minutes, it will automatically set to False.
    :type fast_plot: bool, optional
    :param interactive: If False, suppress plot display and use the 'Agg' backend, defaults to False.
    :type interactive: bool, optional
    :return: List of file paths to the saved figures.
    :rtype: list of str

    Examples:
    ---------
    from suncasa.utils import ovsas_spectrogram as ovsp
    # Example 1: Plotting OVRO-LWA and EOVSA spectrograms along with STIX and GOES light curves for a specific timestamp
    ovsp.plot(datetime(2024, 7, 31), figdir='/data1/workdir/', combine=True)

    # Example 2: Plotting OVRO-LWA and EOVSA spectrograms along with STIX and GOES light curves for a specific time range
    ovsp.plot(timerange=[datetime(2024, 7, 31, 17, 20), datetime(2024, 7, 31, 20, 40)],
        figdir='/data1/workdir/', combine=True, fast_plot=True, clip=[5, 99.995])
    """
    import time
    t0 = time.time()

    import astropy.units as u
    import numpy as np
    import pandas as pd
    import sunpy.timeseries as ts
    from astropy.time import Time
    from stixdcpy.quicklook import LightCurves
    from sunpy.net import Fido
    from sunpy.net import attrs as a
    from sunpy.time import parse_time
    from suncasa.dspec import dspec as ds
    import matplotlib
    import matplotlib.colors as mcolors
    from matplotlib import pyplot as plt
    from matplotlib.dates import AutoDateFormatter, MinuteLocator
    from matplotlib.ticker import FuncFormatter
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if add_logo:
        from ovrolwasolar.visualization import njit_logo_str, nsf_logo
        import base64
        import io
        import matplotlib.image as mpimg

    if not interactive:
        try:
            matplotlib.use('Agg')  # Use the Agg backend
        except Exception as e:
            print(f'Error: {e}')
            print('Failed to set the backend to Agg. Use default backend.')
        plt.ioff()

    if timestamp is None:
        if timerange is not None:
            timestamp = timerange[0]
            if timerange[1] - timerange[0] < timedelta(minutes=30):
                fast_plot = False
        else:
            # Set the current timestamp if none is provided
            timestamp = datetime.utcnow()

    # Ensure the output directory exists
    os.makedirs(figdir, exist_ok=True)

    if combine:
        figname = os.path.join(figdir, f'fig-OVSAs_spec_{timestamp.strftime("%Y%m%d")}.png')
    else:
        # Define file names for each figure
        figname_eovsa = os.path.join(figdir, f'fig-eovsa_spec_{timestamp.strftime("%Y%m%d")}.png')
        figname_ovrolwa = os.path.join(figdir, f'fig-ovrolwa_spec_{timestamp.strftime("%Y%m%d")}.png')
        figname_stix = os.path.join(figdir, f'fig-stix_lc_{timestamp.strftime("%Y%m%d")}.png')
        figname_goes = os.path.join(figdir, f'fig-goes_lc_{timestamp.strftime("%Y%m%d")}.png')

    # Define paths for the data files
    ovrolwapath = '/common/lwa/spec_v2/fits/'
    ovrolwa_specfile = os.path.join(ovrolwapath, timestamp.strftime("%Y%m%d.fits"))
    eovsapath = f'/data1/eovsa/fits/synoptic/{timestamp.strftime("%Y/%m/%d")}/'
    eovsa_specfile = os.path.join(eovsapath, timestamp.strftime("EOVSA_TPall_%Y%m%d.fts"))

    # Define default time range if no data is available
    default_start_time = Time(timestamp.replace(hour=14, minute=0, second=0))
    default_end_time = Time((timestamp + timedelta(days=1)).replace(hour=2, minute=0, second=0))

    # Function to format the y-axis as integer frequencies
    def int_formatter(x, pos):
        return f"{int(x)}"

    # Function to set up the plot formatting
    def setup_time_axis(ax, start, end):
        locator = MinuteLocator(byminute=range(0, 60, 60))
        ax.xaxis.set_major_locator(locator)
        formatter = AutoDateFormatter(locator)
        formatter.scaled[1.0] = '%H:%M'  # For intervals of 1 day
        formatter.scaled[1 / 24] = '%H:%M'  # For intervals of 1 hour
        formatter.scaled[1 / (24 * 60)] = '%H:%M'  # For intervals of 1 minute
        formatter.scaled[1 / (24 * 60 * 60)] = '%H:%M'  # For intervals of 1 second
        ax.xaxis.set_major_formatter(formatter)
        minor_locator = MinuteLocator(byminute=range(0, 60, 10))
        ax.xaxis.set_minor_locator(minor_locator)
        ax.set_xlim(start.plot_date, end.plot_date)

    if combine:
        fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True, height_ratios=[2, 2, 3, 4])
        ax_goes, ax_stix, ax_eovsa, ax_ovrolwa = axs
    else:
        # Initialize figures and axes
        fig_ovrolwa, ax_ovrolwa = plt.subplots(1, 1, figsize=(15, 4))
        fig_eovsa, ax_eovsa = plt.subplots(1, 1, figsize=(15, 4))
        fig_stix, ax_stix = plt.subplots(1, 1, figsize=(15, 3))
        fig_goes, ax_goes = plt.subplots(1, 1, figsize=(15, 3))

    # Load OVRO-LWA spectrogram data
    print(f'processing OVRO-LWA spectrogram data for {timestamp.strftime("%Y-%m-%d")}')
    if os.path.exists(ovrolwa_specfile):
        cmap = 'viridis'
        ovrolwa_bkgtim = None
        d_ovrolwa = ds.Dspec()
        d_ovrolwa.read(ovrolwa_specfile, source='lwa')
        d_ovrolwa.pol = 'IV'
        timerange_ovrolwa = None
        if timerange is not None:
            timerange_ovrolwa = list(Time(timerange).iso)
        vmax = np.nanpercentile(d_ovrolwa.data, 99.995)
        vmin = np.nanpercentile(d_ovrolwa.data, 5)
        norm_I_ovrolwa = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        d_ovrolwa.plot(pol='I', timerange=timerange_ovrolwa, bkgtim=ovrolwa_bkgtim, plot_fast=fast_plot,
                       norm=norm_I_ovrolwa,
                       percentile=clip, minmaxpercentile=True, freq_unit='MHz', cmap=cmap, axes=ax_ovrolwa)
        ovrolwa_tim = d_ovrolwa.time_axis
        ovro_lwa_start, ovro_lwa_end = ovrolwa_tim[0], ovrolwa_tim[-1]
    else:
        ax_ovrolwa.plot([], [])
        ax_ovrolwa.text(0.5, 0.5, 'No OVRO-LWA data available', transform=ax_ovrolwa.transAxes,
                        ha='center', va='center', fontsize=12, color='gray')
        ax_ovrolwa.set_ylabel('Frequency [MHz]')
        ax_ovrolwa.set_ylim(29.033934, 83.871824)
        ovro_lwa_start, ovro_lwa_end = default_start_time, default_end_time
        divider = make_axes_locatable(ax_ovrolwa)
        cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
        cax_spec.set_visible(False)

    ax_ovrolwa.set_yscale('log')
    ax_ovrolwa.yaxis.set_minor_formatter(FuncFormatter(int_formatter))
    ax_ovrolwa.set_title(f'OVRO-LWA Dynamic Spectrum (Stokes I)')

    print(f'processing EOVSA spectrogram data for {timestamp.strftime("%Y-%m-%d")}')
    # Load EOVSA spectrogram data
    ax = ax_eovsa
    if os.path.exists(eovsa_specfile):
        cmap = 'viridis'
        eovsa_bkgtim = None
        timerange_eovsa = None
        if timerange is not None:
            timerange_eovsa = list(Time(timerange).iso)

        # norm_I_eovsa = mcolors.Normalize(vmin=None, vmax=None)
        d_eovsa = ds.Dspec()
        d_eovsa.read(eovsa_specfile, source='eovsa')
        data = d_eovsa.data
        bkgspec = np.nanpercentile(data, 1, axis=1)
        d_eovsa.data = data - bkgspec[:, np.newaxis]
        d_eovsa.data[d_eovsa.data < 0] = 0
        vmax = np.nanpercentile(d_eovsa.data, 99.995)
        vmin = np.nanpercentile(d_eovsa.data, 5)
        norm_I_eovsa = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        d_eovsa.plot(pol='I', timerange=timerange_eovsa, bkgtim=eovsa_bkgtim, plot_fast=False, norm=norm_I_eovsa,
                     percentile=clip, minmaxpercentile=True, freq_unit='GHz', cmap=cmap, axes=ax_eovsa)
        eovsa_tim = d_eovsa.time_axis
        eovsa_start, eovsa_end = eovsa_tim[0], eovsa_tim[-1]
    else:
        ax_eovsa.plot([], [])
        ax_eovsa.text(0.5, 0.5, 'No EOVSA data available', transform=ax_eovsa.transAxes,
                      ha='center', va='center', fontsize=12, color='gray')
        ax_eovsa.set_ylabel('Frequency [GHz]')
        ax_eovsa.set_ylim(1.1053711, 17.979687)
        eovsa_start, eovsa_end = default_start_time, default_end_time
        divider = make_axes_locatable(ax_eovsa)
        cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
        cax_spec.set_visible(False)

    ax_eovsa.set_yscale('log')
    ax_eovsa.yaxis.set_major_formatter(FuncFormatter(int_formatter))
    ax_eovsa.yaxis.set_minor_formatter(FuncFormatter(int_formatter))
    ax_eovsa.set_title(f'EOVSA Dynamic Spectrum (Stokes I)')

    if timerange is None:
        # Determine the overall time range based on available data
        overall_start = min(ovro_lwa_start, eovsa_start)
        overall_end = max(ovro_lwa_end, eovsa_end)
    else:
        ## if timerange is provided, override the overall time range
        overall_start = Time(timerange[0])
        overall_end = Time(timerange[1])
        figname_eovsa = os.path.join(figdir,
                                     f'fig-eovsa_spec_{timerange[0].strftime("%Y%m%dT%H%M%S")}-{timerange[1].strftime("%Y%m%dT%H%M%S")}.png')
        figname_ovrolwa = os.path.join(figdir,
                                       f'fig-ovrolwa_spec_{timerange[0].strftime("%Y%m%dT%H%M%S")}-{timerange[1].strftime("%Y%m%dT%H%M%S")}.png')
        figname_stix = os.path.join(figdir,
                                    f'fig-stix_lc_{timerange[0].strftime("%Y%m%dT%H%M%S")}-{timerange[1].strftime("%Y%m%dT%H%M%S")}.png')
        figname_goes = os.path.join(figdir,
                                    f'fig-goes_lc_{timerange[0].strftime("%Y%m%dT%H%M%S")}-{timerange[1].strftime("%Y%m%dT%H%M%S")}.png')

    print(f'processing STIX light curves for {timestamp.strftime("%Y-%m-%d")}')
    # Load STIX light curves
    lc = LightCurves.from_sdc(start_utc=overall_start.iso, end_utc=overall_end.iso, ltc=True)
    if len(lc.counts) > 0:
        lc.peek(ax=ax_stix)
    else:
        ax_stix.plot([], [])
        ax_stix.set_ylabel('Counts')
        ax_stix.set_ylim(20, 2e6)
        ax_stix.text(0.5, 0.5, 'No STIX data available', transform=ax_stix.transAxes,
                     ha='center', va='center', fontsize=12, color='gray')
        ax_stix.set_yscale('log')

    ax_stix.set_title(f'STIX Quick-look Light Curves')
    divider = make_axes_locatable(ax_stix)
    cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
    cax_spec.set_visible(False)
    if combine:
        pass
    else:
        fig_stix.subplots_adjust(bottom=0.15)

    print(f'processing GOES X-ray light curves for {timestamp.strftime("%Y-%m-%d")}')
    # Download GOES X-ray light curves
    goes_query = Fido.search(a.Time(overall_start.datetime, overall_end.datetime), a.Instrument('XRS'))
    goes_files = Fido.fetch(goes_query)
    if goes_files:
        goes = ts.TimeSeries(goes_files)
        goes.plot(axes=ax_goes, title=False)
        ax_goes.set_yscale('log')
        ax_goes.set_ylabel('Flux [W/m²]')
    else:
        goes_json_data = pd.read_json("https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json")
        # This will get us the short wavelength data.
        goes_short = goes_json_data[goes_json_data["energy"] == "0.05-0.4nm"]
        # This will get us the long wavelength data.
        goes_long = goes_json_data[goes_json_data["energy"] == "0.1-0.8nm"]

        time_array = parse_time(goes_short["time_tag"])
        filtered_indices = (time_array >= overall_start) & (time_array <= overall_end)
        filtered_short = goes_short[filtered_indices]
        filtered_long = goes_long[filtered_indices]

        # Create a DataFrame with the filtered data
        filtered_time_array = time_array[filtered_indices]
        goes_data = pd.DataFrame({
            "xrsa": filtered_short["flux"].values,
            "xrsb": filtered_long["flux"].values
        }, index=filtered_time_array.datetime)
        units = dict([("xrsa", u.W / u.m ** 2), ("xrsb", u.W / u.m ** 2)])
        meta = dict({"instrument": "GOES X-ray sensor", "measurements": "primary", "type": "quicklook"})
        # goes_data = pd.DataFrame({"xrsa": goes_short["flux"].values, "xrsb": goes_long["flux"].values},
        #                          index=time_array.datetime)

        goes_ts = ts.TimeSeries(goes_data, meta, units, source="xrs")
        goes_ts.plot(axes=ax_goes)
        # ax_goes.plot([], [])
        # ax_goes.text(0.5, 0.5, 'No GOES data available', transform=ax_goes.transAxes,
        #              ha='center', va='center', fontsize=12, color='gray')
        ax_goes.set_ylabel('Flux [W/m²]')
        # ax_goes.set_ylim(1e-9, 1e-3)
        divider = make_axes_locatable(ax_goes)
        cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
        cax_spec.set_visible(False)
    ax_goes.set_title(f'GOES X-ray Light Curves')
    if combine:
        pass
    else:
        fig_goes.subplots_adjust(bottom=0.15)

    # Set up the time axis for each plot
    setup_time_axis(ax_ovrolwa, overall_start, overall_end)
    setup_time_axis(ax_eovsa, overall_start, overall_end)
    setup_time_axis(ax_stix, overall_start, overall_end)
    setup_time_axis(ax_goes, overall_start, overall_end)

    for ax in [ax_ovrolwa, ax_eovsa, ax_stix, ax_goes]:
        ax.set_xlabel(f'Time [UTC, {timestamp.strftime("%Y %b %d")}]')


    if combine:
        for ax in axs[:-1]:
            ax.tick_params(axis='x', which='both', labelbottom=False)
            ax.set_xlabel('')
        fig.tight_layout()
        if add_logo:
            ax_logo1 = fig.add_axes([0.89, 0.91, 0.15, 0.08])
            img1 = base64.b64decode(njit_logo_str)
            img1 = io.BytesIO(img1)
            img1 = mpimg.imread(img1, format='png')
            ax_logo1.imshow(img1)
            ax_logo1.axis('off')

            ax_logo2 = fig.add_axes([0.81, 0.91, 0.15, 0.09])
            img2 = base64.b64decode(nsf_logo)
            img2 = io.BytesIO(img2)
            img2 = mpimg.imread(img2, format='png')
            ax_logo2.imshow(img2)
            ax_logo2.axis('off')
        fig.savefig(figname, dpi=300)
        print(f'Saved combined figure to {figname}')
        figout=[figname]
    else:
        if add_logo:
            for fig in [fig_ovrolwa, fig_eovsa, fig_stix, fig_goes]:
                ax_logo1 = fig.add_axes([0.89, 0.91, 0.15, 0.08])
                img1 = base64.b64decode(njit_logo_str)
                img1 = io.BytesIO(img1)
                img1 = mpimg.imread(img1, format='png')
                ax_logo1.imshow(img1)
                ax_logo1.axis('off')

                ax_logo2 = fig.add_axes([0.81, 0.91, 0.15, 0.09])
                img2 = base64.b64decode(nsf_logo)
                img2 = io.BytesIO(img2)
                img2 = mpimg.imread(img2, format='png')
                ax_logo2.imshow(img2)
                ax_logo2.axis('off')
        fig_eovsa.savefig(figname_eovsa, dpi=300)
        fig_ovrolwa.savefig(figname_ovrolwa, dpi=300)
        fig_stix.savefig(figname_stix, dpi=300)
        fig_goes.savefig(figname_goes, dpi=300)
        print(f'Saved EOVSA figure to {figname_eovsa}')
        print(f'Saved OVRO-LWA figure to {figname_ovrolwa}')
        print(f'Saved STIX figure to {figname_stix}')
        print(f'Saved GOES figure to {figname_goes}')
        figout=[figname_eovsa, figname_ovrolwa, figname_stix, figname_goes]
    #
    # plt.show()  # Optionally show all figures
    if not interactive:
        plt.close('all')

    t1 = time.time()
    print(f'Plotting took {t1 - t0:.2f} seconds.')
    return figout
