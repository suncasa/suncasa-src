import argparse
import inspect
import logging
import os
import shutil
import socket
from datetime import datetime, time, timedelta
from glob import glob

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from functools import wraps
import numpy.ma as ma
import pandas as pd
import astropy.units as u
from astropy.time import Time
from scipy import ndimage
from sunpy import map as smap
from tqdm import tqdm

from suncasa.casa_compat import import_casatasks
from suncasa.io import ndfits
from suncasa.utils import helioimage2fits as hf
from suncasa.utils import mstools as mstl
from suncasa.eovsa import wrap_wsclean as ww
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.coordinates import Helioprojective, propagate_with_solar_surface

hostname = socket.gethostname()
is_on_server = hostname in ['pipeline', 'inti.hpcnet.campus.njit.edu']
is_on_inti = hostname == 'inti.hpcnet.campus.njit.edu'
script_mode = True

logging.basicConfig(level=logging.INFO, format='EOVSA pipeline: [%(levelname)s] - %(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

tasks = import_casatasks('gaincal', 'applycal', 'clearcal', 'delmod', 'ft', 'uvsub', 'split', 'concat', 'flagmanager',
                         'flagdata', 'tclean', 'hanningsmooth', 'imhead')
gaincal = tasks.get('gaincal')
applycal = tasks.get('applycal')
clearcal = tasks.get('clearcal')
delmod = tasks.get('delmod')
ft = tasks.get('ft')
uvsub = tasks.get('uvsub')
split = tasks.get('split')
concat = tasks.get('concat')
flagmanager = tasks.get('flagmanager')
flagdata = tasks.get('flagdata')
tclean = tasks.get('tclean')
hanningsmooth = tasks.get('hanningsmooth')
imhead = tasks.get('imhead')

from suncasa.casa_compat import import_casatools

tools = import_casatools(['qatool', 'iatool', 'cltool', 'mstool', 'tbtool'])
qatool = tools['qatool']
iatool = tools['iatool']
cltool = tools['cltool']
mstool = tools['mstool']
tbtool = tools['tbtool']
qa = qatool()
ia = iatool()
ms = mstool()
tb = tbtool()



def log_print(level, message):
    # Get the caller's stack frame and extract information
    frame = inspect.stack()[1]
    module_name = inspect.getmodulename(frame[1])
    function_name = frame[3]

    # Format the custom context information
    context = f"{module_name}::{function_name}"

    # Prepare the full log message
    full_message = f"{context} - {message}"

    # Fetch the logger
    logger = logging.getLogger()

    # Log the message with the specified level
    if level.upper() == 'WARNING':
        logger.warning(full_message)
    elif level.upper() == 'INFO':
        logger.info(full_message)
    elif level.upper() == 'DEBUG':
        logger.debug(full_message)
    elif level.upper() == 'ERROR':
        logger.error(full_message)
    elif level.upper() == 'CRITICAL':
        logger.critical(full_message)
    else:
        logger.info(full_message)  # Default to INFO if an unsupported level is given


def detect_noisy_images(data_stack, rsun_pix, crpix1, crpix2, rsun_ratio=[1.2, 1.3], rms_threshold=None,
                        snr_threshold=None, showplt=False):
    """
    Detect noisy images based on the RMS noise calculated from pixels outside a circular mask.

    This function computes the RMS noise for each image in the stack using only the pixels
    outside the circle defined by the center (crpix1, crpix2) and radius (rsun_pix * rsun_ratio).
    An image is flagged as noisy if its RMS exceeds the specified threshold. If no threshold is provided,
    it is set to the median RMS value plus one standard deviation.

    :param data_stack: 3D array of images (n_images, height, width).
    :type data_stack: numpy.ndarray
    :param rsun_pix: Radius (in pixels) of the reference circle.
    :type rsun_pix: int or float
    :param crpix1: X-coordinate of the circle center.
    :type crpix1: int or float
    :param crpix2: Y-coordinate of the circle center.
    :type crpix2: int or float
    :param rsun_ratio: Ratio to scale the reference circle radius, defaults to 1.2.
    :type rsun_ratio: float, optional
    :param rms_threshold: RMS threshold for flagging an image as noisy. If None, it is set to (median + std)
                          of the computed RMS values.
    :type rms_threshold: float, optional
    :param snr_threshold: SNR threshold for flagging an image as noisy. If None, it is set to (median - std)

    :param showplt: If True, plot the RMS values for each image.
    :type showplt: bool, optional

    :return: A tuple (noisy_indices, rms_values) where:
             - noisy_indices is a boolean array indicating which images are noisy.
             - rms_values is an array of computed RMS noise values for each image.
    :rtype: tuple(numpy.ndarray, numpy.ndarray)

    Example:
        noisy, rms = detect_noisy_images(data_stack, rsun_pix=100, crpix1=256, crpix2=256, rsun_ratio=1.2, rms_threshold=5.0)
    """
    # Transpose the data stack if necessary to have shape (n_images, height, width)
    images = np.squeeze(data_stack)
    ny, nx = images.shape[1], images.shape[2]
    y, x = np.ogrid[:ny, :nx]
    pix_dist = np.sqrt((x - crpix1) ** 2 + (y - crpix2) ** 2)
    mask_ring = (pix_dist > (rsun_pix * rsun_ratio[0])) & (pix_dist < (rsun_pix * rsun_ratio[1]))
    mask_disk = pix_dist <= (rsun_pix * 1.02)
    # Compute RMS noise for each image using only the pixels outside the mask
    rms_values = []
    snr_values = []
    for img in images:
        noise_region = img[mask_ring]
        # signal = np.nanmax(img[mask_disk])
        signal = np.nanpercentile(img[mask_disk], 99.999)
        rms = np.sqrt(np.nanmean(noise_region ** 2))
        rms_values.append(rms)
        snr = signal / rms
        snr_values.append(snr)
    rms_values = np.array(rms_values)
    snr_values = np.array(snr_values)

    # Determine threshold if not provided
    if rms_threshold is None:
        # rms_threshold = np.nanpercentile(rms_values, 80)
        median_rms = np.nanmedian(rms_values)
        std_rms = np.nanstd(rms_values)
        rms_threshold = median_rms + std_rms

    if snr_threshold is None:
        # snr_threshold = np.nanpercentile(snr_values, 20)
        median_snr = np.nanmedian(snr_values)
        std_snr = np.sqrt(np.nanmean((snr_values - median_snr) ** 2))
        snr_threshold = median_snr - 1 * std_snr
        print(f'std_snr={std_snr}, median_snr={median_snr}, snr_threshold={snr_threshold}')

    # Flag the images with RMS noise above the threshold or SNR below the threshold
    noisy_indices = (snr_values < snr_threshold) | (rms_values > rms_threshold)

    # prevent all images from being flagged as noisy
    select_indices = ~noisy_indices
    if np.sum(select_indices) == 0:
        # Get the indices of the top two SNR values
        top_two_indices = np.argsort(snr_values)[-2:]
        # Set the corresponding indices in select_indices to True
        select_indices[top_two_indices] = True

    # --- END NEW ---

    if showplt:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(8, 8))
        axs[0].plot(np.arange(len(rms_values)), rms_values, 'bo-', label='RMS Noise')
        axs[0].axhline(rms_threshold, color='r', linestyle='--', label=f'Threshold ({rms_threshold:.2f})')
        axs[0].set_xlabel('Image Index')
        axs[0].set_ylabel('RMS Noise')
        axs[0].set_title('RMS Noise Outside Masked Region')
        axs[0].legend()

        axs[1].plot(np.arange(len(snr_values)), snr_values, 'go-', label='SNR')
        axs[1].axhline(snr_threshold, color='r', linestyle='--', label=f'Threshold ({snr_threshold:.2f})')
        axs[1].set_xlabel('Image Index')
        axs[1].set_ylabel('SNR')
        axs[1].set_title('SNR (Signal = 99th percentile inside mask)')
        axs[1].legend()
        plt.tight_layout()
        plt.show()

    return select_indices, rms_values, snr_values


def is_factor_of_60_minutes(tdt):
    """
    Check if tdt is a factor of 60 minutes or a harmonic of 60 minutes.

    :param tdt: Time duration to check.
    :type tdt: timedelta
    :return: True if tdt is a factor or harmonic of 60 minutes, False otherwise.
    :rtype: bool
    """
    minutes = tdt.total_seconds() / 60
    return 60 % minutes == 0 or minutes % 60 == 0


def split_ms_by_scan(msfile, timerange='', outdir='./', antenna='0~12&&0~12', datacolumn='data', verbose=False):
    """
    Splits a Measurement Set (MS) file into sub-MS files based on scan numbers.

    This function reads the scan numbers from the MS file and creates a sub-MS file
    for each scan. The output sub-MS files are named using the scan number.

    :param msfile: Path to the Measurement Set (MS) file.
    :type msfile: str
    :param timerange: Time range to consider for splitting the MS, defaults to an empty string.
    :type timerange: str, optional
    :param antenna: Antenna selection string, defaults to '0~12&&0~12'.
    :type antenna: str, optional
    :param datacolumn: Data column to use for splitting, defaults to 'data'.
    :type datacolumn: str, optional
    :param outdir: Directory to save the sub-MS files, defaults to the current directory.
    :type outdir: str, optional
    :param verbose: If True, prints additional information during the process, defaults to False.
    :type verbose: bool, optional
    :return: List of paths to the sub-MS files created.
    :rtype: list[str]
    """
    ms.open(msfile)
    scans = ms.getscansummary()
    ms.close()

    sub_ms_files = []
    for scan in scans:
        run_start_time = datetime.now()
        scan_number = scan
        sub_ms_filename = os.path.basename(f"{msfile.replace('.ms', '')}.scan{scan_number}.ms")
        sub_ms_file = os.path.join(outdir, sub_ms_filename)
        sub_ms_files.append(sub_ms_file)
        if verbose:
            log_print('INFO', f"Splitting scan {scan_number} to {sub_ms_filename}")
        if os.path.exists(sub_ms_file):
            log_print('WARNING', f"Sub-MS file {sub_ms_filename} already exists. Skipping.")
            continue
        split(vis=msfile, outputvis=sub_ms_file, timerange=timerange, scan=str(scan_number), antenna=antenna,
              datacolumn=datacolumn)
        run_end_time = datetime.now()
        elapsed_time = run_end_time - run_start_time
        elapsed_time = elapsed_time.total_seconds() / 60
        if verbose:
            log_print('INFO', f"Splitting scan {scan_number} took {elapsed_time:.1f} minutes.")

    return sub_ms_files


def extract_time_flags_from_ms(msfile, threshold=0.75, plotflag=False, return_flag=False):
    """
    Extracts time steps and their corresponding flag status from a Measurement Set (MS).

    This function processes the flag data from all spectral windows (SPWs) in the MS file
    to compute the average flag status across baselines and channels. A time step is
    considered flagged if the mean flag value across all data exceeds a specified threshold.

    :param msfile: Path to the Measurement Set (MS) file.
    :type msfile: str
    :param threshold: Threshold for considering a time step as flagged, defaults to 0.75.
                      A time step is flagged if the mean flag fraction > threshold.
    :type threshold: float, optional
    :param plotflag: If True, generates a plot of the flag data and mean flag fraction, defaults to False.
    :type plotflag: bool, optional
    :return: A tuple containing:
             - `timedata` (numpy.ndarray): Array of time steps from the MS.
             - `timeflag` (numpy.ndarray): Boolean array where `True` indicates
                                           flagged time steps based on the threshold.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    t0 = datetime.now()

    ms.open(msfile)

    # Get spectral window information
    spwinfo = ms.getspectralwindowinfo()

    # Initialize storage for flag data
    flagdata_list = []
    timedata = None  # Initialize timedata to None for clarity

    # Loop through all spectral windows
    for spw_index, spw in enumerate(spwinfo.keys()):
        ms.selectinit(datadescid=0, reset=True)
        ms.selectinit(datadescid=spw_index)
        # Extract time data once (assume it's the same for all SPWs)
        if timedata is None:
            timedata = ms.getdata(['time'], ifraxis=True)['time']

        if spw_index == 0 and not return_flag:
            ms.close()
            t1 = datetime.now()
            log_print('INFO',
                      f"Extracting time flags from {os.path.basename(msfile)} took {(t1 - t0).total_seconds():.1f} seconds")
            return timedata, None
        data = ms.getdata(['flag'], ifraxis=True)
        flagdata_spw = data['flag']  # Shape: (npol, nchan, nbl, ntim)

        # Use only the first polarization (XX) and average across channels
        flagdata_mean = np.nanmean(flagdata_spw[0], axis=1)  # Shape: (nbl, ntim)
        flagdata_list.append(flagdata_mean)

    # Combine flag data from all SPWs
    flagdata_combined = np.vstack(flagdata_list)

    # Compute time flags based on the threshold
    flagdata_combined_mean = np.nanmean(flagdata_combined, axis=0)
    timeflag = flagdata_combined_mean > threshold

    # Plot results if required
    if plotflag:
        fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
        axs[0].pcolormesh(flagdata_combined, cmap='viridis', vmin=0, vmax=1)
        axs[0].set_title('Flag spectrogram')
        axs[0].set_ylabel('Frequency chan #')
        axs[1].plot(flagdata_combined_mean, label='Mean Flag Fraction', color='blue')
        axs[1].axhline(threshold, color='red', linestyle='--', label=f'Threshold ({threshold})')
        axs[1].set_title('Mean Flag Fraction Per Time Step')
        axs[1].set_xlabel('Time Index')
        axs[1].set_ylabel('Mean Flag Fraction')
        axs[1].legend()
        plt.tight_layout()
        plt.show()

    ms.close()
    t1 = datetime.now()
    log_print('INFO',
              f"Extracting time flags from {os.path.basename(msfile)} took {(t1 - t0).total_seconds():.1f} seconds")
    return timedata, timeflag


def get_timedelta(tim):
    """
    Calculates the most frequent time difference (delta) in minutes from a sequence of time values.

    This function computes the differences between consecutive Modified Julian Dates (MJDs)
    in the provided `tim` object, converts them to minutes, and identifies the most frequently
    occurring time difference.

    :param tim: An object with an `.mjd` attribute containing time values in Modified Julian Date (MJD) format.
    :type tim: object
    :return: The most frequently occurring time difference(s) in minutes.
    :rtype: numpy.ndarray
    """
    # Calculate time differences in minutes
    time_deltas = np.round(np.diff(tim.mjd) * 24 * 60)

    # Get unique values and their counts
    unique_values, counts = np.unique(time_deltas, return_counts=True)

    # Identify the most frequent time difference(s)
    most_frequent_deltas = unique_values[counts == counts.max()]

    return most_frequent_deltas[0]





def plot_style(func):
    """
    A decorator to apply common plotting styles to the plotting functions.
    """

    @wraps(func)
    def wrapper(self, filename=None, *args, **kwargs):
        # Create a figure with two rows
        fig, axs = plt.subplots(nrows=2, figsize=(6, 5), sharex=True)

        # Call the original plotting function
        func(self, axs, *args, **kwargs)

        # Apply common styling
        ax_major = axs[0]
        ax_minor = axs[1]

        # Customize the major intervals plot
        ax_major.set_ylim(0 - 0.3, 1 + 0.3)
        ax_major.set_yticks([0, 1])
        ax_major.set_yticklabels(['False', 'True'])
        ax_major.set_title('Time Intervals for WSClean Imaging')
        ax_major.legend(loc='best')

        # Customize the minor intervals plot
        ax_minor.set_title('Time Intervals for Solar Differential Rotation')
        ax_minor.set_ylim(0, self.nintervals_minor + 1)
        ax_minor.yaxis.set_major_locator(plt.MultipleLocator(1))
        ax_minor.set_ylabel('Sub-Interval Index')
        ax_minor.set_xlabel('Time [UT]')
        ax_minor.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

        plt.tight_layout()

        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    return wrapper


class WSCleanTimeIntervals:
    def __init__(self, msfile, combine_scans=True):
        """
        A class to manage and compute time intervals for WSClean imaging from a Measurement Set (MS).

        This class can process each scan individually or treat all scans as a single entity.

        :param msfile: Path to the Measurement Set (MS) file.
        :type msfile: str
        :param combine_scans: If True, treat all scans as a single entity. If False, process each scan individually.
        :type combine_scans: bool
        """
        self.msfile = msfile
        self.combine_scans = combine_scans
        self.scans = self._get_scans()
        self.scan_results = []
        self.tim = None
        self.tim_flag = None
        self.tdur = None
        self.tdt = None
        self.intervals_indices = None
        self.time_intervals = None
        self.nintervals_major = None
        self.nintervals_minor = 6  # Default value for minor intervals
        self.time_intervals_major_avg = None
        self.time_intervals_minor_avg = None

        if self.combine_scans:
            # Extract time data and flags for the entire MS
            self.tim, self.tim_flag = extract_time_flags_from_ms(self.msfile)
            self.tim = Time(self.tim / 24 / 3600, format='mjd')
            self.tdt = get_timedelta(self.tim)  # The most frequent time delta
            self.ntim = len(self.tim)
            self.tdur = self.tdt * self.ntim  # Total duration in minutes

    def _get_scans(self):
        """
        Extract scan information from the MS file.

        Returns:
            dict: Dictionary containing scan information.
        """
        ms.open(self.msfile)
        scans = ms.getscansummary()
        ms.close()
        return scans

    def _get_nintervals_major(self, tdur, interval_length=50):
        """
        Calculate the number of major intervals for WSClean.

        Args:
            tdur (float): Total duration of the valid time range for WSClean, in minutes.
            interval_length (int, optional): Length of each major interval in minutes. Default is 50.

        Returns:
            int: Number of major intervals.
        """
        thour = tdur / interval_length
        nintervals_major = int(np.round(thour))
        if nintervals_major == 0:
            nintervals_major = 1
            print(
                f'WARNING: The time duration for WSClean is less than {interval_length:.0f} min. Setting nintervals_major to 1.')
        return nintervals_major

    def compute_intervals(self, nintervals_minor=6):
        """
        Compute time intervals for WSClean imaging.

        Args:
            nintervals_minor (int, optional): Number of minor intervals within each major interval. Default is 6.
        """
        self.nintervals_minor = nintervals_minor
        if self.combine_scans:
            # Process all scans as a single entity
            self._compute_intervals_combined(nintervals_minor)
        else:
            # Process each scan individually
            self._compute_intervals_individual(nintervals_minor)

    def _compute_intervals_combined(self, nintervals_minor):
        """
        Compute time intervals for WSClean imaging, treating all scans as a single entity.
        """
        # Calculate number of major intervals
        self.nintervals_major = self._get_nintervals_major(self.tdur)

        # Divide the valid time range into equal intervals
        tim_indices = np.arange(self.ntim)
        self.intervals_indices = np.array_split(tim_indices, self.nintervals_major)
        self.time_intervals = [self.tim[seg] for seg in self.intervals_indices]
        self.time_intervals_major_avg = Time([np.nanmean(seg.mjd) for seg in self.time_intervals], format='mjd')

        # Compute minor intervals
        nintervals_all = self.nintervals_major * nintervals_minor
        if nintervals_all > self.ntim:
            print(
                f'WARNING: The full time span only has {self.ntim} time steps. The number of intervals for solar differential rotation correction is {nintervals_all}. Setting nintervals_all to {self.ntim}.')
            nintervals_all = self.ntim

        intervals_all_indices = np.array_split(tim_indices, nintervals_all)
        self.time_intervals_minor = []
        self.time_intervals_minor_avg = []
        for interval_idx in range(self.nintervals_major):
            start_idx = interval_idx * nintervals_minor
            end_idx = min((interval_idx + 1) * nintervals_minor, len(intervals_all_indices))
            interval_minor = [self.tim[seg] for seg in intervals_all_indices[start_idx:end_idx]]
            self.time_intervals_minor.append(interval_minor)
            self.time_intervals_minor_avg.append(Time([np.nanmean(seg.mjd) for seg in interval_minor], format='mjd'))

    def _compute_intervals_individual(self, nintervals_minor):
        """
        Compute time intervals for WSClean imaging for each scan individually.
        """
        for scan_id, scan_info in self.scans.items():
            for subscan_id, subscan_data in scan_info.items():
                # Extract scan information
                begin_time = subscan_data['BeginTime']
                end_time = subscan_data['EndTime']
                integration_time = subscan_data['IntegrationTime']

                # Calculate time indices for the scan
                tim = np.arange(begin_time, end_time, integration_time / 86400)  # Convert seconds to days
                tim = Time(tim, format='mjd')
                ntim = len(tim)
                tdur = (end_time - begin_time) * 1440  # Convert days to minutes

                # Calculate number of major intervals
                nintervals_major = self._get_nintervals_major(tdur)

                # Divide the valid time range into equal intervals
                tim_indices = np.arange(ntim)
                intervals_indices = np.array_split(tim_indices, nintervals_major)
                time_intervals = [tim[seg] for seg in intervals_indices]
                time_intervals_major_avg = Time([np.nanmean(seg.mjd) for seg in time_intervals], format='mjd')

                # Compute minor intervals
                nintervals_all = nintervals_major * nintervals_minor
                if nintervals_all > ntim:
                    print(
                        f'WARNING: The full time span only has {ntim} time steps. The number of intervals for solar differential rotation correction is {nintervals_all}. Setting nintervals_all to {ntim}.')
                    nintervals_all = ntim

                intervals_all_indices = np.array_split(tim_indices, nintervals_all)
                time_intervals_minor = []
                time_intervals_minor_avg = []
                for interval_idx in range(nintervals_major):
                    start_idx = interval_idx * nintervals_minor
                    end_idx = min((interval_idx + 1) * nintervals_minor, len(intervals_all_indices))
                    interval_minor = [tim[seg] for seg in intervals_all_indices[start_idx:end_idx]]
                    time_intervals_minor.append(interval_minor)
                    time_intervals_minor_avg.append(Time([np.nanmean(seg.mjd) for seg in interval_minor], format='mjd'))

                # Store results for this scan
                self.scan_results.append({
                    'ScanId': scan_id,
                    'SubScanId': subscan_id,
                    'BeginTime': begin_time,
                    'EndTime': end_time,
                    'IntegrationTime': integration_time,
                    'tdur': tdur,
                    'nintervals_major': nintervals_major,
                    'nintervals_minor': nintervals_minor,
                    'intervals_indices': intervals_indices,
                    'time_intervals': time_intervals,
                    'time_intervals_major_avg': time_intervals_major_avg,
                    'time_intervals_minor': time_intervals_minor,
                    'time_intervals_minor_avg': time_intervals_minor_avg
                })

    @plot_style
    def _plot_combined(self, axs):
        """
        Plot the time intervals for the entire MS file.
        """
        # Plot major intervals
        cmap = plt.cm.get_cmap('tab20', self.nintervals_major)
        for sidx, seg_times in enumerate(self.time_intervals):
            color = cmap(sidx)
            axs[0].plot(seg_times.plot_date, [1] * len(seg_times), marker='.', linestyle='', color=color,
                        label=f'Major Interval {sidx + 1}')

        # Plot minor intervals
        cmap2 = plt.cm.get_cmap('tab20', self.nintervals_major)
        y = np.arange(1, self.nintervals_minor + 1)
        for sidx, interval_minor in enumerate(self.time_intervals_minor):
            color = cmap2(sidx)
            for sridx, seg_times in enumerate(interval_minor):
                axs[1].plot(seg_times.plot_date, [y[sridx]] * len(seg_times), marker='.', linestyle='', color=color)

    @plot_style
    def _plot_individual(self, axs):
        """
        Plot the time intervals for each scan individually.
        """
        # Generate a colormap for major intervals
        nintervals_major_total = sum([scan_result['nintervals_major'] for scan_result in self.scan_results])
        cmap = plt.cm.get_cmap('tab20', nintervals_major_total)

        color_idx = 0  # Index to track colors across major and minor intervals
        for scan_result in self.scan_results:
            time_intervals = scan_result['time_intervals']
            time_intervals_minor = scan_result['time_intervals_minor']
            nintervals_major = scan_result['nintervals_major']
            nintervals_minor = scan_result['nintervals_minor']

            # Plot major intervals
            for sidx, seg_times in enumerate(time_intervals):
                color = cmap(color_idx)
                axs[0].plot(seg_times.plot_date, [1] * len(seg_times), marker='.', linestyle='', color=color,
                            label=f'Major Interval {color_idx + 1}')
                interval_minor = time_intervals_minor[sidx]
                # Plot minor intervals
                for sridx, seg_times_minor in enumerate(interval_minor):
                    axs[1].plot(seg_times_minor.plot_date, [sridx + 1] * len(seg_times_minor), marker='.',
                                linestyle='',
                                color=color)
                color_idx += 1

    def plot(self, filename=None):
        """
        Plot the time intervals and optionally save the plot to a file.

        Args:
            filename (str, optional): Path to save the plot. If None, the plot is displayed.
        """
        if self.combine_scans:
            if self.time_intervals is None:
                raise ValueError("Intervals have not been computed. Call `compute_intervals` first.")
            self._plot_combined(filename)
        else:
            if not self.scan_results:
                raise ValueError("Intervals have not been computed. Call `compute_intervals` first.")
            self._plot_individual(filename)

    def results(self):
        """
        Returns the computed intervals and related information.

        :return: If combine_scans=True, a dictionary containing:
                 - `intervals_indices` (list of lists): Indices of time steps in each interval.
                 - `time_intervals` (list of `astropy.Time`): Time values for each interval.
                 - `tdur` (float): Total duration of the valid time range for WSClean, in minutes.
                 - `tdt` (float): The most frequent time delta, in minutes.
                 - `nintervals_major` (int): Number of major intervals.
                 - `nintervals_minor` (int): Number of minor intervals within each major interval.
                 - `time_intervals_major_avg` (`astropy.Time`): Average time for each major interval.
                 - `time_intervals_minor_avg` (list of `astropy.Time`): Average time for each minor interval within major intervals.
                 If combine_scans=False, a list of dictionaries, each containing:
                 - `ScanId`: The ID of the scan.
                 - `SubScanId`: The ID of the sub-scan.
                 - `BeginTime`: The start time of the scan.
                 - `EndTime`: The end time of the scan.
                 - `IntegrationTime`: The integration time of the scan.
                 - `tdur`: Total duration of the valid time range for WSClean, in minutes.
                 - `nintervals_major`: Number of major intervals.
                 - `nintervals_minor`: Number of minor intervals within each major interval.
                 - `intervals_indices`: Indices of time steps in each major interval.
                 - `time_intervals`: Time values for each major interval.
                 - `time_intervals_major_avg`: Average time for each major interval.
                 - `time_intervals_minor`: Time values for each minor interval.
                 - `time_intervals_minor_avg`: Average time for each minor interval within major intervals.
        :rtype: dict or list of dict
        """
        if self.combine_scans:
            return {
                'tdur': self.tdur,
                'tdt': self.tdt,
                'intervals_indices': self.intervals_indices,
                'time_intervals': self.time_intervals,
                'nintervals_major': self.nintervals_major,
                'nintervals_minor': self.nintervals_minor,
                'time_intervals_major_avg': self.time_intervals_major_avg,
                'time_intervals_minor_avg': self.time_intervals_minor_avg
            }
        else:
            return self.scan_results


def trange2timerange(trange):
    """
    Convert a time range tuple in datetime format to a string representation.

    :param trange: A tuple containing start and end times as datetime objects.
    :type trange: tuple
    :return: A string representation of the time range.
    :rtype: str
    """

    sttime, edtime = trange
    if isinstance(sttime, str):
        sttime = pd.to_datetime(sttime)
    if isinstance(edtime, str):
        edtime = pd.to_datetime(edtime)
    timerange = f"{sttime.strftime('%Y/%m/%d/%H:%M:%S')}~{edtime.strftime('%Y/%m/%d/%H:%M:%S')}"
    return timerange


def rotateimage(data, xc_centre, yc_centre, p_angle):
    """
    Rotate an image around a specified point (xc_centre, yc_centre) by a given angle.

    :param data: The image data.
    :type data: numpy.ndarray
    :param xc_centre: The x-coordinate of the rotation center.
    :type xc_centre: int
    :param yc_centre: The y-coordinate of the rotation center.
    :type yc_centre: int
    :param p_angle: The rotation angle in degrees.
    :type p_angle: float
    :return: The rotated image.
    :rtype: numpy.ndarray
    """

    padX = [data.shape[1] - xc_centre, xc_centre]
    padY = [data.shape[0] - yc_centre, yc_centre]
    imgP = np.pad(data, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, p_angle, reshape=False, order=0, prefilter=False)
    return imgR[padY[0]:-padY[1], padX[0]:-padX[1]]


def solar_diff_rot_heliofits(in_fits, newtime, out_fits, in_time=None, template_fits=None, showplt=False,
                             overwrite_prev=True):
    """
    Reproject a FITS file to account for solar differential rotation to a new observation time.

    Parameters
    ----------
    in_fits : str
        Path to the input FITS file to be reprojected
    newtime : astropy.time.Time
        The new time to which the map is reprojected
    out_fits : str
        Path for the output FITS file
    in_time : astropy.time.Time, optional
        The reference time for the input map, defaults to the middle of the exposure time. This is useful when exposure time is not provided in the FITS header.
    template_fits : str, optional
        Path to template FITS file for output. If None, uses in_fits as template
    showplt : bool, optional
        Show plots of the original and reprojected maps, defaults to False

    Returns
    -------
    str
        Path to the output FITS file
    """

    # Convert FITS to SunPy Map
    in_map = smap.Map(in_fits)

    # Calculate reference time and output time
    if in_time is None:
        in_time = in_map.date  # + in_map.exposure_time / 2
    out_time = in_map.date + (Time(newtime) - in_time)
    # out_time = Time(newtime)

    # Set up output frame
    out_frame = Helioprojective(observer=in_map.observer_coordinate,
                                obstime=out_time,
                                rsun=in_map.coordinate_frame.rsun)

    # Define output center and reference pixel
    out_center = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=out_frame)
    out_ref_pixel = [in_map.reference_pixel.x.value,
                     in_map.reference_pixel.y.value] * in_map.reference_pixel.x.unit

    # Create output header
    out_header = smap.make_fitswcs_header(in_map.data.shape,
                                          out_center,
                                          reference_pixel=out_ref_pixel,
                                          scale=u.Quantity(in_map.scale))
    out_wcs = WCS(out_header)

    # Perform reprojection
    with propagate_with_solar_surface():
        out_map, footprint = in_map.reproject_to(out_wcs, return_footprint=True)

    # Fill missing data from original map
    out_data = out_map.data
    out_data[footprint == 0] = in_map.data[footprint == 0]

    # Set NaN values to 0
    out_data[np.isnan(out_data)] = 0

    # Create new map with reprojected data
    out_map = smap.Map(out_data, out_map.meta)
    out_map.meta['p_angle'] = in_map.meta['p_angle']

    # Display plots if requested
    if showplt:
        fig = plt.figure(figsize=(12, 4))

        ax1 = fig.add_subplot(121, projection=in_map)
        in_map.plot(axes=ax1, title='Original map')
        plt.colorbar()

        ax2 = fig.add_subplot(122, projection=out_map)
        out_map.plot(axes=ax2,
                     title=f"Reprojected to an Earth observer {(Time(newtime) - in_time).to('day')} later")
        plt.colorbar()
        plt.show()

    # Save to FITS file using template if provided
    if template_fits is None:
        template_fits = in_fits

    # Load template header
    template_map = smap.Map(template_fits)
    template_header = template_map.meta.copy()

    # Update necessary header information
    template_header.update({
        'DATE-OBS': out_map.date.iso,
        'DATE': out_map.date.iso,
        'CRVAL1': out_map.reference_coordinate.Tx.value,
        'CRVAL2': out_map.reference_coordinate.Ty.value,
        'CRPIX1': out_map.reference_pixel.x.value,
        'CRPIX2': out_map.reference_pixel.y.value,
        'PC1_1': out_map.rotation_matrix[0, 0],
        'PC1_2': out_map.rotation_matrix[0, 1],
        'PC2_1': out_map.rotation_matrix[1, 0],
        'PC2_2': out_map.rotation_matrix[1, 1]
    })

    # Create new map with template header and save to FITS
    final_map = smap.Map(out_data, template_header)
    outdir = os.path.dirname(out_fits)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    final_map.save(out_fits, overwrite=overwrite_prev)

    return out_fits


def sunpyfits_to_j2000fits(in_fits, out_fits, template_fits=None, overwrite_prev=True):
    """
    Rotate a solar FITS file from helioprojective to RA-DEC coordinates and save to a new FITS file.
    Originally written by Dr. Peijin Zhang
    Parameters
    ----------
    in_fits : str
        Path to input FITS file in helioprojective coordinates
    out_fits : str
        Path for output FITS file in RA-DEC coordinates
    template_fits : str, optional
        Path to template FITS file for output format. If None, uses in_fits as template
    overwrite_prev : bool, optional
        If True, overwrites existing output file. Defaults to True

    Returns
    -------
    str
        Path to the output FITS file
    """
    import astropy.units as u
    from astropy.io import fits

    # Load input map
    in_map = smap.Map(in_fits)
    p_ang = in_map.meta['p_angle'] * u.deg

    # Get reference pixel coordinates
    ref_x = int(in_map.reference_pixel.x.value)
    ref_y = int(in_map.reference_pixel.y.value)

    # Perform rotation
    data_rot = rotateimage(in_map.data, ref_x, ref_y, -p_ang.to('deg').value)

    # Set NaN values to 0
    data_rot[np.isnan(data_rot)] = 0

    # Use template if provided, otherwise use input file
    if template_fits is None:
        template_fits = in_fits

    # Read template header
    with fits.open(template_fits) as hdul:
        template_header = hdul[0].header.copy()

    # Update necessary header information
    template_header.update({
        'DATE-OBS': in_map.date.iso,
        'DATE': in_map.date.iso,
    })

    # Create and write output FITS file
    hdu = fits.PrimaryHDU(data_rot, header=template_header)
    hdu.writeto(out_fits, overwrite=overwrite_prev)

    return out_fits


from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
import numpy as np


def get_beam_kernel(header=None, bmaj=None, bmin=None, bpa=None, pixel_scale_arcsec=None):
    """
    Returns a normalized 2D elliptical Gaussian kernel representing the beam.

    The beam can be specified in two ways:

    1. By providing a FITS header (via the 'header' parameter). In this case, if not explicitly provided,
       the beam parameters are extracted from the header using the following keywords:
         - 'BMAJ' (FWHM of the major axis, in degrees)
         - 'BMIN' (FWHM of the minor axis, in degrees)
         - 'BPA'  (beam position angle in degrees, measured counterclockwise from North)
       In addition, if 'pixel_scale_arcsec' is not provided, it is determined from the header using:
         - 'CDELT1' (pixel scale, in degrees or arcseconds)
         - 'CUNIT1' (unit of CDELT1)

    2. By directly providing the beam parameters and pixel scale:
         - bmaj: FWHM of the major axis in degrees.
         - bmin: FWHM of the minor axis in degrees.
         - bpa:  Beam position angle in degrees (counterclockwise from North).
         - pixel_scale_arcsec: Pixel scale in arcseconds per pixel.

    The function converts the FWHM values (in degrees) to arcseconds and then to pixels using the provided
    pixel scale. It then converts the FWHM to sigma (standard deviation) via sigma = FWHM / 2.35482. The beam
    position angle is converted to a rotation angle (theta) relative to the x-axis via theta = radians(90 - BPA).

    :param header: FITS header containing the beam and pixel scale information, defaults to None.
    :type header: astropy.io.fits.header.Header, optional
    :param bmaj: FWHM of the beam major axis in degrees, defaults to None.
    :type bmaj: float, optional
    :param bmin: FWHM of the beam minor axis in degrees, defaults to None.
    :type bmin: float, optional
    :param bpa: Beam position angle in degrees (counterclockwise from North), defaults to None.
    :type bpa: float, optional
    :param pixel_scale_arcsec: Pixel scale in arcseconds per pixel, defaults to None.
    :type pixel_scale_arcsec: float, optional
    :raises KeyError: If a header is provided but required keywords ('BMAJ', 'BMIN', 'BPA', 'CDELT1', or 'CUNIT1') are missing.
    :raises ValueError: If neither a header nor the required beam parameters and pixel_scale_arcsec are provided, or if 'CUNIT1' is unrecognized.
    :return: A 2D numpy array representing the normalized elliptical Gaussian beam kernel.
    :rtype: numpy.ndarray
    """
    if header is not None:
        # If a header is provided, extract beam parameters and pixel scale first.
        try:
            if header["BMAJ"] > 0 and header["BMIN"] > 0:
                bmaj = header["BMAJ"]
                bmin = header["BMIN"]
                bpa = header["BPA"]
            else:
                print(f"Beam parameters 'BMAJ' and 'BMIN' must be positive. Use the provided values: ")
        except KeyError as e:
            raise KeyError(
                "Beam parameters 'BMAJ', 'BMIN', and 'BPA' must be provided either as inputs or in the header.") from e
        if pixel_scale_arcsec is None:
            try:
                cdelt1 = header["CDELT1"]
                cunit1 = header["CUNIT1"].strip().lower()
            except KeyError as e:
                raise KeyError(
                    "When header is provided, it must contain 'CDELT1' and 'CUNIT1' for pixel scale conversion.") from e
            if cunit1 == "deg":
                pixel_scale_arcsec = abs(cdelt1) * 3600.0
            elif cunit1 in ("arcsec", "arcsec."):
                pixel_scale_arcsec = abs(cdelt1)
            else:
                raise ValueError(f"Unrecognized CUNIT1: expected 'deg' or 'arcsec', got '{cunit1}'.")
    else:
        # If no header is provided, all beam parameters must be explicitly given.
        if bmaj is None or bmin is None or bpa is None or pixel_scale_arcsec is None:
            raise ValueError("Either provide a header or explicitly supply bmaj, bmin, bpa, and pixel_scale_arcsec.")

    # Convert FWHM from degrees to arcseconds.
    bmaj_arcsec = bmaj * 3600.0
    bmin_arcsec = bmin * 3600.0

    # Convert FWHM from arcseconds to pixels.
    fwhm_major_pixels = bmaj_arcsec / pixel_scale_arcsec
    fwhm_minor_pixels = bmin_arcsec / pixel_scale_arcsec

    # Convert FWHM to standard deviation (sigma) in pixels.
    sigma_major = fwhm_major_pixels / 2.35482
    sigma_minor = fwhm_minor_pixels / 2.35482

    # Convert beam position angle (BPA) to rotation angle (theta) relative to the x-axis.
    # BPA is measured counterclockwise from North (the y-axis), so theta = radians(90 - BPA).
    theta = np.deg2rad(90.0 - bpa)

    # Determine an appropriate kernel size (e.g., 8 sigma to capture most of the Gaussian)
    size = int(np.ceil(8 * max(sigma_major, sigma_minor)))
    if size % 2 == 0:
        size += 1  # ensure the kernel size is odd

    # Create the Gaussian kernel using astropy.convolution.Gaussian2DKernel.
    # The kernel is normalized so that its integral (sum) is unity.
    kernel = Gaussian2DKernel(x_stddev=sigma_major, y_stddev=sigma_minor, theta=theta,
                              x_size=size, y_size=size, mode='oversample')
    return kernel.array


def add_convolved_disk_to_fits(in_fits, out_fits, dsz, fdn, ignore_data=False, toTb=False, rfreq=None, add_limb=False, bmaj=None):
    """
    Adds a circular disk to the image from a FITS file, optionally applies a 10%-increased limb,
    convolves the disk with the beam, and writes the modified image to a new FITS file.

    Both dsz and fdn can be either a float or a list of floats. If a list is provided for fdn,
    the final value used is the mean (ignoring NaNs). If a list is provided for dsz, a disk image
    is calculated for each value and then averaged to obtain the final disk image.

    If both in_fits and out_fits are provided as lists, the same computed disk image is used for all
    input files.

    :param in_fits: Path to the input FITS file, or list of paths.
    :param out_fits: Path to the output FITS file, or list of paths.
    :param dsz: Disk radius in arcseconds, or list of radii.
    :param fdn: The uniform value to add within the disk, or list of values.
    :param ignore_data: If True, ignore the image data and only add the convolved disk.
    :param add_limb: If True, add a limb to the disk by increasing the intensity of the annulus by 60%.
    :raises KeyError: If required header keywords are missing.
    :raises ValueError: If image data is not at least 2D or if CUNIT1 is unrecognized.
    :return: None
    """
    from scipy.signal import fftconvolve

    import numbers

    def is_number(x):
        return isinstance(x, numbers.Number)

    # If fdn is a list, compute its mean (ignoring NaNs)
    if not is_number(fdn):
        fdn_val = np.nanmean(fdn)
    else:
        fdn_val = fdn

    # Function to compute disk image given a header, a disk radius value (in arcsec),
    # and the uniform disk value (fdn_val).
    def compute_disk_image(header, dsz_value, fdn_value, add_limb, bmaj=None):
        # Extract pixel scale from header.
        try:
            cdelt1 = header["CDELT1"]
            cunit1 = header["CUNIT1"].strip().lower()
        except KeyError as e:
            raise KeyError("Header must contain 'CDELT1' and 'CUNIT1'.") from e
        if cunit1 == "deg":
            pixel_scale_arcsec = abs(cdelt1) * 3600.0
        elif cunit1 in ("arcsec", "arcsec."):
            pixel_scale_arcsec = abs(cdelt1)
        else:
            raise ValueError(f"Unrecognized CUNIT1: expected 'deg' or 'arcsec', got '{cunit1}'.")

        # Calculate disk radius in pixels.
        disk_radius_pixels = dsz_value / pixel_scale_arcsec

        # Get image dimensions.
        ny, nx = header["NAXIS2"], header["NAXIS1"]

        # Get center coordinates (FITS convention: 1-indexed).
        try:
            crpix1 = header["CRPIX1"]
            crpix2 = header["CRPIX2"]
        except KeyError as e:
            raise KeyError("Header must contain 'CRPIX1' and 'CRPIX2'.") from e
        x_center = crpix1 - 1
        y_center = crpix2 - 1

        # Create coordinate grid.
        y_indices, x_indices = np.ogrid[:ny, :nx]
        distance = np.sqrt((x_indices - x_center) ** 2 + (y_indices - y_center) ** 2)

        # Create disk mask and compute uniform disk value per pixel so that total flux is fdn_value.
        disk_mask = distance <= disk_radius_pixels
        n_disk_pixels = np.count_nonzero(disk_mask)
        disk_value = fdn_value / n_disk_pixels
        disk_img = np.zeros((ny, nx), dtype=np.float64)
        disk_img[disk_mask] = disk_value

        if add_limb:
            # Compute limb width: 6 times the beam FWHM in pixels.
            if 'BMAJ' in header:
                if  header["BMAJ"]>0:
                    bmaj = header["BMAJ"]
                else:
                    if bmaj is None:
                        raise KeyError("Header must contain 'BMAJ' for beam FWHM computation.")
                    else:
                        print(f"Beam parameters 'BMAJ' must be positive. Use the provided value: {bmaj}")
            beam_fwhm_arcsec = bmaj * 3600.0
            beam_fwhm_pixels = beam_fwhm_arcsec / pixel_scale_arcsec
            limb_width = 6 * beam_fwhm_pixels
            # Limb is defined as annulus: (disk_radius_pixels - limb_width) < distance <= disk_radius_pixels.
            limb_mask = (distance > (disk_radius_pixels - limb_width)) & (distance <= disk_radius_pixels)
            disk_img[limb_mask] *= 1.6
        return disk_img

    if isinstance(in_fits, (list, tuple)):
        with fits.open(in_fits[0]) as hdul:
            header = hdul[0].header
    else:
        with fits.open(in_fits) as hdul:
            header = hdul[0].header

    # Compute disk image.
    if not is_number(dsz):
        # dsz is a list: compute a disk image for each and take the mean.
        disk_images = []
        for d in dsz:
            disk_images.append(compute_disk_image(header, d, fdn_val, add_limb, bmaj))
        disk_img_final = np.nanmean(np.array(disk_images), axis=0)
    else:
        # dsz is a single float.
        disk_img_final = compute_disk_image(header, dsz, fdn_val, add_limb, bmaj)

    # Obtain the beam kernel from the header (assumed to be the same for all files).
    beam_kernel = get_beam_kernel(header=header, bmaj=bmaj, bmin=bmaj, bpa=0)

    # Convolve the disk image with the beam.
    convolved_disk = fftconvolve(disk_img_final, beam_kernel, mode='same')
    convolved_disk = convolved_disk / np.nansum(convolved_disk) * fdn_val
    if toTb:
        from scipy import constants
        cell = (header['cdelt1'] * u.Unit(header['cunit1']) + header['cdelt2'] * u.Unit(header['cunit2'])) / 2.0
        k_b = constants.k
        c_l = constants.c
        const = 2. * k_b / c_l ** 2
        pix_area = (cell.to(u.rad).value) ** 2
        jy_to_si = 1e-26
        factor2 = 1.
        freq = float(rfreq.rstrip('GHz')) * 1e9
        factor = const * freq ** 2  # SI unit
        jy2tb = jy_to_si / pix_area / factor * factor2
        convolved_disk = convolved_disk * jy2tb
        print(f'tb disk: {np.nanmax(convolved_disk):.0f} K')

    # print(fdn_val, np.nansum(convolved_disk))
    # Define a helper function to process a single file.
    def process_file(in_file, out_file):
        with fits.open(in_file) as hdul:
            data = hdul[0].data.copy()
            hdr = hdul[0].header
        if ignore_data:
            data_modified = convolved_disk
        else:
            # print(f"Adding disk to {in_file}. max disk {np.nanmax(convolved_disk)}")
            data_modified = data + convolved_disk
        fits.writeto(out_file, data_modified, hdr, overwrite=True)
        # print(f"Disk data write to {out_file}.")

    # If in_fits and out_fits are lists, apply the same disk image to all.
    if isinstance(in_fits, (list, tuple)) and isinstance(out_fits, (list, tuple)):
        for infile, outfile in zip(in_fits, out_fits):
            process_file(infile, outfile)
    else:
        process_file(in_fits, out_fits)


def calc_diskmodel(slashdate, nbands, freq, defaultfreq):
    from astropy.time import Time
    # Default disk size measured for 2019/09/03
    # todo add monthly fitting procedure for the disk size and flux density
    defaultsize = np.array([990.6, 989.4, 988.2, 987.1, 986.0, 984.9, 983.8, 982.7, 981.7, 980.7,
                            979.7, 978.8, 977.8, 976.9, 976.0, 975.2, 974.3, 973.5, 972.7, 972.0,
                            971.2, 970.5, 969.8, 969.1, 968.5, 967.8, 967.2, 966.7, 966.1, 965.6,
                            965.1, 964.6, 964.1, 963.7, 963.3, 962.9, 962.5, 962.1, 961.8, 961.5,
                            961.3, 961.0, 960.8, 960.6, 960.4, 960.2, 960.1, 960.0, 959.9, 959.8])

    # Get current solar distance and modify the default size accordingly
    try:
        from sunpy.coordinates.sun import earth_distance
        fac = earth_distance('2019/09/03') / earth_distance(slashdate)
    except:
        import sunpy.coordinates.ephemeris as eph
        fac = eph.get_sunearth_distance('2019/09/03') / eph.get_sunearth_distance(slashdate)

    newsize = defaultsize * fac.to_value()
    if nbands == 34:
        if Time(slashdate.replace('/', '-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectral windows (bands 4-34 -> spw 0-30)
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[[0] + list(range(2, 34))])
    dsize = np.array([str(i)[:5] + 'arcsec' for i in newsize], dtype='U12')

    # These are nominal flux densities * 2, determined on 2019/09/03
    defaultfdens = np.array([891282, 954570, 1173229, 1245433, 1373730, 1506802,
                             1613253, 1702751, 1800721, 1946756, 2096020, 2243951,
                             2367362, 2525968, 2699795, 2861604, 3054829, 3220450,
                             3404182, 3602625, 3794312, 3962926, 4164667, 4360683,
                             4575677, 4767210, 4972824, 5211717, 5444632, 5648266,
                             5926634, 6144249, 6339863, 6598018, 6802707, 7016012,
                             7258929, 7454951, 7742816, 7948976, 8203206, 8411834,
                             8656720, 8908130, 9087766, 9410760, 9571365, 9827078,
                             10023598, 8896671])
    fdens = defaultfdens
    if nbands == 34:
        if Time(slashdate.replace('/', '-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectal windows (bands 4-34 -> spw 0-30)
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[[0] + list(range(2, 34))])
    return dsize, fdens


def ant_trange(vis):
    ''' Figure out nominal times for tracking of old EOVSA antennas, and return time
        range in CASA format
    '''
    from eovsapy import eovsa_array as ea
    from astropy.time import Time

    # Get timerange from the visibility file
    # msinfo = dict.fromkeys(['vis', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs', 'observatory'])
    ms.open(vis)
    # metadata = ms.metadata()
    scans = ms.getscansummary()
    ms.close()
    sk = sorted(list(scans.keys()))
    vistrange = np.array([scans[sk[0]]['0']['BeginTime'], scans[sk[-1]]['0']['EndTime']])

    # Get the Sun transit time, based on the date in the vis file name (must have UDByyyymmdd in the name)
    aa = ea.eovsa_array()
    date = vis.split('UDB')[-1][:8]
    slashdate = date[:4] + '/' + date[4:6] + '/' + date[6:8]
    aa.date = slashdate
    sun = aa.cat['Sun']
    mjd_transit = Time(aa.next_transit(sun).datetime(), format='datetime').mjd
    # Construct timerange limits based on +/- 3h55m from transit time (when all dishes are nominally tracking)
    # and clip the visibility range not to exceed those limits
    mjdrange = np.clip(vistrange, mjd_transit - 0.1632, mjd_transit + 0.1632)
    trange = Time(mjdrange[0], format='mjd').iso[:19] + '~' + Time(mjdrange[1], format='mjd').iso[:19]
    trange = trange.replace('-', '/').replace(' ', '/')
    return trange


def format_spw(spw):
    """
    Formats the spectral window (spw) string for file naming, ensuring start and end values are separated by a dash and zero-padded to two digits.

    :param str spw: The spectral window in the format "start~end", where start and end are integers.
    :return: A formatted string with start and end values zero-padded to two digits, separated by a dash.
    :rtype: str
    """
    return '-'.join(['{:02d}'.format(int(sp_)) for sp_ in spw.split('~')])



class FrequencySetup:
    """
    Manages frequency setup based on observation date for radio astronomy imaging.

    This class is initialized with an observation time and calculates
    essential frequency parameters such as effective observing frequencies (eofreq)
    and spectral windows (spws) based on the observation date. It provides methods
    to calculate reference frequency and bandwidth for given spectral windows.

    :param tim: Observation time used to determine the frequency setup.
    :type tim: astropy.time.Time

    Attributes:
    - tim (astropy.time.Time): The observation time.
    - spw2band (numpy.ndarray): An array mapping spectral window indices to band numbers.
    - bandwidth (float): The bandwidth in GHz.
    - defaultfreq (numpy.ndarray): The default effective observing frequencies in GHz.
    - nbands (int): Number of bands.
    - eofreq (numpy.ndarray): Effective observing frequencies based on the observation time.
    - spws (list of str): Spectral window selections based on the observation time.

    :Example:

    >>> from astropy.time import Time
    >>> tim = Time('2022-01-01T00:00:00', format='isot')
    >>> freq_setup = FrequencySetup(tim)
    >>> crval, cdelt = freq_setup.get_reffreq_and_cdelt('5~10')
    >>> print(crval, cdelt)
    """

    def __init__(self, tim=None):
        if tim is None:
            tim = Time.now()
        self.tim = tim
        self.spw2band = np.array([0, 1] + list(range(4, 52)))
        self.bandwidth = 0.325  # 325 MHz
        self.defaultfreq = 1.1 + self.bandwidth * (self.spw2band + 0.5)

        if self.tim.mjd > 58536:
            self.eofreq = self.defaultfreq
            # self.spws = ['0~1', '2~4', '5~10', '11~20', '21~30', '31~40', '41~49']
            self.spws = ['0~1', '2~4', '5~10', '11~20', '21~30', '31~43', '44~49']
            self.nbands = len(self.spw2band)
        else:
            self.bandwidth = 0.5  # 500 MHz
            self.nbands = 34
            self.eofreq = 1.419 + np.arange(self.nbands) * self.bandwidth
            self.spws = ['1~3', '4~9', '10~16', '17~24', '25~30']

        self.spws_indices = []
        for sp in self.spws:
            s1, s2 = sp.split('~')
            self.spws_indices.append(','.join(map(str, range(int(s1), int(s2) + 1))))

    def get_reffreq_and_cdelt(self, spw, return_bmsize=False):
        """
        Calculates the reference frequency (CRVAL) and the frequency delta (CDELT)
        for a given spectral window range.

        This method takes a spectral window selection and computes the mean of the effective
        observing frequencies (eofreq) within that range as the reference frequency. It also
        calculates the bandwidth covered by the spectral window range as the frequency delta.

        :param spw: Spectral window selection, specified as a range 'start~end' or a single value.
        :type spw: str
        :param return_bmsize: If True, return the beam size based on the reference frequency. Defaults to False.
        :type return_bmsize: bool

        :return: A tuple containing the reference frequency and frequency delta, both in GHz.
        :rtype: (str, str)

        :Example:

        >>> crval, cdelt = freq_setup.get_reffreq_and_cdelt('5~10')
        >>> print(crval, cdelt)
        """
        sp_st, sp_ed = (int(s) for s in spw.split('~')) if '~' in spw else (int(spw), int(spw))
        crval = f'{np.mean(self.eofreq[sp_st:sp_ed + 1]):.4f}GHz'
        nband = (self.eofreq[sp_ed] - self.eofreq[sp_st]) / self.bandwidth + 1
        cdelt = f'{nband * self.bandwidth:.4f}GHz'

        if return_bmsize:
            bmsz = self.get_bmsize(float(crval.rstrip('GHz')), minbmsize=5.0)
            return crval, cdelt, bmsz
        else:
            return crval, cdelt

    def get_bmsize(self, cfreq, refbmsize=60.0, reffreq=1.0, minbmsize=4.0):
        """
        Calculate the beam size at given frequencies based on a reference beam size at a reference frequency.
        This function supports both single frequency values and lists of frequencies.

        :param cfreq: Input frequencies in GHz, can be a float or a list of floats.
        :type cfreq: float or list
        :param refbmsize: Reference beam size in arcsec, defaults to 60.0.
        :type refbmsize: float, optional
        :param reffreq: Reference frequency in GHz, defaults to 1.0.
        :type reffreq: float, optional
        :param minbmsize: Minimum beam size in arcsec, defaults to 4.0.
        :type minbmsize: float, optional
        :return: Beam size at the given frequencies, same type as input (float or numpy array).
        :rtype: float or numpy.ndarray
        """
        # Ensure cfreq is an array for uniform processing
        cfreq = np.array(cfreq, dtype=float)

        # Calculate beam size
        bmsize = refbmsize * reffreq / cfreq

        # Enforce minimum beam size
        bmsize = np.maximum(bmsize, minbmsize)

        # If the original input was a single float, return a single float
        if bmsize.size == 1:
            return bmsize.item()  # Convert numpy scalar to Python float
        else:
            return bmsize

def merge_FITSfiles(fitsfilesin, outfits, snr_weight=None, deselect_index=None,
                    overwrite=True, snr_threshold=None):
    """
    Merges multiple FITS files into a single output file by calculating the mean of stacked data.

    This function stacks the data from a list of input FITS files along a new axis, optionally applying exposure time
    weighting and suppression of on-disk residuals based on a threshold relative to the disk brightness temperature.
    The mean of the stacked data is then written to a specified output FITS file.

    Parameters
    ----------
    fitsfilesin : list of str
        List of file paths to the input FITS files to be merged.
    outfits : str
        File path for the output FITS file where the merged data will be saved.
    snr_weight : list of float, optional
        List of signal-to-noise ratios (SNRs) to use as weights for calculating the mean of the stacked data.
        If both exposure time and SNR weights are provided, the weights will be the product of the two. Defaults to None.
    deselect_index : list of int, optional
        index of the images to be deselected. Defaults to None.
    overwrite : bool, optional
        If True, the output file is overwritten if it already exists. Defaults to True.
    snr_threshold : float, optional
        SNR threshold for selecting images to be deselected. Defaults to 10.

    Raises
    ------
    ValueError
        If any of the input FITS files cannot be read.

    Returns
    -------
    None

    Notes
    -----
    The suppression of on-disk residuals is achieved by applying a sigmoid function to pixels within a specified range
    of the disk brightness temperature, effectively reducing the contribution of residuals while preserving the overall
    structure.

    Examples
    --------
    Merge three FITS files without exposure time weighting and with on-disk residual suppression:

    >>> merge_FITSfiles(['sun_01.fits', 'sun_02.fits', 'sun_03.fits'], 'sun_merged.fits',
    ...                 overwrite=True)
    """
    from astropy.io import fits
    from astropy.time import Time
    from datetime import datetime, timedelta


    exptimes = []
    data = []
    for fidx, file in enumerate(fitsfilesin):
        with fits.open(file) as hdulist:
            for idx, hdu in enumerate(hdulist):
                if 'CDELT1' in hdu.header:
                    header = hdu.header
                    data.append(np.squeeze(hdu.data))
                    exptimes.append(header['EXPTIME'])
                    break
    data_stack = np.array(data)
    if deselect_index is None:
        rsun_pix = header['RSUN_OBS'] / header['CDELT1']
        crpix1, crpix2 = header['CRPIX1'], header['CRPIX2']
        select_indices, rms, snr = detect_noisy_images(data_stack, rsun_pix=rsun_pix, crpix1=crpix1, crpix2=crpix2,
                                                       rsun_ratio=[1.2, 1.3], showplt=False,snr_threshold=snr_threshold)
        deselect_index = ~select_indices
        snr_weight = snr
        print(snr)
    if snr_weight is not None:
        weights = np.array(snr_weight)
    else:
        weights = np.ones(data_stack.shape[-1])
    weights[deselect_index] = 0
    weights = weights / np.sum(weights)
    weights = weights[:, np.newaxis, np.newaxis]
    date_merged = np.nansum(data_stack * weights, axis=0)
    newheader = header.copy()
    exptime = np.nansum(exptimes)
    date_obs = Time(header['date-obs']).datetime.replace(hour=20, minute=0, second=0) - timedelta(
        seconds=exptime / 2.0)
    newheader.update({'EXPTIME': exptime, 'DATE-OBS': date_obs.strftime('%Y-%m-%dT%H:%M:%S')})
    newheader['HISTORY'] = 'Merged from multiple images'
    hdu = fits.PrimaryHDU(date_merged, header=newheader)
    hdu.writeto(outfits, overwrite=overwrite)
    log_print('INFO', f'{np.count_nonzero(deselect_index)} out of {len(fitsfilesin)} images are not selected for merging due to low SNR.')
    return outfits




class MSselfcal:
    """
    A class for imaging using WSClean and selfcal visibility data.

    This class facilitates the creation of images by imaging major
    intervals for WSClean imaging and applying solar differential rotation corrections to minor intervals (new model).
    the model will be written to the visibility data for selfcalibration.

    :param msfile: Path to the Measurement Set (MS) file.
    :type msfile: str
    :param time_intervals: List of major time intervals.
    :type time_intervals: list
    :param time_intervals_major_avg: Average times for major intervals.
    :type time_intervals_major_avg: list
    :param time_intervals_minor_avg: Average times for minor intervals within major intervals.
    :type time_intervals_minor_avg: list
    :param spw: Spectral window index or range.
    :type spw: str
    :param sp_index: String representation of the spectral window indices.
    :type sp_index: str
    :param N1: Number of major intervals.
    :type N1: int
    :param N2: Number of minor intervals within each major interval.
    :type N2: int
    :param workdir: Directory for storing intermediate and output files.
    :type workdir: str
    :param image_marker: Optional marker to add to output image names, defaults to ''
    :type image_marker: str, optional
    :param niter: Maximum number of iterations for WSClean, defaults to 5000.
    :type niter: int, optional
    :param data_column: Column name in the MS file to process, defaults to "CORRECTED_DATA".
    :type data_column: str, optional
    :param beam_size: Beam size for imaging, defaults to None.
    :type beam_size: float, optional
    :param reftime_daily: Reference time for daily rotation corrections, defaults to None.
    :type reftime_daily: `astropy.Time`, optional

    :raises OSError: If any file operation fails.
    :raises Exception: For errors during heliocentric frame registration or rotation corrections.

    :return: Path to the final model name string.
    :rtype: str
    """

    def __init__(self,
                 msfile,
                 time_intervals,
                 time_intervals_major_avg,
                 time_intervals_minor_avg,
                 spw,
                 N1,
                 N2,
                 workdir,
                 image_marker='',
                 niter=5000,
                 gain = 0.1,
                 briggs=0.0,
                 auto_mask=6,
                 auto_threshold=3,
                 maxuv_l=None,
                 minuv_l=None,
                 maxuvw_m=None,
                 minuvw_m=None,
                 data_column="CORRECTED_DATA",
                 pols='XX',
                 beam_size=None,
                 reftime_daily=None):
        self.msfile = msfile
        self.time_intervals = time_intervals
        self.time_intervals_major_avg = time_intervals_major_avg
        self.time_intervals_minor_avg = time_intervals_minor_avg
        self.spw = str(spw)
        if '~' in spw:
            try:
                start_str, end_str = spw.split('~')
                start = int(start_str)
                end = int(end_str)
                self.sp_index = ','.join(str(i) for i in range(start, end + 1))
            except Exception as e:
                raise ValueError(f"Invalid spw range format: '{spw}'. Expected format 'start~end'.") from e
        else:
            self.sp_index = spw
        self.N1 = N1
        self.N2 = N2
        self.workdir = workdir
        self.image_marker = image_marker
        self.niter = niter
        self.gain = gain
        self.pols=pols
        self.data_column = data_column
        self.beam_size = beam_size
        self.reftime_daily = reftime_daily
        self.model_minor_name_str = None
        self.model_ref_name_str = None
        self.briggs = briggs
        self.maxuv_l = maxuv_l
        self.minuv_l = minuv_l
        self.maxuvw_m = maxuvw_m
        self.minuvw_m = minuvw_m
        self.auto_mask = auto_mask
        self.auto_threshold = auto_threshold
        self.succeeded =True

        # self.intervals_out = None

    def run(self, imaging_only=False, predict=True, gen_minor_model=True, clearcache=True):
        """
        Processes the visibility data and generates model images using WSClean.

        This method performs the following steps:
        1. Initial imaging to generate rough model images.
        2. Rotate each model image to the corresponding minor time intervals. This step converts N1 images to N1*N2 images.
           This model will be used for self-calibration and will be subtracted from the visibility data after self-calibration.
        3. Rotate each model image to 20 UT of the day (if the reference time daily is provided; otherwise, skip this step). This model will be added back to the visibility data.
        4. Predict visibilities based on the rotated models.

        :param imaging_only: If True, only perform the initial imaging step and skip the rotation and prediction steps. Default is False.
        :type imaging_only: bool
        :param predict: If True, predict visibilities based on the rotated models. Default is True.
        :type predict: bool
        :param gen_minor_model: If True, generate model images for minor time intervals. Default is True. If False, skip the rotation and prediction steps.
        :type gen_minor_model: bool

        :raises OSError: If any file operation fails.
        :raises Exception: For errors during heliocentric frame registration or rotation corrections.

        :return: Path to the final model name string.
        :rtype: str
        """
        if not gen_minor_model:
            predict = False
        spwstr = format_spw(self.spw)
        msfilename = os.path.basename(self.msfile)
        imname_strlist = ["eovsa", "major", f"{msfilename}", f"sp{spwstr}"]
        if self.image_marker:
            imname_strlist.append(self.image_marker)
        imname = '-'.join(imname_strlist)
        imname_minor = imname.replace('major', 'minor')


        os.system(f'rm -rf {os.path.join(self.workdir, imname)}*.fits')
        clean_obj = ww.WSClean(self.msfile)
        clean_obj.setup(size=1024, scale="2.5asec", weight_briggs=self.briggs, pol=self.pols,
                        niter=self.niter,
                        mgain=0.85,
                        gain=self.gain,
                        data_column=self.data_column,
                        name=os.path.join(self.workdir, imname),
                        multiscale=True,
                        multiscale_scales=[0, 5, 12.5],
                        # multiscale_gain=0.3,
                        maxuv_l=self.maxuv_l, minuv_l=self.minuv_l,
                        maxuvw_m=self.maxuvw_m, minuvw_m=self.minuvw_m,
                        auto_mask=self.auto_mask, auto_threshold=self.auto_threshold,
                        intervals_out=self.N1,
                        no_negative=True, quiet=True,
                        spws=self.sp_index,
                        beam_size=self.beam_size,
                        circular_beam=True)
        clean_obj.run(dryrun=False)

        if clearcache:
            extensions = ['dirty', 'psf', 'residual']
            for ext in extensions:
                for f in glob(os.path.join(self.workdir, f"{imname}-*{ext}.fits")):
                    if os.path.exists(f):
                        os.remove(f)
        if imaging_only:
            return
        # Step 2: Rotate each model image to the corresponding minor time intervals
        if self.N1 == 1:
            fitsname = self.workdir + f"{imname}-t0000-model.fits"
            os.rename(fitsname.replace('-t0000', ''), fitsname)
            modelfitsfiles = [fitsname]
        else:
            modelfitsfiles = sorted(glob(self.workdir + f"{imname}-t*-model.fits"))

        model_dir = os.path.join(self.workdir,
                                 f"modelrot_{msfilename}_{self.image_marker}") if self.image_marker else os.path.join(
            self.workdir, f"modelrot_{msfilename}")

        if not os.path.exists(model_dir): os.makedirs(model_dir)

        self.model_minor_name_str = os.path.join(model_dir, imname_minor)
        os.system(f'rm -rf {self.model_minor_name_str}*.fits')
        if self.reftime_daily is not None:
            self.model_ref_name_str = os.path.join(model_dir, f"{imname}-ref_daily")
        for idx_major, fitsname in enumerate(tqdm(modelfitsfiles)):
            timerange_ = trange2timerange([self.time_intervals[idx_major][0], self.time_intervals[idx_major][-1]])
            heliofitsname = fitsname.replace("model.fits", "model.helio.fits")
            try:
                hf.imreg(vis=self.msfile, imagefile=fitsname, fitsfile=heliofitsname, timerange=timerange_, )
            except Exception as e:
                log_print('ERROR', f'Failed to register {fitsname} to heliocentric frame due to {e}')
                self.succeeded = False
                continue

            if self.reftime_daily is not None:
                log_print('INFO',
                          f"Rotating {fitsname} to daily reftime at {self.reftime_daily.datetime.strftime('%H:%M:%S')} UT...")
                heliorotname_ref_daily = fitsname.replace("model.fits", "model.helio.ref_daily.fits")
                solar_diff_rot_heliofits(heliofitsname, self.reftime_daily,
                                         heliorotname_ref_daily, in_time=self.time_intervals_major_avg[idx_major])
                new_model_file_ref_daily = f"{self.model_ref_name_str}-t{idx_major:04d}-model.fits"
                sunpyfits_to_j2000fits(heliorotname_ref_daily, new_model_file_ref_daily,
                                       template_fits=fitsname, overwrite_prev=True)

            if gen_minor_model:
                log_print('INFO', f"Rotating {fitsname} to {self.N2} minor time intervals ...")
                for idx_minor, time_interval_minor_avg in enumerate(self.time_intervals_minor_avg[idx_major]):
                    heliorotname = heliofitsname.replace('major', 'minor')
                    heliorotname = heliorotname.replace("model.helio.fits",
                                                        f"t{idx_major * self.N2 + idx_minor:04d}-model.helio.rot.fits")
                    heliorotname = os.path.join(model_dir, os.path.basename(heliorotname))
                    solar_diff_rot_heliofits(heliofitsname, time_interval_minor_avg,
                                             heliorotname, in_time=self.time_intervals_major_avg[idx_major])

                    new_model_file = f"{self.model_minor_name_str}-t{idx_major * self.N2 + idx_minor:04d}-model.fits"
                    sunpyfits_to_j2000fits(heliorotname, new_model_file,
                                           template_fits=fitsname, overwrite_prev=True)

        if predict:
            # Step 3: Predict visibilities for each model image
            cmd = f"wsclean -predict -reorder -spws {self.sp_index} -name {self.model_minor_name_str} -quiet -intervals-out {self.N1 * self.N2} {self.msfile}"
            log_print('INFO', f"Running wsclean predict: {cmd}")
            os.system(cmd)
            # os.system(f'rm -rf {self.workdir}/{imname}*model.helio.fits')
            # os.system(f'rm -rf {model_dir}/{imname}*model.helio.rot.fits')
        return

def imaging(vis, workdir=None, imgoutdir=None, pols='XX', data_column='DATA'):

    if isinstance(vis, list):
        msfile = vis[0]
    else:
        msfile=vis

    viz_timerange = ant_trange(msfile)
    (tstart, tend) = viz_timerange.split('~')
    tbg_msfile = Time(qa.quantity(tstart, 's')['value'], format='mjd').to_datetime()
    ted_msfile = Time(qa.quantity(tend, 's')['value'], format='mjd').to_datetime()
    reftime_daily = Time(datetime.combine(tbg_msfile.date(), time(20, 0)))
    date_str = tbg_msfile.strftime('%Y%m%d')
    freq_setup = FrequencySetup(Time(tbg_msfile))
    spws_indices = freq_setup.spws_indices
    spws = freq_setup.spws
    spws_imaging = spws
    defaultfreq = freq_setup.defaultfreq
    freq = defaultfreq
    nbands = freq_setup.nbands
    slashdate = tbg_msfile.strftime('%Y/%m/%d')
    dsize, fdens = calc_diskmodel(slashdate, nbands, freq, defaultfreq)
    fdens = fdens / 2.0  ## convert from I to XX

    briggs = 0.0
    gain = 0.3
    wsclean_intervals = WSCleanTimeIntervals(msfile, combine_scans=True)
    ## this step use a single model from the entire MS file and rotate it to each of the minor intervals.
    wsclean_intervals.compute_intervals(nintervals_minor=3)
    wscln_tim_info_combscan = wsclean_intervals.results()

    N1 = wscln_tim_info_combscan['nintervals_major']
    time_intervals = wscln_tim_info_combscan['time_intervals']
    time_intervals_major_avg = wscln_tim_info_combscan['time_intervals_major_avg']

    timeranges = []

    for tidx, trange in enumerate(time_intervals):
        timeranges.append(
            trange[0].datetime.strftime('%Y/%m/%d/%H:%M:%S') + '~' + trange[-1].datetime.strftime('%Y/%m/%d/%H:%M:%S'))

    for sidx, sp_index in enumerate(spws_indices):
        spwstr = format_spw(spws[sidx])
        # reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        # imname_strlist = ["eovsa", "major", f"{msfilename}", f"sp{spwstr}", 'final','disk']
        # briggs = 0.5 if sidx in [0] else 0.0
        # flagdata(vis=msfile, mode="tfcrop", spw=spws[sidx], action='apply', display='',
        #          timecutoff=3.0, freqcutoff=3.0, maxnpieces=2, flagbackup=False)
        if isinstance(vis, list):
            msfile = vis[sidx]
            sp_index = []
        msname, _ = os.path.splitext(msfile)
        msfilename = os.path.basename(msfile)
        imname_strlist = ["eovsa", "major", f"{msfilename}", f"sp{spwstr}", 'final']
        imname = '-'.join(imname_strlist)
        clean_obj = ww.WSClean(msfile)
        clean_obj.setup(size=1024, scale="2.5asec", pol=pols,
                        weight_briggs=briggs,
                        niter=2000,
                        mgain=0.85,
                        gain=gain,
                        # gain=0.2,
                        data_column=data_column.upper(),
                        name=os.path.join(workdir, imname),
                        multiscale=True,
                        multiscale_gain=0.3,
                        # multiscale_gain=0.5,
                        multiscale_scales=[0, 5, 12.5],
                        auto_mask=2, auto_threshold=1,
                        # auto_mask=1.5, auto_threshold=1,
                        # minuv_l=500,
                        minuv_l=None,
                        minuvw_m=8,
                        maxuvw_m=1500,
                        intervals_out=N1,
                        no_negative=True, quiet=True,
                        circular_beam=True,
                        spws=sp_index)
        clean_obj.run(dryrun=False)

        fitsname = sorted(glob(os.path.join(workdir, imname + '*image.fits')))
        fitsname_helio = [f.replace('image.fits', 'image.helio.fits') for f in fitsname]
        fitsname_helio_ref_daily = [f.replace('image.fits', 'image.helio.ref_daily.fits') for f in fitsname]
        hf.imreg(vis=msfile, imagefile=fitsname, fitsfile=fitsname_helio, timerange=timeranges, toTb=True)
        for eoidx, (fits_helio, fits_helio_ref_daily) in enumerate(zip(fitsname_helio, fitsname_helio_ref_daily)):
            solar_diff_rot_heliofits(fits_helio, reftime_daily, fits_helio_ref_daily,
                                     in_time=time_intervals_major_avg[eoidx])
        fitsfilefinal = fitsname_helio_ref_daily

        log_print('INFO', f"Final imaging for SPW {spws[sidx]}: completed")

        reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        sp_st, sp_ed = spws[sidx].split('~')
        dszs = [float(dsize[int(sp)].rstrip('arcsec')) for sp in range(int(sp_st), int(sp_ed) + 1)]
        fdns = [fdens[int(sp)] for sp in range(int(sp_st), int(sp_ed) + 1)]
        # fdns = [fdens[int(sp)]*50 for sp in range(int(sp_st), int(sp_ed)+1)]
        synfitsfiles = []
        for eoidx, eofile in enumerate(fitsname):
            datetimestr = time_intervals_major_avg[eoidx].datetime.strftime('%Y%m%dT%H%M%SZ')
            synfitsfile = os.path.join(imgoutdir,
                                       f"eovsa.synoptic.{datetimestr}.s{spwstr}.tb.disk.fits")
            synfitsfiles.append(synfitsfile)
        add_convolved_disk_to_fits(fitsfilefinal, synfitsfiles, dszs, fdns, ignore_data=False, toTb=True,
                                   rfreq=reffreq, bmaj=bmsize / 3600.)
    outfits_all = []
    for spw in spws:
        spwstr = format_spw(spw)
        outfits = os.path.join(imgoutdir, f'eovsa.synoptic_daily.{date_str}T200000Z.s{spwstr}.tb.disk.fits')
        synfitsfiles = sorted(glob(os.path.join(imgoutdir,
                                                f"eovsa.synoptic.{date_str[:-1]}?T??????Z.s{spwstr}.tb.disk.fits")))
        if len(synfitsfiles) > 0:
            log_print('INFO', f"Merging synoptic images for SPW {spwstr} to {outfits} ...")
            merge_FITSfiles(synfitsfiles, outfits, overwrite=True, snr_threshold=10)
            outfits_all.append(outfits)
        else:
            outfits_all.append(None)
    return outfits_all


def pipeline_run(vis, outputvis='', workdir=None, slfcaltbdir=None, imgoutdir=None,
                 pols='XX', verbose=True, do_diskslfcal=True, hanning=False, do_sbdcal=False, overwrite=False, skip_slfcal=False, segmented_imaging=True):
    """
    Executes the EOVSA data processing pipeline for solar observation data.

    :param vis: Path to the input measurement set (MS) data.
    :type vis: str
    :param outputvis: Output path for the processed MS, defaults to an empty string.
    :type outputvis: str, optional
    :param workdir: Working directory for intermediate data, defaults to None which sets it to '/data1/workdir'.
    :type workdir: str, optional
    :param slfcaltbdir: Directory for storing calibration tables, defaults to None.
    :type slfcaltbdir: str, optional
    :param imgoutdir: Output directory for image files, defaults to None.
    :type imgoutdir: str, optional
    :param figoutdir: Output directory for figures, defaults to None.
    :type figoutdir: str, optional
    :param clearcache: If True, clears all cache files after processing, defaults to False.
    :type clearcache: bool, optional
    :param clearlargecache: If True, clears the large cache after processing, defaults to False. If True, clearcache is ignored.
    :type clearlargecache: bool, optional
    :param pols: Polarization types to process, defaults to 'XX'.
    :type pols: str, optional
    :param mergeFITSonly: If True, skips processing and only merges FITS files, defaults to False.
    :type mergeFITSonly: bool, optional
    :param verbose: Enables verbose output during processing, defaults to True.
    :type verbose: bool, optional
    :param do_diskslfcal: If True, performs disk self-calibration, defaults to True.
    :type do_diskslfcal: bool, optional
    :param overwrite: If True, overwrites existing files, defaults to False.
    :type overwrite: bool, optional
    :param niter_init: Initial number of iterations for imaging, defaults to 200.
    :type niter_init: int, optional
    :param ncpu: Specifies the number of CPUs for parallel processing, defaults to 'auto'.
    :type ncpu: str or int, optional
    :param tr_series_imaging: Time ranges for imaging, defaults to None.
    :type tr_series_imaging: list of tuple, optional
    :param spws_imaging: Spectral windows selected for imaging, defaults to None.
    :type spws_imaging: list of str, optional
    :param hanning: If True, applies Hanning smoothing to the data, defaults to False.
    :type hanning: bool, optional
    :param do_sbdcal: Boolean flag to perform single-band delay calibration, defaults to False.
    :type do_sbdcal: bool
    :param skip_slfcal: If True, skips self-calibration and directly images the data, defaults to False.
    :param segmented_imaging: Determines whether to segment the all-day data into chunks for imaging,
                          defaults to True.
    :type segmented_imaging: bool (optional)
    :return: Path to the processed visibility data.
    :rtype: str

    :example:
    ## if you want to specify the spectral windows for imaging
    >>> spws_imaging = ['5~10', '11~20', '21~30']
    ## if you want to specify the time range for imaging
    >>> from datetime import datetime, timedelta
    >>> from suncasa.eovsa import eovsa_synoptic_imaging_pipeline as esip
    >>> tbg_imaging = datetime(2024, 4, 8, 16, 58, 0)
    >>> ted_imaging = datetime(2024, 4, 8, 19, 00, 0)
    >>> tdt_imaging = timedelta(minutes=2)
    >>> tr_series_imaging = esip.generate_trange_series(tbg_imaging, ted_imaging, tdt_imaging)
    """
    from eovsapy.util import Time
    if os.path.exists(outputvis) and not overwrite:
        log_print('INFO', f"Output MS file {outputvis} already exists. Skipping processing.")
        return outputvis


    if not workdir.endswith('/'):
        workdir += '/'
    outfits_all = []
    debug_mode = False
    if debug_mode:
        workdir = './'
        slfcaltbdir = None
        imgoutdir = './'
        pols = 'XX'
        verbose = True
        do_diskslfcal = True
        overwrite = False
        hanning = False
        do_sbdcal = False
        segmented_imaging = True
        udbdir = '/data1/eovsa/fits/UDB/'
        from suncasa.suncasatasks import calibeovsa
        from suncasa.suncasatasks import importeovsa
        from eovsapy.dump_tsys import findfiles

        # datein = datetime(2024, 12, 15, 20, 0, 0)
        datein = datetime(2025, 2, 6, 20, 0, 0)
        pols = 'XX,YY'
        # datein = datetime(2021, 11, 25, 20, 0, 0)
        trange = Time(datein)
        if trange.mjd == np.fix(trange.mjd):
            # if only date is given, move the time from 00 to 12 UT
            trange = Time(trange.mjd + 0.5, format='mjd')

        tdatetime = trange.to_datetime()
        dhr = trange.LocalTime.utcoffset().total_seconds() / 60 / 60 / 24
        btime = Time(np.fix(trange.mjd + dhr) - dhr, format='mjd')
        etime = Time(btime.mjd + 1, format='mjd')
        trange = Time([btime, etime])
        print('Selected timerange in UTC: ', trange.iso)
        inpath = '{}{}/'.format(udbdir, tdatetime.strftime("%Y"))
        outpath = './'
        msfile_synoptic = os.path.join(outpath, 'UDB' + tdatetime.strftime("%Y%m%d") + '.ms')
        # sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
        sclist = findfiles(trange, projid='NormalObserving', srcid='Sun')
        udbfilelist = sclist['scanlist']
        udbfilelist = [os.path.basename(ll) for ll in udbfilelist]
        udbfilelist_set = set(udbfilelist)
        msfiles = []
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))
        ncpu = 1
        importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1,
                    visprefix=outpath, nocreatms=False,
                    doconcat=False, modelms="", doscaling=False, keep_nsclms=False, udb_corr=True)

        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob('{}UDB*.ms*'.format(outpath)) if ll.endswith('.ms') or ll.endswith('.ms.tar.gz')]
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))

        invis = [outpath + ll + '.ms' for ll in sorted(list(msfiles))]
        vis_out = os.path.join(os.path.dirname(invis[0]), os.path.basename(invis[0])[:11] + '.ms')
        vis = calibeovsa(invis, caltype=['refpha','phacal'], caltbdir='./', interp='nearest',
                         doflag=True,
                         flagant='13~15',
                         doimage=False, doconcat=True,
                         concatvis=vis_out, keep_orig_ms=True)

        vis=vis_out


    # workdir = '/data1/sjyu/eovsa/20241215'
    # if not workdir.endswith('/'):
    #     workdir += '/'
    # vis = 'UDB20241215.ms'
    #
    # workdir = '/data1/sjyu/eovsa/20211124'
    # if not workdir.endswith('/'):
    #     workdir += '/'
    # vis = 'UDB20211124.ms'

    # workdir = '/data1/sjyu/eovsa/20250101'
    # if not workdir.endswith('/'):
    #     workdir += '/'
    # vis = 'UDB20250101.ms'
    #
    # workdir = '/data1/sjyu/eovsa/20250214'
    # if not workdir.endswith('/'):
    #     workdir += '/'
    # vis = 'UDB20250214.ms'


    msfile = vis.rstrip('/')
    if workdir is None:
        workdir = './'
    os.chdir(workdir)
    if slfcaltbdir is None:
        slfcaltbdir = workdir + '/'
    if imgoutdir is None:
        imgoutdir = workdir + '/'

    msname, _ = os.path.splitext(msfile)
    msname = os.path.basename(msname)

    msfile_copy = os.path.join(workdir, f'{msname}.ms')
    if not os.path.exists(msfile_copy):
        if msfile.lower().endswith('.ms'):
            shutil.copytree(msfile, msfile_copy)
        elif msfile.lower().endswith('.ms.tar.gz'):
            os.system(f'tar -zxf {msfile} -C {workdir}')
        else:
            raise ValueError(f"Unsupported file format: {msfile}")
    msfile = msfile_copy

    ## the msfile we use 60 min model image to correct the data in 10 min interval. the model image is shifted to the reftime_master (20:00 UT of each day).
    # tdt_imaging = timedelta(hours=1)
    # tdt_sub = timedelta(minutes=10)

    viz_timerange = ant_trange(msfile)
    (tstart, tend) = viz_timerange.split('~')
    tbg_msfile = Time(qa.quantity(tstart, 's')['value'], format='mjd').to_datetime()
    ted_msfile = Time(qa.quantity(tend, 's')['value'], format='mjd').to_datetime()
    reftime_daily = Time(datetime.combine(tbg_msfile.date(), time(20, 0)))
    date_str = tbg_msfile.strftime('%Y%m%d')
    freq_setup = FrequencySetup(Time(tbg_msfile))
    spws_indices = freq_setup.spws_indices
    spws = freq_setup.spws
    spws_imaging = spws
    defaultfreq = freq_setup.defaultfreq
    freq = defaultfreq
    nbands = freq_setup.nbands
    # bandinfo = mstl.get_bandinfo(msfile, returnbdinfo=True)

    slashdate = tbg_msfile.strftime('%Y/%m/%d')
    dsize, fdens = calc_diskmodel(slashdate, nbands, freq, defaultfreq)
    fdens = fdens / 2.0  ## convert from I to XX

    spwidx2proc = [0,1,2,3,4,5,6]
    # spwidx2proc = [0,1,2]
    # spwidx2proc = [1, 2]
    # spwidx2proc = [1]
    # spwidx2proc = [2]
    # spwidx2proc = [3]
    ## disk slfcal after feature cal works better for spwidx [0,1,2, 3, !4, !5]
    # spwidx2proc = [4]
    # spwidx2proc = [5]

    ## set a maxmium uv range for disk self-calibration. 30 meters is a good choice for EOVSA. A minumum of 200 lambda.
    uvmax_l_list = []
    ## set a minimum uv range for feature self-calibration. 225 meters is a good choice for EOVSA. A minumum of 1500 lambda.
    uvmin_l_list = []

    ## a flag indicates if the disk self-calibration shall be performed before the feature self-calibration.
    diskslfcal_first = []

    for sidx, sp_index in enumerate(spws_indices):
        # diskslfcal_first.append(False)
        if sidx in [0, 1, 2, 3]:
            diskslfcal_first.append(False)
        else:
            diskslfcal_first.append(True)
        if sidx not in spwidx2proc:
            uvmax_l_list.append(0)
            uvmin_l_list.append(0)
        else:
            # uvmax_d = 20
            # uvmin_d = 130
            # if sidx in [0]:
            #     uvmin_d=65
            uvmax_d = 15
            uvmin_d = 110
            if sidx in [0]:
                uvmin_d=55
            sp_st, sp_ed = (int(s) for s in spws[sidx].split('~'))
            uvmax_l_list.append(uvmax_d / ((3e8 / (np.nanmean(freq[sp_st:sp_st + 1]) * 1e9))))
            uvmin_l_list.append(uvmin_d / ((3e8 / (np.nanmean(freq[sp_st:sp_st + 1]) * 1e9))))
    uvmax_l_list = np.array(uvmax_l_list).astype(float)
    uvmin_l_list = np.array(uvmin_l_list).astype(float)
    # uvmax_l_list[uvmax_l_list<200] = 200
    uvmax_l_list[uvmax_l_list<100] = 100
    uvmax_l_list[uvmax_l_list>1200] = 1200
    uvmin_l_list[uvmin_l_list>1800] = 1800
    uvmax_l_str = [f'<{l:.0f}lambda' for l in uvmax_l_list]
    uvmin_l_list = np.array(uvmin_l_list).astype(float)
    uvmin_l_str = [f'>{l:.0f}lambda' for l in uvmin_l_list]


    run_start_time = datetime.now()
    ### --------------------------------------------------------------###
    ## pre-processing block. flagging, hanning smoothing...
    ### --------------------------------------------------------------###
    run_start_time_pre_proc = datetime.now()
    # if not mergeFITSonly:
    if os.path.isdir(msfile + '.flagversions') == True:
        if verbose:
            log_print('INFO', f'Flagversion of {msfile} already exists. Skipped...')
    else:
        if verbose:
            log_print('INFO', f'Flagging any high amplitudes viz from flares or RFIs in {msfile}')
        flagmanager(msfile, mode='save', versionname='pipeline_init')
        flagdata(vis=msfile, mode="tfcrop", spw='', action='apply', display='',
                 timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
        flagmanager(msfile, mode='save', versionname='pipeline_remove_RFI-and-BURSTS')

    if hanning:
        msfile_hanning = msfile + '.hanning'
        if not os.path.exists(msfile_hanning):
            if verbose:
                log_print('INFO', f'Perform hanningsmooth for {msfile}...')
            hanningsmooth(vis=msfile, datacolumn='data', outputvis=msfile_hanning)
        else:
            if verbose:
                log_print('INFO', f'The hanningsmoothed {msfile} already exists. Skipped...')
        msfile = msfile_hanning

    delmod(msfile)

    wsclean_intervals = WSCleanTimeIntervals(msfile, combine_scans=True)
    ## this step use a single model from the entire MS file and rotate it to each of the minor intervals.
    wsclean_intervals.compute_intervals(nintervals_minor=3)
    wscln_tim_info_combscan = wsclean_intervals.results()

    N1 = wscln_tim_info_combscan['nintervals_major']
    N2 = wscln_tim_info_combscan['nintervals_minor']
    tdur = wscln_tim_info_combscan['tdur']
    time_intervals = wscln_tim_info_combscan['time_intervals']
    time_intervals_major_avg = wscln_tim_info_combscan['time_intervals_major_avg']
    time_intervals_minor_avg = wscln_tim_info_combscan['time_intervals_minor_avg']

    ## this is for spw 0~1 due to its low SNR.
    ## this step use a single model from the entire MS file and rotate it to N1 minor intervals.
    wsclean_intervals_comb = WSCleanTimeIntervals(msfile, combine_scans=True)
    wsclean_intervals_comb.compute_intervals(nintervals_minor=3)
    wscln_tim_info_combscan_comb = wsclean_intervals_comb.results()
    N1_comb = 1
    ## change N2 to the total number of minor intervals
    N2_comb = wscln_tim_info_combscan_comb['nintervals_major'] * wscln_tim_info_combscan_comb['nintervals_minor']
    time_intervals_comb = [
        Time([wscln_tim_info_combscan_comb['time_intervals'][0][0], wscln_tim_info_combscan_comb['time_intervals'][-1][-1]])]
    ## combine all the major intervals to a single one
    time_intervals_major_avg_comb = Time(
        [np.nanmean(np.hstack([l.mjd for l in wscln_tim_info_combscan_comb['time_intervals']]))], format='mjd')
    ## combine all the major's minor intervals to a single one.
    time_intervals_minor_avg_comb = [
        Time(np.hstack([l.mjd for l in wscln_tim_info_combscan_comb['time_intervals_minor_avg']]), format='mjd')]
    # N1_comb = 1
    # N2_comb = N1
    # time_intervals_comb = [
    #     Time([wscln_tim_info_combscan_comb['time_intervals'][0][0], wscln_tim_info_combscan_comb['time_intervals'][-1][-1]])]
    # time_intervals_major_avg_comb = Time(
    #     [np.nanmean(np.hstack([l.mjd for l in wscln_tim_info_combscan['time_intervals']]))], format='mjd')
    # time_intervals_minor_avg_comb = [Time([l.mjd],format='mjd') for l in time_intervals_major_avg]

    timeranges = []
    for tidx, trange in enumerate(time_intervals):
        timeranges.append(
            trange[0].datetime.strftime('%Y/%m/%d/%H:%M:%S') + '~' + trange[-1].datetime.strftime('%Y/%m/%d/%H:%M:%S'))

    msname, _ = os.path.splitext(msfile)
    msname = os.path.basename(msname)
    msfilename = os.path.basename(msfile)

    diskxmlfile = msfile + '.SOLDISK.xml'
    run_start_time_shift_corr = datetime.now()

    caltbs_all = []
    slfcal_init_objs = []
    slfcal_rnd1_objs = []
    slfcal_rnd2_objs = []
    imaging_objs = []
    ## alldaymode_spidx is the spw index for using the all-day viz data for imaging in the initial round of self-calibration.
    ## for observation after 2023-01-15, we use spw 0, 5, 6 for all-day mode imaging. The selection of the spw is mostly empirical as it yields the best results.
    # if reftime_daily>Time(datetime(2023, 1, 15, 0, 0, 0)):
    #     alldaymode_spidx = [0,5,6]
    #     # alldaymode_spidx = [0,1,2, 5, 6]
    # else:
    #     alldaymode_spidx = [0]
    alldaymode_spidx = [0]
    for sidx, sp_index in enumerate(spws_indices):
        spwstr = format_spw(spws[sidx])
        msfile_sp = f'{msname}.sp{spwstr}.slfcaled.ms'
        caltbs = []
        ncaltbs = len(caltbs)

        if os.path.exists(msfile_sp):
            if overwrite:
                os.system('rm -rf ' + msfile_sp)
            else:
                log_print('INFO', f"MS file {msfile_sp} already exists. Skipping processing.")
                slfcal_init_objs.append(None)
                slfcal_rnd1_objs.append(None)
                slfcal_rnd2_objs.append(None)
                caltbs_all.append(caltbs)
                continue
        if sidx not in spwidx2proc:
            slfcal_init_objs.append(None)
            slfcal_rnd1_objs.append(None)
            slfcal_rnd2_objs.append(None)
            caltbs_all.append(caltbs)
            continue
        clearcal(msfile, spw=spws[sidx])
        reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        log_print('INFO', f"Processing SPW {spws[sidx]} for msfile {os.path.basename(msfile)} ...")
        if sidx in [0]:
            auto_mask = 2
            auto_threshold = 1.5
        # elif sidx in [5,6]:
        #     auto_mask = 5
        #     auto_threshold = 2.5
        else:
            auto_mask = 6
            auto_threshold = 3
        if sidx in alldaymode_spidx:
            slfcal_obj = MSselfcal(msfile, time_intervals_comb, time_intervals_major_avg_comb,
                                   time_intervals_minor_avg_comb,
                                   spws[sidx], N1_comb, N2_comb, workdir, image_marker='init', niter=100,
                                   data_column="DATA",
                                   pols=pols,
                                   auto_mask=auto_mask, auto_threshold=auto_threshold,)
        else:
            slfcal_obj = MSselfcal(msfile, time_intervals, time_intervals_major_avg, time_intervals_minor_avg,
                                   spws[sidx], N1, N2, workdir, image_marker='init', niter=100,
                                   data_column="DATA",
                                   pols=pols,
                                   auto_mask=auto_mask, auto_threshold=auto_threshold,)
        # data_column="DATA", beam_size=bmsize)
        slfcal_obj.run()

        ## check if all data in the fits files are flagged. If so, skip the self-calibration for this spw.
        # slfcal_obj.run(imaging_only=True)
        if slfcal_obj.succeeded:
            slfcal_init_objs.append(slfcal_obj)
        else:
            slfcal_init_objs.append(None)
            slfcal_rnd1_objs.append(None)
            slfcal_rnd2_objs.append(None)
            caltbs_all.append(caltbs)
            continue

        caltb = os.path.join(workdir, f"caltb_init_inf_sp{spwstr}.pha")
        if os.path.exists(caltb): os.system('rm -rf ' + caltb)
        ## select a subset of the data for the initial phase self-calibration.
        if sidx in [0]:
            # timerange_sub = trange2timerange([time_intervals[min(1, N1 // 2 - 1)][0], time_intervals[N1 // 2 + 1][-1]])
            timerange_sub = viz_timerange
        else:
            timerange_sub = trange2timerange([time_intervals[0][0], time_intervals[N1 // 2 + 1][-1]])
        gaincal(vis=msfile, caltable=caltb, selectdata=True,
                timerange=timerange_sub,
                uvrange=uvmin_l_str[sidx],
                spw=spws[sidx],
                combine="scan", antenna='0~12&0~12', refant='0',
                solint='inf', refantmode="flex",
                gaintype='G',
                gaintable=caltbs,
                minsnr=1.0, calmode='p', append=False)
        if os.path.exists(caltb): caltbs.append(caltb)

        if sidx in alldaymode_spidx:
            caltb = os.path.join(workdir, f"caltb_init_int_sp{spwstr}.pha")
            if os.path.exists(caltb): os.system('rm -rf ' + caltb)
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    timerange=timerange_sub,
                    uvrange=uvmin_l_str[sidx],
                    spw=spws[sidx],
                    combine="scan", antenna='0~12&0~12', refant='0',
                    solint=f'{tdur / N1:.0f}s', refantmode="flex",
                    gaintype='G',
                    gaintable=caltbs,
                    minsnr=1.0, calmode='p', append=False)
            if os.path.exists(caltb): caltbs.append(caltb)


        if do_diskslfcal and diskslfcal_first[sidx]:
            run_start_time_disk_slfcal = datetime.now()
            delmod(vis=msfile)
            imname_strlist = ["eovsa", "disk", f"{msfilename}", f"sp{spwstr}"]
            imname = '-'.join(imname_strlist)
            clean_obj = ww.WSClean(msfile)
            clean_obj.setup(size=1024, scale="2.5asec",
                            pol=pols,
                            weight_briggs=0.0,
                            niter=1,
                            mgain=0.85,
                            interval=[0, 10],
                            data_column='DATA',
                            name=os.path.join(workdir, imname),
                            multiscale=True,
                            multiscale_gain=0.3,
                            # multiscale_gain=0.5,
                            multiscale_scales=[0, 5, 12.5],
                            auto_mask=6, auto_threshold=3,
                            intervals_out=1,
                            no_negative=True, quiet=True,
                            circular_beam=True,
                            spws=sp_index)
            clean_obj.run(dryrun=False)

            in_fits = os.path.join(workdir, f"{imname}-model.fits")
            for sp in sp_index.split(','):
                # run_start_time_disk_slfcal_sp = datetime.now()
                reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(f'{sp}~{sp}', return_bmsize=True)
                sp = int(sp)
                model_imname = imname.replace(format_spw(spws[sidx]), f'{sp:02d}') + '_adddisk'
                out_fits = model_imname + '-model.fits'
                dsz = float(dsize[sp].rstrip('arcsec'))
                fdn = fdens[sp]
                add_convolved_disk_to_fits(in_fits, out_fits, dsz, fdn, ignore_data=True, bmaj=bmsize/3600.)
                cmd = f"wsclean -predict -reorder -spws {sp}" + f" \
                    -name {model_imname} -quiet -intervals-out 1 {msfile}"
                log_print('INFO', f"Running wsclean predict: {cmd}")
                os.system(cmd)

            caltb = os.path.join(workdir, f"caltb_disk_inf_sp{spwstr}.pha")
            if os.path.exists(caltb): os.system('rm -rf ' + caltb)
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    uvrange=uvmax_l_str[sidx],
                    spw=spws[sidx],
                    combine="scan",
                    antenna='0~12&0~12', refant='0', solint='inf', refantmode="flex",
                    gaintype='G',
                    gaintable=caltbs,
                    minsnr=1,
                    minblperant=2,
                    calmode='p', append=False)
            if os.path.exists(caltb): caltbs.append(caltb)

            caltb = os.path.join(workdir,  f"caltb_disk_inf_sp{spwstr}.amp")
            if os.path.exists(caltb): os.system('rm -rf ' + caltb)
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    uvrange=uvmax_l_str[sidx],
                    spw=spws[sidx],
                    combine="scan",
                    antenna='0~12&0~12', refant='0', solint='inf', refantmode="flex",
                    gaintype='G',
                    gaintable=caltbs,
                    minsnr=1,
                    minblperant=2,
                    calmode='a', append=False)
            if os.path.exists(caltb): caltbs.append(caltb)

            run_end_time_disk_slfcal = datetime.now()
            elapsed_time = run_end_time_disk_slfcal - run_start_time_disk_slfcal
            elapsed_time_disk_slfcal = elapsed_time.total_seconds() / 60
            log_print('INFO', f"Disk slfcal processe for spw {spwstr} completed. Elapsed time: {elapsed_time_disk_slfcal:.1f} minutes")

        if len(caltbs) > ncaltbs:
            log_print('INFO', f'Applying initial feature selfcal solution (spw {spwstr})  to the data')
            applycal(vis=msfile, selectdata=True, antenna="0~12",
                     spw=spws[sidx],
                     gaintable=caltbs, interp='nearest',
                     calwt=False, applymode="calonly")
            ncaltbs = len(caltbs)
        else:
            log_print('WARNING', 'Round 1 Feature Slfcal: No new solutions. Skipping applycal.')

        if do_diskslfcal and diskslfcal_first[sidx]:
            # # step 5: subtract the disk model from the viz
            run_start_time_disk_slfcal_subtract = datetime.now()
            # for sidx, sp_index in enumerate(spws_indices):
            #     if not diskslfcal_first[sidx]: continue
            #     if sidx not in spwidx2proc: continue
            #     spwstr = format_spw(spws[sidx])
            #     msfilename = os.path.basename(msfile)
            #     imname_strlist = ["eovsa", "disk", f"{msfilename}", f"sp{spwstr}"]
            #     imname = '-'.join(imname_strlist)
            #     for sp in sp_index.split(','):
            #         sp = int(sp)
            #         model_imname = imname.replace(format_spw(spws[sidx]), f'{sp:02d}') + '_adddisk'
            #         cmd = f"wsclean -predict -reorder -spws {sp}" + f" \
            #             -name {model_imname} -quiet -intervals-out 1 {msfile}"
            #         log_print('INFO', f"Running wsclean predict: {cmd}")
            #         os.system(cmd)
            uvsub(vis=msfile)
            run_end_time_disk_slfcal_subtract = datetime.now()
            elapsed_time = run_end_time_disk_slfcal_subtract - run_start_time_disk_slfcal_subtract
            elapsed_time_disk_slfcal_subtract = elapsed_time.total_seconds() / 60
            log_print('INFO', f"Disk subtraction for spw {spwstr} completed. Elapsed time: {elapsed_time_disk_slfcal_subtract:.1f} minutes")

        run_end_time_pre_proc = datetime.now()
        elapsed_time = run_end_time_pre_proc - run_start_time_pre_proc
        elapsed_time_pre_proc = elapsed_time.total_seconds() / 60
        log_print('INFO', f"Pre-processing for SPW {spwstr}: completed in {elapsed_time_pre_proc:.1f} minutes")

        ### --------------------------------------------------------------###
        ## shift correction block.
        ### --------------------------------------------------------------###


        run_start_time_scan_shift_corr = datetime.now()

        # reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        log_print('INFO', f"Processing SPW {spws[sidx]} for msfile {os.path.basename(msfile)} ...")
        # briggs = 0.5 if sidx in [0] else 0.0
        briggs = 0.0
        if sidx in [0]:
            auto_mask = 2
            auto_threshold = 1.5
        # elif sidx in [5,6]:
        #     auto_mask = 4
        #     auto_threshold = 1.5
        else:
            auto_mask = 4
            auto_threshold = 2
        slfcal_obj = MSselfcal(msfile, time_intervals, time_intervals_major_avg, time_intervals_minor_avg,
                               spws[sidx], N1, N2, workdir, image_marker='round1', niter=500,
                               briggs=briggs,
                               minuvw_m=8,
                               maxuvw_m=1500,
                               auto_mask=auto_mask, auto_threshold=auto_threshold,
                               pols=pols,
                               data_column="CORRECTED_DATA")
        slfcal_obj.run()
        slfcal_rnd1_objs.append(slfcal_obj)


        caltb = os.path.join(workdir, f"caltb_rnd1_int_sp{spwstr}.pha")
        if os.path.exists(caltb): os.system('rm -rf ' + caltb)
        gaincal(vis=msfile, caltable=caltb, selectdata=True,
                uvrange=uvmin_l_str[sidx],
                spw=spws[sidx],
                combine="scan", antenna='0~12&0~12', refant='0',
                # solint=f'{tdur / N2 / N1:.0f}s', refantmode="flex",
                solint=f'{tdur / N1:.0f}s', refantmode="flex",
                gaintype='G',
                gaintable=caltbs,
                minsnr=1.0, calmode='p', append=False)
        if os.path.exists(caltb): caltbs.append(caltb)

        if do_sbdcal:
            caltb_k = os.path.join(workdir, f"caltb_inf_sp{spwstr}.sbd")
            if os.path.isdir(caltb_k): os.system('rm -rf ' + caltb_k)
            gaincal(vis=msfile, caltable=caltb_k, solint='inf', combine='scan',
                    # uvrange='>1.5klambda',
                    uvrange=uvmin_l_str[sidx],
                    refant='0', gaintype='K',
                    spw=spws[sidx],
                    gaintable=caltbs,
                    calmode='p', refantmode='flex', minsnr=2.0, minblperant=4,
                    append=False)
            if os.path.isdir(caltb_k):
                caltbs.append(caltb_k)

        if len(caltbs) > ncaltbs:
            applycal(vis=msfile, selectdata=True, antenna="0~12",
                     spw=spws[sidx],
                     gaintable=caltbs, interp='linear',
                     calwt=False, applymode="calonly")
            ncaltbs = len(caltbs)
            log_print('INFO', f"Round 1 Feature Selfcal for SPW {spws[sidx]}: solution applied")
        else:
            log_print('WARNING', f"Round 1 Feature Selfcal for SPW {spws[sidx]}: no new solutions, skipping applycal")

        log_print('INFO', f"Processing SPW {spws[sidx]}  ...")
        # briggs = 0.5 if sidx in [0] else 0.0
        briggs = 0.0
        if sidx in [0]:
            auto_mask = 2
            auto_threshold = 1
        # elif sidx in [5,6]:
        #     auto_mask = 2.5
        #     auto_threshold = 1
        else:
            auto_mask = 3
            auto_threshold = 1
        # reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        slfcal_rnd2_obj = MSselfcal(msfile, time_intervals, time_intervals_major_avg, time_intervals_minor_avg,
                                    spws[sidx], N1, N2, workdir, image_marker='round2', niter=1000,
                                    briggs=briggs,
                                    minuv_l=uvmin_l_list[sidx]/3.0,
                                    minuvw_m=8,
                                    maxuvw_m=1500,
                                    auto_mask=auto_mask, auto_threshold=auto_threshold,
                                    pols=pols,
                                    data_column="CORRECTED_DATA"
                                    )
        # beam_size=bmsize)
        slfcal_rnd2_obj.run()
        slfcal_rnd2_objs.append(slfcal_rnd2_obj)

        caltb = os.path.join(workdir, f"caltb_rnd2_int_sp{spwstr}.pha")
        if os.path.exists(caltb): os.system('rm -rf ' + caltb)
        gaincal(vis=msfile, caltable=caltb, selectdata=True,
                uvrange=uvmin_l_str[sidx],
                spw=spws[sidx],
                combine="scan", antenna='0~12&0~12', refant='0', solint=f'{tdur / N2 / N1:.0f}s', refantmode="flex",
                gaintype='G',
                gaintable=caltbs,
                minsnr=1.0, calmode='p', append=False)
        if os.path.exists(caltb): caltbs.append(caltb)

        caltb = os.path.join(workdir, f"caltb_rnd2_inf_sp{spwstr}.amp")
        if os.path.exists(caltb): os.system('rm -rf ' + caltb)
        gaincal(vis=msfile, caltable=caltb, selectdata=True,
                uvrange=uvmin_l_str[sidx],
                spw=spws[sidx],
                # combine="scan", antenna='0~12&0~12', refant='0', solint=f'{tdur / N2 / N1:.0f}s', refantmode="flex",
                combine="scan", antenna='0~12&0~12', refant='0', solint='inf', refantmode="flex",
                gaintype='G',
                gaintable=caltbs,
                minsnr=1.0, calmode='a', append=False)
        if os.path.exists(caltb): caltbs.append(caltb)

        if do_diskslfcal and (not diskslfcal_first[sidx]):
            delmod(vis=msfile)
            run_start_time_disk_slfcal = datetime.now()
            log_print('INFO', f"Disk self-calibration for SPW {spws[sidx]} ...")
            in_fits_files = sorted(glob(slfcal_init_objs[sidx].model_minor_name_str + '-t????-model.fits'))
            for sp in sp_index.split(','):
                run_start_time_disk_slfcal_sp = datetime.now()
                sp = int(sp)
                model_dir_disk_name_str = slfcal_init_objs[sidx].model_minor_name_str.replace(format_spw(spws[sidx]),
                                                                                              f'{sp:02d}') + '_adddisk'
                dsz = float(dsize[sp].rstrip('arcsec'))
                fdn = fdens[sp]
                if len(in_fits_files) > 0:
                    for in_fits in in_fits_files:
                        out_fits = in_fits.replace(slfcal_init_objs[sidx].model_minor_name_str, model_dir_disk_name_str)
                        add_convolved_disk_to_fits(in_fits, out_fits, dsz, fdn, ignore_data=True, bmaj=bmsize/3600.)
                    cmd = f"wsclean -predict -reorder -spws {sp}" + f" \
                        -name {model_dir_disk_name_str} -quiet -intervals-out {N1 * N2} {msfile}"
                    log_print('INFO', f"Running wsclean predict: {cmd}")
                    os.system(cmd)
                else:
                    log_print('WARNING',
                              f"No model files found for SPW {spws[sidx]}. Skipping disk self-calibration...")
            run_end_time_disk_slfcal = datetime.now()
            elapsed_time = run_end_time_disk_slfcal - run_start_time_disk_slfcal
            elapsed_time_disk_slfcal = elapsed_time.total_seconds() / 60
            log_print('INFO',
                      f"Disk self-calibration for SPW {spwstr}: completed in {elapsed_time_disk_slfcal:.1f} minutes")

            caltb = os.path.join(workdir, f"caltb_disk_inf_sp{spwstr}.pha")
            if os.path.exists(caltb): os.system('rm -rf ' + caltb)
            gaincal(vis=msfile, caltable=caltb, selectdata=True,
                    uvrange=uvmax_l_str[sidx],
                    spw=spws[sidx],
                    combine="scan",
                    antenna='0~12&0~12', refant='0', solint='inf', refantmode="flex",
                    gaintype='G',
                    gaintable=caltbs,
                    minsnr=1,
                    minblperant=2,
                    calmode='p', append=False)
            if os.path.exists(caltb):
                caltbs.append(caltb)

        if len(caltbs) > ncaltbs:
            applycal(vis=msfile, selectdata=True, antenna="0~12",
                     spw=spws[sidx],
                     gaintable=caltbs, interp='linear',
                     calwt=False, applymode="calonly")
            # ncaltbs = len(caltbs)
            log_print('INFO', f"Round 2 Feature Selfcal for SPW {spws[sidx]}: solution applied to the data")
        else:
            log_print('WARNING', f"Round 2 Feature Selfcal for SPW {spws[sidx]}: no new solutions, skipping applycal")

        if do_diskslfcal and (not diskslfcal_first[sidx]):
            # # step 5: subtract the disk model from the viz
            run_start_time_disk_slfcal_subtract = datetime.now()
            uvsub(vis=msfile)
            run_end_time_disk_slfcal_subtract = datetime.now()
            elapsed_time = run_end_time_disk_slfcal_subtract - run_start_time_disk_slfcal_subtract
            elapsed_time_disk_slfcal_subtract = elapsed_time.total_seconds() / 60
            log_print('INFO',
                      f"Disk subtraction for SPW {spwstr}: completed in {elapsed_time_disk_slfcal_subtract:.1f} minutes")

        caltbs_all.append(caltbs)

        briggs = 0.0
        gain = 0.3
        # reffreq, cdelt4_real, bmsize = freq_setup.get_reffreq_and_cdelt(spws[sidx], return_bmsize=True)
        # imname_strlist = ["eovsa", "major", f"{msfilename}", f"sp{spwstr}", 'final','disk']
        imname_strlist = ["eovsa", "major", f"{msfilename}", f"sp{spwstr}", 'final']
        imname = '-'.join(imname_strlist)
        # briggs = 0.5 if sidx in [0] else 0.0
        flagdata(vis=msfile, mode="tfcrop", spw=spws[sidx], action='apply', display='',
                 timecutoff=3.0, freqcutoff=3.0, maxnpieces=2, flagbackup=False)
        if segmented_imaging:
            clean_obj = ww.WSClean(msfile)
            clean_obj.setup(size=1024, scale="2.5asec", pol=pols,
                            weight_briggs=briggs,
                            niter=2000,
                            mgain=0.85,
                            gain=gain,
                            # gain=0.2,
                            data_column='CORRECTED_DATA',
                            name=os.path.join(workdir, imname),
                            multiscale=True,
                            multiscale_gain=0.3,
                            # multiscale_gain=0.5,
                            multiscale_scales=[0, 5, 12.5],
                            auto_mask=2, auto_threshold=1,
                            # auto_mask=1.5, auto_threshold=1,
                            # minuv_l=500,
                            minuv_l=uvmin_l_list[sidx]/4.0,
                            # minuv_l=None,
                            minuvw_m=8,
                            maxuvw_m=1500,
                            intervals_out=N1,
                            no_negative=True, quiet=True,
                            circular_beam=True,
                            spws=sp_index)
            clean_obj.run(dryrun=False)

            fitsname = sorted(glob(os.path.join(workdir, imname + '*image.fits')))
            fitsname_helio = [f.replace('image.fits', 'image.helio.fits') for f in fitsname]
            fitsname_helio_ref_daily = [f.replace('image.fits', 'image.helio.ref_daily.fits') for f in fitsname]
            hf.imreg(vis=msfile, imagefile=fitsname, fitsfile=fitsname_helio, timerange=timeranges, toTb=True)
            for eoidx, (fits_helio, fits_helio_ref_daily) in enumerate(zip(fitsname_helio, fitsname_helio_ref_daily)):
                solar_diff_rot_heliofits(fits_helio, reftime_daily, fits_helio_ref_daily,
                                         in_time=time_intervals_major_avg[eoidx])
            fitsfilefinal = fitsname_helio_ref_daily

            log_print('INFO', f"Final imaging for SPW {spws[sidx]}: completed")
        else:
            clean_obj = ww.WSClean(msfile)
            imaging_obj = MSselfcal(msfile, time_intervals, time_intervals_major_avg, time_intervals_minor_avg,
                                    spws[sidx], N1, N2, workdir, image_marker='temp', niter=1000,
                                    briggs=briggs,
                                    gain=gain,
                                    minuv_l=uvmin_l_list[sidx] / 3.0,
                                    # minuv_l=None,
                                    minuvw_m=8,
                                    maxuvw_m=1500,
                                    auto_mask=2, auto_threshold=1,
                                    data_column="CORRECTED_DATA",
                                    pols=pols,
                                    reftime_daily=reftime_daily
                                    )
            imaging_obj.run()
            imaging_objs.append(imaging_obj)

            uvsub(vis=msfile)
            delmod(vis=msfile)
            model_dir_name_str = imaging_objs[sidx].model_ref_name_str
            cmd = "wsclean -predict -reorder -spws " + sp_index + f" \
                -name {model_dir_name_str} -intervals-out {N1} {msfile}"
            os.system(cmd)
            uvsub(vis=msfile, reverse=True)  # Now the viz data contains residuals + corrected models

            clean_obj.setup(size=1024, scale="2.5asec",
                            pol=pols,
                            weight_briggs=briggs,
                            niter=2000,
                            mgain=0.85,
                            gain=gain,
                            # gain=0.2,
                            data_column='CORRECTED_DATA',
                            name=os.path.join(workdir, imname),
                            multiscale=True,
                            multiscale_gain=0.3,
                            # multiscale_gain=0.5,
                            multiscale_scales=[0, 5, 12.5],
                            auto_mask=2, auto_threshold=1,
                            minuv_l=uvmin_l_list[sidx] / 4.0,
                            # minuv_l=None,
                            minuvw_m=8,
                            maxuvw_m=1500,
                            intervals_out=1,
                            no_negative=True, quiet=True,
                            circular_beam=True,
                            spws=sp_index)
            clean_obj.run(dryrun=False)

            log_print('INFO', f"Final imaging for SPW {spws[sidx]}: completed")

            fitsname = sorted(glob(os.path.join(workdir, imname + '*image.fits')))
            fitsname_helio = [f.replace('image.fits', 'image.helio.fits') for f in fitsname]
            hf.imreg(vis=msfile, imagefile=fitsname, fitsfile=fitsname_helio, timerange=[viz_timerange], toTb=True)
            fitsfilefinal = fitsname_helio

        sp_st, sp_ed = spws[sidx].split('~')
        dszs = [float(dsize[int(sp)].rstrip('arcsec')) for sp in range(int(sp_st), int(sp_ed) + 1)]
        fdns = [fdens[int(sp)] for sp in range(int(sp_st), int(sp_ed) + 1)]
        # fdns = [fdens[int(sp)]*50 for sp in range(int(sp_st), int(sp_ed)+1)]
        synfitsfiles = []
        for eoidx, eofile in enumerate(fitsname):
            datetimestr = time_intervals_major_avg[eoidx].datetime.strftime('%Y%m%dT%H%M%SZ')
            if segmented_imaging:
                synfitsfile = os.path.join(imgoutdir,
                                           f"eovsa.synoptic.{datetimestr}.s{spwstr}.tb.disk.fits")
            else:
                synfitsfile = os.path.join(imgoutdir,
                                           f"eovsa.synoptic_daily.{date_str}T200000Z.s{spwstr}.tb.disk.fits")
            synfitsfiles.append(synfitsfile)
        add_convolved_disk_to_fits(fitsfilefinal, synfitsfiles, dszs, fdns, ignore_data=False, toTb=True,
                                   rfreq=reffreq, bmaj=bmsize / 3600.)
        if not segmented_imaging:
            if len(synfitsfiles) > 0:
                outfits_all.append(synfitsfiles[0])
            else:
                outfits_all.append(None)

        split(vis=msfile, outputvis=msfile_sp, spw=spws[sidx],datacolumn='corrected')

        run_end_time = datetime.now()
        elapsed_time = run_end_time - run_start_time
        elapsed_time_total = elapsed_time.total_seconds() / 60
        log_print('INFO', f"Pipeline for SPW {spws[sidx]}: completed in {elapsed_time_total:.1f} minutes")

    if segmented_imaging:
        for spw in spws:
            spwstr = format_spw(spw)
            outfits= os.path.join(imgoutdir, f'eovsa.synoptic_daily.{date_str}T200000Z.s{spwstr}.tb.disk.fits')
            synfitsfiles = sorted(glob(os.path.join(imgoutdir,
                                                    f"eovsa.synoptic.{date_str[:-1]}?T??????Z.s{spwstr}.tb.disk.fits")))
            if len(synfitsfiles) > 0:
                log_print('INFO', f"Merging synoptic images for SPW {spwstr} to {outfits} ...")
                merge_FITSfiles(synfitsfiles, outfits, overwrite=True, snr_threshold=10)
                outfits_all.append(outfits)
            else:
                outfits_all.append(None)
                log_print('WARNING', f"No synoptic images found for SPW {spwstr}. Skipping merge_FITSfiles.")

    ms2concat = glob(f'{msname}.sp*.slfcaled.ms')
    if os.path.exists(outputvis):
        if overwrite:
            os.system('rm -rf ' + outputvis)
            concat(vis=ms2concat, concatvis=outputvis)
            log_print('INFO', f"Concatenated MS file written to {outputvis}")
        else:
            log_print('WARNING', f"Output MS file {outputvis} already exists. Skipping concatenation.")
    else:
        concat(vis=ms2concat, concatvis=outputvis)
        log_print('INFO', f"Concatenated MS file written to {outputvis}")

    log_print('INFO', f"Moving calibration tables to {slfcaltbdir} ...")
    caltbs_comb = [item for sublist in caltbs_all for item in sublist]
    for caltb in caltbs_comb:
        if os.path.exists(caltb): os.system(f'mv {caltb} {slfcaltbdir}/')
    if os.path.isdir(msfile + '.flagversions') == True:
        os.system(f'mv {msfile}.flagversions {os.path.dirname(outputvis)}/')

    # if clearlargecache:
    #     clearcache = False
    #     shutil.rmtree(os.path.join(subdir, 'images'), ignore_errors=True)
    #     shutil.rmtree(os.path.join(subdir, '*.ms'), ignore_errors=True)
    #     shutil.rmtree(os.path.join(subdir, '*.model'), ignore_errors=True)
    #     print(f'files in {os.path.join(subdir, "images")} removed.')
    #     print(f'files in {os.path.join(subdir, "*.ms")} removed.')
    #     print(f'files in {os.path.join(subdir, "*.model")} removed.')

    return outfits_all

if __name__ == '__main__':
    '''
    this code is trying to address the issue of smearing effect in all-day synthesis images.
    description of the issue: https://www.ovsa.njit.edu/wiki/index.php/All-Day_Synthesis_Issues
    This code is built on python 3.8.
    '''
    description = 'this code is trying to address the issue of smearing effect in all-day synthesis images. Description of the issue: https://www.ovsa.njit.edu/wiki/index.php/All-Day_Synthesis_Issues. This code is built on python 3.8.'
    parser = argparse.ArgumentParser(
        description="Executes the EOVSA data processing pipeline for solar observation data.")
    parser.add_argument('vis', type=str, help='Path to the input measurement set (MS) data.')
    parser.add_argument('--outputvis', type=str, default='', help='Output path for the processed MS.')
    parser.add_argument('--workdir', type=str, help='Working directory for intermediate data.')
    parser.add_argument('--slfcaltbdir', type=str, help='Directory for storing calibration tables.')
    parser.add_argument('--imgoutdir', type=str, help='Output directory for image files.')
    parser.add_argument('--figoutdir', type=str, help='Output directory for figures.')
    parser.add_argument('--clearcache', action='store_true', help='Clears all cache files after processing.')
    parser.add_argument('--clearlargecache', action='store_true',
                        help='Clears the large cache files after processing. If True, clearcache will be set to False.')
    parser.add_argument('--pols', type=str, default='XX', help='Polarization types to process.')
    parser.add_argument('--mergeFITSonly', action='store_true', help='Skips processing and only merges FITS files.')
    # parser.add_argument('--verbose', action='store_true', help='Enables verbose output during processing.')
    # parser.add_argument('--do_diskslfcal', action='store_true', help='Performs disk self-calibration.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrites existing files.')
    parser.add_argument('--niter_init', type=int, default=200, help='Initial number of iterations for imaging.')
    parser.add_argument('--ncpu', type=str, default='auto',
                        help="Specifies the number of CPUs for parallel processing.")
    parser.add_argument('--tr_series_imaging', nargs='*', help='Time ranges for imaging, expects a list of tuples.')
    parser.add_argument('--spws_imaging', nargs='*', help='Spectral windows selected for imaging.')
    parser.add_argument('--hanning', action='store_true', help='Applies Hanning smoothing to the data.')
    parser.add_argument('--do_sbdcal', action='store_true', help='Perform single-band delay calibration.')

    args = parser.parse_args()

    pipeline_run(
        vis=args.vis,
        outputvis=args.outputvis,
        workdir=args.workdir,
        slfcaltbdir=args.slfcaltbdir,
        imgoutdir=args.imgoutdir,
        pols=args.pols,
        hanning=args.hanning,
        do_sbdcal=args.do_sbdcal
    )
