import os
import h5py
import numpy as np
import datetime
from astropy.time import Time
import itertools
from astropy import units as u
import pandas as pd
from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz

def rebin1d(arr, new_len):
    shape = (new_len, len(arr) // new_len)
    return arr.reshape(shape).mean(1)


def rebin2d(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def timestamp_to_mjd(times):
    # This is from Ivey Davis's BeamTools.py
    t_flat = np.array(list(itertools.chain(*times)))
    ts_inds = np.linspace(0, len(times)-1, len(times), dtype  = int)*2
    other_inds = ts_inds + 1
    ts = Time(t_flat[ts_inds],format = 'unix') + t_flat[other_inds]* u.s
    ts = ts.mjd
    return ts


def read_data(filename, stokes='I', timerange=[], freqrange=[], timebin=1, freqbin=1, verbose=True, 
            flux_factor_file=None, bkg_file=None,  do_pb_correction=False, 
            flux_factor_calfac_x = None, flux_factor_calfac_y = None, bkg_flux_arr = None):
    '''
    :param filename: name of the OVRO-LWA hdf5 beamforming file; 
              This can be a string (single file) or a list of strings (multiple files)
    :param stokes: currently supporting 'XX', 'YY', 'I', 'Q', 'U', 'V', 'IV' 
    :param timerange: list of [start_time, end_time], start_time and end_time should be recognized by astropy.time.Time
            e.g., ['2023-09-22T18:00', '2023-09-22T18:10']
    :param freqrange: list of [start_frequency, end_frequency] in MHz. Example: [23, 82]
    :param timebin: number to bin in time
    :param freqbin: number to bin in frequency
    :param verbose: if True, print extra information
    :param flux_factor_file: Path to the csv file that contains the flux correction factors
    :param bkg_file: Path to the csv file that contains the raw background (off-Sun scan) measurements (before scaling); 
            currently it only contains Stokes I.
    :param do_pb_correction: if True, apply primary beam correction. Currently only use the analytical form in Stokes I only.
    :param flux_factor_calfac_x: user input correction factor for the X polarization
    :param flux_factor_calfac_y: user input correction factor for the Y polarization
    :param bkg_flux_arr: user input background flux in Jy
    '''
    # Check the input filename
    if type(filename) == str:
        filename = [filename]

    obs = EarthLocation.of_site('ovro')
    filelist = []
    for ll in filename:
        if not os.path.exists(ll):
            print('{} does not exist. Skip this one.'.format(ll))
        else:
            filelist.append(ll)

    if not filelist:
        print('I cannot find any file. Abort.')
        return False
    else:
        filelist.sort()
    
    n0 = 0 # this is the index for the first file being read successfully
    firstset_read = False 
    for n, file in enumerate(filelist):
        if verbose:
            print('Processing {0:d} of {1:d} files'.format(n+1, len(filelist)))
        try:
            data = h5py.File(file, 'r',swmr=True)
        except:
            print('Cannot read {0:s}. Skip this file.'.format(file))
            continue
        freqs = data['Observation1']['Tuning1']['freq'][:]
        ts = data['Observation1']['time'][:]
        # The following line works the same way as timestamp_to_mjd(), but a bit too slow
        # times_mjd = np.array([(Time(t[0], format='unix') + TimeDelta(t[1], format='sec')).mjd for t in ts])
        times_mjd = timestamp_to_mjd(ts)
        idx0, = np.where(times_mjd > 50000.) # filter out those prior to 1995 (obviously wrong for OVRO-LWA)

        # read the flux factors file if provided
        if not (flux_factor_file is None): 
            try:
                out = pd.read_csv(flux_factor_file)
                calfac_x = np.array(out['calfac_x'])
                calfac_y = np.array(out['calfac_y'])
            except:
                print('Failed in reading the flux factor csv file. Setting correction factors to unity.')
        else:
            print('Flux factor csv file does not exist. Setting correction factors to unity.')
            calfac_x = np.ones_like(freqs)
            calfac_y = np.ones_like(freqs)

        if not (flux_factor_calfac_x is None) and not (flux_factor_calfac_y is None):
            # user input correction factor
            calfac_x = calfac_x*flux_factor_calfac_x
            calfac_y = calfac_y*flux_factor_calfac_y
        
        if not (bkg_flux_arr is None):
            # add the user input background flux
            bkg_flux += bkg_flux_arr

        # read background flux file if provided
        if not (bkg_file is None): 
            try:
                out = pd.read_csv(bkg_file)
                bkg_flux = out['bkg_flux']
                print('Using the provided raw aackground flux csv file.')
            except:
                print('Failed in reading the background flux csv file. Setting background flux to zero.')
        else:
            print('No background csv file provided. Setting background flux to zero.')
            bkg_flux = np.zeros_like(freqs)

        if verbose:
            print('Data time range is from {0:s} to {1:s}'.format(Time(times_mjd[idx0][0], format='mjd').isot, 
                Time(times_mjd[idx0][-1], format='mjd').isot))
            print('Data has {0:d} time stamps and {1:d} frequency channels'.format(len(times_mjd[idx0]), len(freqs)))

        # Select time range
        if len(timerange) > 0:
            try:
                timerange_obj = Time(timerange)
                # Take the larger value of the supplied start time and the first time stamp of the data
                t0 = max(timerange_obj[0].mjd, min(times_mjd[idx0]))
                # Take the smaller value of the supplied end time and the last time stamp of the data
                t1 = min(timerange_obj[1].mjd, max(times_mjd[idx0]))
                ti0 = np.argmin(np.abs(times_mjd - t0))
                ti1 = np.argmin(np.abs(times_mjd - t1)) 
                if ti1 - ti0 < timebin:
                    print('Selected number of time samples {0:d} is less than the timebin {1:d}. Skip this file.'.format(ti1-ti0, timebin))
                    if not firstset_read:
                        n0 += 1
                    continue
                if verbose:
                    print('Selected time range is from {0:s} to {1:s}'.format(Time(times_mjd[ti0], format='mjd').isot, 
                                                                              Time(times_mjd[ti1], format='mjd').isot))
            except:
                print('timerange not parsed correctly. Use the full range in the data.')
                ti0 = 0
                ti1 = len(times_mjd) 
        else:
            ti0=0
            ti1=len(times_mjd)

        firstset_read = True
        times_mjd = times_mjd[ti0:ti1] 

        # Select frequency range
        if len(freqrange) > 0:
            if type(freqrange) == list:
                try:
                    f0 = freqrange[0]
                    f1 = freqrange[1]
                    if f0 > 100. or f1 > 100.:
                        # I am assuming input frequency range is in Hz
                        print('Input frequency range is greater than 100. Assuming unit in Hz.')
                    else:
                        # I am assuming input frequency range is in MHz
                        print('Input frequency range is less than 100. Assuming unit in MHz.')
                        f0 *= 1e6
                        f1 *= 1e6
                    fi0 = np.argmin(np.abs(freqs - f0))
                    fi1 = np.argmin(np.abs(freqs - f1)) 
                except:
                    print('freqrange not parsed correctly. Use the full range.')
                    fi0 = 0
                    fi1 = len(freqs) 
        else:
            fi0=0
            fi1=len(freqs)

        freqs = freqs[fi0:fi1] 

        # select stokes
        stokes_valid = ['XX', 'YY', 'I', 'Q', 'U', 'V', 'IV']
        if verbose:
            print('Reading dynamic spectrum for stokes {0:s}'.format(stokes))

        if stokes not in stokes_valid:
            raise Exception("Provided Stokes {0:s} is not in 'XX, YY, RR, LL, I, Q, U, V'".format(stokes))
        if stokes.upper() == 'XX':
            spec = data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] / calfac_x[None, fi0:fi1]
            stokes_out = ['XX']
        if stokes.upper() == 'YY':
            spec = data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1] / calfac_y[None, fi0:fi1]
            stokes_out = ['YY']
        if stokes.upper() == 'I' or stokes.upper() == 'IV':
            spec_I = (data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] / calfac_x[None, fi0:fi1] + \
                   data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1] / calfac_y[None, fi0:fi1]) / 2. - \
                   bkg_flux[None, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.) 
            if verbose:
                print('Median of the subtracted background flux (Jy)', np.median(bkg_flux[None, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.)))
                print('RMS of the subtracted background flux (Jy)', np.std(bkg_flux[None, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.)))
            if do_pb_correction:
                pbfacs = np.ones_like(times_mjd)
                t0 = times_mjd[0]
                t1 = times_mjd[-1]
                if t1-t0 > 5./1440.:
                    nstep = int((t1-t0)/(5./1440.))
                    ts_ref = np.linspace(t0, t1, nstep)
                    pbfacs_ref = np.ones_like(ts_ref) 
                    print('Duration of the file is {0:.1f} hours, interpolating into {1:d} steps.'.format((t1-t0)*24., nstep))
                    for i, t_ref in enumerate(ts_ref):
                        sun_loc = get_body('sun', Time(t_ref, format='mjd'), location=obs)
                        alt = sun_loc.transform_to(AltAz(obstime=Time(t_ref, format='mjd'), location=obs)).alt.radian
                        if np.degrees(alt) > 5.:
                            pbfacs_ref[i]=np.sin(alt)**1.6
                        else:
                            print('Warning! Calculated solar altitude is lower than 5 degrees. Something is wrong with the data (non-solar)?') 
                            pbfacs_ref[i]=np.sin(np.radians(5.))**1.6

                    pbfacs = np.interp(times_mjd, ts_ref, pbfacs_ref)
                else:
                    print('Duration of the file is {0:.1f} minutes, no interpolation will be done.'.format((t1-t0)*24.*60.))
                    t_ref = (t0+t1)/2.
                    sun_loc = get_body('sun', Time(t_ref, format='mjd'), location=obs)
                    alt = sun_loc.transform_to(AltAz(obstime=Time(t_ref, format='mjd'), location=obs)).alt.radian
                    if np.degrees(alt) > 5.:
                        pbfacs *= np.sin(alt)**1.6
                    else:
                        print('Warning! Calculated solar altitude is lower than 5 degrees. Something is wrong with the data (non-solar)?') 
                        pbfacs *= np.sin(np.radians(5.))**1.6
                spec_I /= pbfacs[:, None]

            if stokes.upper() == 'IV':
                spec_V = data['Observation1']['Tuning1']['XY_imag'][ti0:ti1, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.)
                spec = np.stack((spec_I, spec_V), axis=2)
                stokes_out = ['I', 'V']
            else:
                spec = spec_I
                stokes_out = ['I'] 
        if stokes == 'V':
            spec = data['Observation1']['Tuning1']['XY_imag'][ti0:ti1, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.)
            stokes_out = ['V'] 
        if stokes == 'Q':
            spec = (data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] / calfac_x[None, fi0:fi1] - \
                    data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1] / calfac_y[None, fi0:fi1]) / 2. 
            stokes_out = ['Q'] 
        if stokes == 'U':
            spec = data['Observation1']['Tuning1']['XY_real'][ti0:ti1, fi0:fi1] / ((calfac_x[None, fi0:fi1] + calfac_y[None, fi0:fi1]) / 2.)
            stokes_out = ['U'] 

        idx, = np.where(times_mjd > 50000.) # filter out those prior to 1995 (obviously wrong for OVRO-LWA)
        times_mjd = times_mjd[idx]
        spec = spec[idx]

        if spec.ndim == 2:
            # TODO: for now I have just ignored the rest of the data that falls outside of the whole factor of timebin * nt_new or freqbin * nf_new 
            nt, nf = spec.shape
            npol = 1
            nt_new, nf_new = (nt // timebin, nf // freqbin)
            spec_new = rebin2d(spec[:nt_new*timebin, :nf_new*freqbin], (nt_new, nf_new))
        if spec.ndim == 3:
            nt, nf, npol = spec.shape
            spec_new_ = [] 
            for i in range(npol):
                nt_new, nf_new = (nt // timebin, nf // freqbin)
                spec_ = spec[:, :, i]
                spec_new_.append(rebin2d(spec_[:nt_new*timebin, :nf_new*freqbin], (nt_new, nf_new)))
            spec_new = np.stack(spec_new_, axis=2)
            
        spec_new = np.transpose(spec_new).reshape((npol, 1, nf_new, nt_new)) / 1e4
        times_mjd_new = rebin1d(times_mjd[:nt_new * timebin], nt_new)
        freqs_new = rebin1d(freqs[:nf_new * freqbin], nf_new)
        if firstset_read:
            if n == n0:
                nfreq0 = len(freqs_new)
                freqs_out = freqs_new
                spec_out = spec_new
                times_mjd_out = times_mjd_new
            elif n > n0:
                if len(freqs_new) != nfreq0:
                    print('Something is wrong in concatenating {}'.format(file)) 
                    print('Dimension of the output frequency {0:d} does not match that of the first file {1:d)'.format(len(freqs_new), nfreq0)) 
                    continue
                else:
                    spec_out = np.concatenate((spec_out, spec_new), axis=3)
                    times_mjd_out = np.concatenate((times_mjd_out, times_mjd_new))
        else:
            continue

    if firstset_read:
        if verbose:
            print('Output time range is from {0:s} to {1:s}'.format(Time(times_mjd_out[0], format='mjd').isot, Time(times_mjd_out[-1], format='mjd').isot))
            print('Output data has {0:d} time stamps and {1:d} frequency channels'.format(len(times_mjd_out), len(freqs_out)))
        return spec_out, times_mjd_out, freqs_out, stokes_out, calfac_x, calfac_y, bkg_flux
    else:
        return False
