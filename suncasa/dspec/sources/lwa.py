import os
import h5py
import numpy as np
import datetime
from astropy.time import Time
import itertools
from astropy import units as u

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


def read_data(filename, stokes='I', timerange=[], freqrange=[], timebin=1, freqbin=1, verbose=True):
    '''
    filename: name of the OVRO-LWA hdf5 beamforming file; 
              This can be a string (single file) or a list of strings (multiple files)
    stokes: currently supporting 'XX', 'YY', 'I', 'Q', 'U', 'V', 'IV' 
    timerange: list of [start_time, end_time], start_time and end_time should be recognized by astropy.time.Time
            e.g., ['2023-09-22T18:00', '2023-09-22T18:10']
    freqrange: list of [start_frequency, end_frequency] in MHz. Example: [23, 82]
    timebin: number to bin in time
    freqbin: number to bin in frequency
    verbose: if True, print extra information
    '''
    # Check the input filename
    if type(filename) == str:
        filename = [filename]

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
        data = h5py.File(file, 'r')
        freqs = data['Observation1']['Tuning1']['freq'][:]
        ts = data['Observation1']['time'][:]
        # The following line works the same way as timestamp_to_mjd(), but a bit too slow
        # times_mjd = np.array([(Time(t[0], format='unix') + TimeDelta(t[1], format='sec')).mjd for t in ts])
        times_mjd = timestamp_to_mjd(ts)
        idx0, = np.where(times_mjd > 50000.) # filter out those prior to 1995 (obviously wrong for OVRO-LWA)

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
        if stokes == 'XX':
            spec = data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1]
        if stokes == 'YY':
            spec = data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1]
        if stokes == 'I':
            spec = data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] + \
                   data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1]
        if stokes == 'V':
            spec = 2 * data['Observation1']['Tuning1']['XY_imag'][ti0:ti1, fi0:fi1]
        if stokes == 'Q':
            spec = data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] - \
                    data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1]
        if stokes == 'U':
            spec = 2 * data['Observation1']['Tuning1']['XY_real'][ti0:ti1, fi0:fi1]

        # Multiple polarizations
        if stokes.upper() == 'IV':
            spec_I = data['Observation1']['Tuning1']['XX'][ti0:ti1, fi0:fi1] + \
                   data['Observation1']['Tuning1']['YY'][ti0:ti1, fi0:fi1]
            spec_V = 2 * data['Observation1']['Tuning1']['XY_imag'][ti0:ti1, fi0:fi1]
            spec = np.stack((spec_I, spec_V), axis=2)
            stokes_out = ['I', 'V']

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
        return spec_out, times_mjd_out, freqs_out, stokes_out
    else:
        return False
