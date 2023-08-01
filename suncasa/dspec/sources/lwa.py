import h5py
import numpy as np
import datetime
from astropy.time import Time

def rebin1d(arr, new_len):
    shape = (new_len, len(arr) // new_len)
    return arr.reshape(shape).mean(1)


def rebin2d(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def read_data(filename, stokes='I', verbose=True, timerange=[], freqrange=[], timebin=10, freqbin=1):
    data = h5py.File(filename, 'r')
    freqs = data['Observation1']['Tuning1']['freq'][:]
    ts = data['Observation1']['time'][:]
    times_mjd = Time([datetime.datetime.fromtimestamp(t[0]+t[1]) for t in ts]).mjd
    idx0, = np.where(times_mjd > 50000.) # filter out those prior to 1995 (obviously wrong for OVRO-LWA)
    if verbose:
        print('Data time range is from {0:s} to {1:s}'.format(Time(times_mjd[idx0][0], format='mjd').isot, 
            Time(times_mjd[idx0][-1], format='mjd').isot))
        print('Data has {0:d} time stamps and {1:d} frequency channels'.format(len(times_mjd[idx0]), len(freqs)))

    # Select time range
    if len(timerange) > 0:
        if type(timerange) == list or type(timerange) == str:
            try:
                timerange = Time(timerange)
                t0 = timerange[0].mjd
                t1 = timerange[1].mjd
                ti0 = np.argmin(np.abs(times_mjd - t0))
                ti1 = np.argmin(np.abs(times_mjd - t1)) 
            except:
                print('timerange not parsed correctly. Use the full range.')
                ti0 = 0
                ti1 = len(times_mjd) 
    else:
        ti0=0
        ti1=len(times_mjd)

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
    stokes_valid = ['XX', 'YY', 'RR', 'LL', 'I', 'Q', 'U', 'V']
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

    idx, = np.where(times_mjd > 50000.) # filter out those prior to 1995 (obviously wrong for OVRO-LWA)
    times_mjd = times_mjd[idx]
    spec = spec[idx]

    nt, nf = spec.shape
    nt_new, nf_new = (nt // timebin, nf // freqbin)
    # TODO: for now I have just ignored the rest of the data that falls outside of the whole factor of timebin * nt_new or freqbin * nf_new 
    spec_new = rebin2d(spec[:nt_new*timebin, :nf_new*freqbin], (nt_new, nf_new))
    times_mjd_new = rebin1d(times_mjd[:nt_new * timebin], nt_new)
    freqs_new = rebin1d(freqs[:nf_new * freqbin], nf_new)
    if verbose:
        print('Selected time range is from {0:s} to {1:s}'.format(Time(times_mjd_new[0], format='mjd').isot, Time(times_mjd_new[-1], format='mjd').isot))
        print('Output data has {0:d} time stamps and {1:d} frequency channels'.format(len(times_mjd_new), len(freqs_new)))
    return spec_new.T/1e4, times_mjd_new, freqs_new, stokes
