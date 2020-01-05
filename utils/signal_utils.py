from scipy.signal import butter, freqz, filtfilt
import matplotlib.pyplot as plt
import numpy as np


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low')
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def lowps_filter(data, cutoff, fs, ix):
    x = data[ix]
    y = butter_lowpass_filter(x, cutoff * fs, fs, order=6)
    return {'idx': ix, 'y': y}


def low_pass_filter(t, data, fs=1. / 4, cutoff=1. / 60, order=6, showplot=False):
    # # Filter requirements.
    # order = 6
    # fs = 1. / 4  # sample rate, Hz
    # cutoff = 1. / 60  # desired cutoff frequency of the filter, Hz

    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_lowpass(cutoff, fs, order)

    # Plot the frequency response.
    w, h = freqz(b, a, worN=10000)
    if showplot:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(1 / (0.5 * fs * w / np.pi), np.abs(h), 'b')
        plt.plot(1 / cutoff, 0.5 * np.sqrt(2), 'ko')
        plt.axvline(1 / cutoff, color='k')
        # plt.xlim(0, 0.5 * fs)
        plt.xlim(1. / (0.5 * fs), 100. / (0.5 * fs))
        plt.title("Lowpass Filter Frequency Response")
        # plt.xlabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.grid()

    # Filter the data, and plot both the original and filtered signals.
    y = butter_lowpass_filter(data, cutoff, fs, order)

    if showplot:
        plt.subplot(2, 1, 2)
        plt.plot(t, data, 'b-', label='data')
        plt.plot(t, y, 'g-', linewidth=2, label='filtered data')
        plt.xlabel('Time [sec]')
        plt.grid()
        plt.legend()

        plt.subplots_adjust(hspace=0.35)
        plt.show()
    return y


def bandpass_filter(t, data, fs=1. / 4, cutoff=1. / 60, order=6, showplot=False):
    from scipy.signal import butter, freqz, filtfilt
    import matplotlib.pyplot as plt

    def butter_bandpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='bandpass')
        return b, a

    def butter_bandpass_filter(data, cutoff, fs, order=5):
        b, a = butter_bandpass(cutoff, fs, order=order)
        y = filtfilt(b, a, data)
        return y

    # # Filter requirements.
    # order = 6
    # fs = 1. / 4  # sample rate, Hz
    # cutoff = 1. / 60  # desired cutoff frequency of the filter, Hz

    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_bandpass(cutoff, fs, order)

    # Plot the frequency response.
    w, h = freqz(b, a, worN=10000)
    if showplot:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(1./(0.5 * fs * w / np.pi), np.abs(h), 'b')
        plt.plot(1.0/cutoff[0], 0.5 * np.sqrt(2), 'ko')
        plt.plot(1.0/cutoff[1], 0.5 * np.sqrt(2), 'ko')
        plt.axvline(1.0/cutoff[0], color='k')
        plt.axvline(1.0/cutoff[1], color='k')
        plt.xlim(1. / (0.5 * fs), 100. / (0.5 * fs))
        plt.title("Bandpass Filter Frequency Response")
        plt.xlabel('Time [sec]')
        plt.grid()

    y = butter_bandpass_filter(data, cutoff, fs, order)

    if showplot:
        plt.subplot(2, 1, 2)
        plt.plot(t, data, 'b-', label='data')
        plt.plot(t, y + np.nanmean(data), 'g-', linewidth=2, label='filtered data')
        plt.xlabel('Time [sec]')
        plt.grid()
        plt.legend()

        plt.subplots_adjust(hspace=0.35)
        plt.show()
    return y


def c_correlateX(a, v, returnx=False, returnav=False, s=0):
    '''

    :param a:
    :param v: a and v can be a dict in following format {'x':[],'y':[]}. The length of a and v can be different.
    :param returnx:
    :return:
    '''
    from scipy.interpolate import splev, splrep
    import numpy.ma as ma

    if isinstance(a, dict):
        max_a = np.nanmax(a['x'])
        min_a = np.nanmin(a['x'])
        max_v = np.nanmax(v['x'])
        min_v = np.nanmin(v['x'])
        max_ = min(max_a, max_v)
        min_ = max(min_a, min_v)
        if not max_ > min_:
            print('the x ranges of a and v have no overlap.')
            return None
        a_x = ma.masked_outside(a['x'].copy(), min_, max_)
        if isinstance(a['y'], np.ma.core.MaskedArray):
            a_y = ma.masked_array(a['y'].copy(), a['y'].mask | a_x.mask)
            a_x = ma.masked_array(a_x, a_y.mask)
        else:
            a_y = ma.masked_array(a['y'].copy(), a_x.mask)
        v_x = ma.masked_outside(v['x'].copy(), min_, max_)
        if isinstance(v['y'], np.ma.core.MaskedArray):
            v_y = ma.masked_array(v['y'].copy(), v['y'].mask | v_x.mask)
            v_x = ma.masked_array(v_x, v_y.mask)
        else:
            v_y = ma.masked_array(v['y'].copy(), v_x.mask)

        dx_a = np.abs(np.nanmean(np.diff(a_x)))
        dx_v = np.abs(np.nanmean(np.diff(v_x)))
        if dx_a >= dx_v:
            v_ = v_y.compressed()
            x_ = v_x.compressed()
            tck = splrep(a_x.compressed(), a_y.compressed(), s=s)
            ys = splev(x_, tck)
            a_ = ys

        elif dx_a < dx_v:
            a_ = a_y.compressed()
            x_ = a_x.compressed()
            tck = splrep(v_x.compressed(), v_y.compressed(), s=s)
            ys = splev(x_, tck)
            v_ = ys

    else:
        a_ = a.copy()
        v_ = v.copy()
        x_ = None
    a_ = (a_ - np.nanmean(a_)) / (np.nanstd(a_) * len(a_))
    v_ = (v_ - np.nanmean(v_)) / np.nanstd(v_)
    if returnx:
        if x_ is None:
            return np.arange(len(a_)) - np.floor(len(a_) / 2.0), np.correlate(a_, v_, mode='same')
        else:
            return (np.arange(len(a_)) - np.floor(len(a_) / 2.0)) * np.nanmean(np.diff(x_)), np.correlate(a_, v_,
                                                                                                          mode='same'), x_, a_, v_
    else:
        return np.correlate(a_, v_, mode='same')
