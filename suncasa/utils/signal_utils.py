from scipy.signal import butter, freqz, filtfilt
import matplotlib.pyplot as plt
import numpy as np


def normalize(y, ymax=None, ymin=None, center=None, yerr=None, symgamma=None):
    '''
    :param y:
    :param ymax:
    :param ymin:
    :param center: option ---- None, zero, 0, mean
    :param symgamma:
    :return:
    '''
    if ymax is None:
        ymax = np.nanmax(y)
    if ymin is None:
        ymin = np.nanmin(y)
    if center is None:
        if yerr is None:
            return (y - ymin) / (ymax - ymin)
        else:
            return [(y - ymin) / (ymax - ymin), yerr / (ymax - ymin)]
    else:
        if center == 'mean':
            ycenter = np.nanmedian(y)
        elif center in ['zero', '0', 0]:
            ycenter = 0.0
        if symgamma is not None:
            # ynew = (y - ycenter) / (ymax - ycenter)
            idx_pos = np.where(y > ycenter)
            idx_neg = np.where(y < ycenter)
            ynew = np.empty_like(y)
            ynew[:] = ycenter
            ynew[idx_pos] = (y[idx_pos]) ** symgamma + ycenter
            ynew[idx_neg] = -np.abs(y[idx_neg]) ** symgamma + ycenter
            ymax = np.nanmax(ynew)
            # ymin = np.nanmin(ynew)
            yout = (ynew - ycenter) / (ymax - ycenter) * 0.5 + 0.5
            if yerr is None:
                return yout
            else:
                return [yout, yerr ** symgamma / (ymax - ycenter) * 0.5]
        else:
            yout = (y - ycenter) / (ymax - ycenter) * 0.5 + 0.5
            if yerr is None:
                return yout
            else:
                return [yout, yerr / (ymax - ycenter) * 0.5]


def smooth(x, window_len=11, window='hanning', mode='same'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    import numpy
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = numpy.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode=mode)
    if mode == 'same':
        return y[np.int(window_len) - 1:-np.int(window_len) + 1]
    else:
        return y[np.int(window_len / 2 - 1):-np.int(window_len / 2)]


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
        plt.plot(1. / (0.5 * fs * w / np.pi), np.abs(h), 'b')
        plt.plot(1.0 / cutoff[0], 0.5 * np.sqrt(2), 'ko')
        plt.plot(1.0 / cutoff[1], 0.5 * np.sqrt(2), 'ko')
        plt.axvline(1.0 / cutoff[0], color='k')
        plt.axvline(1.0 / cutoff[1], color='k')
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


def c_correlateX(a, v, returnx=False, returnav=False, s=0, xran=None, coarse=False, interp='spl'):
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
        if xran is not None:
            max_ = min(max_, xran[1])
            min_ = max(min_, xran[0])
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

        if coarse:
            if dx_a < dx_v:
                icase = 0
            elif dx_a >= dx_v:
                icase = 1
        else:
            if dx_a >= dx_v:
                icase = 0
            elif dx_a < dx_v:
                icase = 1

        if icase == 0:
            v_ = v_y.compressed()
            x_ = v_x.compressed()
            if interp == 'spl':
                tck = splrep(a_x.compressed(), a_y.compressed(), s=s)
                ys = splev(x_, tck)
            else:
                ys = np.interp(x_, a_x.compressed(), a_y.compressed())
            a_ = ys
        else:
            a_ = a_y.compressed()
            x_ = a_x.compressed()
            if interp == 'spl':
                tck = splrep(v_x.compressed(), v_y.compressed(), s=s)
                ys = splev(x_, tck)
            else:
                ys = np.interp(x_, v_x.compressed(), v_y.compressed())
            v_ = ys
    else:
        a_ = a.copy()
        v_ = v.copy()
        x_ = None
    a_ = (a_ - np.nanmean(a_)) / (np.nanstd(a_) * len(a_))
    v_ = (v_ - np.nanmean(v_)) / np.nanstd(v_)
    print(a_.shape,v_.shape,x_.shape)
    if returnx:
        if x_ is None:
            return [np.arange(len(a_)) - np.floor(len(a_) / 2.0), np.correlate(a_, v_, mode='same')]
        else:
            return [(np.arange(len(a_)) - np.floor(len(a_) / 2.0)) * np.nanmean(np.diff(x_)), np.correlate(a_, v_,
                                                                                                           mode='same'),
                    x_, a_, v_]
    else:
        return np.correlate(a_, v_, mode='same')


def get_xcorr_info(xcorr, cwidth_guess=2.5 / 24 / 60, showplt=False, verbose=False):
    '''
    calculate the error in time lag using equation (3) in Gaskell & Peterson 1987
    :param xcorr:
    :param cwidth_guess:
    :return:
    '''
    import numpy.ma as ma
    from scipy.optimize import curve_fit
    def gauss_func(x, A, mu, sigma, c):
        return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + c

    xcorr_idxmax = np.nanargmax(xcorr[1])
    lag, peak = xcorr[0][xcorr_idxmax], xcorr[1][xcorr_idxmax]
    mxcorr_peakx = ma.masked_outside(xcorr[0], lag - cwidth_guess, lag + cwidth_guess)
    mxcorr_peaky = ma.masked_array(xcorr[1], mxcorr_peakx.mask)
    # mxcorr_residx = ma.masked_array(xcorr[0],~mxcorr_peakx.mask)
    mxcorr_residy = ma.masked_array(xcorr[1], ~mxcorr_peakx.mask)
    popt, pcov = curve_fit(gauss_func, mxcorr_peakx.compressed(), mxcorr_peaky.compressed(),
                           p0=[peak, lag, 2. / 24 / 60, 0])

    ##
    popt[2] = np.abs(popt[2])
    wc = np.sqrt(2 * np.log(2)) * popt[2]  ## half-width at half-maximum of the peak in the (xcorr function) CCF
    h = peak  ## height of the peak
    # rms = np.nanstd(mxcorr_residy)  ## noise level in the CCF
    rms = np.sqrt(np.mean(mxcorr_residy ** 2))  ## noise level in the CCF

    prof_peak = gauss_func(mxcorr_peakx, *popt)
    if showplt:
        plt.figure()
        plt.plot(xcorr[0], xcorr[1])
        plt.plot(mxcorr_peakx, mxcorr_peaky)
        plt.plot(mxcorr_peakx, prof_peak)
    err = 0.75 * wc / (1 + h / rms)
    if verbose:
        print({'lag': lag, 'peak': peak, 'lagerr': err, 'wc': wc, 'rms': rms, 'popt': popt})
    return {'lag': lag, 'peak': peak, 'lagerr': err}


def plot_wavelet(t, dat, dt, pl, pr, period_pltlim=None, ax=None, ax2=None, stscale=2, siglev=0.95, cmap='viridis',
                 title='', levels=None,
                 label='', units='', tunits='',
                 sav_img=False):
    import pycwt as wavelet
    from pycwt.helpers import find
    import numpy as np
    import matplotlib.pyplot as plt
    from copy import copy
    import numpy.ma as ma

    t_ = copy(t)
    t0 = t[0]
    # print(Time(t[-1:], format='plot_date').iso)
    # We also create a time array in years.
    N = dat.size
    t = np.arange(0, N) * dt + t0
    # print(Time(t[-1:], format='plot_date').iso)
    # We write the following code to detrend and normalize the input data by its
    # standard deviation. Sometimes detrending is not necessary and simply
    # removing the mean value is good enough. However, if your dataset has a well
    # defined trend, such as the Mauna Loa CO\ :sub:`2` dataset available in the
    # above mentioned website, it is strongly advised to perform detrending.
    # Here, we fit a one-degree polynomial function and then subtract it from the
    # original data.
    p = np.polyfit(t - t0, dat, 1)
    dat_notrend = dat - np.polyval(p, t - t0)
    std = dat_notrend.std()  # Standard deviation
    var = std ** 2  # Variance
    dat_norm = dat_notrend / std  # Normalized dataset

    # The next step is to define some parameters of our wavelet analysis. We
    # select the mother wavelet, in this case the Morlet wavelet with
    # :math:`\omega_0=6`.
    mother = wavelet.Morlet(6)
    s0 = stscale * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
    dj = 1 / 12  # Twelve sub-octaves per octaves
    J = -1  # 7 / dj  # Seven powers of two with dj sub-octaves
    alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise

    # The following routines perform the wavelet transform and inverse wavelet
    # transform using the parameters defined above. Since we have normalized our
    # input time-series, we multiply the inverse transform by the standard
    # deviation.
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J, mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

    # We calculate the normalized wavelet and Fourier power spectra, as well as
    # the Fourier equivalent periods for each wavelet scale.
    power = (np.abs(wave)) ** 2
    fft_power = np.abs(fft) ** 2
    period = 1 / freqs

    # We could stop at this point and plot our results. However we are also
    # interested in the power spectra significance test. The power is significant
    # where the ratio ``power / sig95 > 1``.
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha, significance_level=siglev, wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    # Then, we calculate the global wavelet spectrum and determine its
    # significance level.
    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha, significance_level=siglev, dof=dof,
                                            wavelet=mother)

    # We also calculate the scale average between 2 years and 8 years, and its
    # significance level.
    sel = find((period >= pl) & (period < pr))
    Cdelta = mother.cdelta
    scale_avg = (scales * np.ones((N, 1))).transpose()
    scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
    scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
    scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha, significance_level=siglev,
                                                 dof=[scales[sel[0]], scales[sel[-1]]],
                                                 wavelet=mother)

    # levels = [0.25, 0.5, 1, 2, 4, 8, 16,32]
    if levels is None:
        levels = np.linspace(0.0, 128., 256)
    # ax.contourf(t, np.log2(period), np.log2(power), np.log2(levels), extend='both', cmap=plt.cm.viridis)
    im = ax.contourf(t_, np.array(period) * 24 * 60, power, levels, extend='both', cmap=cmap, zorder=-20)
    # for pathcoll in im.collections:
    #     pathcoll.set_rasterized(True)
    ax.set_rasterization_zorder(-10)
    # im = ax.pcolormesh(t_, np.array(period) * 24 * 60, power,vmax=32.,vmin=0, cmap=cmap)
    # im = ax.contourf(t, np.array(period)*24*60, np.log2(power), np.log2(levels), extend='both', cmap=cmap)
    extent = [t_.min(), t_.max(), 0, max(period) * 24 * 60]
    # ax.contour(t, np.log2(period), sig95, [-99, 1], colors='k', linewidths=1, extent=extent)
    CS = ax.contour(t_, np.array(period) * 24 * 60, sig95 * siglev, [-99, 1.0 * siglev], colors='k', linewidths=1,
                    extent=extent)
    ax.clabel(CS, inline=1, fmt='%1.3f')
    ax.fill(np.concatenate([t_, t_[-1:] + dt, t_[-1:] + dt, t_[:1] - dt, t_[:1] - dt]),
            np.concatenate(
                [np.array(coi), [2 ** (1e-9)], np.array(period[-1:]), np.array(period[-1:]),
                 [2 ** (1e-9)]]) * 24 * 60,
            color='k', alpha=0.75, edgecolor='None', facecolor='k', hatch='x')
    # ### not Matplotlib does not display hatching when rendering to pdf. Here is a workaround.
    # ax.fill(np.concatenate([t_, t_[-1:] + dt, t_[-1:] + dt, t_[:1] - dt, t_[:1] - dt]),
    #         np.concatenate(
    #             [np.array(coi), [2 ** (1e-9)], np.array(period[-1:]), np.array(period[-1:]),
    #              [2 ** (1e-9)]]) * 24 * 60,
    #         color='None', alpha=1.0, edgecolor='k', hatch='x')
    # ax.set_title('b) {} Wavelet Power Spectrum ({})'.format(label, mother.name))
    #
    # ax.set_rasterization_zorder(20)
    # Yticks = np.arange(np.ceil(np.array(period.min()*24*60)), np.ceil(np.array(period.max()*24*60)))
    # ax.set_yticks(np.array(Yticks))
    # ax.set_yticklabels(Yticks)

    ax2.plot(glbl_signif, np.array(period) * 24 * 60, 'k--')
    # ax2.plot(var * fft_theor, np.array(period) * 24 * 60, '--', color='#cccccc')
    # ax2.plot(var * fft_power, np.array(1. / fftfreqs) * 24 * 60, '-', color='#cccccc',
    #          linewidth=1.)
    ax2.plot(var * glbl_power, np.array(period) * 24 * 60, 'k-', linewidth=1)
    mperiod = ma.masked_outside(np.array(period), period_pltlim[0], period_pltlim[1])
    mpower = ma.masked_array(var * glbl_power, mask=mperiod.mask)
    # ax2.set_title('c) Global Wavelet Spectrum')
    ax2.set_xlabel(r'Power'.format(units))
    ax2.set_xlim([0, mpower.compressed().max() + var])
    # print(glbl_power)
    # ax2.set_ylim(np.array([period.min(), period.max()]))
    # ax2.set_yticks(np.array(Yticks))
    # ax2.set_yticklabels(Yticks)
    plt.setp(ax2.get_yticklabels(), visible=False)

    if period_pltlim:
        ax.set_ylim(np.array(period_pltlim) * 24 * 60)
    else:
        ax.set_ylim(np.array([period.min(), period.max()]) * 24 * 60)

    return im
