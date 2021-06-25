import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import numpy as np
import os
import datetime
import jdutil
import pdb
import signalsmooth
import struct
from scipy.io.idl import readsav
from datetime import datetime
from taskinit import *


def get_dspec(vis=None, specfile=None, bl=None, spw=None, timeran=None):
    """
    Note: antennas specified in "bl" is antennas INDEX but not antenna NAME.
    REQUIRED INPUTS:
        vis: name of the measurement set file
        bl: baseline pair to use for the dynamic spectrum in CASA msselect format. e.g., bl='10&11'
    OPTIONAL INPUTS:
        specfile: name of the dynamic spectrum file
        spw: in CASA msselect format, e.g., '0~7:10~30', '0~3' (other formats not supported yet)
        timeran: in CASA msselect format, e.g., '18:00:00~18:30:00'
    EXAMPLE:
        import dspec
        dspec.get_dspec(mspath='./',msfile='yourcasadata.ms',spw='0~7',
                        bl='2&3',timeran='18:00:00~18:30:00')
    """
    # Open the ms and plot dynamic spectrum
    if not spw:
        spw = ''
    if not timeran:
        timeran = ''
    if not bl:
        bl = ''
    if not specfile:
        specfile = vis + '.spec.npz'
    print 'Splitting selected dynamic spectral data...'
    ms.open(vis, nomodify=True)
    vis_spl = './tmpms.splitted'
    ms.split(outputms=vis_spl, whichcol='DATA', time=timeran, spw=spw, baseline=bl)
    ms.close()
    ms.open(vis_spl, nomodify=False)
    print 'Regridding into a single spectral window...'
    ms.cvel(outframe='LSRK', mode='frequency', interp='nearest')
    ms.selectinit(datadescid=0, reset=True)
    specdata = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
    ms.close()
    os.system('rm -rf ' + vis_spl)
    npol = specdata['amplitude'].shape[0]
    nfreq = specdata['amplitude'].shape[1]
    nbl = specdata['amplitude'].shape[2]
    ntim = specdata['amplitude'].shape[3]
    print 'shape of the spectral data: ', specdata['amplitude'].shape
    print 'number of selected polarizations: ', npol
    print 'number of selected baselines: ', nbl
    print 'number of selected frequency channels: ', nfreq
    print 'number of selected time pixels: ', ntim
    spec = np.swapaxes(specdata['amplitude'], 2, 1)
    # CASA 4.3.1
    freq = specdata['axis_info']['freq_axis']['chan_freq'].reshape(nfreq)
    tim = specdata['time']
    # Save variables
    np.savez(specfile, spec=spec, tim=tim, freq=freq,
             timeran=timeran, bl=bl, spw=spw,
             npol=npol, nbl=nbl, nfreq=nfreq, ntim=ntim)
    print 'Dynamic spectrum saved as: ' + specfile


def get_med_dspec(mspath=None, msfile=None, specfile=None, \
                  smoothmode='base', tsmoothwidth=40, tsmoothtype='hanning', \
                  fsmoothwidth=4, fsmoothtype='hannning', \
                  fsmooth=1, timeran=None, chanran=None, dmin=0.5, dmax=2.0):
    # Open the ms and plot dynamic spectrum
    print 'retrieving spectral data...'
    ms.open(mspath + '/' + msfile)
    ms.selectinit(datadescid=0)
    if timeran and (type(timeran) == str):
        ms.msselect({'time': timeran})
    if chanran and (type(chanran) == str):
        [bchan, echan] = chanran.split('~')
        nchan = int(echan) - int(bchan) + 1
        ms.selectchannel(nchan, int(bchan), 1, 1)
    specdata = ms.getdata(['data', 'time', 'axis_info'], ifraxis=True)
    npol = specdata['data'].shape[0]
    nchan = specdata['data'].shape[1]
    nbl = specdata['data'].shape[2]
    ntime = specdata['data'].shape[3]
    print 'shape of the spectral data: ', specdata['data'].shape
    print 'number of time pixels: ', specdata['time'].shape[0]
    print 'number of frequency channels selected: ', specdata['axis_info']['freq_axis']['chan_freq'].shape[0]
    print 'number of baselines: ', specdata['axis_info']['ifr_axis']['baseline'].shape[0]
    # ratio over time-average
    if smoothmode == 'hipass':
        specsmooth = np.copy(specdata['data'])
        # time smooth
        for i in range(npol):
            for j in range(nchan):
                for k in range(nbl):
                    specsmooth[i, j, k, :] = signalsmooth.smooth(specdata['data'][i, j, k, :], tsmoothwidth,
                                                                 tsmoothtype)
        # frequency smooth
        for i in range(npol):
            for j in range(nbl):
                for k in range(ntime):
                    specsmooth[i, :, j, k] = signalsmooth.smooth(specsmooth[i, :, j, k], fsmoothwidth, fsmoothtype)
        specratio = specdata['data'] / specsmooth
    if smoothmode == 'base':
        specsmooth = np.mean(specdata['data'], axis=3)
        specratio = specdata['data'] / specsmooth[:, :, :, None]

    spec_med = np.median(specratio, axis=2)
    if fsmooth:
        for i in range(npol):
            for j in range(ntime):
                spec_med[i, :, j] = signalsmooth.smooth(spec_med[i, :, j], fsmoothwidth, fsmoothtype)

    tim = specdata['time']
    freq = specdata['axis_info']['freq_axis']['chan_freq'].reshape(nchan)
    tim0 = tim[0]
    tim0str = qa.time(qa.quantity(tim0, 's'), prec=8)[0]
    tim_ = tim - tim[0]
    freqghz = freq / 1e9
    f = plt.figure(figsize=(8, 8), dpi=100)
    ax1 = f.add_subplot(211)
    f.subplots_adjust(hspace=0.4)
    ax1.pcolormesh(tim_, freqghz, np.abs(spec_med[0, :, :]), cmap='jet', vmin=dmin, vmax=dmax)
    ax1.set_xlim([tim_[0], tim_[-1]])
    ax1.set_ylim([freqghz[0], freqghz[-1]])
    ax1.set_title('Median-Filtered Dynamic Spectrum (RCP)')
    ax1.set_xlabel('Time (seconds) since ' + tim0str)
    ax1.set_ylabel('Frequency (GHz)')
    ax2 = f.add_subplot(212)
    ax2.pcolormesh(tim_, freqghz, np.abs(spec_med[1, :, :]), cmap='jet', vmin=dmin, vmax=dmax)
    ax2.set_xlim([tim_[0], tim_[-1]])
    ax2.set_ylim([freqghz[0], freqghz[-1]])
    ax2.set_title('Median-Filtered Dynamic Spectrum (LCP)')
    ax2.set_xlabel('Time (seconds) since ' + tim0str)
    ax2.set_ylabel('Frequency (GHz)')
    if specfile:
        np.savez(mspath + '/' + specfile, spec_med=spec_med, tim=tim, freq=freq)


def plt_dspec(specfile=None, pol='I', dmin=None, dmax=None,
              timerange=None, freqrange=None, timestr=True,
              movie=False, framedur=60., dtframe=10.,
              goessav=None, goes_trange=None,
              savepng=True, savepdf=False):
    """
    timerange: format: ['2012/03/10/18:00:00','2012/03/10/19:00:00']
    freqrange: format: [1000.,1500.] in MHz
    movie: do a movie of dynamic spectrum?
    framedur: time range of each frame
    dtframe: time difference of consecutive frames
    goessav: provide an IDL save file from the sswidl GOES widget output
    goes_trange: plot only the specified time range for goes
    timestr: display time as strings on X-axis -- currently the times do not update themselves when zooming in
    """
    # Set up variables 
    if pol != 'RR' and pol != 'LL' and pol != 'RRLL' and pol != 'I' and pol != 'V' and pol != 'IV':
        print "Please enter 'RR', 'LL', 'RRLL', 'I', 'V', 'IV' for pol"
        return 0

    specdata = np.load(specfile)
    spec = specdata['spec']
    tim = specdata['tim']
    freq = specdata['freq']
    bl = specdata['bl'].item()
    #npol = specdata['npol']
    #nbl = specdata['nbl']
    (npol, nbl, nfreq, ntim)=spec.shape

    if timerange:
        if type(timerange[0]) is str:
            timerange = [qa.convert(qa.quantity(t), 's')['value'] for t in timerange]
        tidx = np.where((tim >= timerange[0]) & (tim <= timerange[1]))[0]
    else:
        tidx = range(ntim)
    if freqrange:
        fidx = np.where((freq >= freqrange[0] * 1e6) & (freq <= freqrange[1] * 1e6))[0]
    else:
        fidx = range(nfreq)

    # setup plot parameters
    print 'ploting dynamic spectrum...'
    spec_med = np.median(np.absolute(spec))
    if not dmin:
        dmin = spec_med / 20.
    if not dmax:
        dmax = spec_med * 5.
    # do the plot
    for b in range(nbl):
        if pol != 'RRLL' and pol != 'IV':
            if pol == 'RR':
                spec_plt = spec[0, b, :, :]
            elif pol == 'LL':
                spec_plt = spec[1, b, :, :]
            elif pol == 'I':
                spec_plt = (spec[0, b, :, :] + spec[1, b, :, :]) / 2.
            elif pol == 'V':
                spec_plt = (spec[0, b, :, :] - spec[1, b, :, :]) / 2.
            if movie:
                f = plt.figure(figsize=(16, 8), dpi=100)
                if goessav:
                    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
                    gs.update(left=0.06, right=0.97, top=0.95, bottom=0.06)
                    ax1 = f.add_subplot(gs[0])
                    ax2 = f.add_subplot(gs[1])
                    if os.path.exists(goessav):
                        goes = readsav(goessav)
                        # IDL anytim 0 sec correspond to 1979 Jan 01, convert to mjd time
                        anytimbase = qa.convert(qa.quantity('1979/01/01/00:00:00'), 's')['value']
                        mjdbase = goes['utbase'] + anytimbase
                        ts = goes['tarray'] + mjdbase
                        lc0 = goes['yclean'][0, :]
                        lc1 = goes['yclean'][1, :]
                else:
                    ax1 = f.add_subplot(211)
                tstart = tim[tidx[0]]
                tend = tim[tidx[-1]]
                tstartstr = qa.time(qa.quantity(tstart, 's'))[0]
                tendstr = qa.time(qa.quantity(tend, 's'))[0]
                nfrm = int((tend - tstart) / dtframe) + 1
                print 'Movie mode set. ' + str(nfrm) + ' frames to plot from ' + tstartstr + ' to ' + tendstr
                for i in range(nfrm):
                    if (i != 0) and (i % 10 == 0):
                        print str(i) + ' frames done'
                    timeran = [tstart + i * dtframe, tstart + i * dtframe + framedur]
                    tidx1 = np.where((tim >= timeran[0]) & (tim <= timeran[1]))[0]
                    tim1 = tim[tidx1]
                    freq1 = freq[fidx] / 1e9
                    # the following is wrong
                    spec_plt1 = spec_plt[fidx, :][:, tidx1]
                    ax1.pcolormesh(tim1, freq1, spec_plt1, cmap='jet', vmin=dmin, vmax=dmax)
                    ax1.set_xlim(tim1[0], tim1[-1])
                    ax1.set_ylim(freq1[0], freq1[-1])
                    ax1.set_ylabel('Frequency (GHz)')
                    ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                    if timestr:
                        labels = ax1.get_xticks().tolist()
                        newlabels = [qa.time(qa.quantity(lb, 's'))[0] for lb in labels]
                        ax1.set_xticklabels(newlabels)
                    ax1.set_autoscale_on(False)
                    if goessav:
                        if goes_trange:
                            if type(goes_trange[0]) is str:
                                goes_trange = [qa.convert(qa.quantity(t), 's')['value'] for t in goes_trange]
                                idx = np.where((ts >= goes_trange[0]) & (ts <= goes_trange[1]))[0]
                        else:
                            idx = range(len(ts))
                        ts_plt = ts[idx]
                        lc0_plt = lc0[idx]
                        ts_plt_dt = [datetime.strptime(qa.time(qa.quantity(t, 's'), form='ymd')[0], '%Y/%m/%d/%H:%M:%S')
                                     for t in ts_plt]
                        utbase = qa.convert(qa.quantity('0001/01/01/00:00:00'), 'd')['value'] + 1
                        ts_plt_d = ts_plt / 3600. / 24. - utbase
                        # ts_plt_d = dates.date2num(ts_plt_dt)
                        ax2.plot_date(ts_plt_d, lc0_plt, 'b-')
                        ax2.axvspan(tim1[0] / 3600. / 24. - utbase, tim1[-1] / 3600. / 24. - utbase, color='red',
                                    alpha=0.5)
                        ax2.set_yscale('log')
                        # ax2.set_xlim(goes_trange)
                        # labels=ax2.get_xticks().tolist()
                        # newlabels=[qa.time(qa.quantity(lb,'s'))[0] for lb in labels]
                        # ax2.set_xticklabels(newlabels)
                        ax2.set_title('GOES 1-8 A')

                    tstartstr_ = qa.time(qa.quantity(tim1[0], 's'))[0]
                    tendstr_ = qa.time(qa.quantity(tim1[-1], 's'))[0]
                    timstr = tstartstr_.replace(':', '') + '-' + tendstr_.replace(':', '')
                    figfile = 'dspec_t' + timstr + '.png'
                    if not os.path.isdir('dspec'):
                        os.makedirs('dspec')
                    f.savefig('dspec/' + figfile)
                    plt.cla()
            else:
                f = plt.figure(figsize=(8, 4), dpi=100)
                ax = f.add_subplot(111)
                freqghz = freq / 1e9
                ax.pcolormesh(tim, freqghz, spec_plt, cmap='jet', vmin=dmin, vmax=dmax)
                ax.set_xlim(tim[tidx[0]], tim[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

                def format_coord(x, y):
                    col = np.argmin(np.absolute(tim - x))
                    row = np.argmin(np.absolute(freqghz - y))
                    if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                        timv = tim[col]
                        timstr = qa.time(qa.quantity(timv, 's'), form='clean', prec=9)[0]
                        flux = spec_plt[row, col]
                        return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy/beam'.format(col, timstr, y, flux)
                    else:
                        return 'x = {0}, y = {1:.3f}'.format(x, y)

                ax.format_coord = format_coord
                ax.set_ylabel('Frequency (GHz)')
                ax.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                if timestr:
                    labels = ax.get_xticks().tolist()
                    newlabels = [qa.time(qa.quantity(lb, 's'))[0] for lb in labels]
                    ax.set_xticklabels(newlabels)
                ax.set_autoscale_on(False)

        else:
            f = plt.figure(figsize=(8, 6), dpi=100)
            R_plot = np.absolute(spec[0, b, :, :])
            L_plot = np.absolute(spec[1, b, :, :])
            I_plot = (R_plot + L_plot) / 2.
            V_plot = (R_plot - L_plot) / 2.
            if pol == 'RRLL':
                spec_plt_1 = R_plot
                spec_plt_2 = L_plot
                polstr = ['RR', 'LL']
            if pol == 'IV':
                spec_plt_1 = I_plot
                spec_plt_2 = V_plot
                polstr = ['I', 'V']

            ax1 = f.add_subplot(211)
            freqghz = freq / 1e9
            ax1.pcolormesh(tim, freqghz, spec_plt_1, cmap='jet', vmin=dmin, vmax=dmax)
            ax1.set_xlim(tim[tidx[0]], tim[tidx[-1]])
            ax1.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timv = tim[col]
                    timstr = qa.time(qa.quantity(timv, 's'), form='clean', prec=9)[0]
                    flux = spec_plt[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy/beam'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax1.format_coord = format_coord
            ax1.set_ylabel('Frequency (GHz)')
            if timestr:
                labels = ax1.get_xticks().tolist()
                newlabels = [qa.time(qa.quantity(lb, 's'))[0] for lb in labels]
                ax1.set_xticklabels(newlabels)
            ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[0])
            ax1.set_autoscale_on(False)
            ax2 = f.add_subplot(212)
            ax2.pcolormesh(tim, freqghz, spec_plt_2, cmap='jet', vmin=dmin, vmax=dmax)
            ax2.set_xlim(tim[tidx[0]], tim[tidx[-1]])
            ax2.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
            if timestr:
                labels = ax2.get_xticks().tolist()
                newlabels = [qa.time(qa.quantity(lb, 's'))[0] for lb in labels]
                ax2.set_xticklabels(newlabels)

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timv = tim[col]
                    timstr = qa.time(qa.quantity(timv, 's'), form='clean', prec=9)[0]
                    flux = spec_plt[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy/beam'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax2.format_coord = format_coord
            ax2.set_ylabel('Frequency (GHz)')
            ax2.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[1])
            ax2.set_autoscale_on(False)
        if savepng:
            figfile = specfile[:(specfile.find('npz'))] + 'png'
            f.savefig(figfile)
        if savepdf:
            figfile = specfile[:(specfile.find('npz'))] + 'pdf'
            f.savefig(figfile)
            # plt.close()
            # return (f,ax)


def wrt_dspec(specfile=None, specdat=None):
    try:
        specfile
    except NameError:
        print 'No input centfile specified for reading. Aborting...'
    if not specdat:
        print 'Output file name is not specified, use the default convention'
        specdat = specfile.replace('npz', 'dat')
    specdata = np.load(specfile)
    spec = np.copy(specdata['spec'][:, :, :, :])
    npl, nbl, nf, nt = spec.shape
    print 'Dimension of the data cube -- # of pol, # of baseline, # of frequency, # of time:'
    print npl, nbl, nf, nt
    nele = npl * nbl * nf * nt
    # need to transpose to nf, nt, npl, nbl
    spec = np.transpose(spec, (2, 3, 0, 1))
    # Turn fmx into a 1d array for writing
    spec = spec.flatten()
    # Pack the data into a byte buffer.  Note the asterisks in the next three lines!
    # If your data includes integers or other data types, use 'i' or other instead of 'f'
    buf = struct.pack(str(nf) + 'f', *specdata['freq'])
    buf += struct.pack(str(nt) + 'd', *specdata['tim'])
    buf += struct.pack(str(nele) + 'f', *spec)
    with open(specdat, 'wb') as f:
        f.write(buf)
    f.close()


def plt_dspec_old(mspath=None, specfile=None, pol='I', dmin=None, dmax=None, fig=None):
    # Set up variables 
    if pol != 'RR' and pol != 'LL' and pol != 'RRLL' and pol != 'I' and pol != 'V' and pol != 'IV':
        print "Please enter 'RR', 'LL', 'RRLL', 'I', 'V', 'IV' for pol"
        return 0

    if not mspath:
        mspath = '.'
    specdata = np.load(mspath + '/' + specfile)
    spec = specdata['spec']
    tim = specdata['tim']
    freq = specdata['freq']

    # setup plot parameters
    print 'ploting dynamic spectrum...'
    # mask the channels from 512 to 519 (no observation)
    spec = np.ma.array(spec)
    spec[:, 512:519, :] = np.ma.masked
    spec_med = np.median(np.absolute(spec))
    # set the time axis
    ntim = tim.shape
    if ntim[0] < 20:
        xticks = np.arange((ntim[0] - 1) / 2 + 1) * 2
    elif ntim[0] < 100:
        xticks = np.arange((ntim[0] - 1) / 20 + 1) * 20
    elif ntim[0] < 600:
        xticks = np.arange((ntim[0] - 1) / 100 + 1) * 100
    elif ntim[0] < 2400:
        xticks = np.arange((ntim[0] - 1) / 400 + 1) * 400
    elif ntim[0] < 12000:
        # 1 min per step
        tstart = np.fix(tim[0] / 60.) + 1
        xstart = np.abs(tim - tstart * 60.).argmin()
        xticks = np.arange(ntim[0] / 1200) * 1200 + xstart
    elif ntim[0] > 12000:
        xticks = np.arange(ntim[0] / 6000 + 1) * 6000
    nticks = xticks.shape
    xticktims = []
    for i in range(nticks[0]):
        # xticktim0=qa.time(qa.quantity(tim[xticks[i]],'s'))
        tims = tim[xticks[i]]  # in seconds
        tims_jd = jdutil.mjd_to_jd(tims / 3600. / 24.)  # to julian date
        tims_dt = jdutil.jd_to_datetime(tims_jd)
        tims_dt2 = tims_dt + datetime.timedelta(seconds=round(tims_dt.microsecond / 1e6))
        tims_char = tims_dt2.strftime('%H:%M:%S')
        xticktims.append(tims_char)
    xticks = list(xticks)
    # do the plot
    f = plt.figure(figsize=(10, 6), dpi=100)
    if not dmin:
        dmin = spec_med / 20.
    if not dmax:
        dmax = spec_med * 5.
    if pol != 'RRLL' and pol != 'IV':
        ax = f.add_subplot(111)
        if pol == 'RR':
            spec_plt = np.absolute(spec[0, :, :])
        elif pol == 'LL':
            spec_plt = np.absolute(spec[1, :, :])
        elif pol == 'I':
            spec_plt = (np.absolute(spec[0, :, :]) + np.absolute(spec[1, :, :])) / 2.
        elif pol == 'V':
            spec_plt = (np.absolute(spec[0, :, :]) - np.absolute(spec[1, :, :])) / 2.
        # ax.imshow(spec_plt,aspect='auto',origin='lower',extent=[0,ntim[0]-1,np.min(freq)/1e6,np.max(freq)/1e6],\
        #          vmin=dmin,vmax=dmax,interpolation='none')
        # ax.set_xticks(xticks)
        # ax.set_xticklabels(xticktims)
        tim_ = tim - tim[0]
        freqghz = freq / 1e9
        ax.pcolormesh(tim_, freqghz, spec_plt, cmap='jet', vmin=dmin, vmax=dmax)
        ax.set_xlim([tim_[0], tim_[-1]])
        ax.set_ylim([freqghz[0], freqghz[-1]])
        ax.set_xlabel('Time (pixel)')
        ax.set_ylabel('Frequency (GHz)')
        ax.set_title('VLA dynamic spectrum for pol ' + pol)
        ax.set_autoscale_on(False)
    else:
        R_plot = np.absolute(spec[0, :, :])
        L_plot = np.absolute(spec[1, :, :])
        I_plot = (Rplot + Lplot) / 2.
        V_plot = (Rplot - Lplot) / 2.
        if pol == 'RRLL':
            spec_plt_1 = R_plot
            spec_plt_2 = L_plot
        if pol == 'IV':
            spec_plt_1 = I_plot
            spec_plt_2 = V_plot

        ax1 = f.add_subplot(211)
        ax1.imshow(spec_plt_1, aspect='auto', origin='lower',
                   extent=[0, ntim[0] - 1, np.min(freq) / 1e6, np.max(freq) / 1e6], \
                   vmin=dmin, vmax=dmax, interpolation='none')
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticktims)
        ax1.set_xlabel('Universal Time (50 ms/pixel)')
        ax1.set_ylabel('Frequency (MHz)')
        ax2.set_title('VLA dynm spec for pol ' + pol[0])
        ax1.set_autoscale_on(False)
        ax2 = f.add_subplot(212)
        ax2.imshow(spec_plt_2, aspect='auto', origin='lower',
                   extent=[0, ntim[0] - 1, np.min(freq) / 1e6, np.max(freq) / 1e6], \
                   vmin=dmin, vmax=dmax, interpolation='none')
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticktims)
        ax2.set_xlabel('Universal Time (50 ms/pixel)')
        ax2.set_ylabel('Frequency (MHz)')
        ax2.set_title('VLA dynm spec for pol ' + pol[-1])
        ax2.set_autoscale_on(False)
    if fig:
        figfile = specfile[:(specfile.find('spec'))] + pol + '.pdf'
        f.savefig(mspath + '/' + figfile)
    # plt.close()
    return (f, ax)
