import matplotlib.gridspec as gridspec
import numpy as np
import os
import datetime
import struct
from scipy.io.idl import readsav
from datetime import datetime
import sys

pversion = sys.version_info.major
if pversion<3:
    ## CASA version < 6
    from taskinit import ms, tb, qa
else:
    ## CASA version >= 6
    from casatools import table as tbtool
    from casatools import ms as mstool
    from casatools import quanta as qatool

    tb = tbtool()
    ms = mstool()
    qa = qatool()

import matplotlib.dates as mdates
from matplotlib.dates import date2num, AutoDateFormatter, AutoDateLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

stokestype = [
    'Undefined',
    # undefined value = 0
    'I',
    'Q',
    'U',
    # standard stokes parameters
    'V',
    #
    'RR',
    'RL',
    'LR',
    # circular correlation products
    'LL',
    #
    'XX',
    'XY',
    'YX',
    # linear correlation products
    'YY',
    #
    'RX',
    'RY',
    'LX',
    'LY',
    'XR',
    'XL',
    'YR',
    # mixed correlation products
    'YL',
    #
    'PP',
    'PQ',
    'QP',
    # general quasi-orthogonal correlation products
    'QQ',
    #
    'RCircular',
    'LCircular',
    # single dish polarization types
    'Linear',
    # Polarized intensity ((Q^2+U^2+V^2)^(1/2))
    'Ptotal',
    # Linearly Polarized intensity ((Q^2+U^2)^(1/2))
    'Plinear',
    # Polarization Fraction (Ptotal/I)
    'PFtotal',
    # Linear Polarization Fraction (Plinear/I)
    'PFlinear',
    # Linear Polarization Angle (0.5 arctan(U/Q)) (in radians)
    'Pangle']

stokesenum = {}
for k, v in zip(range(len(stokestype)), stokestype):
    stokesenum[k] = v


def get_dspec(vis=None, savespec=True, specfile=None, bl='', uvrange='', field='', scan='', datacolumn='data',
              domedian=False, timeran=None, spw=None, timebin='0s', regridfreq=False, fillnan=None, verbose=False):
    # from split_cli import split_cli as split
    if vis.endswith('/'):
        vis = vis[:-1]
    msfile = vis
    if not spw:
        spw = ''
    if not timeran:
        timeran = ''
    if not bl:
        bl = ''
    if domedian:
        if not uvrange:
            uvrange = '0.2~0.8km'
        # bl = ''
    else:
        uvrange = ''
    # Open the ms and plot dynamic spectrum
    if verbose:
        print('Splitting selected data...')
    vis_spl = './tmpms.splitted'
    if os.path.exists(vis_spl):
        os.system('rm -rf ' + vis_spl)

    # split(vis=msfile, outputvis=vis_spl, timerange=timeran, antenna=bl, field=field, scan=scan, spw=spw,
    #       uvrange=uvrange, timebin=timebin, datacolumn=datacolumn)

    try:
        from split_cli import split_cli as split
        split(vis=msfile, outputvis=vis_spl, datacolumn=datacolumn, timerange=timeran, spw=spw, antenna=bl, field=field,
              scan=scan, uvrange=uvrange, timebin=timebin)
    except:
        ms.open(msfile, nomodify=True)
        ms.split(outputms=vis_spl, whichcol=datacolumn, time=timeran, spw=spw, baseline=bl, field=field, scan=scan,
                 uvrange=uvrange, timebin=timebin)
        ms.close()

    if verbose:
        print('Regridding into a single spectral window...')
        # print('Reading data spw by spw')

    try:
        tb.open(vis_spl + '/POLARIZATION')
        corrtype = tb.getcell('CORR_TYPE', 0)
        pols = [stokesenum[p] for p in corrtype]
        tb.close()
    except:
        pols = []

    if regridfreq:
        ms.open(vis_spl, nomodify=False)
        ms.cvel(outframe='LSRK', mode='frequency', interp='nearest')
        ms.selectinit(datadescid=0, reset=True)
        data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
        specamp = data['amplitude']
        freq = data['axis_info']['freq_axis']['chan_freq']

    else:
        ms.open(vis_spl)
        ms.selectinit(datadescid=0, reset=True)
        spwinfo = ms.getspectralwindowinfo()
        specamp = []
        freq = []
        time = []
        for descid in range(len(spwinfo.keys())):
            ms.selectinit(datadescid=0, reset=True)
            ms.selectinit(datadescid=descid)
            data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
            specamp_ = data['amplitude']
            freq_ = data['axis_info']['freq_axis']['chan_freq']
            time_ = data['time']
            if fillnan is not None:
                flag_ = ms.getdata(['flag', 'time', 'axis_info'], ifraxis=True)['flag']
                if type(fillnan) in [int, float, long]:
                    specamp_[flag_] = float(fillnan)
                else:
                    specamp_[flag_] = 0.0
            specamp.append(specamp_)
            freq.append(freq_)
            time.append(time_)
        specamp = np.concatenate(specamp, axis=1)
        freq = np.concatenate(freq, axis=0)
        ms.selectinit(datadescid=0, reset=True)
    ms.close()
    os.system('rm -rf ' + vis_spl)
    (npol, nfreq, nbl, ntim) = specamp.shape
    freq = freq.reshape(nfreq)

    if verbose:
        print('npol, nfreq, nbl, ntime:', (npol, nfreq, nbl, ntim))
    spec = np.swapaxes(specamp, 2, 1)

    tim = data['time']

    if domedian:
        if verbose:
            print('doing median of all the baselines')
        # mask zero values before median
        spec_masked = np.ma.masked_where(spec < 1e-9, spec)
        spec_masked2 = np.ma.masked_invalid(spec)
        spec_masked = np.ma.masked_array(spec,mask=np.logical_or(spec_masked.mask, spec_masked2.mask))
        spec_med = np.ma.filled(np.ma.median(spec_masked, axis=1), fill_value=0.)
        nbl = 1
        ospec = spec_med.reshape((npol, nbl, nfreq, ntim))
    else:
        ospec = spec
    # Save the dynamic spectral data
    if savespec:
        if not specfile:
            specfile = msfile + '.dspec.npz'
        if os.path.exists(specfile):
            os.system('rm -rf ' + specfile)
        np.savez(specfile, spec=ospec, tim=tim, freq=freq,
                 timeran=timeran, spw=spw, bl=bl, uvrange=uvrange, pol=pols)
        if verbose:
            print('Median dynamic spectrum saved as: ' + specfile)

    return {'spec': ospec, 'tim': tim, 'freq': freq, 'timeran': timeran, 'spw': spw, 'bl': bl, 'uvrange': uvrange,
            'pol': pols}


def plt_dspec(specdata, pol='I', dmin=None, dmax=None, norm=None,
              timerange=None, freqrange=None, timestr=True,
              movie=False, framedur=60., dtframe=10.,
              goessav=None, goes_trange=None, cmap='jet',
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
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    import numpy
    from numpy import log10
    from astropy.time import Time
    if pol not in ['RR', 'LL', 'RRLL', 'XX', 'YY', 'XXYY', 'I', 'V', 'IV']:
        print("Please enter 'RR', 'LL', 'RRLL','XX', 'YY', 'XXYY', 'I', 'V', 'IV' for pol")
        return 0

    if type(specdata) is str:
        try:
            specdata = np.load(specdata)
        except:
            specdata = np.load(specdata, encoding='latin1')
        bl = specdata['bl'].item()
    try:
        (npol, nbl, nfreq, ntim) = specdata['spec'].shape
        spec = specdata['spec'] / 1.e4
        tim = specdata['tim']
        tim_ = Time(tim / 3600. / 24., format='mjd')
        tim_plt = tim_.plot_date
        freq = specdata['freq']
        if not 'bl' in vars():
            bl = specdata['bl']
        import matplotlib
        from pkg_resources import parse_version
        if parse_version(matplotlib.__version__) < parse_version('1.5.0'):
            spec[np.isnan(spec)] = 0.0

    except:
        print('format of specdata not recognized. Check your input')
        return -1

    if 'pol' in specdata.keys():
        pols = specdata['pol']
    else:
        pols = []
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
    print('ploting dynamic spectrum...')
    spec_med = np.median(np.absolute(spec))
    # if not dmin:
    #    dmin = spec_med / 20.
    # if not dmax:
    #    dmax = spec_med * 5.
    # do the plot
    if norm is None:
        norm = colors.Normalize(vmax=dmax, vmin=dmin)
    for b in range(nbl):
        if pol not in ['RRLL', 'IV', 'XXYY']:
            if pol in ['RR', 'XX']:
                spec_plt = spec[0, b, :, :]
            elif pol in ['LL', 'YY']:
                spec_plt = spec[1, b, :, :]
            elif pol == 'I':
                if ('XX' in pols) or ('YY' in pols):
                    spec_plt = spec[0, b, :, :] + spec[1, b, :, :]
                else:
                    spec_plt = (spec[0, b, :, :] + spec[1, b, :, :]) / 2.
            elif pol == 'V':
                if ('XX' in pols) or ('YY' in pols):
                    spec_plt = spec[0, b, :, :] - spec[1, b, :, :]
                else:
                    spec_plt = (spec[0, b, :, :] - spec[1, b, :, :]) / 2.
            if movie:
                fig = plt.figure(figsize=(16, 8), dpi=100)
                if goessav:
                    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
                    gs.update(left=0.06, right=0.97, top=0.95, bottom=0.06)
                    ax1 = fig.add_subplot(gs[0])
                    ax2 = fig.add_subplot(gs[1])
                    if os.path.exists(goessav):
                        goes = readsav(goessav)
                        # IDL anytim 0 sec correspond to 1979 Jan 01, convert to mjd time
                        anytimbase = qa.convert(qa.quantity('1979/01/01/00:00:00'), 's')['value']
                        mjdbase = goes['utbase'] + anytimbase
                        ts = goes['tarray'] + mjdbase
                        lc0 = goes['yclean'][0, :]
                        lc1 = goes['yclean'][1, :]
                else:
                    ax1 = fig.add_subplot(211)
                tstart = tim[tidx[0]]
                tend = tim[tidx[-1]]
                tstartstr = qa.time(qa.quantity(tstart, 's'))[0]
                tendstr = qa.time(qa.quantity(tend, 's'))[0]
                nfrm = int((tend - tstart) / dtframe) + 1
                print('Movie mode set. ' + str(nfrm) + ' frames to plot from ' + tstartstr + ' to ' + tendstr)
                for i in range(nfrm):
                    if (i != 0) and (i % 10 == 0):
                        print(str(i) + ' frames done')
                    timeran = [tstart + i * dtframe, tstart + i * dtframe + framedur]
                    tidx1 = np.where((tim >= timeran[0]) & (tim <= timeran[1]))[0]
                    tim1 = tim_[tidx1]
                    freq1 = freq[fidx] / 1e9
                    spec_plt1 = spec_plt[fidx, :][:, tidx1]
                    ax1.pcolormesh(tim1.plot_date, freq1, spec_plt1, cmap=cmap, norm=norm)
                    ax1.set_xlim(tim1[0].plot_date, tim1[-1].plot_date)
                    ax1.set_ylim(freq1[0], freq1[-1])
                    ax1.set_ylabel('Frequency (GHz)')
                    ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                    if timestr:
                        # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                        # ax1.xaxis_date()
                        # ax1.xaxis.set_major_formatter(date_format)
                        locator = AutoDateLocator(minticks=2)
                        ax1.xaxis.set_major_locator(locator)
                        ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
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
                        utbase = qa.convert(qa.quantity('0001/01/01/00:00:00'), 'd')['value'] + 1
                        ts_plt_d = ts_plt / 3600. / 24. - utbase
                        ax2.plot_date(ts_plt_d, lc0_plt, 'b-')
                        ax2.axvspan(tim1[0].mjd - utbase, tim1[-1].mjd - utbase, color='red',
                                    alpha=0.5)
                        ax2.set_yscale('log')
                        ax2.set_title('GOES 1-8 A')

                    tstartstr_ = tim1[0].datetime.strftime('%Y-%m-%dT%H%M%S.%f')[:-3]
                    tendstr_ = tim1[1].datetime.strftime('%H%M%S.%f')[:-3]
                    timstr = tstartstr_ + '-' + tendstr_
                    figfile = 'dspec_t' + timstr + '.png'
                    if not os.path.isdir('dspec'):
                        os.makedirs('dspec')
                    fig.savefig('dspec/' + figfile)
                    plt.cla()
            else:
                fig = plt.figure(figsize=(8, 4), dpi=100)
                ax = fig.add_subplot(111)
                freqghz = freq / 1e9
                im = ax.pcolormesh(tim_plt, freqghz, spec_plt, cmap=cmap, norm=norm)
                divider = make_axes_locatable(ax)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax, cax=cax_spec)
                clb_spec.set_label('Flux [sfu]')
                ax.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
                try:
                    from sunpy import lightcurve
                    from sunpy.time import TimeRange, parse_time
                    t1 = tim_[tidx[0]]
                    t2 = tim_[tidx[-1]]
                    tr = TimeRange(t1.iso, t2.iso)
                    goes = lightcurve.GOESLightCurve.create(tr)
                    goes.data['xrsb'] = 2 * (np.log10(goes.data['xrsb'])) + 26
                    xx = [str(ll) for ll in np.array(goes.data.index)]
                    yy = np.array(goes.data['xrsb'])
                    ax.plot(Time(xx).mjd * 24 * 3600, yy, c='yellow')
                    rightaxis_label_time = Time(xx[-1]).mjd * 24 * 3600
                    ax.text(rightaxis_label_time, 9.6, 'A', fontsize='15')
                    ax.text(rightaxis_label_time, 11.6, 'B', fontsize='15')
                    ax.text(rightaxis_label_time, 13.6, 'C', fontsize='15')
                    ax.text(rightaxis_label_time, 15.6, 'M', fontsize='15')
                    ax.text(rightaxis_label_time, 17.6, 'X', fontsize='15')
                except:
                    pass

                def format_coord(x, y):
                    col = np.argmin(np.absolute(tim_plt - x))
                    row = np.argmin(np.absolute(freqghz - y))
                    if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                        timstr = tim_[col].isot
                        flux = spec_plt[row, col]
                        return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} sfu'.format(col, timstr, y, flux)
                    else:
                        return 'x = {0}, y = {1:.3f}'.format(x, y)

                ax.format_coord = format_coord
                ax.set_ylabel('Frequency (GHz)')
                if bl:
                    ax.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                else:
                    ax.set_title('Median dynamic spectrum')
                if timestr:
                    # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                    # ax.xaxis_date()
                    # ax.xaxis.set_major_formatter(date_format)
                    locator = AutoDateLocator(minticks=2)
                    ax.xaxis.set_major_locator(locator)
                    ax.xaxis.set_major_formatter(AutoDateFormatter(locator))
                ax.set_autoscale_on(False)

        else:
            fig = plt.figure(figsize=(8, 6), dpi=100)
            R_plot = np.absolute(spec[0, b, :, :])
            L_plot = np.absolute(spec[1, b, :, :])
            I_plot = (R_plot + L_plot) / 2.
            V_plot = (R_plot - L_plot) / 2.
            if pol in ['RRLL', 'XXYY']:
                spec_plt_1 = R_plot
                spec_plt_2 = L_plot
                polstr = [pol[:2], pol[2:]]
            if pol == 'IV':
                spec_plt_1 = I_plot
                spec_plt_2 = V_plot
                polstr = ['I', 'V']

            ax1 = fig.add_subplot(211)
            freqghz = freq / 1e9
            im = ax1.pcolormesh(tim_plt, freqghz, spec_plt_1, cmap=cmap, norm=norm)
            divider = make_axes_locatable(ax1)
            cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
            clb_spec = plt.colorbar(im, ax=ax1, cax=cax_spec)
            clb_spec.set_label('Flux [sfu]')

            ax1.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
            ax1.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim_plt - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timstr = tim_[col].isot
                    flux = spec_plt_1[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} sfu'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax1.format_coord = format_coord
            ax1.set_ylabel('Frequency (GHz)')
            if timestr:
                # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                # ax1.xaxis_date()
                # ax1.xaxis.set_major_formatter(date_format)
                locator = AutoDateLocator(minticks=2)
                ax1.xaxis.set_major_locator(locator)
                ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
            ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[0])
            ax1.set_autoscale_on(False)
            ax2 = fig.add_subplot(212)
            im = ax2.pcolormesh(tim_plt, freqghz, spec_plt_2, cmap=cmap, norm=norm)
            divider = make_axes_locatable(ax2)
            cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
            clb_spec = plt.colorbar(im, ax=ax2, cax=cax_spec)
            clb_spec.set_label('Flux [sfu]')
            ax2.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
            ax2.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
            if timestr:
                # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                # ax2.xaxis_date()
                # ax2.xaxis.set_major_formatter(date_format)
                locator = AutoDateLocator(minticks=2)
                ax2.xaxis.set_major_locator(locator)
                ax2.xaxis.set_major_formatter(AutoDateFormatter(locator))

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim_plt - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timstr = tim_[col].isot
                    flux = spec_plt_2[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} sfu'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax2.format_coord = format_coord
            ax2.set_ylabel('Frequency (GHz)')
            ax2.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[1])
            ax2.set_autoscale_on(False)

        fig.tight_layout()


def wrt_dspec(specfile=None, specdat=None):
    try:
        specfile
    except NameError:
        print('No input centfile specified for reading. Aborting...')
    if not specdat:
        print('Output file name is not specified, use the default convention')
        specdat = specfile.replace('npz', 'dat')
    specdata = np.load(specfile)
    spec = np.copy(specdata['spec'][:, :, :, :])
    npl, nbl, nf, nt = spec.shape
    print('Dimension of the data cube -- # of pol, # of baseline, # of frequency, # of time:')
    print(npl, nbl, nf, nt)
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


def common_elements(a, b):
    a_set = set(a)
    b_set = set(b)

    # check length
    if len(a_set.intersection(b_set)) > 0:
        return (a_set.intersection(b_set))
    else:
        return ("no common elements")


def concat_dspec(specfiles, outfile=None, savespec=False):
    '''
    concatenate a list of specfiles in time axis
    :param specfiles: a list of specfile to concatenate
    :return: concatenated specdata
    '''
    from tqdm import tqdm
    if isinstance(specfiles, list):
        if len(specfiles) > 1:
            specfiles = sorted(specfiles)
        else:
            print('Abort. Only one specfile is provided.')
            return -1
    else:
        print('Please provide a list of specfiles')
        return -1

    specdata = {}

    specdata_ = np.load(specfiles[0])
    for k, v in specdata_.iteritems():
        specdata[k] = v

    for spfile in tqdm(specfiles[1:]):
        specdata_ = np.load(spfile)
        specdata['spec'] = np.concatenate((specdata['spec'], specdata_['spec']), axis=-1)
        specdata['tim'] = np.hstack((specdata['tim'], specdata_['tim']))

    if savespec:
        if not outfile:
            specfile = 'dspec.npz'
        if os.path.exists(specfile):
            os.system('rm -rf ' + specfile)
        np.savez(specfile, **specdata)
    return specdata
