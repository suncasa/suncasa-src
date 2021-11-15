def lightcurves(timerange, outdir='./', specfile=None, goes=True, hessifile=None, fermifile=None, ylog=False, hessi_smoth=0, dspec_cmap='cubehelix',
                vmax=None, vmin=None):
    from sunpy.lightcurve import GOESLightCurve
    from sunpy.time import TimeRange, parse_time
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import numpy as np
    import numpy.ma as ma
    from scipy.signal import medfilt
    import matplotlib.colors as colors
    from scipy.io import readsav
    from suncasa.utils import DButil
    import os

    timerange = Time(timerange)
    if hessifile:
        if not os.path.exists(hessifile):
            hessi_script = 'HESSI_lc.pro'
            print('Run the script {} in SSWIDL to download RHESSI summary file first!'.format(hessi_script))
            fi = open(hessi_script, 'wb')
            fi.write(b'search_network, /enable \n')
            text = "time_range = ['{}','{}'] \n".format(timerange[0].datetime.strftime('%Y-%b-%d %H:%M:%S'),
                                                          timerange[1].datetime.strftime('%Y-%b-%d %H:%M:%S'))

            fi.write(text.encode())
            fi.write(b'obs_obj = hsi_obs_summary(obs_time_interval=time_range) \n')
            fi.write(b'data = obs_obj->getdata() \n')
            fi.write(b'info = obs_obj->get( / info) \n')
            fi.write(b'obs_data = data.countrate \n')
            fi.write(b'obs_times = obs_obj->getdata(/time) \n')
            fi.write(b'obs_times_str  = anytim(obs_times,/CCSDS) \n')
            fi.write(b'obs_energies  = fltarr(2,n_elements(info.energy_edges)-1) \n')
            fi.write(b'for ll=0,n_elements(info.energy_edges)-2 do begin \n')
            fi.write(b'    obs_energies[0,ll]=info.energy_edges[ll] \n')
            fi.write(b'    obs_energies[1,ll]=info.energy_edges[ll+1] \n')
            fi.write(b'endfor \n')
            text = 'save,filename="{}",OBS_DATA,OBS_ENERGIES,obs_times_str \n'.format(hessifile)
            fi.write(text.encode())
            fi.write(b'end \n')
            fi.close()
            return -1
        hessi = readsav(hessifile)
        hessi_tim = Time(list(hessi['obs_times_str']))
        hessi_tim_plt = hessi_tim.plot_date

    specdata = np.load(specfile)
    spec = specdata['spec']
    if len(spec.shape) == 4:
        (npol, nbl, nfreq, ntim) = spec.shape
    else:
        (nfreq, ntim) = spec.shape
    spec = np.mean(np.mean(spec, axis=0), axis=0)
    freq = specdata['freq']
    freqghz = freq / 1e9
    spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')

    fidx_plt = np.linspace(0, nfreq - nfreq / 8.0 / 2.0, 8).astype(np.int) + nfreq / 8.0 / 2.0

    try:
        plt.style.use('seaborn-bright')
        params = {'font.size': 8, 'axes.grid': False, 'axes.facecolor': 'w', 'xtick.color': '#555555', 'ytick.color': '#555555',
                  'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'xtick.minor.size': 1.0, 'ytick.minor.size': 1.0, 'axes.axisbelow': False,
                  'axes.xmargin': 0.0, 'axes.ymargin': 0.0, 'axes.linewidth': 0.5, 'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5,
                  'lines.linewidth': 1.0}
        mpl.rcParams.update(params)
    except:
        pass

    if goes:
        tr = TimeRange(timerange.iso)
        goes = GOESLightCurve.create(tr)
        dates = mpl.dates.date2num(parse_time(goes.data.index))

    tr_plt = Time(timerange)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 6))
    ax = axs[0]
    tidx_hessi, = np.where((hessi_tim >= tr_plt[0]) & (hessi_tim <= tr_plt[1]))
    for idx, eg in enumerate(hessi['obs_energies']):
        flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
        if hessi_smoth > 0:
            flux_hessi = DButil.smooth(flux_hessi, 20)
            flux_hessi = flux_hessi / np.nanmax(flux_hessi)
            ax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth), label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
        else:
            ax.step(hessi_tim_plt[tidx_hessi], flux_hessi, label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))

    if ylog:
        ax.set_yscale("log", nonposy='clip')
    else:
        ax.set_yscale("linear", nonposy='clip')
    ax.legend()
    ax.set_ylabel('Count Rate [s$^{-1}$ detector$^{-1}$ ]')

    ax = axs[1]
    tidx_goes, = np.where((dates >= tr_plt[0].plot_date) & (dates <= tr_plt[1].plot_date))
    if goes:
        ax.plot(dates[tidx_goes], goes.data['xrsb'][tidx_goes] / np.nanmax(goes.data['xrsb'][tidx_goes]), label='GOES 1.0--8.0 $\AA$')
        ax.plot(dates[tidx_goes], goes.data['xrsa'][tidx_goes] / np.nanmax(goes.data['xrsa'][tidx_goes]), label='GOES 0.5--4.0 $\AA$')
    tidx_spec, = np.where((spec_tim >= tr_plt[0]) & (spec_tim <= tr_plt[1]))
    if len(tidx_spec) > 1:
        spec_tim_plt = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
        flux_colors = []
        for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
            flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 7)
            p = ax.plot(spec_tim_plt, flux_plt / np.nanmax(flux_plt), label='{:.2f} GHz'.format(freqghz[min(fidx, nfreq - 1)]))
            flux_colors.append(p[0].get_color())
        ax.set_ylabel('Flux (Normalized)')
        ax.set_ylim(0, 1.1)
        ax.legend()

        ax = axs[2]
        spec_plt = spec[:, tidx_spec[0]:tidx_spec[-1]]
        ax.pcolormesh(spec_tim_plt, freqghz, spec_plt, cmap=dspec_cmap, norm=colors.LogNorm(vmax=vmax, vmin=vmin))
        ax.set_ylabel('Frequency [GHz]')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(formatter)
        ax.fmt_xdata = formatter
        ax.set_xlim(tr_plt.plot_date)
        ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
            ax.axhline(freqghz[min(fidx, nfreq - 1)], color=flux_colors[idx], ls=':')
    else:
        print('Warning: No radio data in the timerange. Proceed without dynamic spectrum.')

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.06)
    imgdir = outdir + '/fig01_{}-{}.png'.format(tr_plt[0].datetime.strftime('%H%M%S'), tr_plt[1].datetime.strftime('%H%M%S'))
    fig.savefig(imgdir, dpi=200)
    print('Save image to ' + imgdir)
