import os

__all__ = ['Dspec']

import numpy as np
import os
import struct
import matplotlib.pyplot as plt
from astropy.time import Time

try:
    ## Full Installation of CASA 4, 5 and 6
    from taskinit import ms, tb, qa
except:
    ## Modular Installation of CASA 6
    from casatools import table as tbtool
    from casatools import ms as mstool
    from casatools import quanta as qatool

    tb = tbtool()
    ms = mstool()
    qa = qatool()

from matplotlib.dates import AutoDateFormatter, AutoDateLocator
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


class Dspec:
    data = None
    time_axis = None
    freq_axis = None
    instrument = None
    t_label = None
    f_label = None
    bl = None
    uvrange = None
    pols = None


    def __init__(self, fname=None, specfile=None, bl='', uvrange='', field='', scan='',
                 datacolumn='data', domedian=False, timeran=None, spw=None, timebin='0s', regridfreq=False,
                 fillnan=None, verbose=False, usetbtool=True):
        if fname:
            if isinstance(fname, str):
                if os.path.exists(fname):
                    if os.path.isdir(fname):
                        self.get_dspec(fname, specfile=specfile, bl=bl, uvrange=uvrange, field=field,
                                       scan=scan, datacolumn=datacolumn, domedian=domedian, timeran=timeran, spw=spw,
                                       timebin=timebin, regridfreq=regridfreq, fillnan=fillnan, verbose=verbose,
                                       usetbtool=usetbtool)
                    else:
                        self.read(fname)
            else:
                self.read(fname)

    def read(self, fname):
        if fname.endswith('.fts'):
            from .sources import eovsa
            s = eovsa.get_dspec(fname, doplot=False)
            self.data = s['spectrogram']
            self.time_axis = Time(s['time_axis'], format='mjd')
            self.freq_axis=s['spectrum_axis'] * 1e9
        else:
            spec, bl, tim, freq, pols = self.rd_dspec(fname, zaxis='amp')
            self.data = spec
            self.time_axis = Time(tim / 24 / 3600, format='mjd')
            self.freq_axis=freq
            self.bl = bl
            self.pols = pols

    def get_dspec(self, fname=None, specfile=None, bl='', uvrange='', field='', scan='',
                  datacolumn='data',
                  domedian=False, timeran=None, spw=None, timebin='0s', regridfreq=False, fillnan=None, verbose=False,
                  usetbtool=False):
        if fname.endswith('/'):
            fname = fname[:-1]
        msfile = fname
        if not spw:
            spw = ''
        if not timeran:
            timeran = ''
        if domedian:
            if not uvrange:
                uvrange = '0.2~0.8km'
        else:
            uvrange = ''
        if not bl:
            bl = ''
        else:
            uvrange = ''
        # Open the ms and plot dynamic spectrum
        if verbose:
            print('Splitting selected data...')

        if usetbtool:
            try:
                tb.open(fname + '/POLARIZATION')
                corrtype = tb.getcell('CORR_TYPE', 0)
                pols = [stokesenum[p] for p in corrtype]
                tb.close()
            except:
                pols = []

            antmask = []
            if uvrange is not '' or bl is not '':
                ms.open(fname)
                ms.selectinit(datadescid=0)
                mdata = ms.metadata()
                antlist = mdata.antennaids()
                mdata.done()
                staql = {'uvdist': uvrange, 'baseline': bl, 'spw': spw, 'field': field, 'scan': scan,
                         'timerange': timeran}
                ### todo the selection only works for uvrange and bl. To make the selection of other items works,
                ## I need to make mask for other items.
                a = ms.msselect(staql)
                mdata = ms.metadata()
                baselines = mdata.baselines()
                for lidx, l in enumerate(antlist):
                    antmask.append(baselines[l][antlist[lidx:]])
                antmask = np.hstack(antmask)
                mdata.done()
                ms.close()

            tb.open(fname)
            spwtb = tbtool()
            spwtb.open(fname + '/SPECTRAL_WINDOW')
            ptb = tbtool()
            ptb.open(fname + '/POLARIZATION')

            ms.open(fname)
            spwlist = []
            mdata = ms.metadata()
            nspw = mdata.nspw()
            nbl = mdata.nbaselines() + mdata.nantennas()
            nscans = mdata.nscans()
            spw_nfrq = []  # List of number of frequencies in each spw
            for i in range(nspw):
                spw_nfrq.append(mdata.nchan(i))
            spw_nfrq = np.array(spw_nfrq)
            nf = np.sum(spw_nfrq)
            smry = mdata.summary()
            scan_ntimes = []  # List of number of times in each scan
            for iscan in range(nscans):
                scan_ntimes.append(
                    smry['observationID=0']['arrayID=0']['scan=' + str(iscan)]['fieldID=0']['nrows'] / nspw / nbl)
            scan_ntimes = np.array(scan_ntimes)
            scan_ntimes_integer = scan_ntimes.astype(np.int)
            if len(np.where(scan_ntimes % scan_ntimes_integer != 0)[0]) != 0:
                # if True:
                scan_ntimes = []  # List of number of times in each scan
                for iscan in range(nscans):
                    scan_ntimes.append(
                        len(smry['observationID=0']['arrayID=0']['scan=' + str(iscan)]['fieldID=0'].keys()) - 6)
                scan_ntimes = np.array(scan_ntimes)
            else:
                scan_ntimes = scan_ntimes_integer

            nt = np.sum(scan_ntimes)
            times = tb.getcol('TIME')
            if times[nbl] - times[0] != 0:
                # This is frequency/scan sort order
                order = 'f'
            elif times[nbl * nspw - 1] - times[0] != 0:
                # This is time sort order
                order = 't'
            npol = ptb.getcol('NUM_CORR', 0, 1)[0]
            ptb.close()
            freq = np.zeros(nf, float)
            times = np.zeros(nt, float)
            if order == 't':
                specamp = np.zeros((npol, nf, nbl, nt), np.complex)
                for j in range(nt):
                    fptr = 0
                    # Loop over spw
                    for i in range(nspw):
                        # Get channel frequencies for this spw (annoyingly comes out as shape (nf, 1)
                        cfrq = spwtb.getcol('CHAN_FREQ', i, 1)[:, 0]
                        if j == 0:
                            # Only need this the first time through
                            spwlist += [i] * len(cfrq)
                        if i == 0:
                            times[j] = tb.getcol('TIME', nbl * (i + nspw * j), 1)  # Get the time
                        spec_ = tb.getcol('DATA', nbl * (i + nspw * j), nbl)  # Get complex data for this spw
                        flag = tb.getcol('FLAG', nbl * (i + nspw * j), nbl)  # Get flags for this spw
                        nfrq = len(cfrq)
                        # Apply flags
                        if type(fillnan) in [int, float]:
                            spec_[flag] = float(fillnan)
                        else:
                            spec_[flag] = 0.0
                        # Insert data for this spw into larger array
                        specamp[:, fptr:fptr + nfrq, :, j] = spec_
                        freq[fptr:fptr + nfrq] = cfrq
                        fptr += nfrq
            else:
                specf = np.zeros((npol, nf, nt, nbl), np.complex)  # Array indexes are swapped
                iptr = 0
                for j in range(nscans):
                    # Loop over scans
                    for i in range(nspw):
                        # Loop over spectral windows
                        s = scan_ntimes[j]
                        f = spw_nfrq[i]
                        s1 = np.sum(scan_ntimes[:j])  # Start time index
                        s2 = np.sum(scan_ntimes[:j + 1])  # End time index
                        f1 = np.sum(spw_nfrq[:i])  # Start freq index
                        f2 = np.sum(spw_nfrq[:i + 1])  # End freq index
                        spec_ = tb.getcol('DATA', iptr, nbl * s)
                        flag = tb.getcol('FLAG', iptr, nbl * s)
                        if j == 0:
                            cfrq = spwtb.getcol('CHAN_FREQ', i, 1)[:, 0]
                            freq[f1:f2] = cfrq
                            spwlist += [i] * len(cfrq)
                        times[s1:s2] = tb.getcol('TIME', iptr, nbl * s).reshape(s, nbl)[:, 0]  # Get the times
                        iptr += nbl * s
                        # Apply flags
                        if type(fillnan) in [int, float]:
                            spec_[flag] = float(fillnan)
                        else:
                            spec_[flag] = 0.0
                        # Insert data for this spw into larger array
                        specf[:, f1:f2, s1:s2] = spec_.reshape(npol, f, s, nbl)
                        # Swap the array indexes back to the desired order
                specamp = np.swapaxes(specf, 2, 3)
            tb.close()
            spwtb.close()
            ms.close()
            if len(antmask) > 0:
                specamp = specamp[:, :, np.where(antmask)[0], :]
            (npol, nfreq, nbl, ntim) = specamp.shape
            tim = times
        else:
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
                split(vis=msfile, outputvis=vis_spl, datacolumn=datacolumn, timerange=timeran, spw=spw, antenna=bl,
                      field=field,
                      scan=scan, uvrange=uvrange, timebin=timebin)
            except:
                ms.open(msfile, nomodify=True)
                ms.split(outputms=vis_spl, whichcol=datacolumn, time=timeran, spw=spw, baseline=bl, field=field,
                         scan=scan,
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

            tim = data['time']

        if verbose:
            print('npol, nfreq, nbl, ntime:', (npol, nfreq, nbl, ntim))
        spec = np.swapaxes(specamp, 2, 1)

        if domedian:
            if verbose:
                print('doing median of all the baselines')
            # mask zero values before median
            # spec_masked = np.ma.masked_where(spec < 1e-9, spec)
            # spec_masked2 = np.ma.masked_invalid(spec)
            # spec_masked = np.ma.masked_array(spec, mask=np.logical_or(spec_masked.mask, spec_masked2.mask))
            # spec_med = np.ma.filled(np.ma.median(spec_masked, axis=1), fill_value=0.)
            spec = np.abs(spec)
            spec_med = np.nanmedian(spec, axis=1)
            nbl = 1
            ospec = spec_med.reshape((npol, nbl, nfreq, ntim))
        else:
            ospec = spec
        # Save the dynamic spectral data
        if not specfile:
            specfile = msfile + '.dspec.npz'
        if os.path.exists(specfile):
            os.system('rm -rf ' + specfile)
        np.savez(specfile, spec=ospec, tim=tim, freq=freq,
                 timeran=timeran, spw=spw, bl=bl, uvrange=uvrange, pol=pols)
        if verbose:
            print('Median dynamic spectrum saved as: ' + specfile)

        self.read(specfile)
        return specfile

    def wrt_dspec(self, specfile=None, specdat=None):
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

    def rd_dspec(self, specdata, zaxis='amp'):
        zaxis = zaxis.lower()
        if type(specdata) is str:
            try:
                specdata = np.load(specdata)
            except:
                specdata = np.load(specdata, encoding='latin1')
        try:
            if np.iscomplexobj(specdata['spec']):
                if zaxis == 'amp':
                    spec = np.abs(specdata['spec']) / 1.e4
                elif zaxis == 'phase':
                    spec = np.angle(specdata['spec'])
                else:
                    raise ValueError('zaxis must be amp or phase!')
            else:
                zaxis = 'amp'
                spec = specdata['spec'] / 1.e4
            try:
                bl = specdata['bl'].item()
            except:
                bl = specdata['bl']

            tim = specdata['tim']
            freq = specdata['freq']
            if 'pol' in specdata.keys():
                pols = specdata['pol']
            else:
                pols = []

            return spec, bl, tim, freq, pols
        except:
            raise ValueError('format of specdata not recognized. Check your input')

    def concat_dspec(self, specfiles, outfile=None, savespec=False):
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


    def peek(self, *args, **kwargs):
        """
        Plot dynamaic spectrum onto current axes.

        Parameters
        ----------
        *args : dict

        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """
        plt.figure()
        ret = self.plot(*args, **kwargs)
        plt.show()
        return ret

    def plot(self, pol='I', vmin=None, vmax=None, norm=None, cmap='turbo', timerange=None, freqrange=None, ignore_gaps=True):
        """
        timerange: format: ['2021-05-07T18:00:00','2021-05-07T19:00:00']
        freqrange: format: [1.,18.] in GHz
        """

        # Set up variables
        import matplotlib.colors as colors
        import matplotlib.pyplot as plt
        from astropy.time import Time
        if pol not in ['RR', 'LL', 'RRLL', 'XX', 'YY', 'XXYY', 'I', 'V', 'IV']:
            print("Please enter 'RR', 'LL', 'RRLL','XX', 'YY', 'XXYY', 'I', 'V', 'IV' for pol")
            return 0

        try:
            cmap = plt.get_cmap(cmap)
        except:
            cmap = plt.get_cmap('viridis')

        if norm is None:
            norm = colors.Normalize(vmax=vmax, vmin=vmin)
        spec = self.data
        bl = self.bl
        freq = self.freq_axis
        pols = self.pols

        if spec.ndim == 2:
            nfreq, ntim = len(self.freq_axis), len(self.time_axis)
            npol = 1
            nbl = 1
            pol = 'I'
        else:
            (npol, nbl, nfreq, ntim) = spec.shape

        tim_ = self.time_axis
        tim_plt = tim_.plot_date

        if timerange:
            if isinstance(timerange[0], str):
                timerange = Time(timerange)
            tidx = np.where((tim_ >= timerange[0]) & (tim_ <= timerange[1]))[0]
        else:
            tidx = range(ntim)

        if ignore_gaps:
            df = np.median(freq[1:] - freq[:-1])
            fbreak, = np.where(freq[1:] - freq[:-1] > 2 * df)
            for n, fb in enumerate(fbreak):
                loc = fb + n + 1
                freq = np.concatenate((freq[:loc], np.array([freq[loc - 1] + df]), freq[loc:]))

        if freqrange:
            fidx = np.where((freq >= freqrange[0] * 1e9) & (freq <= freqrange[1] * 1e9))[0]
        else:
            fidx = range(nfreq)

        # setup plot parameters
        print('ploting dynamic spectrum...')

        for b in range(nbl):
            if pol not in ['RRLL', 'IV', 'XXYY']:
                if spec.ndim == 2:
                    spec_plt = spec
                else:
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
                if ignore_gaps:
                    for n, fb in enumerate(fbreak):
                        loc = fb + n + 1
                        spec_plt = np.concatenate((spec_plt[:loc], np.zeros((1, ntim)) + np.nan, spec_plt[loc:]), 0)

                fig = plt.figure(figsize=(8, 4), dpi=100)
                ax = fig.add_subplot(111)
                freqghz = freq / 1e9
                im = ax.pcolormesh(tim_plt, freqghz, spec_plt, cmap=cmap, norm=norm, shading='auto')
                divider = make_axes_locatable(ax)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax, cax=cax_spec)
                clb_spec.set_label('Flux [sfu]')
                ax.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

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
                # elif pol in ['IV']:
                else:
                    spec_plt_1 = I_plot
                    spec_plt_2 = V_plot
                    polstr = ['I', 'V']

                if ignore_gaps:
                    for n, fb in enumerate(fbreak):
                        loc = fb + n + 1
                        spec_plt_1 = np.concatenate((spec_plt_1[:loc], np.zeros((1, ntim)) + np.nan, spec_plt_1[loc:]),
                                                    0)
                        spec_plt_2 = np.concatenate((spec_plt_2[:loc], np.zeros((1, ntim)) + np.nan, spec_plt_2[loc:]),
                                                    0)
                ax1 = fig.add_subplot(211)
                freqghz = freq / 1e9
                im = ax1.pcolormesh(tim_plt, freqghz, spec_plt_1, cmap=cmap, norm=norm, shading='auto')
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

                locator = AutoDateLocator(minticks=2)
                ax1.xaxis.set_major_locator(locator)
                ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
                ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[0])
                ax1.set_autoscale_on(False)
                ax2 = fig.add_subplot(212)
                im = ax2.pcolormesh(tim_plt, freqghz, spec_plt_2, cmap=cmap, norm=norm, shading='auto')
                divider = make_axes_locatable(ax2)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax2, cax=cax_spec)
                clb_spec.set_label('Flux [sfu]')
                ax2.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
                ax2.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

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
        return fig