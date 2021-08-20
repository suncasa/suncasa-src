import numpy as np
# import sys
import math
import os, sys, platform
import astropy.units as u
from sunpy import map as smap
from astropy.coordinates import SkyCoord
from suncasa.io import ndfits
from . import gstools  # initialization library - located either in the current directory or in the system path
from suncasa.utils import mstools
import lmfit
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from suncasa.utils import mstools
from suncasa.utils import qlookplot as ql
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from astropy.io import fits
import numpy.ma as ma


# name of the fast gyrosynchrotron codes shared library
if platform.system() == 'Linux' or platform.system() == 'Darwin':
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'binaries/MWTransferArr.so')
if platform.system() == 'Windows':
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'binaries/MWTransferArr64.dll')


def kev2k(eng):
    return 11604525.00617 * eng


def ff_emission(em, T=1.e7, Z=1., mu=1.e10):
    from astropy import constants as const
    import astropy.units as u

    T = T * u.k
    mu = mu * u.Hz
    esu = const.e.esu
    k_B = const.k_B.cgs
    m_e = const.m_e.cgs
    c = const.c.cgs
    bmax = (3 * k_B * T * u.k / m_e) ** 0.5 / 2.0 / np.pi / (mu * u.Hz)
    bmin = Z * esu ** 2 / 3. / k_B / T
    lnbb = np.log((bmax / bmin).value)
    ka_mu = 1. / mu ** 2 / T ** 1.5 * (
            Z ** 2 * esu ** 6 / c / np.sqrt(2. * np.pi * (m_e * k_B) ** 3)) * np.pi ** 2 / 4.0 * lnbb
    # print(ka_mu, em)
    opc = ka_mu * em
    return T.value * (1 - np.exp(-opc.value))


def sfu2tb(freq, flux, area):
    # frequency in Hz
    # flux in sfu
    # area: area of the radio source in arcsec^2
    sfu2cgs = 1e-19
    vc = 2.998e10
    kb = 1.38065e-16
    # sr = np.pi * (size[0] / 206265. / 2.) * (size[1] / 206265. / 2.)
    sr = area / 206265. ** 2
    Tb = flux * sfu2cgs * vc ** 2. / (2. * kb * freq ** 2. * sr)
    return Tb


def tb2sfu(freq, tb, area):
    # frequency in Hz
    # brightness temperature in K
    # area: area of the radio source in arcsec^2
    sfu2cgs = 1e-19
    vc = 2.998e10
    kb = 1.38065e-16
    # sr = np.pi * (size[0] / 206265. / 2.) * (size[1] / 206265. / 2.)
    sr = area / 206265. ** 2
    flux = tb / (sfu2cgs * vc ** 2. / (2. * kb * freq ** 2. * sr))
    return flux


def initspecplot(axes, cplts):
    errobjs = []
    for cpltidx, cplt in enumerate(cplts):
        errobjs.append(axes.errorbar([], [], yerr=[], linestyle='', marker='o', mfc='none', mec=cplt, alpha=1.0))

    axes.set_yscale("log")
    axes.set_xscale("log")
    axes.set_xlim([1, 20])
    axes.set_ylim([0.1, 1000])
    axes.set_xticks([1, 5, 10, 20])
    axes.set_xticklabels([1, 5, 10, 20])
    axes.set_xticks([1, 5, 10, 20])
    axes.set_yticks([])
    axes.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
    axes.set_ylabel('T$_b$ [MK]')
    axes.set_xlabel('Frequency [GHz]')

    x = np.linspace(1, 20, 10)
    for ll in [-1, 0, 1, 2, 3, 4]:
        y = 10. ** (-2 * np.log10(x) + ll)
        axes.plot(x, y, 'k--', alpha=0.1)
        # y2 = 10. ** (-4 * np.log10(x) + ll)
        # y3 = 10. ** (-8 * np.log10(x) + ll)
        # ax_eospec.plot(x, y, 'k--', x, y2, 'k:', x, y3, 'k-.', alpha=0.1)
    return errobjs


def set_errorobj(xout, yout, errobj, yerr=None):
    eospec, dummy, (errbar_eospec,) = errobj
    eospec.set_data(xout, yout)
    if yerr is not None:
        yerr_top = yout + yerr
        yerr_bot = yout - yerr
        new_segments_y = [np.array([[x, yt], [x, yb]]) for x, yt, yb in zip(xout, yerr_top, yerr_bot)]
        errbar_eospec.set_segments(new_segments_y)


def mwspec2min_1src(params, freqghz, tb=None, tb_err=None, arcsec2cm=0.725e8, showplt=False):
    # params are defined by lmfit.Paramters()
    '''
    params: parameters defined by lmfit.Paramters()
    freqghz: frequencies in GHz
    ssz: pixel size in arcsec
    tb: reference brightness temperature in K
    tb_err: uncertainties of reference brightness temperature in K
    '''

    from scipy import interpolate
    GET_MW = gstools.initGET_MW(libname)  # load the library

    ssz = float(params['ssz'].value)  # # source area in arcsec^2
    depth = float(params['depth'].value)  # total source depth in arcsec
    Bmag = float(params['Bmag'].value)  # magnetic field strength in G
    Tth = float(params['Tth'].value)  # thermal temperature in MK
    nth = float(params['nth'].value)  # thermal density in 1e10 cm^{-3}
    nrlh = 10. ** float(params['lognrlh'].value)  # total nonthermal density above 0.1 MeV
    delta = float(params['delta'].value)  # powerlaw index
    theta = float(params['theta'].value)  # viewing angle in degrees
    Emin = float(params['Emin'].value)  # low energy cutoff of nonthermal electrons in MeV
    Emax = float(params['Emax'].value)  # high energy cutoff of nonthermal electrons in MeV
    E_hi = 0.1
    nrl = nrlh * (Emin ** (1. - delta) - Emax * (1. - delta)) / (E_hi ** (1. - delta) - Emax ** (1. - delta))

    Nf = 100  # number of frequencies
    NSteps = 1  # number of nodes along the line-of-sight

    N_E = 15  # number of energy nodes
    N_mu = 15  # number of pitch-angle nodes

    Lparms = np.zeros(11, dtype='int32')  # array of dimensions etc.
    Lparms[0] = NSteps
    Lparms[1] = Nf
    Lparms[2] = N_E
    Lparms[3] = N_mu

    Rparms = np.zeros(5, dtype='double')  # array of global floating-point parameters
    Rparms[0] = ssz * arcsec2cm ** 2  # Area, cm^2
    # Rparms[0] = 1e20  # area, cm^2
    Rparms[1] = 1e9  # starting frequency to calculate spectrum, Hz
    Rparms[2] = 0.02  # logarithmic step in frequency
    Rparms[3] = 12  # f^C
    Rparms[4] = 12  # f^WH

    ParmLocal = np.zeros(24, dtype='double')  # array of voxel parameters - for a single voxel
    ParmLocal[0] = depth * arcsec2cm / NSteps  # voxel depth, cm
    ParmLocal[1] = Tth * 1e6  # T_0, K
    ParmLocal[2] = nth * 1e10  # n_0 - thermal electron density, cm^{-3}
    ParmLocal[3] = Bmag  # B - magnetic field, G

    Parms = np.zeros((24, NSteps), dtype='double', order='F')  # 2D array of input parameters - for multiple voxels
    for i in range(NSteps):
        Parms[:, i] = ParmLocal  # most of the parameters are the same in all voxels
        # if NSteps > 1:
        #     Parms[4, i] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS
        # else:
        #     Parms[4, i] = 50.0  # the viewing angle varies from 50 to 80 degrees along the LOS
        Parms[4, i] = theta

    # parameters of the electron distribution function
    n_b = nrl  # n_b - nonthermal electron density, cm^{-3}
    mu_c = np.cos(np.pi * 70 / 180)  # loss-cone boundary
    dmu_c = 0.2  # Delta_mu

    E_arr = np.logspace(np.log10(Emin), np.log10(Emax), N_E, dtype='double')  # energy grid (logarithmically spaced)
    mu_arr = np.linspace(-1.0, 1.0, N_mu, dtype='double')  # pitch-angle grid

    f0 = np.zeros((N_E, N_mu), dtype='double')  # 2D distribution function array - for a single voxel

    # computing the distribution function (equivalent to PLW & GLC)
    A = n_b / (2.0 * np.pi) * (delta - 1.0) / (Emin ** (1.0 - delta) - Emax ** (1.0 - delta))
    B = 0.5 / (mu_c + dmu_c * np.sqrt(np.pi) / 2 * math.erf((1.0 - mu_c) / dmu_c))
    for i in range(N_E):
        for j in range(N_mu):
            amu = abs(mu_arr[j])
            f0[i, j] = A * B * E_arr[i] ** (-delta) * (1.0 if amu < mu_c else np.exp(-((amu - mu_c) / dmu_c) ** 2))

    f_arr = np.zeros((N_E, N_mu, NSteps), dtype='double',
                     order='F')  # 3D distribution function array - for multiple voxels
    for k in range(NSteps):
        f_arr[:, :, k] = f0  # electron distribution function is the same in all voxels

    RL = np.zeros((7, Nf), dtype='double', order='F')  # input/output array

    # calculating the emission for array distribution (array -> on)
    res = GET_MW(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL)

    if res:
        # retrieving the results
        f = RL[0]
        I_L = RL[5]
        I_R = RL[6]

        if showplt:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1)
            ax.plot(f, I_L + I_R)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title('Total intensity (array)')
            ax.set_xlabel('Frequency, GHz')
            ax.set_ylabel('Intensity, sfu')

        flx_model = I_L + I_R
        flx_model = np.nan_to_num(flx_model) + 1e-11
        logf = np.log10(f)
        logflx_model = np.log10(flx_model)
        logfreqghz = np.log10(freqghz)
        interpfunc = interpolate.interp1d(logf, logflx_model, kind='linear')
        logmflx = interpfunc(logfreqghz)
        mflx = 10. ** logmflx
        mtb = sfu2tb(np.array(freqghz) * 1.e9, mflx, ssz)
    else:
        print("Calculation error!")

    if tb is None:
        return mtb
    if tb_err is None:
        # return mTb - Tb
        return mtb - tb
    # wt = 1./flx_err
    # wt = 1./(Tb_err/Tb/np.log(10.))
    # residual = np.abs((logmTb - np.log10(Tb))) * wt
    # residual = np.abs((mflx - flx)) * wt
    residual = (mtb - tb) / tb_err
    return residual


class RegionSelector:

    # def set_errorobj(self, xout, yout, errobj, yerr):
    #     eospec, dummy, (errbar_eospec,) = errobj
    #     eospec.set_data(xout, yout)
    #     if yerr is not None:
    #         yerr_top = yout + yerr
    #         yerr_bot = yout - yerr
    #         new_segments_y = [np.array([[x, yt], [x, yb]]) for x, yt, yb in zip(xout, yerr_top, yerr_bot)]
    #         errbar_eospec.set_segments(new_segments_y)
    #     return 1

    def subdata(self, xs, ys, rfile):
        rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rdelts = ndfits.read(rfile)
        ny, nx = rmap.data.shape
        tr_coord = rmap.top_right_coord
        bl_coord = rmap.bottom_left_coord
        x0 = bl_coord.Tx.to(u.arcsec).value
        y0 = bl_coord.Ty.to(u.arcsec).value
        x1 = tr_coord.Tx.to(u.arcsec).value
        y1 = tr_coord.Ty.to(u.arcsec).value
        dx = rmap.scale.axis1.to(u.arcsec / u.pix).value
        dy = rmap.scale.axis2.to(u.arcsec / u.pix).value
        mapx, mapy = np.linspace(x0, x1, nx) - dx / 2.0, np.linspace(y0, y1, ny) - dy / 2.0
        xsmin = np.nanmin(xs)
        xsmax = np.nanmax(xs)
        ysmin = np.nanmin(ys)
        ysmax = np.nanmax(ys)
        if np.abs(xsmax - xsmin) < dx:
            xsmax = xsmin + dx
        if np.abs(ysmax - ysmin) < dy:
            ysmax = ysmin + dy
        xmask = np.logical_and(mapx >= xsmin, mapx <= xsmax)
        nxnew = np.count_nonzero(xmask)
        ymask = np.logical_and(mapy >= ysmin, mapy <= ysmax)
        nynew = np.count_nonzero(ymask)
        xmask = np.tile(xmask, ny).reshape(ny, nx)
        ymask = np.tile(ymask, nx).reshape(nx, ny).transpose()
        mask = xmask & ymask
        # print(np.count_nonzero(mask))
        self.npix = np.count_nonzero(mask)
        self.area = self.npix * dx * dy
        data = rdata[:, mask]
        # print(rdata[:, :, mask])
        # print(mask.shape, rdata.shape, data.shape)
        data = np.squeeze(data)
        # print(data.shape)
        return data

    def __init__(self, clkpnts, boxlines, eofiles, errobjs, cfreqs=None, rms=None, eofile_ref=None, errobj_ref=None,
                 wTmap=None, outspec_ff=None, scatter_gsfit=None,
                 get_peak=False, get_sum=False):
        self.boxline = []
        self.clkpnt = []
        self.xs = list(clkpnts[0].get_xdata())
        self.ys = list(clkpnts[0].get_ydata())
        self.npix = None
        self.area = None
        self.xout = []
        self.yout = []
        self.xouterr = []
        self.youterr = []
        for errobj in errobjs:
            eospec, dummy, (errbar_eospec,) = errobj
            self.xout.append(eospec.get_xdata())
            self.yout.append(eospec.get_ydata())
        self.errobjs = errobjs
        self.errobj_ref = errobj_ref
        self.outspec_ff = outspec_ff
        self.scatter_gsfit = scatter_gsfit
        self.cfreqs = cfreqs
        self.rms = rms
        self.eofiles = eofiles
        self.eofile_ref = eofile_ref
        self.wTmap = wTmap
        self.wT = None
        self.em = None
        self.get_peak = get_peak
        self.get_sum = get_sum
        self.tps = []
        self.params = None
        for idx, s in enumerate(clkpnts):
            self.boxline.append(boxlines[idx])
            self.clkpnt.append(s)
            self.cid = s.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        axes = [clkpnt.axes for clkpnt in self.clkpnt]
        if self.clkpnt[0].figure.canvas.toolbar.mode == '':
            if event.inaxes not in axes:
                return
            nxs = len(self.xs)
            if event.button == 1:
                if nxs < 2:
                    self.xs.append(event.xdata)
                    self.ys.append(event.ydata)
                else:
                    self.xs = [event.xdata]
                    self.ys = [event.ydata]
            elif event.button == 3:
                if len(self.xs) > 0:
                    self.xs.pop()
                    self.ys.pop()

            self.get_flux()

    def get_flux(self):
        if len(self.xs) > 0:
            xs = np.array(self.xs, dtype=np.float64)
            ys = np.array(self.ys, dtype=np.float64)
            for clkpnt in self.clkpnt:
                clkpnt.set_data(xs, ys)
        else:
            for clkpnt in self.clkpnt:
                clkpnt.set_data([], [])
        nxs = len(self.xs)
        if nxs <= 1:
            for line in self.boxline:
                line.set_data([], [])
        elif nxs == 2:
            datas = []
            # eofile = self.eofiles[0]
            # rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rdelts = ndfits.read(eofile)
            # data = self.subdata(xs, ys, eofile)
            # datas.append(data)
            for tidx, eofile in enumerate(self.eofiles):
                data = self.subdata(xs, ys, eofile)
                datas.append(data)

            if self.eofile_ref is not None:
                data_ref = self.subdata(xs, ys, self.eofile_ref)
            if self.wTmap is not None:
                datawT = self.subdata(xs, ys, self.wTmap)

            if self.get_peak:
                youts_outspec = []
                for data in datas:
                    if data.ndim > 1:
                        youts_outspec.append(np.nanmax(data, axis=-1) / 1e6)
                    else:
                        youts_outspec.append(data / 1e6)
                if self.eofile_ref is not None:
                    youts_outspec_ref = np.nanmax(data_ref[0, dd, :, :]) / 1e6
            else:
                youts_outspec = []
                for data in datas:
                    if data.ndim > 1:
                        youts_outspec.append(np.nanmean(data, axis=-1) / 1e6)
                    else:
                        youts_outspec.append(data / 1e6)
                if self.eofile_ref is not None:
                    if data.ndim > 1:
                        youts_outspec_ref = np.nanmean(data_ref, axis=-1) / 1e6
                    else:
                        youts_outspec_ref = data_ref / 1e6
            self.tps = []
            for data in datas:
                if data.ndim > 1:
                    self.tps.append(np.nansum(data, axis=-1) / 1e6)
                else:
                    self.tps.append(data / 1e6)
            xout = self.cfreqs
            for tidx, errobj in enumerate(self.errobjs):
                set_errorobj(xout, youts_outspec[tidx], errobj, self.rms)
            if self.eofile_ref is not None:
                set_errorobj(xout, youts_outspec_ref, self.errobj_ref, self.rms)
            if self.wTmap is not None:
                print(datawT.shape)
                wT = np.nanmean(datawT[..., 1]) * 1e6
                em = np.nanmean(datawT[..., 0])
                arcsec2cm = (self.wTmap[0].rsun_meters / self.wTmap[0].rsun_obs).to(u.cm / u.arcsec).value
                # nele = 4.0e10
                # depth = em / nele ** 2 / arcsec2cm
                # print('Temperature: {:.1f} MK, EM: {:.2e} cm-5, depth: {:.1f} arcsec if nele is {:.2e} cm-3'.format(wT / 1e6, em, depth, nele))
                depth = 20.  ## arcsec
                nele = np.sqrt(em / (depth * arcsec2cm))
                print('Temperature: {:.1f} MK, EM: {:.2e} cm-5, nele: {:.2e} cm-3 if depth is {:.1f} arcsec'.format(
                    wT / 1e6, em, nele, depth))
                self.wT = wT
                self.em = em
                yout_ff = np.array([ff_emission(em, T=wT, Z=1., mu=ll) for ll in xout * 1e9]) / 1.e6

                self.outspec_ff.set_data(xout, yout_ff)

            self.errobjs[0][0].figure.canvas.draw_idle()
            for line in self.boxline:
                line.set_data([xs[0], xs[1], xs[1], xs[0], xs[0]], [ys[0], ys[0], ys[1], ys[1], ys[0]])

        clkpnt.figure.canvas.draw_idle()


class GStool:
    # def get_showaia(self):
    #     return self._showaia
    #
    # def set_showaia(self, value):
    #     self._showaia = value
    #
    # showaia = property(fget=get_showaia, fset=set_showaia, doc="`Boolean`-like: Display AIA image or not")

    def __init__(self, eofiles, aiafile=None, xycen=None, fov=None, freqghz_bound=[-1, 100], calpha=0.5,
                 clevels=np.array([0.3, 1.0]), opencontour=None):
        self.aiafile = aiafile
        self.eofiles = eofiles
        self.xycen = xycen
        self.fov = fov
        self.calpha = calpha
        self.clevels = clevels
        self.freqghz_bound = freqghz_bound
        self.opencontour = opencontour
        self._showaia = False

        rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rdelts = ndfits.read(eofiles[0])
        self.bdinfo = bdinfo = ndfits.get_bdinfo(rfreqs, rdelts)
        self.cfreqs = cfreqs = bdinfo['cfreqs']
        self.cfreqs_all = cfreqs_all = bdinfo['cfreqs_all']
        self.freq_dist = lambda fq: (fq - cfreqs_all[0]) / (cfreqs_all[-1] - cfreqs_all[0])

        self.ntim = ntim = len(eofiles)
        self.xlim = xlim = xycen[0] + np.array([-1, 1]) * 0.5 * fov[0]
        self.ylim = ylim = xycen[1] + np.array([-1, 1]) * 0.5 * fov[1]

        nspw = len(rfreqs)
        eodate = Time(rmap.date.mjd + rmap.exposure_time.value / 2. / 24 / 3600, format='mjd')
        ny, nx = rmap.data.shape
        x0, x1 = (np.array([1, rmap.meta['NAXIS1']]) - rmap.meta['CRPIX1']) * rmap.meta['CDELT1'] + \
                 rmap.meta['CRVAL1']
        y0, y1 = (np.array([1, rmap.meta['NAXIS2']]) - rmap.meta['CRPIX2']) * rmap.meta['CDELT2'] + \
                 rmap.meta['CRVAL2']
        dx = rmap.meta['CDELT1']
        dy = rmap.meta['CDELT2']
        mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)

        fig = plt.figure(figsize=(15, 6))
        self.fig = fig

        grids = fig.add_gridspec(ncols=3, nrows=1, width_ratios=[1, 1, 0.6])
        self.grids = grids
        axs = []
        axs.append(fig.add_subplot(grids[0, 0]))
        axs.append(fig.add_subplot(grids[0, 1], sharex=axs[-1], sharey=axs[-1]))
        axs.append(fig.add_subplot(grids[0, 2]))
        if aiafile:
            if os.path.exists(aiafile):
                try:
                    aiacmap = plt.get_cmap('gray_r')
                    aiamap = smap.Map(aiafile)
                    ax = axs[0]
                    aiamap.plot(axes=ax, cmap=aiacmap)
                    ax = axs[1]
                    aiamap.plot(axes=ax, cmap=aiacmap)
                    self._showaia = True
                except:
                    self._showaia = False

        if self._showaia:
            if self.opencontour is None:
                self.opencontour = False
        else:
            if self.opencontour is None:
                self.opencontour = True
        ## Plot EOVSA images as filled contour on top of the AIA image
        icmap = plt.get_cmap('RdYlBu')
        cts = []
        ## color map for spectra from the image series
        tcmap = plt.get_cmap('turbo')

        for s, sp in enumerate(rfreqs):
            data = rdata[s, ...]
            clvls = clevels * np.nanmax(data)
            rcmap = [icmap(self.freq_dist(self.cfreqs[s]))] * len(clvls)
            if self.opencontour:
                cts.append(ax.contour(mapx, mapy, data, levels=clvls,
                                      colors=rcmap,
                                      alpha=calpha))
            else:
                cts.append(ax.contourf(mapx, mapy, data, levels=clvls,
                                       colors=rcmap,
                                       alpha=calpha))
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        for ax in axs[:2]:
            ax.set_xlabel('Solar-X [arcsec]')
            ax.set_ylabel('Solar-y [arcsec]')
            ax.set_title('')
        ax.text(0.02, 0.01,
                ' '.join(['AIA {:.0f} Å'.format(aiamap.wavelength.value),
                          aiamap.date.datetime.strftime('%Y-%m-%dT%H:%M:%S')]),
                ha='left',
                va='bottom',
                color='k', transform=ax.transAxes)
        ax.text(0.02, 0.05, ' '.join(['EOVSA     ', eodate.datetime.strftime('%Y-%m-%dT%H:%M:%S')]), ha='left',
                va='bottom',
                color='k', transform=ax.transAxes)

        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes("right", size="8%", pad=0.08)
        cax.set_visible(False)

        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="8%", pad=0.08)

        ticks, bounds, vmax, vmin, freqmask = ql.get_colorbar_params(bdinfo)

        cb = colorbar.ColorbarBase(cax, norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=icmap,
                                   orientation='vertical', boundaries=bounds, spacing='proportional',
                                   ticks=ticks, format='%4.1f', alpha=calpha)

        for fbd_lo, fbd_hi in freqmask:
            if fbd_hi is not None:
                cax.axhspan(fbd_lo, fbd_hi, hatch='//', edgecolor='k', facecolor='#BBBBBB')
        plt.text(0.5, 1.05, 'MW', ha='center', va='bottom', transform=cax.transAxes, color='k', fontweight='normal')
        plt.text(0.5, 1.01, '[GHz]', ha='center', va='bottom', transform=cax.transAxes, color='k',
                 fontweight='normal')

        cax.xaxis.set_visible(False)
        cax.tick_params(axis="y", pad=-20., length=0, colors='k', labelsize=7)
        cax.axhline(vmin, xmin=1.0, xmax=1.2, color='k', clip_on=False)
        cax.axhline(vmax, xmin=1.0, xmax=1.2, color='k', clip_on=False)
        cax.text(1.25, 0.0, '{:.1f}'.format(vmin), fontsize=9, transform=cax.transAxes, va='center', ha='left')
        cax.text(1.25, 1.0, '{:.1f}'.format(vmax), fontsize=9, transform=cax.transAxes, va='center', ha='left')

        boxlines = []
        clkpnts = []

        for idx, ax in enumerate(axs[:2]):
            if idx == 0:
                c = 'g'
            elif idx == 1:
                c = 'b'
            else:
                c = 'k'
            line, = ax.plot([], [], '-', c=c, alpha=1.0)  # empty line
            boxlines.append(line)
            clkpnt, = ax.plot([], [], '+', c='white', alpha=0.7)  # empty line
            clkpnts.append(clkpnt)

        if ntim < 2:
            cplts = ['k']
        else:
            cplts = tcmap(np.linspace(0, 1, ntim))

        self.cplts = cplts
        self.ax_eospec = axs[-1]

        errobjs = initspecplot(self.ax_eospec, cplts)
        grids.tight_layout(fig)

        self.region = RegionSelector(clkpnts, boxlines, eofiles, errobjs, cfreqs=cfreqs, rms=None, wTmap=None)
        self.scatter_eospecs_fit = []
        self.scatter_eospecs = []

    def set_params(self, params):
        ssz = self.region.area  # source area in arcsec^2
        params.add('ssz', value=ssz, vary=False)  # pixel size in arcsec^2
        self.params = params

    def plot_components(self):
        ti = 0
        tb = self.region.errobjs[ti][0].get_ydata() * 1e6
        tb_ma = ma.masked_less_equal(tb, 0)
        freqghz = self.region.errobjs[0][0].get_xdata()
        # freqghz_ma = ma.masked_outside(freqghz, 1.0, 15.0)
        freqghz_ma = ma.masked_outside(freqghz, self.freqghz_bound[0], self.freqghz_bound[1])
        mask_fit = np.logical_or(freqghz_ma.mask, tb_ma.mask)
        freqghz_ma = ma.masked_array(freqghz, mask_fit)
        tb_ma = ma.masked_array(tb, mask_fit)

        # scatter_eospecs_fit.append(
        #     ax_spec.plot(freqghz_ma, tb_ma / 1.e6, marker='o', linestyle='', c=cplts[ti]))

        # flx_rms = rms
        tb_err = tb * 0.0
        tb_err[:] = 1.e6
        tb_err_ma = ma.masked_array(tb_err, tb_ma.mask)
        if len(self.scatter_eospecs_fit) == 0:
            for ti, cplt in enumerate(self.cplts):
                self.scatter_eospecs_fit.append(
                    self.ax_eospec.errorbar(freqghz_ma, tb_ma / 1.e6, yerr=tb_err_ma / 1.e6, marker='.', ms=1,
                                            linestyle='',
                                            c=cplt))
        else:
            for ti, cplt in enumerate(self.cplts):
                set_errorobj(freqghz_ma, tb_ma / 1.e6, self.scatter_eospecs_fit[ti], yerr=tb_err_ma / 1.e6)

    def fit(self):
        ti = 0
        tb = self.region.errobjs[ti][0].get_ydata() * 1e6
        tb_ma = ma.masked_less_equal(tb, 0)
        freqghz = self.region.errobjs[0][0].get_xdata()
        # freqghz_ma = ma.masked_outside(freqghz, 1.0, 15.0)
        freqghz_ma = ma.masked_outside(freqghz, self.freqghz_bound[0], self.freqghz_bound[1])
        mask_fit = np.logical_or(freqghz_ma.mask, tb_ma.mask)
        freqghz_ma = ma.masked_array(freqghz, mask_fit)
        tb_ma = ma.masked_array(tb, mask_fit)

        # scatter_eospecs_fit.append(
        #     ax_spec.plot(freqghz_ma, tb_ma / 1.e6, marker='o', linestyle='', c=cplts[ti]))

        # flx_rms = rms
        tb_err = tb * 0.1
        # tb_err[:] = 0.2e6
        tb_err_ma = ma.masked_array(tb_err, tb_ma.mask)
        if len(self.scatter_eospecs_fit) == 0:
            for ti, cplt in enumerate(self.cplts):
                self.scatter_eospecs_fit.append(
                    self.ax_eospec.errorbar(freqghz_ma, tb_ma / 1.e6, yerr=tb_err_ma / 1.e6, marker='.', ms=1,
                                            linestyle='', c=cplt))
        else:
            for ti, cplt in enumerate(self.cplts):
                set_errorobj(freqghz_ma, tb_ma / 1.e6, self.scatter_eospecs_fit[ti], yerr=tb_err_ma / 1.e6)

        mini = lmfit.Minimizer(mwspec2min_1src, self.params, fcn_args=(freqghz_ma.compressed(),),
                               fcn_kws={'tb': tb_ma.compressed(), 'tb_err': tb_err_ma.compressed()},
                               nan_policy='omit')
        method = 'nelder'
        # # method = 'differential_evolution'
        mi = mini.minimize(method=method)
        print(method + ' minimization results')
        print(lmfit.fit_report(mi.params))

        tb_fit = mwspec2min_1src(mi.params, freqghz)
        if len(self.scatter_eospecs) == 0:
            for ti, cplt in enumerate(self.cplts):
                self.scatter_eospecs.append(self.ax_eospec.plot(freqghz, tb_fit / 1.e6, linestyle='-', c=cplt))
        else:
            for ti, cplt in enumerate(self.cplts):
                self.scatter_eospecs[ti][0].set_data(freqghz, tb_fit / 1.e6)
