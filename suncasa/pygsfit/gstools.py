import numpy as np
import os, platform
import astropy.units as u
from astropy import constants as const
import ctypes
from numpy.ctypeslib import ndpointer
from scipy import interpolate
import warnings

warnings.simplefilter("default")


def initGET_MW(libname):
    """
    Python wrapper for fast gyrosynchrotron codes.
    Identical to GScodes.py in https://github.com/kuznetsov-radio/gyrosynchrotron
    This is for the single thread version
    @param libname: path for locating compiled shared library
    @return: An executable for calling the GS codes in single thread
    """
    _intp = ndpointer(dtype=ctypes.c_int32, flags='F')
    _doublep = ndpointer(dtype=ctypes.c_double, flags='F')

    libc_mw = ctypes.CDLL(libname)
    mwfunc = libc_mw.pyGET_MW
    mwfunc.argtypes = [_intp, _doublep, _doublep, _doublep, _doublep, _doublep, _doublep]
    mwfunc.restype = ctypes.c_int

    return mwfunc


def sfu2tb(frequency, flux, area=None, size=None, square=True, reverse=False, verbose=False):
    """
        frequency: single element or array, in Hz
        flux: single element or array of flux, in sfu; if reverse, it is brightness temperature in K
        area: area in arcsec^2
        size: Two-dimensional width of the radio source, [major, minor], in arcsec.
              Ignored if both area and size are provided
        reverse: if True, convert brightness temperature in K to flux in sfu integrated uniformly withing the size
    """

    c = const.c.cgs
    k_B = const.k_B.cgs
    sfu = u.jansky * 1e4

    if (not 'area' in vars()) and (not 'size' in vars()):
        print('Neither area nor size is provided. Abort...')

    if not hasattr(frequency, 'unit'):
        # assume frequency is in Hz
        frequency = frequency * u.Hz

    if not hasattr(flux, 'unit'):
        # assume flux is in sfu
        if reverse:
            flux = flux * u.K
        else:
            flux = flux * sfu

    has_area = area is not None
    has_size = size is not None
    if has_area:
        if not hasattr(area, 'unit'):
            # assume area is in arcsec^2
            area = area * u.arcsec ** 2

    if has_size and (not has_area):
        if not isinstance(size, list):
            size = [size]

        if len(size) > 2:
            print('size needs to have 1 or 2 elements.')
        elif len(size) < 2:
            if verbose:
                print('Only one element is provided for source size. Assume symmetric source')
            if not hasattr(size[0], 'unit'):
                # assume size in arcsec
                size[0] = size[0] * u.arcsec
            # define half size
            a = b = size[0] / 2.
        else:
            if not hasattr(size[0], 'unit'):
                # assume size in arcsec
                size[0] = size[0] * u.arcsec
            if not hasattr(size[1], 'unit'):
                # assume size in arcsec
                size[1] = size[1] * u.arcsec
            # define half size
            a = size[0] / 2.
            b = size[1] / 2.
        if square:
            if verbose:
                print('Assume square-shaped source.')
            area = 4. * a * b
        else:
            if verbose:
                print('Assume elliptical-shaped source.')
            area = np.pi * a * b

    sr = area.to(u.radian ** 2)
    factor = c ** 2. / (2. * k_B * frequency ** 2. * sr)

    if reverse:
        # returned value is flux in sfu
        if verbose:
            print('converting input brightness temperature in K to flux density in sfu.')
        return (flux / factor).to(sfu, equivalencies=u.dimensionless_angles())
    else:
        # returned value is brightness temperature in K
        if verbose:
            print('converting input flux density in sfu to brightness temperature in K.')
        return (flux * factor).to(u.K, equivalencies=u.dimensionless_angles())


def ff_emission(em, T=1.e7, Z=1., mu=1.e10):
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

    opc = ka_mu * em
    return T.value * (1 - np.exp(-opc.value))


class GSCostFunctions:
    def SinglePowerLawMinimizerOneSrc(fit_params, freqghz, spec=None, spec_err=None,
                                      spec_in_tb=True, pgplot_widget=None, show_plot=False, debug=False, verbose=False):
        """
        params: parameters defined by lmfit.Paramters()
        freqghz: frequencies in GHz
        spec: input spectrum, can be brightness temperature in K, or flux density in sfu
        spec_err: uncertainties of spectrum in K or sfu
        spec_in_tb: if True, input is brightness temperature in K, otherwise is flux density in sfu
        calc_flux: Default (False) is to return brightness temperature.
                    True if return the calculated flux density. Note one needs to provide src_area/src_size for this
                        option. Otherwise assumes src_size = 2 arcsec (typical EOVSA pixel size).
        @rtype: 1. If no tb/tb_err or flux/flux_err is provided, return the calculated
                    brightness temperature or flux for each input frequency.
                2. If tb/tb_err or flux/flux_err are provided, return the
                    (scaled) residual for each input frequency
        """
        if platform.system() == 'Linux' or platform.system() == 'Darwin':
            libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'binaries/MWTransferArr.so')
        if platform.system() == 'Windows':  ##TODO: not yet tested on Windows platform
            libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'binaries/MWTransferArr64.dll')
        GET_MW = initGET_MW(libname)  # load the library

        asec2cm = 0.725e8
        if 'area_asec2' in fit_params.keys():
            src_area = float(fit_params['area_asec2'].value)  # source area in arcsec^2
        else:
            src_area = 4.  # arcsec^2. default area for bright temperature spectral fitting. Will be divided out.
            if not spec_in_tb:
                print('=======Warning: no source area is provided for flux density calculation. '
                      'Use area = 4 arcsec^2 as the place default (1 EOVSA pixel).======')
        src_area_cm2 = src_area * asec2cm ** 2.  # source area in cm^2
        depth_cm = float(fit_params['depth_asec'].value) * asec2cm  # total source depth in cm
        Bmag = float(fit_params['Bx100G'].value) * 100.  # magnetic field strength in G
        Tth = float(fit_params['T_MK'].value) * 1e6  # thermal temperature in K
        nth = 10. ** float(fit_params['log_nth'].value)  # thermal density
        nrl = 10. ** float(fit_params['log_nnth'].value)  # total nonthermal density above E_min
        delta = float(fit_params['delta'].value)  # powerlaw index
        theta = float(fit_params['theta'].value)  # viewing angle in degrees
        Emin = float(fit_params['Emin_keV'].value) / 1e3  # low energy cutoff of nonthermal electrons in MeV
        Emax = float(fit_params['Emax_MeV'].value)  # high energy cutoff of nonthermal electrons in MeV
        if debug:
            # debug against previous codes
            print('depth, Bmag, Tth, nth/1e10, lognrl, delta, theta, Emin, Emax: '
                  '{0:.1f}, {1:.1f}, {2:.1f}, {3:.1f}, {4:.1f}, '
                  '{5:.1f}, {6:.1f}, {7:.2f}, {8:.1f}'.format(depth_cm / 0.725e8, Bmag, Tth / 1e6, nth / 1e10,
                                                              np.log10(nrl), delta, theta, Emin, Emax))
        # E_hi = 0.1
        # nrl = nrlh * (Emin ** (1. - delta) - Emax * (1. - delta)) / (E_hi ** (1. - delta) - Emax ** (1. - delta))

        Nf = 100  # number of frequencies
        NSteps = 1  # number of nodes along the line-of-sight

        Lparms = np.zeros(11, dtype='int32')  # array of dimensions etc.
        Lparms[0] = NSteps
        Lparms[1] = Nf

        Rparms = np.zeros(5, dtype='double')  # array of global floating-point parameters
        Rparms[0] = src_area_cm2  # Area, cm^2
        Rparms[1] = 0.8e9  # starting frequency to calculate spectrum, Hz
        Rparms[2] = 0.02  # logarithmic step in frequency
        Rparms[3] = 0  # f^C
        Rparms[4] = 0  # f^WH

        ParmLocal = np.zeros(24, dtype='double')  # array of voxel parameters - for a single voxel
        ParmLocal[0] = depth_cm / NSteps  # voxel depth, cm
        ParmLocal[1] = Tth  # T_0, K
        ParmLocal[2] = nth  # n_0 - thermal electron density, cm^{-3}
        ParmLocal[3] = Bmag  # B - magnetic field, G
        ParmLocal[6] = 3  # distribution over energy (PLW is chosen)
        ParmLocal[7] = nrl  # n_b - nonthermal electron density, cm^{-3}
        ParmLocal[9] = Emin  # E_min, MeV
        ParmLocal[10] = Emax  # E_max, MeV
        ParmLocal[12] = delta  # \delta_1
        ParmLocal[14] = 0  # distribution over pitch-angle (isotropic is chosen)
        ParmLocal[15] = 90  # loss-cone boundary, degrees
        ParmLocal[16] = 0.2  # \Delta\mu

        Parms = np.zeros((24, NSteps), dtype='double', order='F')  # 2D array of input parameters - for multiple voxels
        for i in range(NSteps):
            Parms[:, i] = ParmLocal  # most of the parameters are the same in all voxels
            Parms[4, i] = theta

        RL = np.zeros((7, Nf), dtype='double', order='F')  # input/output array
        dummy = np.array(0, dtype='double')

        # calculating the emission for array distribution (array -> on)
        res = GET_MW(Lparms, Rparms, Parms, dummy, dummy, dummy, RL)
        if res:
            # retrieving the results
            f = RL[0]
            I_L = RL[5]
            I_R = RL[6]
            all_zeros = not RL.any()
            if not all_zeros:
                flux_model = I_L + I_R
                flux_model = np.nan_to_num(flux_model) + 1e-11
                logf = np.log10(f)
                logflux_model = np.log10(flux_model)
                logfreqghz = np.log10(freqghz)
                interpfunc = interpolate.interp1d(logf, logflux_model, kind='linear')
                logmflux = interpfunc(logfreqghz)
                mflux = 10. ** logmflux
                mtb = sfu2tb(np.array(freqghz) * 1.e9, mflux, area=src_area).value

                if pgplot_widget:
                    ##todo: figure out a way to update the main widget
                    import pyqtgraph as pg
                    all_items = pgplot_widget.getPlotItem().listDataItems()
                    if len(all_items) > 0:
                        pgplot_widget.removeItem(all_items[-1])
                    spec_fitplot = pgplot_widget.plot(x=np.log10(freqghz), y=np.log10(mtb),
                                                      pen=dict(color=pg.mkColor(0), width=3),
                                                      symbol=None, symbolBrush=None)
                    pgplot_widget.addItem(spec_fitplot)

                if show_plot:
                    import matplotlib.pyplot as plt
                    fig, (ax1, ax2) = plt.subplots(1, 2)
                    ax1.plot(freqghz, mflux, 'k')
                    # ax1.set_xlim([1, 20])
                    ax1.set_xlabel('Frequency (GHz)')
                    ax1.set_ylabel('Flux (sfu)')
                    ax1.set_title('Flux Spectrum')
                    ax1.set_xscale('log')
                    ax1.set_yscale('log')
                    ax2.plot(freqghz, mtb, 'k')
                    # ax2.set_xlim([1, 20])
                    ax2.set_xlabel('Frequency (GHz)')
                    ax2.set_ylabel('Brightness Temperature (K)')
                    ax2.set_title('Brightness Temperature Spectrum')
                    ax2.set_xscale('log')
                    ax2.set_yscale('log')
                    ax2.legend()
                    plt.show()
            else:
                print("Calculation error! Assign an unrealistically huge number")
                mflux = np.ones_like(freqghz) * 1e4
                mtb = sfu2tb(np.array(freqghz) * 1.e9, mflux, area=src_area).value
        else:
            print("Calculation error! Assign an unrealistically huge number")
            mflux = np.ones_like(freqghz) * 1e9
            mtb = sfu2tb(np.array(freqghz) * 1.e9, mflux, area=src_area).value

        # Return values
        if spec_in_tb:
            if spec is None:
                # nothing is provided, return the model spectrum
                return mtb
            if spec_err is None:
                # no uncertainty provided, return absolute residual
                return mtb - spec
            # Return scaled residual
            return (mtb - spec) / spec_err
        else:
            if spec is None:
                # nothing is provided, return the model spectrum
                return mflux
            if spec_err is None:
                # no uncertainty provided, return absolute residual
                return mflux - spec
            # Return scaled residual
            return (mflux - spec) / spec_err
