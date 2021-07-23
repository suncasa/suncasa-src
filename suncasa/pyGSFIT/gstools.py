def ff_emission(em, T=1.e7, Z=1., mu=1.e10):
    from astropy import constants as const
    import astropy.units as u
    import numpy as np

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

class GSCostFunctions:
    def SinglePowerLawMinimizerOneSrc(fit_params, freqghz, tb=None, tb_err=None, showplt=False):
        from suncasa.pyGSFIT import GScodes
        import os, sys
        import numpy as np
        import math
        libname = os.path.join(os.path.dirname(GScodes.__file__),
                       'binaries/MWTransferArr.so')
        # params are defined by lmfit.Paramters()
        '''
        params: parameters defined by lmfit.Paramters()
        freqghz: frequencies in GHz
        ssz: pixel size in arcsec
        tb: reference brightness temperature in K
        tb_err: uncertainties of reference brightness temperature in K
        '''

        from scipy import interpolate
        GET_MW = GScodes.initGET_MW(libname)  # load the library

        #ssz = float(fit_params['ssz'].value)  # source area in arcsec^2
        ssz = 10. ** 2 ##TODO currently assuming 10 arcsec^2
        ssz_cm2 = ssz * (0.725e8) ** 2. # now in cm2
        depth = float(fit_params['depth_Mm'].value) * 1e8  # total source depth in cm
        Bmag = float(fit_params['Bx100G'].value) * 100.  # magnetic field strength in G
        Tth = float(fit_params['T_MK'].value) * 1e6  # thermal temperature in K
        nth = 10. ** float(fit_params['log_nth'].value)  # thermal density
        nrl = 10. ** float(fit_params['log_nnth'].value)  # total nonthermal density above E_min
        delta = float(fit_params['delta'].value)  # powerlaw index
        theta = float(fit_params['theta'].value)  # viewing angle in degrees
        Emin = float(fit_params['Emin_keV'].value) / 1e3  # low energy cutoff of nonthermal electrons in MeV
        Emax = float(fit_params['Emax_MeV'].value)  # high energy cutoff of nonthermal electrons in MeV
        #E_hi = 0.1
        #nrl = nrlh * (Emin ** (1. - delta) - Emax * (1. - delta)) / (E_hi ** (1. - delta) - Emax ** (1. - delta))

        Nf = 100  # number of frequencies
        NSteps = 1  # number of nodes along the line-of-sight

        Lparms = np.zeros(11, dtype='int32')  # array of dimensions etc.
        Lparms[0] = NSteps
        Lparms[1] = Nf

        Rparms = np.zeros(5, dtype='double')  # array of global floating-point parameters
        Rparms[0] = ssz_cm2  # Area, cm^2
        Rparms[1] = 0.8e9  # starting frequency to calculate spectrum, Hz
        Rparms[2] = 0.02  # logarithmic step in frequency
        Rparms[3] = 0  # f^C
        Rparms[4] = 0  # f^WH

        ParmLocal = np.zeros(24, dtype='double')  # array of voxel parameters - for a single voxel
        ParmLocal[0] = depth / NSteps  # voxel depth, cm
        ParmLocal[1] = Tth  # T_0, K
        ParmLocal[2] = nth  # n_0 - thermal electron density, cm^{-3}
        ParmLocal[3] = Bmag  # B - magnetic field, G
        ParmLocal[6] = 3     # distribution over energy (PLW is chosen)
        ParmLocal[7] = nrl   # n_b - nonthermal electron density, cm^{-3}
        ParmLocal[9] = Emin   # E_min, MeV
        ParmLocal[10] = Emax # E_max, MeV
        ParmLocal[12] = delta  # \delta_1
        ParmLocal[14] = 0    # distribution over pitch-angle (isotropic is chosen)
        ParmLocal[15] = 90   # loss-cone boundary, degrees
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