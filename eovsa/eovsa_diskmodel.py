import numpy as np
from tqdm import tqdm
import os


# from taskinit import ms, tb, qa
# from taskinit import iatool as ia
# from taskinit import cltool as cl
# from ft_cli import ft_cli as ft
# from delmod_cli import delmod_cli as delmod
# from uvsub_cli import uvsub_cli as uvsub


def read_ms(vis):
    ''' Read a CASA ms file and return a dictionary of amplitude, phase, uvdistance,
        uvangle, frequency (GHz) and time (MJD).  Currently only returns the XX IF channel.

        vis     Name of the visibility (ms) folder
    '''
    ms.open(vis)
    spwinfo = ms.getspectralwindowinfo()
    nspw = len(spwinfo.keys())
    for i in range(nspw):
        print('Working on spw', i)
        ms.selectinit(datadescid=0, reset=True)
        ms.selectinit(datadescid=i)
        if i == 0:
            spw = ms.getdata(['amplitude', 'phase', 'u', 'v', 'axis_info'], ifraxis=True)
            xxamp = spw['amplitude']
            xxpha = spw['phase']
            fghz = spw['axis_info']['freq_axis']['chan_freq'][:, 0] / 1e9
            band = np.ones_like(fghz) * i
            mjd = spw['axis_info']['time_axis']['MJDseconds'] / 86400.
            uvdist = np.sqrt(spw['u'] ** 2 + spw['v'] ** 2)
            uvang = np.angle(spw['u'] + 1j * spw['v'])
        else:
            spw = ms.getdata(['amplitude', 'phase', 'axis_info'], ifraxis=True)
            xxamp = np.concatenate((xxamp, spw['amplitude']), 1)
            xxpha = np.concatenate((xxpha, spw['phase']), 1)
            fg = spw['axis_info']['freq_axis']['chan_freq'][:, 0] / 1e9
            fghz = np.concatenate((fghz, fg))
            band = np.concatenate((band, np.ones_like(fg) * i))
    ms.close()
    return {'amp': xxamp, 'phase': xxpha, 'fghz': fghz, 'band': band, 'mjd': mjd, 'uvdist': uvdist, 'uvangle': uvang}


def fit_diskmodel(out, bidx, rstn_flux, uvfitrange=[1, 150], angle_tolerance=np.pi / 2, doplot=True):
    ''' Given the result returned by read_ms(), plots the amplitude vs. uvdistance
        separately for polar and equatorial directions rotated for P-angle, then overplots
        a disk model for a disk enlarged by eqfac in the equatorial direction, and polfac
        in the polar direction.  Also requires the RSTN flux spectrum for the date of the ms,
        determined from (example for 2019-09-01):
           import rstn
           frq, flux = rstn.rd_rstnflux(t=Time('2019-09-01'))
           rstn_flux = rstn.rstn2ant(frq, flux, out['fghz']*1000, t=Time('2019-09-01'))

    '''
    from util import bl2ord, lobe
    import matplotlib.pylab as plt
    import sun_pos
    from scipy.special import j1
    import scipy.constants
    mperns = scipy.constants.c / 1e9  # speed of light in m/ns
    # Rotate uv angle for P-angle
    pa, b0, r = sun_pos.get_pb0r(out['mjd'][0], arcsec=True)
    uvangle = lobe(out['uvangle'] - pa * np.pi / 180.)
    a = 2 * r * np.pi ** 2 / (180. * 3600.)  # Initial scale for z, uses photospheric radius of the Sun
    if doplot: f, ax = plt.subplots(3, 1)
    uvmin, uvmax = uvfitrange
    uvdeq = []
    uvdpol = []
    ampeq = []
    amppol = []
    zeq = []
    zpol = []
    # Loop over antennas 1-4
    antmax = 7
    at = angle_tolerance
    for i in range(4):
        fidx, = np.where(out['band'] == bidx)  # Array of frequency indexes for channels in this band
        for j, fi in enumerate(fidx):
            amp = out['amp'][0, fi, bl2ord[i, i + 1:antmax]].flatten() / 10000.  # Convert to sfu
            # Use only non-zero amplitudes
            good, = np.where(amp != 0)
            amp = amp[good]
            uva = uvangle[bl2ord[i, i + 1:antmax]].flatten()[good]
            # Equatorial points are within +/- pi/8 of solar equator
            eq, = np.where(np.logical_or(np.abs(uva) < at / 2, np.abs(uva) >= np.pi - at / 2))
            # Polar points are within +/- pi/8 of solar pole
            pol, = np.where(np.logical_and(np.abs(uva) >= np.pi / 2 - at / 2, np.abs(uva) < np.pi / 2 + at / 2))
            uvd = out['uvdist'][bl2ord[i, i + 1:antmax]].flatten()[good] * out['fghz'][fi] / mperns  # Wavelengths
            # Add data for this set of baselines to global arrays
            uvdeq.append(uvd[eq])
            uvdpol.append(uvd[pol])
            ampeq.append(amp[eq])
            amppol.append(amp[pol])
            zeq.append(uvd[eq])
            zpol.append(uvd[pol])
    uvdeq = np.concatenate(uvdeq)
    uvdpol = np.concatenate(uvdpol)
    uvdall = np.concatenate((uvdeq, uvdpol))
    ampeq = np.concatenate(ampeq)
    amppol = np.concatenate(amppol)
    ampall = np.concatenate((ampeq, amppol))
    zeq = np.concatenate(zeq)
    zpol = np.concatenate(zpol)
    zall = np.concatenate((zeq, zpol))
    # These indexes are for a restricted uv-range to be fitted
    ieq, = np.where(np.logical_and(uvdeq > uvmin, uvdeq <= uvmax))
    ipol, = np.where(np.logical_and(uvdpol > uvmin, uvdpol <= uvmax))
    iall, = np.where(np.logical_and(uvdall > uvmin, uvdall <= uvmax))
    if doplot:
        # Plot all of the data points
        ax[0].plot(uvdeq, ampeq, 'k+')
        ax[1].plot(uvdpol, amppol, 'k+')
        ax[2].plot(uvdall, ampall, 'k+')
        # Overplot the fitted data points in a different color
        ax[0].plot(uvdeq[ieq], ampeq[ieq], 'b+')
        ax[1].plot(uvdpol[ipol], amppol[ipol], 'b+')
        ax[2].plot(uvdall[iall], ampall[iall], 'b+')
    # Minimize ratio of points to model
    ntries = 300
    solfac = np.linspace(1.0, 1.3, ntries)
    d2m_eq = np.zeros(ntries, np.float)
    d2m_pol = np.zeros(ntries, np.float)
    d2m_all = np.zeros(ntries, np.float)
    sfac = np.zeros(ntries, np.float)
    sfacall = np.zeros(ntries, np.float)
    # Loop over ntries (300) models of solar disk size factor ranging from 1.0 to 1.3 r_Sun
    for k, sizfac in enumerate(solfac):
        eqpts = rstn_flux[fidx][0] * 2 * np.abs(j1(a * sizfac * zeq[ieq]) / (a * sizfac * zeq[ieq]))
        polpts = rstn_flux[fidx[0]] * 2 * np.abs(j1(a * sizfac * zpol[ipol]) / (a * sizfac * zpol[ipol]))
        sfac[k] = (np.nanmedian(ampeq[ieq] / eqpts) + np.nanmedian(amppol[ipol] / polpts)) / 2
        eqpts = rstn_flux[fidx[0]] * (2 * sfac[k]) * np.abs(j1(a * sizfac * zeq[ieq]) / (a * sizfac * zeq[ieq]))
        polpts = rstn_flux[fidx[0]] * (2 * sfac[k]) * np.abs(j1(a * sizfac * zpol[ipol]) / (a * sizfac * zpol[ipol]))
        allpts = rstn_flux[fidx[0]] * (2 * sfac[k]) * np.abs(j1(a * sizfac * zall[iall]) / (a * sizfac * zall[iall]))
        sfacall[k] = np.nanmedian(ampall[iall] / allpts)
        d2m_eq[k] = np.nanmedian(abs(ampeq[ieq] / eqpts - 1))
        d2m_pol[k] = np.nanmedian(abs(amppol[ipol] / polpts - 1))
        d2m_all[k] = np.nanmedian(abs(ampall[iall] / allpts - 1))
    keq = np.argmin(d2m_eq)
    kpol = np.argmin(d2m_pol)
    kall = np.argmin(d2m_all)
    eqradius = solfac[keq] * r
    polradius = solfac[kpol] * r
    allradius = solfac[kall] * r
    sfactor = sfac[keq]
    sfall = sfacall[kall]
    sflux = sfall * rstn_flux[fidx[0]]
    if doplot:
        z = np.linspace(1.0, 1000.0, 10000)
        # Overplot the best fit
        ax[0].plot(z, rstn_flux[fidx[0]] * (2 * sfactor) * np.abs(j1(a * solfac[keq] * z) / (a * solfac[keq] * z)))
        ax[1].plot(z, rstn_flux[fidx[0]] * (2 * sfactor) * np.abs(j1(a * solfac[kpol] * z) / (a * solfac[kpol] * z)))
        ax[2].plot(z, rstn_flux[fidx[0]] * (2 * sfall) * np.abs(j1(a * solfac[kall] * z) / (a * solfac[kall] * z)))
        # ax[1].plot(zpol,polpts,'y.')
        ax[0].set_title(
            str(out['fghz'][fidx][0])[:4] + 'GHz. R_eq:' + str(eqradius)[:6] + '". R_pol' + str(polradius)[:6]
            + '". R_all' + str(allradius)[:6] + '". Flux scl fac:' + str(sfall)[:4])
        # ax[0].plot(uvdeq,ampeq/eqpts,'k+')
        # ax[0].plot([0,1000],np.array([1,1])*np.nanmedian(ampeq/eqpts))
        # ax[1].plot(uvdpol,amppol/polpts,'k+')
        # ax[1].plot([0,1000],np.array([1,1])*np.nanmedian(amppol/polpts))
        for i in range(3):
            ax[i].set_xlim(0, 1000)
            ax[i].set_ylim(0.01, rstn_flux[fidx[0]] * 2 * sfactor)
            ax[i].set_yscale('log')
            ax[2].set_xlabel('UV Distance (wavelengths)')
            ax[i].set_ylabel('Amplitude (sfu)')
            ax[i].text(850, 125, ['Equator', 'Pole', 'All'][i])
    return bidx, out['fghz'][fidx[0]], eqradius, polradius, allradius, sfall, sflux


def fit_vs_freq(out):
    import matplotlib.pylab as plt
    import rstn
    from astropy.time import Time
    t = Time(out['mjd'][0], format='mjd')
    frq, flux = rstn.rd_rstnflux(t=t)
    rstn_flux = rstn.rstn2ant(frq, flux, out['fghz'] * 1000, t=t)
    band = []
    fghz = []
    eqrad = []
    polrad = []
    allrad = []
    sfac = []
    sflux = []
    for i in range(50):
        uvfitrange = np.array([10, 150]) + np.array([1, 18]) * i
        a, b, c, d, e, f, g = fit_diskmodel(out, i, rstn_flux, uvfitrange=uvfitrange, angle_tolerance=np.pi / 2,
                                            doplot=False)
        band.append(a)
        fghz.append(b)
        eqrad.append(c)
        polrad.append(d)
        allrad.append(e)
        sfac.append(f)
        sflux.append(g)
        if (i % 10) == 0: print(i)
    result = {'band': np.array(band), 'fghz': np.array(fghz), 'eqradius': np.array(eqrad),
              'polradius': np.array(polrad),
              'radius': np.array(allrad), 'flux_correction_factor': np.array(sfac), 'disk_flux': np.array(sflux) * 2.}
    plt.figure()
    plt.plot(result['fghz'], result['eqradius'], 'o', label='Equatorial Radius')
    plt.plot(result['fghz'], result['polradius'], 'o', label='Polar Radius')
    plt.plot(result['fghz'], result['radius'], 'o', label='Circular Radius')
    plt.legend()
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Radius [arcsec]')
    plt.title('Frequency-dependent Solar Disk Size for 2019-Sep-01')
    return result


def diskmodel(outname='disk', bdwidth='325MHz', direction='J2000 10h00m00.0s 20d00m00.0s',
              reffreq='2.8GHz', flux=660000.0, eqradius='16.166arcmin', polradius='16.166arcmin',
              pangle='21.1deg', index=None, cell='2.0arcsec'):
    ''' Create a blank solar disk model image (or optionally a data cube)

        outname       String to use for part of the image and fits file names (default 'disk')
        direction     String specifying the position of the Sun in RA and Dec.  Default
                        means use the standard string "J2000 10h00m00.0s 20d00m00.0s"
        reffreq       The reference frequency to use for the disk model (the frequency at which
                        the flux level applies). Default is '2.8GHz'.
        flux          The flux density, in Jy, for the entire disk. Default is 66 sfu.
        eqradius      The equatorial radius of the disk.  Default is
                        16 arcmin + 10" (for typical extension of the radio limb)
        polradius     The polar radius of the disk.  Default is
                        16 arcmin + 10" (for typical extension of the radio limb)
        pangle        The solar P-angle (geographic position of the N-pole of the Sun) in
                        degrees E of N.  This only matters if eqradius != polradius
        index         The spectral index to use at other frequencies.  Default None means
                        use a constant flux density for all frequencies.
        cell          The cell size (assumed square) to use for the image.  The image size
                        is determined from a standard radius of 960" for the Sun, divided by
                        cell size, increased to nearest power of 512 pixels. The default is '2.0arcsec',
                        which results in an image size of 1024 x 1024.
        Note that the frequency increment used is '325MHz', which is the width of EOVSA bands
          (not the width of individual science channels)
    '''
    cl.done()
    try:
        aspect = 1.01  # Enlarge the equatorial disk by 1%
        eqradius = qa.quantity(eqradius)
        diamajor = qa.quantity(2 * aspect * eqradius['value'], eqradius['unit'])
        polradius = qa.quantity(polradius)
        diaminor = qa.quantity(2 * polradius['value'], polradius['unit'])
        solrad = qa.convert(polradius, 'arcsec')
    except:
        print('Radius', eqradius, polradius,
              'does not have the expected format, number + unit where unit is arcmin or arcsec')
        return
    try:
        cell = qa.convert(qa.quantity(cell), 'arcsec')
        cellsize = float(cell['value'])
        diskpix = solrad['value'] * 2 / cellsize
        cell_rad = qa.convert(cell, 'rad')
    except:
        print('Cell size', cell, 'does not have the expected format, number + unit where unit is arcmin or arcsec')
        return

    # Add 90 degrees to pangle, due to angle definition in addcomponent() -- it puts the majoraxis vertical
    pangle = qa.add(qa.quantity(pangle), qa.quantity('90deg'))
    mapsize = ((int(diskpix) / 512) + 1) * 512
    # Flux density is doubled because it is split between XX and YY
    cl.addcomponent(dir=direction, flux=flux * 2, fluxunit='Jy', freq=reffreq, shape='disk',
                    majoraxis=diamajor, minoraxis=diaminor, positionangle=pangle)
    cl.setrefdirframe(0, 'J2000')
    diskim = outname + reffreq + '.im'
    if os.path.exists(diskim):
        os.system('rm -rf {}'.format(diskim))
    ia.fromshape(diskim, [mapsize, mapsize, 1, 1], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad_val = cell_rad['value']
    cs.setincrement([-cell_rad_val, cell_rad_val], 'direction')
    epoch, ra, dec = direction.split()
    cs.setreferencevalue([qa.convert(ra, 'rad')['value'], qa.convert(dec, 'rad')['value']], type="direction")
    cs.setreferencevalue(reffreq, 'spectral')
    cs.setincrement(bdwidth, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()
    ia.done()
    cl.close()
    cl.done()
    return diskim


def insertdiskmodel(vis):
    fdens = np.array([891282, 954570, 1173229, 1245433, 1373730, 1506802,
                      1613253, 1702751, 1800721, 1946756, 2096020, 2243951,
                      2367362, 2525968, 2699795, 2861604, 3054829, 3220450,
                      3404182, 3602625, 3794312, 3962926, 4164667, 4360683,
                      4575677, 4767210, 4972824, 5211717, 5444632, 5648266,
                      5926634, 6144249, 6339863, 6598018, 6802707, 7016012,
                      7258929, 7454951, 7742816, 7948976, 8203206, 8411834,
                      8656720, 8908130, 9087766, 9410760, 9571365, 9827078,
                      10023598, 8896671])
    dsize = np.array(['1228.0arcsec', '1194.0arcsec', '1165.0arcsec', '1139.0arcsec', '1117.0arcsec',
                      '1097.0arcsec', '1080.0arcsec', '1065.0arcsec', '1053.0arcsec', '1042.0arcsec', '1033.0arcsec',
                      '1025.0arcsec', '1018.0arcsec', '1012.0arcsec',
                      '1008.0arcsec', '1003.0arcsec', '1000.0arcsec', '997.0arcsec', '994.0arcsec', '992.0arcsec',
                      '990.0arcsec', '988.0arcsec', '986.0arcsec', '985.0arcsec', '983.0arcsec', '982.0arcsec',
                      '980.0arcsec', '979.0arcsec', '978.0arcsec', '976.0arcsec', '975.0arcsec', '974.0arcsec',
                      '972.0arcsec', '971.0arcsec', '970.0arcsec', '969.0arcsec', '968.0arcsec', '967.0arcsec',
                      '966.0arcsec', '965.0arcsec', '964.0arcsec', '964.0arcsec', '963.0arcsec', '962.0arcsec',
                      '962.0arcsec', '961.0arcsec', '960.0arcsec', '959.0arcsec', '957.0arcsec', '956.0arcsec'])

    msfile = vis
    ms.open(msfile)
    diskim = []
    ms.open(msfile)
    spwinfo = ms.getspectralwindowinfo()
    nspw = len(spwinfo.keys())
    ms.close()
    diskimdir = 'diskim/'
    if not os.path.exists(diskimdir):
        os.makedirs(diskimdir)
    frq = []
    for sp in range(nspw):
        spw = spwinfo[str(sp)]
        frq.append('{:.4f}GHz'.format((spw['RefFreq'] + spw['TotalWidth'] / 2.0) / 1e9))
    frq = np.array(frq)
    tb.open(msfile + '/FIELD')
    phadir = tb.getcol('PHASE_DIR').flatten()
    tb.close()
    ra = phadir[0]
    dec = phadir[1]
    direction = 'J2000 ' + str(ra) + 'rad ' + str(dec) + 'rad'

    for sp in tqdm(range(nspw),desc='Generating {} disk models'.format(nspw)):
        diskim.append(
            diskmodel(outname=diskimdir+'disk{:02d}_'.format(sp), bdwidth=spwinfo[str(sp)], direction=direction, reffreq=frq[sp],
                      flux=fdens[sp], eqradius=dsize[sp], polradius=dsize[sp]))

    delmod(msfile, otf=True, scr=True)
    for sp in tqdm(range(nspw),desc='Inserting disk model'):
        ft(vis=msfile, spw=str(sp), field='', model=str(diskim[sp]), nterms=1,
           reffreq="", complist="", incremental=False, usescratch=True)

    uvsub(vis=msfile)
    return msfile
