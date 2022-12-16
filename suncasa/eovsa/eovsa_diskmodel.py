from tqdm import tqdm
from taskinit import ms, tb, qa
from taskinit import iatool
from taskinit import cltool
from delmod_cli import delmod_cli as delmod
from clearcal_cli import clearcal_cli as clearcal
from suncasa.utils import mstools as mstl
from suncasa.utils import helioimage2fits as hf
import shutil, os
import sunpy.coordinates.ephemeris as eph
import numpy as np
from gaincal_cli import gaincal_cli as gaincal
from applycal_cli import applycal_cli as applycal
from flagdata_cli import flagdata_cli as flagdata
from flagmanager_cli import flagmanager_cli as flagmanager
from uvsub_cli import uvsub_cli as uvsub
from split_cli import split_cli as split
from tclean_cli import tclean_cli as tclean
from ft_cli import ft_cli as ft
from suncasa.utils import mstools as mstl


# def ant_trange(vis):
#     ''' Figure out nominal times for tracking of old EOVSA antennas, and return time
#         range in CASA format
#     '''
#     import eovsa_array as ea
#     from astropy.time import Time
#     # Get the Sun transit time, based on the date in the vis file name (must have UDByyyymmdd in the name)
#     aa = ea.eovsa_array()
#     date = vis.split('UDB')[-1][:8]
#     slashdate = date[:4] + '/' + date[4:6] + '/' + date[6:8]
#     aa.date = slashdate
#     sun = aa.cat['Sun']
#     mjd_transit = Time(aa.next_transit(sun).datetime(), format='datetime').mjd
#     # Construct timerange based on +/- 3h55m from transit time (when all dishes are nominally tracking)
#     trange = Time(mjd_transit - 0.1632, format='mjd').iso[:19] + '~' + Time(mjd_transit + 0.1632, format='mjd').iso[:19]
#     trange = trange.replace('-', '/').replace(' ', '/')
#     return trange

def ant_trange(vis):
    ''' Figure out nominal times for tracking of old EOVSA antennas, and return time
        range in CASA format
    '''
    import eovsa_array as ea
    from astropy.time import Time
    from taskinit import ms

    # Get timerange from the visibility file
    # msinfo = dict.fromkeys(['vis', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs', 'observatory'])
    ms.open(vis)
    # metadata = ms.metadata()
    scans = ms.getscansummary()
    sk = np.sort(scans.keys())
    vistrange = np.array([scans[sk[0]]['0']['BeginTime'], scans[sk[-1]]['0']['EndTime']])

    # Get the Sun transit time, based on the date in the vis file name (must have UDByyyymmdd in the name)
    aa = ea.eovsa_array()
    date = vis.split('UDB')[-1][:8]
    slashdate = date[:4] + '/' + date[4:6] + '/' + date[6:8]
    aa.date = slashdate
    sun = aa.cat['Sun']
    mjd_transit = Time(aa.next_transit(sun).datetime(), format='datetime').mjd
    # Construct timerange limits based on +/- 3h55m from transit time (when all dishes are nominally tracking)
    # and clip the visibility range not to exceed those limits
    mjdrange = np.clip(vistrange, mjd_transit - 0.1632, mjd_transit + 0.1632)
    trange = Time(mjdrange[0], format='mjd').iso[:19] + '~' + Time(mjdrange[1], format='mjd').iso[:19]
    trange = trange.replace('-', '/').replace(' ', '/')
    return trange


def gaussian2d(x, y, amplitude, x0, y0, sigma_x, sigma_y, theta):
    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = amplitude * np.exp(- (a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) ** 2)))
    return g


def writediskxml(dsize, fdens, freq, xmlfile='SOLDISK.xml'):
    import xml.etree.ElementTree as ET
    # create the file structure
    sdk = ET.Element('SOLDISK')
    sdk_dsize = ET.SubElement(sdk, 'item')
    sdk_fdens = ET.SubElement(sdk, 'item')
    sdk_freqs = ET.SubElement(sdk, 'item')
    sdk_dsize.set('disk_size', ','.join(dsize))
    sdk_fdens.set('flux_dens', ','.join(['{:.1f}Jy'.format(s) for s in fdens]))
    sdk_freqs.set('freq', ','.join(freq))

    # create a new XML file with the results
    mydata = ET.tostring(sdk)
    if os.path.exists(xmlfile):
        os.system('rm -rf {}'.format(xmlfile))
    with open(xmlfile, 'w') as sf:
        sf.write(mydata)
    return xmlfile


def readdiskxml(xmlfile):
    import astropy.units as u
    import xml.etree.ElementTree as ET
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    diskinfo = {}
    for elem in root:
        d = elem.attrib
        for k, v in d.items():
            v_ = v.split(',')
            v_ = [u.Unit(f).to_string().split(' ') for f in v_]
            diskinfo[k] = []
            for val, uni in v_:
                diskinfo[k].append(float(val))
            diskinfo[k] = np.array(diskinfo[k]) * u.Unit(uni)

    return diskinfo


def image_adddisk(eofile, diskinfo, edgeconvmode='frommergeddisk', caltbonly=False):
    '''
    :param eofile:
    :param diskxmlfile:
    :param edgeconvmode: available mode: frommergeddisk,frombeam
    :return:
    '''

    from sunpy import map as smap
    from suncasa.utils import plot_mapX as pmX
    from scipy import constants
    import astropy.units as u
    from sunpy import io as sio
    dsize = diskinfo['disk_size']
    fdens = diskinfo['flux_dens']
    freqs = diskinfo['freq']
    eomap = smap.Map(eofile)
    eomap_ = pmX.Sunmap(eomap)
    header = eomap.meta
    bmaj = header['bmaj'] * 3600 * u.arcsec
    bmin = header['bmin'] * 3600 * u.arcsec
    cell = (header['cdelt1'] * u.Unit(header['cunit1']) + header['cdelt2'] * u.Unit(header['cunit2'])) / 2.0
    bmsize = (bmaj + bmin) / 2.0
    data = eomap.data  # remember the data order is reversed due to the FITS convension
    keys = header.keys()
    values = header.values()
    mapx, mapy = eomap_.map2wcsgrids(cell=False)
    mapx = mapx[:-1, :-1]
    mapy = mapy[:-1, :-1]
    rdisk = np.sqrt(mapx ** 2 + mapy ** 2)

    k_b = constants.k
    c_l = constants.c
    const = 2. * k_b / c_l ** 2
    pix_area = (cell.to(u.rad).value) ** 2
    jy_to_si = 1e-26
    factor2 = 1.
    faxis = keys[values.index('FREQ')][-1]

    if caltbonly:
        edgeconvmode = ''
    if edgeconvmode == 'frommergeddisk':

        nul = header['CRVAL' + faxis] + header['CDELT' + faxis] * (1 - header['CRPIX' + faxis])
        nuh = header['CRVAL' + faxis] + header['CDELT' + faxis] * (header['NAXIS' + faxis] - header['CRPIX' + faxis])
        ## get the frequency range of the image
        nu_bound = (np.array([nul, nuh]) + 0.5 * np.array([-1, 1]) * header['CDELT' + faxis]) * u.Unit(
            header['cunit' + faxis])
        nu_bound = nu_bound.to(u.GHz)
        ## get the frequencies of the disk models
        fidxs = np.logical_and(freqs > nu_bound[0], freqs < nu_bound[1])
        ny, nx = rdisk.shape
        freqs_ = freqs[fidxs]
        fdens_ = fdens[fidxs] / 2.0  # divide by 2 because fdens is 2x solar flux density
        dsize_ = dsize[fidxs]
        fdisk_ = np.empty((len(freqs_), ny, nx))
        fdisk_[:] = np.nan
        for fidx, freq in enumerate(freqs_):
            fdisk_[fidx, ...][rdisk <= dsize_[fidx].value] = 1.0
            # nu = header['CRVAL' + faxis] + header['CDELT' + faxis] * (1 - header['CRPIX' + faxis])
            factor = const * freq.to(u.Hz).value ** 2  # SI unit
            jy2tb = jy_to_si / pix_area / factor * factor2
            fdisk_[fidx, ...] = fdisk_[fidx, ...] / np.nansum(fdisk_[fidx, ...]) * fdens_[fidx].value
            fdisk_[fidx, ...] = fdisk_[fidx, ...] * jy2tb
        #         # fdisk_[np.isnan(fdisk_)] = 0.0
        tbdisk = np.nanmean(fdisk_, axis=0)
        tbdisk[np.isnan(tbdisk)] = 0.0

        sig2fwhm = 2.0 * np.sqrt(2 * np.log(2))
        x0, y0 = 0, 0
        sigx, sigy = bmaj.value / sig2fwhm, bmin.value / sig2fwhm
        theta = -(90.0 - header['bpa']) * u.deg
        x = (np.arange(31) - 15) * cell.value
        y = (np.arange(31) - 15) * cell.value
        x, y = np.meshgrid(x, y)
        kernel = gaussian2d(x, y, 1.0, x0, y0, sigx, sigy, theta.to(u.radian).value)
        kernel = kernel / np.nansum(kernel)
        from scipy import signal
        tbdisk = signal.fftconvolve(tbdisk, kernel, mode='same')
    else:
        nu = header['CRVAL' + faxis] + header['CDELT' + faxis] * (1 - header['CRPIX' + faxis])
        freqghz = nu / 1.0e9
        factor = const * nu ** 2  # SI unit
        jy2tb = jy_to_si / pix_area / factor * factor2
        p_dsize = np.poly1d(np.polyfit(freqs.value, dsize.value, 15))
        p_fdens = np.poly1d(
            np.polyfit(freqs.value, fdens.value, 15)) / 2.  # divide by 2 because fdens is 2x solar flux density
        if edgeconvmode == 'frombeam':
            from scipy.special import erfc
            factor_erfc = 2.0  ## erfc function ranges from 0 to 2
            fdisk = erfc((rdisk - p_dsize(freqghz)) / bmsize.value) / factor_erfc
        else:
            fdisk = np.zeros_like(rdisk)
            fdisk[rdisk <= p_dsize(freqghz)] = 1.0
        fdisk = fdisk / np.nansum(fdisk) * p_fdens(freqghz)
        tbdisk = fdisk * jy2tb

    tb_disk = np.nanmax(tbdisk)
    if caltbonly:
        return tb_disk
    else:
        datanew = data + tbdisk
        # datanew[np.isnan(data)] = 0.0
        header['TBDISK'] = tb_disk
        header['TBUNIT'] = 'K'
        eomap_disk = smap.Map(datanew, header)
        nametmp = eofile.split('.')
        nametmp.insert(-1, 'disk')
        outfits = '.'.join(nametmp)
        datanew = datanew.astype(np.float32)
        if os.path.exists(outfits):
            os.system('rm -rf {}'.format(outfits))
        sio.write_file(outfits, datanew, header)
        return eomap_disk, tb_disk, outfits


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


def im2cl(imname, clname, convol=True, verbose=False):
    if os.path.exists(clname):
        os.system('rm -rf {}'.format(clname))
    ia = iatool()
    ia.open(imname)
    ia2 = iatool()
    ia2.open(imname.replace('.model', '.image'))
    bm = ia2.restoringbeam()
    bmsize = (qa.convert(qa.quantity(bm['major']), 'arcsec')['value'] +
              qa.convert(qa.quantity(bm['minor']), 'arcsec')['value']) / 2.0
    if convol:
        im2 = ia.sepconvolve(types=['gaussian', 'gaussian'], widths="{0:}arcsec {0:}arcsec".format(2.5*bmsize),
                             overwrite=True)
        ia2.done()
    else:
        im2 = ia
    cl = cltool()
    srcs = im2.findsources(point=False, cutoff=0.3, width=int(np.ceil(bmsize/2.5)))
    # srcs = ia.findsources(point=False, cutoff=0.1, width=5)
    if verbose:
        for k, v in srcs.iteritems():
            if k.startswith('comp'):
                ## note: Stokes I to XX
                print(srcs[k]['flux']['value'])
                # srcs[k]['flux']['value'] = srcs[k]['flux']['value'] / 2.0
    cl.fromrecord(srcs)
    cl.rename(clname)
    cl.done()
    ia.done()
    im2.done()


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


def calc_diskmodel(slashdate, nbands, freq, defaultfreq):
    from astropy.time import Time
    # Default disk size measured for 2019/09/03
    # todo add monthly fitting procedure for the disk size and flux density
    defaultsize = np.array([990.6, 989.4, 988.2, 987.1, 986.0, 984.9, 983.8, 982.7, 981.7, 980.7,
                            979.7, 978.8, 977.8, 976.9, 976.0, 975.2, 974.3, 973.5, 972.7, 972.0,
                            971.2, 970.5, 969.8, 969.1, 968.5, 967.8, 967.2, 966.7, 966.1, 965.6,
                            965.1, 964.6, 964.1, 963.7, 963.3, 962.9, 962.5, 962.1, 961.8, 961.5,
                            961.3, 961.0, 960.8, 960.6, 960.4, 960.2, 960.1, 960.0, 959.9, 959.8])

    # Get current solar distance and modify the default size accordingly
    fac = eph.get_sunearth_distance('2019/09/03') / eph.get_sunearth_distance(slashdate)
    newsize = defaultsize * fac.to_value()
    if nbands == 34:
        if Time(slashdate.replace('/','-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectral windows (bands 4-34 -> spw 0-30)
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            newsize = np.polyval(np.polyfit(defaultfreq, newsize, 5), freq[[0]+range(2,34)])
    dsize = np.array([str(i)[:5] + 'arcsec' for i in newsize], dtype='S12')

    # These are nominal flux densities * 2, determined on 2019/09/03
    defaultfdens = np.array([891282, 954570, 1173229, 1245433, 1373730, 1506802,
                             1613253, 1702751, 1800721, 1946756, 2096020, 2243951,
                             2367362, 2525968, 2699795, 2861604, 3054829, 3220450,
                             3404182, 3602625, 3794312, 3962926, 4164667, 4360683,
                             4575677, 4767210, 4972824, 5211717, 5444632, 5648266,
                             5926634, 6144249, 6339863, 6598018, 6802707, 7016012,
                             7258929, 7454951, 7742816, 7948976, 8203206, 8411834,
                             8656720, 8908130, 9087766, 9410760, 9571365, 9827078,
                             10023598, 8896671])
    fdens = defaultfdens
    if nbands == 34:
        if Time(slashdate.replace('/','-')).mjd < Time('2018-03-13').mjd:
            # Interpolate size to 31 spectal windows (bands 4-34 -> spw 0-30)
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[3:])
        else:
            # Dates between 2018-03-13 have 33 spectral windows
            fdens = np.polyval(np.polyfit(defaultfreq, fdens, 5), freq[[0]+range(2,34)])
    return dsize, fdens


def mk_diskmodel(outname='disk', direction='J2000 10h00m00.0s 20d00m00.0s',
                 reffreq='2.8GHz', flux=660000.0, eqradius='16.166arcmin', polradius='16.166arcmin',
                 pangle='21.1deg', overwrite=True):
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

    diskcl = outname + reffreq + '.cl'
    if os.path.exists(diskcl):
        if overwrite:
            os.system('rm -rf {}'.format(diskcl))
        else:
            return diskcl

    cl = cltool()

    try:
        aspect = 1.01  # Enlarge the equatorial disk by 1%
        eqradius = qa.quantity(eqradius)
        diamajor = qa.quantity(2 * aspect * eqradius['value'], eqradius['unit'])
        polradius = qa.quantity(polradius)
        diaminor = qa.quantity(2 * polradius['value'], polradius['unit'])
    except:
        print('Radius', eqradius, polradius,
              'does not have the expected format, number + unit where unit is arcmin or arcsec')
        return

    # Add 90 degrees to pangle, due to angle definition in addcomponent() -- it puts the majoraxis vertical
    pangle = qa.add(qa.quantity(pangle), qa.quantity('90deg'))
    # Flux density is split between XX and YY
    cl.addcomponent(dir=direction, flux=flux / 2.0, fluxunit='Jy', freq=reffreq, shape='disk',
                    majoraxis=diamajor, minoraxis=diaminor, positionangle=pangle)
    cl.setrefdirframe(0, 'J2000')
    cl.rename(diskcl)
    cl.done()
    return diskcl


def insertdiskmodel(vis, sizescale=1.0, fdens=None, dsize=None, xmlfile='SOLDISK.xml', writediskinfoonly=False,
                    active=False, overwrite=True):
    # Apply size scale adjustment (default is no adjustment)
    for i in range(len(dsize)):
        num, unit = dsize[i].split('arc')
        dsize[i] = str(float(num) * sizescale)[:6] + 'arc' + unit

    msfile = vis

    ms.open(msfile)
    spwinfo = ms.getspectralwindowinfo()
    nspw = len(spwinfo.keys())
    ms.close()
    diskcldir = 'diskcl/'
    if not os.path.exists(diskcldir):
        os.makedirs(diskcldir)
    frq = []
    spws = range(nspw)
    for sp in spws:
        spw = spwinfo[str(sp)]
        frq.append('{:.4f}GHz'.format((spw['RefFreq'] + spw['TotalWidth'] / 2.0) / 1e9))
    frq = np.array(frq)
    writediskxml(dsize, fdens, frq, xmlfile=xmlfile)

    if not writediskinfoonly:
        tb.open(msfile + '/FIELD')
        phadir = tb.getcol('PHASE_DIR').flatten()
        tb.close()
        ra = phadir[0]
        dec = phadir[1]
        direction = 'J2000 ' + str(ra) + 'rad ' + str(dec) + 'rad'

        diskcl = []
        for sp in tqdm(spws, desc='Generating {} disk models'.format(nspw), ascii=True):
            diskcl.append(
                mk_diskmodel(outname=diskcldir + 'disk{:02d}_'.format(sp),
                             direction=direction, reffreq=frq[sp],
                             flux=fdens[sp], eqradius=dsize[sp], polradius=dsize[sp], overwrite=overwrite))

        if not active:
            delmod(msfile, otf=True, scr=True)
            for sp in tqdm(spws, desc='Inserting disk model', ascii=True):
                ft(vis=msfile, spw=str(sp), field='', model="", nterms=1,
                   reffreq="", complist=str(diskcl[sp]), incremental=False, usescratch=True)
        else:
            for sp in tqdm(spws, desc='Inserting disk model', ascii=True):
                model_ft = mstl.getmodel(msfile, spw=str(sp))
                ft(vis=msfile, spw=str(sp), field='', model="", nterms=1,
                   reffreq="", complist=str(diskcl[sp]), incremental=False, usescratch=True)
                model_disk = mstl.getmodel(msfile, spw=str(sp))
                mstl.putmodel(msfile, spw=str(sp), model=model_ft + model_disk)

        return msfile, diskcl


def disk_slfcal(vis, slfcaltbdir='./', active=False, clearcache=False, pols='XX'):
    ''' Starting with the name of a calibrated ms (vis, which must have 'UDByyyymmdd' in the name)
        add a model disk based on the solar disk size for that date and perform multiple selfcal
        adjustments (two phase and one amplitude), and write out a final selfcaled database with
        the disk subtracted.  Returns the name of the final database.
    '''
    trange = ant_trange(vis)

    if vis.endswith('/'):
        vis = vis[:-1]
    # Use vis name to determine date, and hence number of bands
    spw2band = np.array([0, 1] + range(4, 52))
    defaultfreq = 1.1 + 0.325 * (spw2band + 0.5)
    # Calculate the center frequency of each spectral window
    if mstl.get_trange(vis)[0].mjd > 58536:
        # After 2019 Feb 22, the band numbers changed to 1-52, and spw from 0-49
        nbands = 52
        freq = defaultfreq
    else:
        # Before 2019 Feb 22, the band numbers were 1-34, and spw from 0-30
        nbands = 34
        freq = 1.419 + np.arange(nbands) / 2.

    caltbs = []
    slashdate = trange[:10]
    # Verify that the vis is not in the current working directory
    '''
    if os.getcwd() == os.path.dirname(vis):
        print('Cannot copy vis file onto itself.')
        print('Please change to a different working directory')
        return None

    # Copy original ms to local directory
    if os.path.exists(os.path.basename(vis)):
        shutil.rmtree(os.path.basename(vis))
    print('Copy {} to working directory {}.'.format(vis, os.getcwd()))
    shutil.copytree(vis, os.path.basename(vis))
    vis = os.path.basename(vis)
    '''

    if not active:
        clearcal(vis)


    dsize, fdens = calc_diskmodel(slashdate, nbands, freq, defaultfreq)
    diskxmlfile = vis + '.SOLDISK.xml'
    # Insert the disk model (msfile is the same as vis, and will be used as the "original" vis file name)
    msfile, diskcl = insertdiskmodel(vis, dsize=dsize, fdens=fdens, xmlfile=diskxmlfile, active=active)

    flagmanager(msfile, mode='save', versionname='diskslfcal_init')
    ## automaticaly flag any high amplitudes from flares or RFI
    flagdata(vis=msfile, mode="tfcrop", spw='', action='apply', display='',
             timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
    flagmanager(msfile, mode='save', versionname='diskslfcal_remove_RFI-and-BURSTS')

    tdate = mstl.get_trange(msfile)[0].datetime.strftime('%Y%m%d')
    caltb = os.path.join(slfcaltbdir, tdate + '_1.pha')
    if os.path.exists(caltb):
        os.system('rm -rf {}*'.format(caltb))

    if pols == 'XXYY':
        caltbs_ = {'XX': [], 'YY': []}
        pols_ = ['XX', 'YY']
        msfileXY = {}
        for pol in pols_:
            msfileXY[pol] = '.'.join([msfile, pol])
            if os.path.exists(msfileXY[pol]):
                os.system('rm -rf {}'.format(msfileXY[pol]))
            mstl.splitX(vis=msfile, outputvis=msfileXY[pol], correlation=pol, datacolumn='data',
                        datacolumn2='MODEL_DATA')
        mstl.gaincalXY(vis=msfile, caltable=caltb, pols=pols, msfileXY=msfileXY, selectdata=True, uvrange="",
                       antenna="0~12&0~12", solint="inf",
                       combine="scan", refant="0", refantmode="strict", minsnr=1.0, gaintype="G", calmode="p",
                       append=False)
        for pol in pols_:
            caltb_ = '.'.join([caltb, pol])
            caltbs_[pol].append(caltb_)
    else:
        gaincal(vis=msfile, caltable=caltb, selectdata=True, uvrange="", antenna="0~12&0~12", solint="inf",
                combine="scan", refant="0", refantmode="strict", minsnr=1.0, gaintype="G", calmode="p", append=False)
    caltbs.append(caltb)

    caltb = os.path.join(slfcaltbdir, tdate + '_2.pha')
    if os.path.exists(caltb):
        os.system('rm -rf {}*'.format(caltb))
    # Second round of phase selfcal on the disk using solution interval "1min"
    if pols == 'XXYY':
        mstl.gaincalXY(vis=msfile, caltable=caltb, pols=pols, msfileXY=msfileXY, gaintableXY=caltbs_,
                       selectdata=True, uvrange="", antenna="0~12&0~12",
                       solint="10min",
                       combine="scan", interp="linear",
                       refant="0", refantmode="strict", minsnr=1.0, gaintype="G", calmode="p", append=False)
        for pol in pols_:
            caltb_ = '.'.join([caltb, pol])
            caltbs_[pol].append(caltb_)
    else:
        gaincal(vis=msfile, caltable=caltb, selectdata=True, uvrange="", antenna="0~12&0~12", solint="10min",
                combine="scan", gaintable=caltbs, interp="linear",
                refant="0", refantmode="strict", minsnr=1.0, gaintype="G", calmode="p", append=False)
    caltbs.append(caltb)

    caltb = os.path.join(slfcaltbdir, tdate + '_3.amp')
    if os.path.exists(caltb):
        os.system('rm -rf {}*'.format(caltb))
    # Final round of amplitude selfcal with 1-h solution interval (restrict to 16-24 UT)
    if pols == 'XXYY':
        mstl.gaincalXY(vis=msfile, caltable=caltb, pols=pols, msfileXY=msfileXY, gaintableXY=caltbs_,
                       selectdata=True, uvrange="", antenna="0~12&0~12",
                       timerange=trange, interp="linear",
                       solint="60min", combine="scan", refant="10", refantmode="flex", minsnr=1.0, gaintype="G",
                       calmode="a",
                       append=False)
        for pol in pols_:
            caltb_ = '.'.join([caltb, pol])
            caltbs_[pol].append(caltb_)
    else:
        gaincal(vis=msfile, caltable=caltb, selectdata=True, uvrange="", antenna="0~12&0~12",
                timerange=trange, gaintable=caltbs, interp="linear",
                solint="60min", combine="scan", refant="10", refantmode="flex", minsnr=1.0, gaintype="G", calmode="a",
                append=False)
    mstl.flagcaltboutliers(caltb, limit=[0.125, 8.0])
    # mstl.flagcaltboutliers(caltb, limit=[0.5, 2.0])
    caltbs.append(caltb)
    # Split out corrected data and model and do uvsub
    vis2 = 'slf3_' + msfile
    if os.path.exists(vis2):
        os.system('rm -rf {}'.format(vis2))
    if os.path.exists(vis2 + '.flagversions'):
        os.system('rm -rf {}'.format(vis2 + '.flagversions'))

    flagmanager(msfile, mode='restore', versionname='diskslfcal_init')
    clearcal(msfile)
    applycal(vis=msfile, selectdata=True, antenna="0~12", gaintable=caltbs, interp="linear", calwt=False,
             applymode="calonly")
    split(msfile, outputvis=vis2, datacolumn="corrected")
    for sp, dkcl in tqdm(enumerate(diskcl), desc='Inserting disk model', ascii=True):
        ft(vis=vis2, spw=str(sp), field='', model="", nterms=1,
           reffreq="", complist=str(dkcl), incremental=False, usescratch=True)
    #     # mstl.modeltransfer(msfile, spw='{}'.format(sp))
    uvsub(vis=vis2, reverse=False)

    # Final split to
    final = 'final_' + msfile
    if os.path.exists(final):
        os.system('rm -rf {}'.format(final))
    if os.path.exists(final + '.flagversions'):
        os.system('rm -rf {}'.format(final + '.flagversions'))
    split(vis2, outputvis=final, datacolumn='corrected')
    os.system('mv {} {}'.format(msfile + '.flagversions', final + '.flagversions'))
    # Remove the interim ms files
    if clearcache:
        if os.path.exists(msfile):
            os.system('rm -rf {}'.format(msfile))
        if os.path.exists(msfile + '.flagversions'):
            os.system('rm -rf {}'.format(msfile + '.flagversions'))
        if os.path.exists(vis2):
            os.system('rm -rf {}'.format(vis2))
        if os.path.exists(vis2 + '.flagversions'):
            os.system('rm -rf {}'.format(vis2 + '.flagversions'))

    # Return the name of the selfcaled ms
    return final, diskxmlfile


def fd_images(vis, cleanup=False, niter=None, spws=['0~1', '2~5', '6~10', '11~20', '21~30', '31~43'], imgoutdir='./',
              bright=None, stokes="XX"):
    ''' Create standard full-disk images in "images" subdirectory of the current directory.
        If cleanup is True, delete those images after completion, leaving only the fits images.
    '''
    # Check if "images" directory exists (if not, create it and mark it for later deletion)
    try:
        if os.stat('images'):
            rm_images = False  # Mark as not removeable
    except:
        os.mkdir('images')
        if cleanup:
            rm_images = True  # Mark as removeable
        else:
            rm_images = False  # Mark as not removeable

    trange = ant_trange(vis)
    tdate = trange.replace('/', '')[:8]
    if niter is None:
        niter = 5000
    if bright is None:
        bright = [True] * len(spws)
    imagefile = []
    fitsfile = []
    for s, sp in enumerate(spws):
        if bright[s]:
            spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
            imname = "images/briggs" + spwstr
            # tclean(vis=vis, selectdata=True, spw=sp, timerange=trange,
            #        antenna="0~12", datacolumn="corrected", imagename=imname, imsize=[1024], cell=['2.5arcsec'],
            #        stokes="XX", projection="SIN", specmode="mfs", interpolation="linear", deconvolver="multiscale",
            #        scales=[0, 5, 15, 30], nterms=2, smallscalebias=0.6, restoration=True, weighting="briggs", robust=0,
            #        niter=niter, gain=0.05, savemodel="none")
            os.system('rm -rf {}.*'.format(imname))

            tclean(vis=vis, selectdata=True, spw=sp, timerange=trange,
                   antenna="0~12", datacolumn="data", imagename=imname, imsize=[1024], cell=['2.5arcsec'],
                   stokes=stokes, projection="SIN", specmode="mfs", interpolation="linear", deconvolver="multiscale",
                   scales=[0, 5, 15, 30], nterms=2, smallscalebias=0.6, restoration=True, weighting="briggs",
                   robust=0.0,
                   niter=niter, gain=0.05, savemodel="none", usemask='auto-multithresh', pbmask=0.0,
                   sidelobethreshold=1.0, noisethreshold=2.5, lownoisethreshold=1.5, negativethreshold=5.0,
                   smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True,
                   minpercentchange=-1.0)
            outfits = os.path.join(imgoutdir, 'eovsa_' + tdate + '.spw' + spwstr + '.tb.fits')

            if os.path.exists(outfits):
                os.system('rm -rf {}'.format(outfits))

            imagefile.append(imname + '.image')
            fitsfile.append(outfits)
    hf.imreg(vis=vis, imagefile=imagefile, fitsfile=fitsfile, timerange=[trange] * len(fitsfile), toTb=True,
             usephacenter=False, overwrite=True)
    if rm_images:
        shutil.rmtree('images')  # Remove all images and the folder named images

    # To add disk model image to the images, I can try scipy.ndimage routines gaussian_filter() and zoom()
    return fitsfile


def feature_slfcal(vis, niter=200, spws=['0~1', '2~5', '6~10', '11~20', '21~30', '31~49'], slfcaltbdir='./',
                   bright=None, pols='XX'):
    ''' Uses images from disk-selfcaled data as model for further self-calibration of outer antennas.
        This is only a good idea if there are bright active regions that provide strong signal on the
        long baselines.
    '''
    trange = ant_trange(vis)

    if bright is None:
        bright = [True] * len(spws)
    # Insert model into ms and do "inf" gaincal, appending to table each subsequent time

    if os.path.exists('images_init'):
        os.system('rm -rf images_init')
    os.system('mv images images_init')
    clearcal(vis, addmodel=True)

    flagmanager(vis, mode='save', versionname='featureslfcal_init')
    ## automaticaly flag any high amplitudes from flares or RFI
    flagdata(vis=vis, mode="tfcrop", spw='', action='apply', display='',
             timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
    flagmanager(vis, mode='save', versionname='featureslfcal_remove_RFI-and-BURSTS')

    fd_images(vis, cleanup=False, niter=niter, spws=spws, bright=bright)  # Does shallow clean for selfcal purposes
    tdate = mstl.get_trange(vis)[0].datetime.strftime('%Y%m%d')
    caltb = os.path.join(slfcaltbdir, tdate + '_d1.pha')
    if os.path.exists(caltb):
        os.system('rm -rf {}*'.format(caltb))
    for s, sp in enumerate(spws):
        if bright[s]:
            spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
            imname = "images/briggs" + spwstr + '.model'
            if sp == '31~49':
                # The high-band image is only made to band 43, so adjust the name
                imname = 'images/briggs31-43.model'
            imcl = imname.replace('.model', '.cl')
            im2cl(imname, imcl)
            ## Note: ft does not work with complist if incremental is True.
            ## Likely, set incremental to False is ok if we make visibility model for spw by spw.
            ## Make model for specified spws will not affect other spws.
            ft(vis=vis, spw=sp, model="", complist=imcl, usescratch=True, incremental=False)

            ## Note: modeltransfer is commented because ft generates model for both XX and YY
            # if pols == 'XXYY':
            #     mstl.modeltransfer(vis, spw=sp)
    if pols == 'XXYY':
        mstl.gaincalXY(vis=vis, caltable=caltb, pols=pols, selectdata=True, timerange=trange, uvrange='>1.5Klambda',
                       combine="scan", antenna='0~12&0~12', refant='0', refantmode="strict", solint='inf', gaintype='G',
                       minsnr=1.0, calmode='p', append=False)
    else:
        gaincal(vis=vis, caltable=caltb, selectdata=True, timerange=trange, uvrange='>1.5Klambda',
                combine="scan", antenna='0~12&0~12', refant='0', refantmode="strict", solint='inf', gaintype='G',
                minsnr=1.0, calmode='p', append=False)
    # Apply the corrections to the data and split to a new ms
    applycal(vis=vis, selectdata=True, antenna="0~12", gaintable=caltb, interp="linear", calwt=False,
             applymode="calonly")
    vis1 = 'dslf1_' + vis
    if os.path.exists(vis1):
        os.system('rm -rf {}'.format(vis1))
    flagmanager(vis, mode='restore', versionname='featureslfcal_init')
    split(vis, outputvis=vis1, datacolumn="corrected")

    caltb = os.path.join(slfcaltbdir, tdate + '_d2.pha')
    if os.path.exists(caltb):
        os.system('rm -rf {}*'.format(caltb))
    # Move the existing images directory so that a new one will be created
    if os.path.exists('images_ftcal_rnd1'):
        os.system('rm -rf images_ftcal_rnd1')
    # shutil.move('images', 'old_images2')
    os.system('mv images images_ftcal_rnd1')
    # Make new model images for another round of selfcal
    flagmanager(vis1, mode='save', versionname='featureslfcal_init')
    ## automaticaly flag any high amplitudes from flares or RFI
    flagdata(vis=vis1, mode="tfcrop", spw='', action='apply', display='',
             timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
    flagmanager(vis1, mode='save', versionname='featureslfcal_remove_RFI-and-BURSTS')
    fd_images(vis1, cleanup=False, niter=niter, spws=spws, bright=bright)
    for s, sp in enumerate(spws):
        if bright[s]:
            spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
            imname = "images/briggs" + spwstr + '.model'
            if sp == '31~49':
                # The high-band image is only made to band 43, so adjust the name
                imname = 'images/briggs31-43.model'
            imcl = imname.replace('.model', '.cl')
            im2cl(imname, imcl)
            ft(vis=vis1, spw=sp, model="", complist=imcl, usescratch=True, incremental=False)
            ## Note: modeltransfer is commented because ft generates model for both XX and YY
            # if pols == 'XXYY':
            #     mstl.modeltransfer(vis1, spw=sp)
    if pols == 'XXYY':
        mstl.gaincalXY(vis=vis1, caltable=caltb, pols=pols, selectdata=True, timerange=trange, uvrange='>1.5Klambda',
                       combine="scan", antenna='0~12&0~12', refant='0', solint='10min', refantmode="strict",
                       gaintype='G', minsnr=1.0, calmode='p', append=False)
    else:
        gaincal(vis=vis1, caltable=caltb, selectdata=True, timerange=trange, uvrange='>1.5Klambda',
                combine="scan", antenna='0~12&0~12', refant='0', solint='10min', refantmode="strict", gaintype='G',
                minsnr=1.0, calmode='p', append=False)
    # Apply the corrections to the data and split to a new ms
    applycal(vis=vis1, selectdata=True, antenna="0~12", gaintable=caltb, interp="linear", calwt=False,
             applymode="calonly")
    vis2 = 'dslf2_' + vis
    if os.path.exists(vis2):
        os.system('rm -rf {}'.format(vis2))
    flagmanager(vis1, mode='restore', versionname='featureslfcal_init')
    mstl.splitX(vis1, outputvis=vis2, datacolumn="corrected", datacolumn2="model_data")
    if os.path.exists('images_ftcal_rnd2'):
        os.system('rm -rf images_ftcal_rnd2')
    os.system('mv images images_ftcal_rnd2')
    return vis2


def plt_eovsa_image(eofiles, figoutdir='./'):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from astropy.visualization.stretch import AsinhStretch, LinearStretch, SqrtStretch
    from astropy.visualization import ImageNormalize
    from suncasa.utils import plot_mapX as pmX
    from sunpy import map as smap
    import astropy.units as u
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.colorbar as colorbar

    # It is expected that nfiles will be either 4 (for older 34-band data) or 6 (for newer 52-band data)
    nfiles = len(eofiles)
    plt.ioff()
    fig = plt.figure(figsize=(5 * nfiles // 2, 9))

    axs = []
    cmap = 'gist_heat'
    for idx, eofile in enumerate(eofiles):
        ax = fig.add_subplot(2, nfiles // 2, idx + 1)
        axs.append(ax)
        # ax = axs[idx]
        eomap = smap.Map(eofile)
        tb_disk = eomap.meta['TBDISK']
        vmaxs = [70.0e4, 30e4, 18e4, 13e4, 8e4, 6e4]
        vmins = [-18.0e3, -8e3, -4.8e3, -3.4e3, -2.1e3, -1.6e3]
        # norm = colors.Normalize(vmin=vmins[idx], vmax=vmaxs[idx])
        stretch = AsinhStretch(a=0.15)
        norm = ImageNormalize(vmin=vmins[idx], vmax=vmaxs[idx], stretch=stretch)

        # norm = colors.Normalize(vmin=tb_disk * (-0.2), vmax=0.5*np.nanmax(eomap.data))
        eomap_ = pmX.Sunmap(eomap)
        eomap_.imshow(axes=ax, cmap=cmap, norm=norm)
        eomap_.draw_limb(axes=ax, lw=0.5, alpha=0.5)
        eomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.5)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2.0%', pad=0.08)
        cax.tick_params(direction='in')
        clb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        clb.set_label(r'T$_b$ [$\times$10$^3$K]')
        if idx != nfiles / 2:
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        ax.tick_params(direction="out")
        ax.text(0.02, 0.98,
                'EOVSA {:.1f} GHz   {}'.format(eomap.meta['CRVAL3'] / 1e9, eomap.date.strftime('%d-%b-%Y 20:00 UT')),
                transform=ax.transAxes, color='w', ha='left', va='top', fontsize=8, fontweight='bold')
        ax.text(0.02, 0.02, 'Max Tb {:.0f} K'.format(np.nanmax(eomap.data)),
                transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=8, fontweight='bold')
        ax.set_xlim(-1200, 1200)
        ax.set_ylim(-1200, 1200)
    fig.tight_layout()
    figname = os.path.join(figoutdir, 'eovsa_qlimg_{}.png'.format(eomap.date.strftime('%Y%m%d')))
    fig.savefig(figname, dpi=150)
    plt.close(fig)
    plt.ion()
    return figname


def pipeline_run(vis, outputvis='', workdir=None, slfcaltbdir=None, imgoutdir=None, figoutdir=None, clearcache=False,
                 pols='XX'):
    from astropy.io import fits

    spw2band = np.array([0, 1] + range(4, 52))
    defaultfreq = 1.1 + 0.325 * (spw2band + 0.5)
    # Use vis name to determine date, and hence number of bands
    if mstl.get_trange(vis)[0].mjd > 58536:
        # After 2019 Feb 22, the band numbers changed to 1-52, and spw from 0-49
        nbands = 52
        freq = defaultfreq
    else:
        # Before 2019 Feb 22, the band numbers were 1-34, and spw from 0-30
        nbands = 34
        freq = 1.419 + np.arange(nbands) / 2.

    slashdate = ant_trange(vis)[:10]

    spws = ['0~1', '2~5', '6~10', '11~20', '21~30', '31~43', '44~49']
    # bright_thresh = [6, 6, 5, 5, 5, 5, 5]
    bright_thresh = [6, 5, 4, 3, 2, 2, 2]
    active = False
    if nbands == 34:
        # These spectral window ranges correspond to the frequency ranges
        # of the last 4 band-ranges of the 52-band case.
        spws = ['1~3', '4~9', '10~16', '17~24', '25~30']
        bright_thresh = [6, 4.5, 3.5, 3, 2]

    if workdir is None:
        workdir = '/data1/workdir'
    os.chdir(workdir)
    if slfcaltbdir is None:
        slfcaltbdir = workdir + '/'
    if imgoutdir is None:
        imgoutdir = workdir + '/'
    if figoutdir is None:
        figoutdir = workdir + '/'
    if outputvis[-1] == '/':
        outputvis = outputvis[:-1]
    if vis[-1] == '/':
        vis = vis[:-1]

    if not os.path.exists(slfcaltbdir):
        os.makedirs(slfcaltbdir)
    # Verify that the vis is not in the current working directory
    if os.getcwd() == os.path.dirname(vis):
        print('Cannot copy vis file onto itself.')
        print('Please change to a different working directory')
        return None

    # Copy original ms to local directory
    if os.path.exists(os.path.basename(vis)):
        shutil.rmtree(os.path.basename(vis))
    print('Copy {} to working directory {}.'.format(vis, os.getcwd()))
    shutil.copytree(vis, os.path.basename(vis))
    vis = os.path.basename(vis)
    # Generate calibrated visibility by self calibrating on the solar disk
    ##ms_slfcaled, diskxmlfile = disk_slfcal(vis, slfcaltbdir=slfcaltbdir)
    flagmanager(vis, mode='save', versionname='pipeline_init')
    ## automaticaly flag any high amplitudes from flares or RFI
    flagdata(vis=vis, mode="tfcrop", spw='', action='apply', display='',
             timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
    flagmanager(vis, mode='save', versionname='pipeline_remove_RFI-and-bursts')

    # Make initial images from self-calibrated visibility file, and check T_b max
    if os.path.exists('images'):
        shutil.rmtree('images')
    outputfits = fd_images(vis, imgoutdir=imgoutdir, spws=spws)

    flagmanager(vis, mode='restore', versionname='pipeline_init')
    # outputfits = fd_images(vis, imgoutdir=imgoutdir, spws=spws, cleanup=True) change cleanup here?
    ####### outputfits is with model
    # Check if any of the images has a bright source (T_b > 300,000 K), and if so, remake images
    # with fewer components and execute feature_slfcal
    diskxmlfile = vis + '.SOLDISK.xml'
    dsize, fdens = calc_diskmodel(slashdate, nbands, freq, defaultfreq)
    insertdiskmodel(vis, dsize=dsize, fdens=fdens, xmlfile=diskxmlfile, writediskinfoonly=True)
    files = outputfits
    diskinfo = readdiskxml(diskxmlfile)
    bright = np.zeros((len(files)), dtype=np.bool)
    for idx, file in enumerate(files):
        if os.path.exists(file):
            tb_disk = image_adddisk(file, diskinfo, caltbonly=True)
            data = fits.getdata(file)
            data.shape = data.shape[-2:]  # gets rid of any leading axes of size 1
            if np.nanmax(data) > bright_thresh[idx] * tb_disk: bright[idx] = True

    if any(bright):
        print('spw {} have bright features on disk.'.format(';'.join(np.array(spws)[np.where(bright)[0]])))
        active = True
        # A bright source exists, so do feature self-calibration
        ms_slfcaled2 = feature_slfcal(vis, niter=200, slfcaltbdir=slfcaltbdir, spws=spws, bright=bright, pols=pols)
        vis = ms_slfcaled2
    else:
        if os.path.exists('images_init'):
            os.system('rm -rf images_init')
        os.system('mv images images_init')

    ms_slfcaled, diskxmlfile = disk_slfcal(vis, slfcaltbdir=slfcaltbdir, active=active, clearcache=clearcache,
                                           pols=pols)

    flagmanager(ms_slfcaled, mode='save', versionname='pipeline_final_init')
    ## automaticaly flag any high amplitudes from flares or RFI
    flagdata(vis=ms_slfcaled, mode="tfcrop", spw='', action='apply', display='',
             timecutoff=3.0, freqcutoff=2.0, maxnpieces=2, flagbackup=False)
    flagmanager(ms_slfcaled, mode='save', versionname='pipeline_final_remove_RFI-and-BURSTS')

    outputfits = fd_images(ms_slfcaled, imgoutdir=imgoutdir, spws=spws)

    flagmanager(ms_slfcaled, mode='restore', versionname='pipeline_final_init')

    if outputvis:
        if os.path.exists(outputvis):
            os.system('rm -rf {}'.format(outputvis))
        os.system('mv {} {}'.format(ms_slfcaled, outputvis))
        os.system('mv {}.flagversions {}.flagversions'.format(ms_slfcaled, outputvis))
        ms_slfcaled = outputvis

        newdiskxmlfile = '{}.SOLDISK.xml'.format(outputvis)
        if os.path.exists(newdiskxmlfile):
            os.system('rm -rf {}'.format(newdiskxmlfile))
        os.system('mv {} {}'.format(diskxmlfile, newdiskxmlfile))
        diskxmlfile = newdiskxmlfile

    eofiles = []
    datestr = mstl.get_trange(ms_slfcaled)[0].datetime.strftime('%Y%m%d')
    for s, sp in enumerate(spws):
        spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
        eofiles.append(imgoutdir + '/eovsa_{}.spw{}.tb.fits'.format(datestr, spwstr))
    eofiles = sorted(eofiles)
    eofiles_new = []
    diskinfo = readdiskxml(diskxmlfile)
    for idx, eofile in enumerate(eofiles):
        eomap_disk, tb_disk, eofile_new = image_adddisk(eofile, diskinfo)
        eofiles_new.append(eofile_new)

    ### obsolete module --- replaced by eovsa_pltQlookImage.py
    # plt_eovsa_image(eofiles_new[:-1], figoutdir)  # skip plotting the image of the highest bands
    return ms_slfcaled
