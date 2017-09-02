import os
import numpy as np
from datetime import datetime
from math import *
# import jdutil
import bisect
import pdb
from taskinit import *

try:
    from astropy.io import fits as pyfits
except:
    try:
        import pyfits
    except ImportError:
        raise ImportError('Neither astropy nor pyfits exists in this CASA installation')


# from astropy.constants import R_sun, au

def read_horizons(ephemfile=None):
    # inputs:
    #   ephemfile:
    #       OBSERVER output from JPL Horizons for topocentric coordinates with for example
    #       target=Sun, observer=VLA=-5
    #       extra precision, quantities 1,17,20, REFRACTION
    #       routine goes through file to find $$SOE which is start of ephemeris and ends with $$EOE
    # outputs: a Python dictionary containing the following:
    #   timestr: date and time as a string
    #   time: modified Julian date
    #   ra: right ascention, in rad
    #   dec: declination, in rad
    #   rastr: ra in string
    #   decstr: dec in string
    #   p0: solar p angle, CCW with respect to the celestial north pole
    #   delta: distance from the disk center to the observer, in AU
    #   delta_dot: time derivative of delta, in the light of sight direction. Negative means it is moving toward the observer
    #
    if not ephemfile or ephemfile.isspace():
        raise ValueError, 'Please specify input ephem file!'
    if not os.path.isfile(ephemfile):
        raise ValueError, 'The specified input ephem file does not exist!'
    # initialize the return dictionary
    ephem0 = dict.fromkeys(['time', 'timestr', 'ra', 'dec', 'rastr', 'decstr', 'delta', 'delta_dot', 'p0'])
    f = open(ephemfile, 'rU')  # force the universal newline mode
    lines = f.readlines()
    f.close()
    nline = len(lines)
    for i in range(nline):
        line = lines[i]
        if line[0:5] == '$$SOE':  # start recording
            istart = i + 1
        if line[0:5] == '$$EOE':  # end recording
            iend = i
    newlines = lines[istart:iend]
    nrec = len(newlines)
    ephem_ = []
    for j in range(nrec):
        line = newlines[j]
        dt = datetime.strptime(line[1:18], '%Y-%b-%d %H:%M')
        newdtstr = dt.strftime('%Y-%m-%dT%H:%M')
        ephem0['timestr'] = newdtstr
        mjd = qa.quantity(newdtstr)
        # jd=jdutil.datetime_to_jd(dt)
        # mjd=jdutil.jd_to_mjd(jd)
        ephem0['time'] = mjd
        ephem0['rastr'] = line[23:36]
        ephem0['decstr'] = line[37:50]
        ephem0['ra'] = {'unit': 'rad', 'value': radians(
            (long(line[23:25]) + long(line[26:28]) / 60. + float(line[29:36]) / 3600.) * 15.)}  # in rad
        ephem0['dec'] = {'unit': 'rad', 'value': radians(
            long(line[38:40]) + long(line[41:43]) / 60. + float(line[44:50]) / 3600.)}  # in rad
        ephem0['p0'] = {'unit': 'deg', 'value': float(line[51:59])}
        ephem0['delta'] = {'unit': 'au', 'value': float(line[70:86])}
        ephem0['delta_dot'] = {'unit': 'km/s', 'value': float(line[88:98])}
        if line[37] == '-':
            ephem0['dec']['value'] = -ephem0['dec']['value']
        ephem_.append(ephem0.copy())
    # convert list of dictionary to a dictionary of arrays
    times = [ep['time']['value'] for ep in ephem_]
    ras = [ep['ra']['value'] for ep in ephem_]
    decs = [ep['dec']['value'] for ep in ephem_]
    p0s = [ep['p0']['value'] for ep in ephem_]
    deltas = [ep['delta']['value'] for ep in ephem_]
    ephem = {'times': times, 'ras': ras, 'decs': decs, 'p0s': p0s, 'deltas': deltas}
    return ephem


def read_msinfo(msfile=None, msinfofile=None):
    # read MS information #
    msinfo = dict.fromkeys(['msfile', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs'])
    ms.open(msfile)
    scans = ms.getscansummary()
    scanids = sorted(scans.keys(), key=lambda x: int(x))
    nscanid = len(scanids)
    btimes = []
    btimestr = []
    etimes = []
    fieldids = []
    inttimes = []
    dirs = []
    ras = []
    decs = []
    for i in range(nscanid):
        btimes.append(scans[scanids[i]]['0']['BeginTime'])
        etimes.append(scans[scanids[i]]['0']['EndTime'])
        fieldid = scans[scanids[i]]['0']['FieldId']
        fieldids.append(fieldid)
        dir = ms.getfielddirmeas('PHASE_DIR', fieldid)
        dirs.append(dir)
        ras.append(dir['m0'])
        decs.append(dir['m1'])
        inttimes.append(scans[scanids[i]]['0']['IntegrationTime'])
    ms.close()
    btimestr = [qa.time(qa.quantity(btimes[i], 'd'), form='fits', prec=10)[0] for i in range(nscanid)]
    msinfo['msfile'] = msfile
    msinfo['scans'] = scans
    msinfo['fieldids'] = fieldids
    msinfo['btimes'] = btimes
    msinfo['btimestr'] = btimestr
    msinfo['inttimes'] = inttimes
    msinfo['ras'] = ras
    msinfo['decs'] = decs
    if msinfofile:
        np.savez(msinfofile, msfile=msfile, scans=scans, fieldids=fieldids, btimes=btimes, btimestr=btimestr,
                 inttimes=inttimes, ras=ras, decs=decs)
    return msinfo


def ephem_to_helio(msinfo=None, ephem=None, reftime=None, polyfit=None):
    ## 1. Take a solar ms database, read the scan and field information, find out the pointings (in RA and DEC)
    ## 2. Compare with the ephemeris of the solar disk center (in RA and DEC)
    ## 3. Generate VLA pointings in heliocentric coordinates
    ## inputs:
    ##      msinfo: CASA MS information, output from read_msinfo
    ##      ephem: solar ephem, output from read_horizons
    ##      reftime: list of reference times (e.g., used for imaging)
    #                CASA standard time format, either a single time (e.g., '2012/03/03/12:00:00'
    ##               or a time range (e.g., '2012/03/03/12:00:00~2012/03/03/13:00:00'. If the latter,
    ##               take the midpoint of the timerange for reference. If no date specified, take
    ##               the date of the first scan
    #       polyfit: ONLY works for MS database with only one source with continously tracking; not recommanded unless scan length is too long and want to have very high accuracy
    ## return value:
    ##      helio: a list of VLA pointing information
    ##              reftimestr: reference time, in FITS format string
    ##              reftime: reference time, in mjd format
    ##              ra: actual RA of VLA pointing at the reference time (interpolated)
    ##              dec: actual DEC of VLA pointing at the reference time (interpolated)
    ##              # CASA knows only RA and DEC of the closest field (e.g. in clean) #
    ##              ra_fld: right ascention of the CASA reference pointing direction
    ##              dec_fld: declination of the CASA reference pointing direction
    ##              raoff: RA offset of the actual VLA pointing to solar center
    ##              decoff: DEC offset of the actual VLA pointing to solar center
    ##              refx: heliocentric X offset of the actual VLA pointing to solar center
    ##              refy: heliocentric Y offset of the actual VLA pointing to solar center
    ######## Example #########
    #    msfile='sun_C_20140910T221952-222952.10s.cal.ms'
    #    ephemfile='horizons_sun_20140910.radecp'
    #    ephem=vla_prep.read_horizons(ephemfile=ephemfile)
    #    msinfo=vla_prep.read_msinfo(msfile=msfile)
    #    polyfit=0
    #    reftime = '22:25:20~22:25:40'
    if not ephem:
        raise ValueError, 'Please provide information of the MS database!'
    if not msinfo:
        raise ValueError, 'Please provide information of the MS database!'
    if isinstance(msinfo, str):
        try:
            msinfo0 = np.load(msinfo)
        except:
            raise ValueError, 'The specified input msinfo file does not exist!'
    elif isinstance(msinfo, dict):
        msinfo0 = msinfo
    else:
        raise ValueError, 'msinfo should be either a numpy npz or a dictionary'
    print 'msinfo is derived from: ', msinfo0['msfile']
    scans = msinfo0['scans']
    fieldids = msinfo0['fieldids']
    btimes = msinfo0['btimes']
    inttimes = msinfo0['inttimes']
    ras = msinfo0['ras']
    decs = msinfo0['decs']
    ra_rads = [ra['value'] for ra in ras]
    dec_rads = [dec['value'] for dec in decs]
    # fit 2nd order polynomial fits to the RAs and DECs #
    if polyfit:
        cra = np.polyfit(btimes, ra_rads, 2)
        cdec = np.polyfit(btimes, dec_rads, 2)

    # find out pointing direction according to the input time or timerange #
    if not reftime:
        raise ValueError, 'Please specify a reference time for pointing/imaging!'
    if isinstance(reftime, str):
        reftime = [reftime]
    if (not isinstance(reftime, list)):
        print 'input "reftime" is not a valid list. Abort...'

    nreftime = len(reftime)
    helio = []
    for reftime0 in reftime:
        helio0 = dict.fromkeys(
            ['reftimestr', 'reftime', 'ra', 'dec', 'ra_fld', 'dec_fld', 'raoff', 'decoff', 'refx', 'refy', 'p0'])
        helio0['reftimestr'] = reftime0
        if '~' in reftime0:
            # if reftime0 is specified as a timerange
            [tbg0, tend0] = reftime0.split('~')
            tbg_d = qa.getvalue(qa.convert(qa.totime(tbg0), 'd'))[0]
            tend_d = qa.getvalue(qa.convert(qa.totime(tend0), 'd'))[0]
            tdur_s = (tend_d - tbg_d) * 3600. * 24.
            # if no date is specified, add up the date of the first scan
            if tend_d < 1.:
                if tend_d >= tbg_d:
                    tend_d += int(btimes[0])
                else:
                    tend_d += int(btimes[0]) + 1
            if tbg_d < 1.:
                tbg_d += int(btimes[0])
            tref_d = (tbg_d + tend_d) / 2.
        else:
            # if reftime0 is specified as a single value
            tref_d = qa.getvalue(qa.convert(qa.totime(reftime0), 'd'))
            # if no date is specified, add up the date of the first scan
            if tref_d < 1.:
                tref_d += int(btimes[0])
            tbg_d = tref_d
            # use the intergration time
            ind = bisect.bisect_left(btimes, tref_d)
            tdur_s = inttims[ind - 1]
        helio0['reftime'] = tref_d
        helio0['date-obs'] = qa.time(qa.quantity(tbg_d, 'd'), form='fits', prec=10)[0]
        helio0['exptime'] = tdur_s

        # find out RA and DEC coords according to the reference time
        # if polyfit, then use the 2nd order polynomial coeffs
        ind = bisect.bisect_left(btimes, tref_d)
        if ind > 1:
            dt = tref_d - btimes[ind - 1]
            if ind < len(btimes):
                scanlen = btimes[ind] - btimes[ind - 1]
                (ra_b, ra_e) = (ras[ind - 1]['value'], ras[ind]['value'])
                (dec_b, dec_e) = (decs[ind - 1]['value'], decs[ind]['value'])
            if ind >= len(btimes):
                scanlen = btimes[ind - 1] - btimes[ind - 2]
                (ra_b, ra_e) = (ras[ind - 2]['value'], ras[ind - 1]['value'])
                (dec_b, dec_e) = (decs[ind - 2]['value'], decs[ind - 1]['value'])
        if ind == 1:  # only one scan exists (e.g., imported from AIPS)
            ra_b = ras[ind - 1]['value']
            ra_e = ra_b
            dec_b = decs[ind - 1]['value']
            dec_e = dec_b
            scanlen = 10.  # radom value
            dt = 0.
        if ind < 1:
            raise ValueError, 'Reference time does not fall into the scan list!'
        if polyfit:
            ra = cra[0] * tref_d ** 2. + cra[1] * tref_d + cra[2]
            dec = cdec[0] * tref_d ** 2. + cdec[1] * tref_d + cdec[2]
        # if not, use linearly interpolated RA and DEC at the beginning of this scan and next scan
        else:
            ra = ra_b + (ra_e - ra_b) / scanlen * dt
            dec = dec_b + (dec_e - dec_b) / scanlen * dt
        if ra < 0:
            ra += 2. * np.pi
        if ra_b < 0:
            ra_b += 2. * np.pi

        # compare with ephemeris from JPL Horizons
        time0s = ephem['times']
        ra0s = ephem['ras']
        dec0s = ephem['decs']
        p0s = ephem['p0s']
        delta0s = ephem['deltas']
        ind = bisect.bisect_left(time0s, tref_d)
        dt0 = time0s[ind] - time0s[ind - 1]
        dt_ref = tref_d - time0s[ind - 1]
        dra0 = ra0s[ind] - ra0s[ind - 1]
        ddec0 = dec0s[ind] - dec0s[ind - 1]
        dp0 = p0s[ind] - p0s[ind - 1]
        ddelta0 = delta0s[ind] - delta0s[ind - 1]
        ra0 = ra0s[ind - 1] + dra0 / dt0 * dt_ref
        dec0 = dec0s[ind - 1] + ddec0 / dt0 * dt_ref
        p0 = p0s[ind - 1] + dp0 / dt0 * dt_ref
        delta0 = delta0s[ind - 1] + ddelta0 / dt0 * dt_ref
        if ra0 < 0:
            ra0 += 2. * np.pi

        # RA and DEC offset in arcseconds
        decoff = degrees((dec - dec0)) * 3600.
        raoff = degrees((ra - ra0) * cos(dec)) * 3600.
        # Convert into heliocentric offsets
        prad = -radians(p0)
        refx = (-raoff) * cos(prad) - decoff * sin(prad)
        refy = (-raoff) * sin(prad) + decoff * cos(prad)
        helio0['ra'] = ra  # ra of the actual pointing
        helio0['dec'] = dec  # dec of the actual pointing
        helio0['ra_fld'] = ra_b  # ra of the field, used as the reference in e.g., clean
        helio0['dec_fld'] = dec_b  # dec of the field, used as the refenrence in e.g., clean
        helio0['raoff'] = raoff
        helio0['decoff'] = decoff
        helio0['refx'] = refx
        helio0['refy'] = refy
        helio0['p0'] = p0
        # helio['r_sun']=np.degrees(R_sun.value/(au.value*delta0))*3600. #in arcsecs
        helio.append(helio0)
    return helio


def getbeam(imagefile=None, beamfile=None):
    if not imagefile:
        raise ValueError, 'Please specify input images'
    bmaj = []
    bmin = []
    bpa = []
    beamunit = []
    bpaunit = []
    chans = []
    nimg = len(imagefile)
    for n in range(nimg):
        img = imagefile[n]
        if not os.path.exists(img):
            raise ValueError, 'The input image does not exist!'
        ia.open(img)
        sum = ia.summary()
        bmaj_ = []
        bmin_ = []
        bpa_ = []
        if sum.has_key('perplanebeams'):  # beam vary with frequency
            nbeams = sum['perplanebeams']['nChannels']
            beams = sum['perplanebeams']['beams']
            chans_ = [key[1:] for key in beams.keys()]
            chans_.sort(key=float)
            for chan in chans_:
                bmaj0 = beams['*' + chan]['*0']['major']['value']
                bmaj_.append(bmaj0)
                bmin0 = beams['*' + chan]['*0']['minor']['value']
                bmin_.append(bmin0)
                bpa0 = beams['*' + chan]['*0']['positionangle']['value']
                bpa_.append(bpa0)
            beamunit_ = beams['*' + chans_[0]]['*0']['major']['unit']
            bpaunit_ = beams['*' + chans_[0]]['*0']['positionangle']['unit']
        if sum.has_key('restoringbeam'):  # only one beam
            bmaj_.append(sum['restoringbeam']['major']['value'])
            bmin_.append(sum['restoringbeam']['minor']['value'])
            bpa_.append(sum['restoringbeam']['positionangle']['value'])
            beamunit_ = sum['restoringbeam']['major']['unit']
            bpaunit_ = sum['restoringbeam']['positionangle']['unit']
            nbeams = 1
            chans_ = [0]

        bmaj.append(bmaj_)
        bmin.append(bmin_)
        bpa.append(bpa_)
        beamunit.append(beamunit_)
        bpaunit.append(bpaunit_)
        chans.append(chans_)
    if beamfile:  # write beams to ascii file
        print 'Writing beam info to ascii file...'
        f = open(beamfile, 'w')
        f.write('CHANNEL No., BMAJ (' + beamunit[0] + '), BMIN (' + beamunit[0] + '), BPA (' + bpaunit[0] + ')')
        f.write("\n")
        for n in range(nimg):
            f.write('----For image: ' + imagefile[n] + '----')
            f.write('\n')
            chans_ = chans[n]
            for i in range(len(chans_)):
                f.write(str(chans_[i]) + ', ' + str(bmaj[n][i]) + ', ' + str(bmin[n][i]) + ', ' + str(bpa[n][i]))
                f.write("\n")
        f.close()
    return bmaj, bmin, bpa, beamunit, bpaunit


def imreg(imagefile=None, fitsfile=None, beamfile=None, helio=None, offsetfile=None, toTb=None, scl100=None,
          verbose=False):
    if not imagefile:
        raise ValueError, 'Please specify input image'
    if not helio:
        raise ValueError, 'Please specify input coordinate info for image registration. Use ephem_to_helio to derive that info'
    if not fitsfile:
        fitsfile = [img + '.fits' for img in imagefile]
    # if len(imagefile) != len(helio):
    #    raise ValueError, 'Number of input images does not equal to number of helio coord headers!'
    if len(imagefile) != len(fitsfile):
        raise ValueError, 'Number of input images does not equal to number of output fits files!'
    # get restoring beam info
    (bmajs, bmins, bpas, beamunits, bpaunits) = getbeam(imagefile=imagefile, beamfile=beamfile)
    nimg = len(imagefile)
    if verbose:
        print str(nimg) + ' images to process...'
    for n in range(nimg):
        if verbose:
            print 'processing image #' + str(n)
        img = imagefile[n]
        fitsf = fitsfile[n]
        hel = helio[n]
        bmaj = bmajs[n]
        bmin = bmins[n]
        beamunit = beamunits[n]
        if not os.path.exists(img):
            raise ValueError, 'Please specify input image'
        if os.path.exists(fitsf):
            raise ValueError, 'Specified fits file already exists!'
        else:
            p0 = hel['p0']
            ia.open(img)
            imr = ia.rotate(pa=str(-p0) + 'deg')
            imr.tofits(fitsf, history=False)
            imr.close()
            sum = ia.summary()
            ia.close()
        # construct the standard fits header
        # RA and DEC of the reference pixel crpix1 and crpix2
        (imra, imdec) = (sum['refval'][0], sum['refval'][1])
        # find out the difference of the image center to the CASA reference center
        # RA and DEC difference in arcseconds
        ddec = degrees((imdec - hel['dec_fld'])) * 3600.
        dra = degrees((imra - hel['ra_fld']) * cos(hel['dec_fld'])) * 3600.
        # Convert into image heliocentric offsets
        prad = -radians(hel['p0'])
        dx = (-dra) * cos(prad) - ddec * sin(prad)
        dy = (-dra) * sin(prad) + ddec * cos(prad)
        if offsetfile:
            try:
                offset = np.load(offsetfile)
            except:
                raise ValueError, 'The specified offsetfile does not exist!'
            reftimes_d = offset['reftimes_d']
            xoffs = offset['xoffs']
            yoffs = offset['yoffs']
            timg_d = hel['reftime']
            ind = bisect.bisect_left(reftimes_d, timg_d)
            xoff = xoffs[ind - 1]
            yoff = yoffs[ind - 1]
        else:
            xoff = hel['refx']
            yoff = hel['refy']
        if verbose:
            print 'offset of image phase center to visibility phase center (arcsec): ', dx, dy
            print 'offset of visibility phase center to solar disk center (arcsec): ', xoff, yoff
        (crval1, crval2) = (xoff + dx, yoff + dy)
        # update the fits header to heliocentric coordinates

        hdu = pyfits.open(fitsf, mode='update')
        header = hdu[0].header
        (cdelt1, cdelt2) = (
            -header['cdelt1'] * 3600., header['cdelt2'] * 3600.)  # Original CDELT1, 2 are for RA and DEC in degrees
        header['cdelt1'] = cdelt1
        header['cdelt2'] = cdelt2
        header['cunit1'] = 'arcsec'
        header['cunit2'] = 'arcsec'
        header['crval1'] = crval1
        header['crval2'] = crval2
        header['ctype1'] = 'HPLN-TAN'
        header['ctype2'] = 'HPLT-TAN'
        header['date-obs'] = hel['date-obs']
        try:
            # this works for pyfits version of CASA 4.7.0 but not CASA 4.6.0
            header.update('exptime', hel['exptime'])
            header.update('p_angle', hel['p0'])
        except:
            # this works for astropy.io.fits
            header.append(('exptime', hel['exptime']))
            header.append(('p_angle', hel['p0']))

        # header.update('comment', 'Fits header updated to heliocentric coordinates by Bin Chen')

        # update intensity units, i.e. to brightness temperature?
        if toTb:
            data = hdu[0].data  # remember the data order is reversed due to the FITS convension
            dim = data.ndim
            sz = data.shape
            keys = header.keys()
            values = header.values()
            # which axis is frequency?
            faxis = keys[values.index('FREQ')][-1]
            faxis_ind = dim - int(faxis)
            if header['BUNIT'].lower() == 'jy/beam':
                header['BUNIT'] = 'K'
                for i in range(sz[faxis_ind]):
                    nu = header['CRVAL' + faxis] + header['CDELT' + faxis] * (i + 1 - header['CRPIX' + faxis])
                    if header['CUNIT' + faxis] == 'KHz':
                        nu *= 1e3
                    if header['CUNIT' + faxis] == 'MHz':
                        nu *= 1e6
                    if header['CUNIT' + faxis] == 'GHz':
                        nu *= 1e9
                    if len(bmaj) > 1:  # multiple (per-plane) beams
                        bmajtmp = bmaj[i]
                        bmintmp = bmin[i]
                    else:  # one single beam
                        bmajtmp = bmaj[0]
                        bmintmp = bmin[0]
                    if beamunit == 'arcsec':
                        bmaj0 = np.radians(bmajtmp / 3600.)
                        bmin0 = np.radians(bmajtmp / 3600.)
                    if beamunit == 'arcmin':
                        bmaj0 = np.radians(bmajtmp / 60.)
                        bmin0 = np.radians(bmintmp / 60.)
                    if beamunit == 'deg':
                        bmaj0 = np.radians(bmajtmp)
                        bmin0 = np.radians(bmintmp)
                    if beamunit == 'rad':
                        bmaj0 = bmajtmp
                        bmin0 = bmintmp
                    beam_area = bmaj0 * bmin0 * np.pi / (4. * log(2.))
                    k_b = qa.constants('k')['value']
                    c_l = qa.constants('c')['value']
                    factor = 2. * k_b * nu ** 2 / c_l ** 2  # SI unit
                    jy_to_si = 1e-26
                    # print nu/1e9, beam_area, factor
                    factor2 = 1.
                    if scl100:
                        factor2 = 100.
                    if faxis == '3':
                        data[:, i, :, :] *= jy_to_si / beam_area / factor * factor2
                    if faxis == '4':
                        data[i, :, :, :] *= jy_to_si / beam_area / factor * factor2

        hdu.flush()
        hdu.close()
