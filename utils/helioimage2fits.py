import os
import numpy as np
from datetime import datetime
from math import *
# import jdutil
import bisect
import pdb
from taskinit import ms, tb, qa, iatool
from astropy.time import Time
from sunpy import sun
import astropy.units as u

try:
    from astropy.io import fits as pyfits
except:
    try:
        import pyfits
    except ImportError:
        raise ImportError('Neither astropy nor pyfits exists in this CASA installation')


# from astropy.constants import R_sun, au

def msclearhistory(msfile):
    from taskinit import tb
    tb.open(msfile + '/HISTORY', nomodify=False)
    nrows = tb.nrows()
    if nrows > 0:
        tb.removerows(range(nrows))
    tb.close()


def read_horizons(t0=None, dur=None, vis=None, observatory=None, verbose=False):
    import urllib2
    import ssl
    if not t0 and not vis:
        t0 = Time.now()
    if not dur:
        dur = 1. / 60. / 24.  # default to 2 minutes
    if t0:
        try:
            btime = Time(t0)
        except:
            print('input time ' + str(t0) + ' not recognized')
            return -1
    if vis:
        if not os.path.exists(vis):
            print 'Input ms data ' + vis + ' does not exist! '
            return -1
        try:
            # ms.open(vis)
            # summary = ms.summary()
            # ms.close()
            # btime = Time(summary['BeginTime'], format='mjd')
            # etime = Time(summary['EndTime'], format='mjd')
            ## alternative way to avoid conflicts with importeovsa, if needed -- more time consuming
            if observatory=='geocentric':
                observatory='500'
            else:
                ms.open(vis)
                metadata = ms.metadata()
                if metadata.observatorynames()[0] == 'EVLA':
                    observatory = '-5'
                elif metadata.observatorynames()[0] == 'EOVSA':
                    observatory = '-81'
                elif metadata.observatorynames()[0] == 'ALMA':
                    observatory = '-7'
                ms.close()
            tb.open(vis)
            btime_vis = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
            etime_vis = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
            tb.close()
            if verbose:
                print "Beginning time of this scan " + btime_vis.iso
                print "End time of this scan " + etime_vis.iso

            # extend the start and end time for jpl horizons by 0.5 hr on each end
            btime = Time(btime_vis.mjd - 0.5 / 24., format='mjd')
            dur = etime_vis.mjd - btime_vis.mjd + 1.0 / 24.
        except:
            print 'error in reading ms file: ' + vis + ' to obtain the ephemeris!'
            return -1

    # default the observatory to VLA, if none provided
    if not observatory:
        observatory = '-5'

    etime = Time(btime.mjd + dur, format='mjd')

    try:
        cmdstr = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE='OBSERVER'&QUANTITIES='1,17,20'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'&CAL_FORMAT='BOTH'&SOLAR_ELONG='0,180'&CENTER='{}@399'&COMMAND='10'&START_TIME='".format(
            observatory) + btime.iso.replace(' ', ',') + "'&STOP_TIME='" + etime.iso[:-4].replace(' ',
                                                                                                  ',') + "'&STEP_SIZE='1m'&SKIP_DAYLT='NO'&EXTRA_PREC='YES'&APPARENT='REFRACTED'"
        cmdstr = cmdstr.replace("'", "%27")
        try:
            context = ssl._create_unverified_context()
            f = urllib2.urlopen(cmdstr, context=context)
        except:
            f = urllib2.urlopen(cmdstr)
        lines = f.readlines()
        f.close()
    except:
        #todo use geocentric coordinate for the new VLA data
        import requests, collections
        params = collections.OrderedDict()
        params['batch'] = '1'
        params['TABLE_TYPE'] = "'OBSERVER'"
        params['QUANTITIES'] = "'1,17,20'"
        params['CSV_FORMAT'] = "'YES'"
        params['ANG_FORMAT'] = "'DEG'"
        params['CAL_FORMAT'] = "'BOTH'"
        params['SOLAR_ELONG'] = "'0,180'"
        if observatory=='500':
            params['CENTER'] = "'500'"
        else:
            params['CENTER'] = "'{}@399'".format(observatory)
        params['COMMAND'] = "'10'"
        params['START_TIME'] = "'{}'".format(btime.iso[:-4].replace(' ', ','))
        params['STOP_TIME'] = "'{}'".format(etime.iso[:-4].replace(' ', ','))
        params['STEP_SIZE'] = "'1m'"
        params['SKIP_DAYLT'] = "'NO'"
        params['EXTRA_PREC'] = "'YES'"
        params['APPAENT'] = "'REFRACTED'"
        results = requests.get("https://ssd.jpl.nasa.gov/horizons_batch.cgi", params=params)
        lines = [ll for ll in results.iter_lines()]

    nline = len(lines)
    istart = 0
    for i in range(nline):
        line = lines[i]
        if line[0:5] == '$$SOE':  # start recording
            istart = i + 1
        if line[0:5] == '$$EOE':  # end recording
            iend = i
    newlines = lines[istart:iend]
    nrec = len(newlines)
    ephem_ = []
    t = []
    ra = []
    dec = []
    p0 = []
    delta = []
    for line in newlines:
        items = line.split(',')
        t.append(Time(float(items[1]), format='jd').mjd)
        ra.append(np.radians(float(items[4])))
        dec.append(np.radians(float(items[5])))
        p0.append(float(items[6]))
        delta.append(float(items[8]))
    # convert list of dictionary to a dictionary of arrays
    ephem = {'time': t, 'ra': ra, 'dec': dec, 'p0': p0, 'delta': delta}
    return ephem


def read_msinfo(vis=None, msinfofile=None, use_scan_time=True):
    import glob
    # read MS information #
    msinfo = dict.fromkeys(['vis', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs', 'observatory'])
    ms.open(vis)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
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
    ephem_file = glob.glob(vis + '/FIELD/EPHEM*SUN.tab')
    if ephem_file:
        print('Loading ephemeris info from {}'.format(ephem_file[0]))
        tb.open(ephem_file[0])
        col_ra = tb.getcol('RA')
        col_dec = tb.getcol('DEC')
        col_mjd = tb.getcol('MJD')
        if use_scan_time:
            from scipy.interpolate import interp1d
            f_ra = interp1d(col_mjd, col_ra)
            f_dec = interp1d(col_mjd, col_dec)
            for idx, scanid in enumerate(scanids):
                btimes.append(scans[scanid]['0']['BeginTime'])
                etimes.append(scans[scanid]['0']['EndTime'])
                fieldid = scans[scanid]['0']['FieldId']
                fieldids.append(fieldid)
                inttimes.append(scans[scanid]['0']['IntegrationTime'])
            ras = f_ra(np.array(btimes))
            decs = f_dec(np.array(btimes))
            ras = qa.convert(qa.quantity(ras, 'deg'), 'rad')
            decs = qa.convert(qa.quantity(decs, 'deg'), 'rad')
        else:
            ras = qa.convert(qa.quantity(col_ra, 'deg'), 'rad')
            decs = qa.convert(qa.quantity(col_dec, 'deg'), 'rad')
    else:
        for idx, scanid in enumerate(scanids):
            btimes.append(scans[scanid]['0']['BeginTime'])
            etimes.append(scans[scanid]['0']['EndTime'])
            fieldid = scans[scanid]['0']['FieldId']
            fieldids.append(fieldid)
            inttimes.append(scans[scanid]['0']['IntegrationTime'])
            dir = ms.getfielddirmeas('PHASE_DIR', fieldid)
            dirs.append(dir)
            ras.append(dir['m0'])
            decs.append(dir['m1'])
    ms.close()
    btimestr = [qa.time(qa.quantity(btimes[idx], 'd'), form='fits', prec=10)[0] for idx in range(nscanid)]
    msinfo['vis'] = vis
    msinfo['scans'] = scans
    msinfo['fieldids'] = fieldids
    msinfo['btimes'] = btimes
    msinfo['btimestr'] = btimestr
    msinfo['inttimes'] = inttimes
    msinfo['ras'] = ras
    msinfo['decs'] = decs
    msinfo['observatory'] = observatory
    if msinfofile:
        np.savez(msinfofile, vis=vis, scans=scans, fieldids=fieldids, btimes=btimes, btimestr=btimestr,
                 inttimes=inttimes, ras=ras, decs=decs,
                 observatory=observatory)
    return msinfo


def ephem_to_helio(vis=None, ephem=None, msinfo=None, reftime=None, polyfit=None, usephacenter=False):
    '''1. Take a solar ms database, read the scan and field information, find out the pointings (in RA and DEC)
       2. Compare with the ephemeris of the solar disk center (in RA and DEC)
       3. Generate VLA pointings in heliocentric coordinates
       inputs:
           msinfo: CASA MS information, output from read_msinfo
           ephem: solar ephem, output from read_horizons
           reftime: list of reference times (e.g., used for imaging)
                    CASA standard time format, either a single time (e.g., '2012/03/03/12:00:00'
                    or a time range (e.g., '2012/03/03/12:00:00~2012/03/03/13:00:00'. If the latter,
                    take the midpoint of the timerange for reference. If no date specified, take
                    the date of the first scan
           polyfit: ONLY works for MS database with only one source with continously tracking;
                    not recommanded unless scan length is too long and want to have very high accuracy
           usephacenter: Bool -- if True, correct for the RA and DEC in the ms file based on solar empheris.
                                 Otherwise assume the phasecenter is correctly pointed to the solar disk center
                                 (EOVSA case)

       return value:
           helio: a list of VLA pointing information
                   reftimestr: reference time, in FITS format string
                   reftime: reference time, in mjd format
                   ra: actual RA of phasecenter in the ms file at the reference time (interpolated)
                   dec: actual DEC of phasecenter in the ms file at the reference time (interpolated)
                   # CASA uses only RA and DEC of the closest field (e.g. in clean) #
                   ra_fld: right ascention of the CASA reference pointing direction
                   dec_fld: declination of the CASA reference pointing direction
                   raoff: RA offset of the phasecenter in the ms file to solar center
                   decoff: DEC offset of the phasecenter in the ms file to solar center
                   refx: heliocentric X offset of the phasecenter in the ms file to solar center
                   refy: heliocentric Y offset of the phasecenter in the ms file to solar center
    ######## Example #########
         msfile='sun_C_20140910T221952-222952.10s.cal.ms'
         ephemfile='horizons_sun_20140910.radecp'
         ephem=vla_prep.read_horizons(ephemfile=ephemfile)
         msinfo=vla_prep.read_msinfo(msfile=msfile)
         polyfit=0
         reftime = '22:25:20~22:25:40'
    '''
    if not vis or not os.path.exists(vis):
        raise ValueError, 'Please provide information of the MS database!'
    if not ephem:
        ephem = read_horizons(vis=vis)
    if not msinfo:
        msinfo0 = read_msinfo(vis)
    else:
        if isinstance(msinfo, str):
            try:
                msinfo0 = np.load(msinfo)
            except:
                raise ValueError, 'The specified input msinfo file does not exist!'
        elif isinstance(msinfo, dict):
            msinfo0 = msinfo
        else:
            raise ValueError, 'msinfo should be either a numpy npz or a dictionary'
    print 'msinfo is derived from: ', msinfo0['vis']
    scans = msinfo0['scans']
    fieldids = msinfo0['fieldids']
    btimes = msinfo0['btimes']
    inttimes = msinfo0['inttimes']
    ras = msinfo0['ras']
    decs = msinfo0['decs']
    if 'observatory' in msinfo0.keys():
        if msinfo0['observatory'] == 'EOVSA':
            usephacenter = False
    if type(ras) is list:
        ra_rads = [ra['value'] for ra in ras]
    elif type(ras) is dict:
        ra_rads = ras['value']
    else:
        print('Type of msinfo0["ras"] unrecognized.')
        return 0
    if type(decs) is list:
        dec_rads = [dec['value'] for dec in decs]
    elif type(decs) is dict:
        dec_rads = decs['value']
    else:
        print('Type of msinfo0["decs"] unrecognized.')
    # fit 2nd order polynomial fits to the RAs and DECs #
    if polyfit:
        cra = np.polyfit(btimes, ra_rads, 2)
        cdec = np.polyfit(btimes, dec_rads, 2)

    # find out phase center infomation in ms according to the input time or timerange #
    if not reftime:
        raise ValueError, 'Please specify a reference time of the image'
    if type(reftime) == str:
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
            try:
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
            except:
                print 'Error in converting the input reftime: ' + str(reftime0) + '. Aborting...'
        else:
            # if reftime0 is specified as a single value
            try:
                tref_d = qa.getvalue(qa.convert(qa.totime(reftime0), 'd'))[0]
                # if no date is specified, add up the date of the first scan
                if tref_d < 1.:
                    tref_d += int(btimes[0])
                tbg_d = tref_d
                # use the intergration time
                # ind = bisect.bisect_left(btimes, tref_d)
                # if msinfo0['etimes']:
                #    tdur_s = inttimes[ind - 1]
                # else:
                #    tdur_s = np.mean(inttimes)
                tdur_s = 1.
            except:
                print 'Error in converting the input reftime: ' + str(reftime0) + '. Aborting...'
        helio0['reftime'] = tref_d
        # helio0['date-obs'] = qa.time(qa.quantity(tbg_d, 'd'), form='fits', prec=10)[0]
        # helio0['exptime'] = tdur_s

        # find out phase center RA and DEC in the measurement set according to the reference time
        # if polyfit, then use the 2nd order polynomial coeffs
        inttime = np.nanmean(inttimes)
        ind = bisect.bisect_left(btimes, tref_d)
        if ind > 1:
            dt = tref_d - btimes[ind - 1]
            if ind < len(btimes):
                scanlen = btimes[ind] - btimes[ind - 1]
                (ra_b, ra_e) = (ra_rads[ind - 1], ra_rads[ind])
                (dec_b, dec_e) = (dec_rads[ind - 1], dec_rads[ind])
            if ind >= len(btimes):
                scanlen = btimes[ind - 1] - btimes[ind - 2]
                (ra_b, ra_e) = (ra_rads[ind - 2], ra_rads[ind - 1])
                (dec_b, dec_e) = (dec_rads[ind - 2], dec_rads[ind - 1])
        if ind == 1:  # only one scan exists (e.g., imported from AIPS)
            ra_b = ra_rads[ind - 1]
            ra_e = ra_b
            dec_b = dec_rads[ind - 1]
            dec_e = dec_b
            scanlen = 10.  # radom value
            dt = 0.
        if ind < 1:
            if np.abs((tref_d - btimes[0]) * 24 * 3600) < inttime / 2.0:
                ind = 1
                ra_b = ra_rads[ind - 1]
                ra_e = ra_b
                dec_b = dec_rads[ind - 1]
                dec_e = dec_b
                scanlen = 10.  # radom value
                dt = 0.
            else:
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
        time0 = ephem['time']
        ra0 = ephem['ra']
        dec0 = ephem['dec']
        p0 = ephem['p0']
        delta0 = ephem['delta']
        if len(time0) > 1:
            ind = bisect.bisect_left(time0, tref_d)
            dt0 = time0[ind] - time0[ind - 1]
            dt_ref = tref_d - time0[ind - 1]
            dra0 = ra0[ind] - ra0[ind - 1]
            ddec0 = dec0[ind] - dec0[ind - 1]
            dp0 = p0[ind] - p0[ind - 1]
            ddelta0 = delta0[ind] - delta0[ind - 1]
            ra0 = ra0[ind - 1] + dra0 / dt0 * dt_ref
            dec0 = dec0[ind - 1] + ddec0 / dt0 * dt_ref
            p0 = p0[ind - 1] + dp0 / dt0 * dt_ref
            delta0 = delta0[ind - 1] + ddelta0 / dt0 * dt_ref
        else:
            try:
                ra0 = ra0[0]
                dec0 = dec0[0]
                p0 = p0[0]
                delta0 = delta0[0]
            except:
                print "Error in retrieving info from ephemeris!"
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
        if usephacenter:
            helio0['refx'] = refx
            helio0['refy'] = refy
        else:
            helio0['refx'] = 0.
            helio0['refy'] = 0.
        helio0['p0'] = p0
        # helio['r_sun']=np.degrees(R_sun.value/(au.value*delta0))*3600. #in arcsecs
        helio.append(helio0)
    return helio


def getbeam(imagefile=None, beamfile=None):
    ia = iatool()
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


def imreg(vis=None, ephem=None, msinfo=None, imagefile=None, timerange=None, reftime=None, fitsfile=None, beamfile=None,
          offsetfile=None, toTb=None,
          scl100=None, verbose=False, p_ang=False, overwrite=True, usephacenter=True, deletehistory=False):
    ''' 
    main routine to register CASA images
           Required Inputs:
               vis: STRING. CASA measurement set from which the image is derived
               imagefile: STRING or LIST. name of the input CASA image
               timerange: STRING or LIST. timerange used to generate the CASA image, must have the same length as the input images. 
                          Each element should be in CASA standard time format, e.g., '2012/03/03/12:00:00~2012/03/03/13:00:00'
           Optional Inputs:
               msinfo: DICTIONARY. CASA MS information, output from read_msinfo. If not provided, generate one from the supplied vis
               ephem: DICTIONARY. solar ephem, output from read_horizons. 
                      If not provided, query JPL Horizons based on time info of the vis (internet connection required)
               fitsfile: STRING or LIST. name of the output registered fits files
               reftime: STRING or LIST. Each element should be in CASA standard time format, e.g., '2012/03/03/12:00:00'
               offsetfile: optionally provide an offset with a series of solar x and y offsets with timestamps 
               toTb: Bool. Convert the default Jy/beam to brightness temperature?
               scl100: Bool. If True, scale the image values up by 100 (to compensate VLA 20 dB attenuator)
               verbose: Bool. Show more diagnostic info if True.
               usephacenter: Bool -- if True, correct for the RA and DEC in the ms file based on solar empheris.
                                     Otherwise assume the phasecenter is correctly pointed to the solar disk center
                                     (EOVSA case)
    '''
    ia = iatool()

    if deletehistory:
        msclearhistory(vis)
    if verbose:
        import time
        t0 = time.time()
        prtidx = 1
        print('point {}: {}'.format(prtidx, time.time() - t0))
        prtidx += 1

    if not imagefile:
        raise ValueError, 'Please specify input image'
    if not timerange:
        raise ValueError, 'Please specify timerange of the input image'
    if type(imagefile) == str:
        imagefile = [imagefile]
    if type(timerange) == str:
        timerange = [timerange]
    if not fitsfile:
        fitsfile = [img + '.fits' for img in imagefile]
    if type(fitsfile) == str:
        fitsfile = [fitsfile]
    nimg = len(imagefile)
    if len(timerange) != nimg:
        raise ValueError, 'Number of input images does not equal to number of timeranges!'
    if len(fitsfile) != nimg:
        raise ValueError, 'Number of input images does not equal to number of output fits files!'
    nimg = len(imagefile)
    if verbose:
        print str(nimg) + ' images to process...'

    if verbose:
        print('point {}: {}'.format(prtidx, time.time() - t0))
        prtidx += 1

    if reftime:  # use as reference time to find solar disk RA and DEC to register the image, but not the actual timerange associated with the image
        if type(reftime) == str:
            reftime = [reftime] * nimg
        if len(reftime) != nimg:
            raise ValueError, 'Number of reference times does not match that of input images!'
        helio = ephem_to_helio(vis, ephem=ephem, msinfo=msinfo, reftime=reftime, usephacenter=usephacenter)
    else:
        # use the supplied timerange to register the image
        helio = ephem_to_helio(vis, ephem=ephem, msinfo=msinfo, reftime=timerange, usephacenter=usephacenter)

    if verbose:
        print('point {}: {}'.format(prtidx, time.time() - t0))
        prtidx += 1

    for n, img in enumerate(imagefile):
        if verbose:
            print 'processing image #' + str(n)
        fitsf = fitsfile[n]
        timeran = timerange[n]
        # obtain duration of the image as FITS header exptime
        try:
            [tbg0, tend0] = timeran.split('~')
            tbg_d = qa.getvalue(qa.convert(qa.totime(tbg0), 'd'))[0]
            tend_d = qa.getvalue(qa.convert(qa.totime(tend0), 'd'))[0]
            tdur_s = (tend_d - tbg_d) * 3600. * 24.
            dateobs = qa.time(qa.quantity(tbg_d, 'd'), form='fits', prec=10)[0]
        except:
            print 'Error in converting the input timerange: ' + str(timeran) + '. Proceeding to the next image...'
            continue

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

        hel = helio[n]
        if not os.path.exists(img):
            raise ValueError, 'Please specify input image'
        if os.path.exists(fitsf) and not overwrite:
            raise ValueError, 'Specified fits file already exists and overwrite is set to False. Aborting...'
        else:
            p0 = hel['p0']
            ia.open(img)
            imr = ia.rotate(pa=str(-p0) + 'deg')
            imr.tofits(fitsf, history=False, overwrite=overwrite)
            imr.close()
            imsum = ia.summary()
            ia.close()

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

        # construct the standard fits header
        # RA and DEC of the reference pixel crpix1 and crpix2
        (imra, imdec) = (imsum['refval'][0], imsum['refval'][1])
        # find out the difference of the image center to the CASA phase center
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

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

        hdu = pyfits.open(fitsf, mode='update')

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

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
        header['date-obs'] = dateobs  # begin time of the image
        if not p_ang:
            hel['p0'] = 0
        try:
            # this works for pyfits version of CASA 4.7.0 but not CASA 4.6.0
            if tdur_s:
                header.set('exptime', tdur_s)
            else:
                header.set('exptime', 1.)
            header.set('p_angle', hel['p0'])
            header.set('dsun_obs', sun.sunearth_distance(Time(dateobs)).to(u.meter).value)
            header.set('rsun_obs', sun.solar_semidiameter_angular_size(Time(dateobs)).value)
            header.set('rsun_ref', sun.constants.radius.value)
            header.set('hgln_obs', 0.)
            header.set('hglt_obs', sun.heliographic_solar_center(Time(dateobs))[1].value)
        except:
            # this works for astropy.io.fits
            if tdur_s:
                header.append(('exptime', tdur_s))
            else:
                header.append(('exptime', 1.))
            header.append(('p_angle', hel['p0']))
            header.append(('dsun_obs', sun.sunearth_distance(Time(dateobs)).to(u.meter).value))
            header.append(('rsun_obs', sun.solar_semidiameter_angular_size(Time(dateobs)).value))
            header.append(('rsun_ref', sun.constants.radius.value))
            header.append(('hgln_obs', 0.))
            header.append(('hglt_obs', sun.heliographic_solar_center(Time(dateobs))[1].value))

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

        # update intensity units, i.e. to brightness temperature?
        if toTb:
            # get restoring beam info
            (bmajs, bmins, bpas, beamunits, bpaunits) = getbeam(imagefile=imagefile, beamfile=beamfile)
            bmaj = bmajs[n]
            bmin = bmins[n]
            beamunit = beamunits[n]
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
                header['BTYPE'] = 'Brightness Temperature'
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

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1

        hdu.flush()
        hdu.close()

        if verbose:
            print('point {}: {}'.format(prtidx, time.time() - t0))
            prtidx += 1
