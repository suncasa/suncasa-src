import os, sys
import numpy as np
import sys
from math import *
import bisect
from astropy.time import Time
import astropy.units as u
import warnings
from suncasa.utils import fitsutils as fu
import ssl
from scipy.interpolate import interp1d

py3 = sys.version_info.major >= 3
if py3:
    ## CASA version >= 6
    from urllib.request import urlopen
else:
    ## CASA version < 6
    from urllib2 import urlopen

try:
    ## Full Installation of CASA 4, 5 and 6
    from taskinit import ms, tb, qa, iatool

    ia = iatool()
except:
    ## Modular Installation of CASA 6
    from casatools import ms as mstool
    from casatools import table, quanta, image

    ms = mstool()
    tb = table()
    qa = quanta()
    ia = image()

import sunpy

# check sunpy version
sunpy1 = sunpy.version.major >= 1
if sunpy1:
    from sunpy.coordinates import sun
else:
    from sunpy import sun

try:
    from astropy.io import fits as pyfits
except:
    try:
        import pyfits
    except ImportError:
        raise ImportError('Neither astropy nor pyfits exists in this CASA installation')


def ms_clearhistory(msfile):
    from taskinit import tb
    tb_history = msfile + '/HISTORY'
    os.system('cp -r {0} {0}_bk'.format(tb_history))
    tb.open(tb_history, nomodify=False)
    nrows = tb.nrows()
    if nrows > 0:
        tb.removerows(range(nrows))
    tb.close()


def normalise(angle):  ### convenience function to normalise angles by 2 pi
    ### angle in radians
    b = int(angle / (2 * np.pi))
    return angle - b * (2 * np.pi)


def ms_restorehistory(msfile):
    tb_history = msfile + '/HISTORY'
    os.system('mv {0}_bk {0}'.format(tb_history))


def read_horizons(t0=None, dur=None, vis=None, observatory=None, verbose=False):
    '''
    This function visits JPL Horizons to retrieve J2000 topocentric RA and DEC of the solar disk center
    as a function of time.

    Keyword arguments:
    t0: Referece time in astropy.Time format
    dur: duration of the returned coordinates in days. Default to 1 minute
    vis: CASA visibility dataset (in measurement set format). If provided, use entire duration from
         the visibility data
    observatory: observatory code (from JPL Horizons). If not provided, use information from visibility.
         if no visibility found, use earth center (code=500)
    verbose: True to provide extra information

    Usage:
    >>> from astropy.time import Time
    >>> out = read_horizons(t0=Time('2017-09-10 16:00:00'), observatory='-81')
    >>> out = read_horizons(vis = 'mydata.ms')

    History:
    BC (sometime in 2014): function was first wrote, followed by a number of edits by BC and SY
    BC (2019-07-16): Added docstring documentation

    '''
    if not t0 and not vis:
        t0 = Time.now()
    if not dur:
        dur = 1. / 60. / 24.  # default to 1 minute
    if t0:
        try:
            btime = Time(t0)
        except:
            print('input time ' + str(t0) + ' not recognized')
            return -1
    if vis:
        if not os.path.exists(vis):
            print('Input ms data ' + vis + ' does not exist! ')
            return -1
        try:
            # ms.open(vis)
            # summary = ms.summary()
            # ms.close()
            # btime = Time(summary['BeginTime'], format='mjd')
            # etime = Time(summary['EndTime'], format='mjd')
            ## alternative way to avoid conflicts with importeovsa, if needed -- more time consuming
            if observatory == 'geocentric' or observatory == '500':
                observatory = '500'
            else:
                ms.open(vis)
                metadata = ms.metadata()
                if metadata.observatorynames()[0] == 'EVLA':
                    observatory = '-5'
                elif metadata.observatorynames()[0] == 'EOVSA' or metadata.observatorynames()[0] == 'FASR':
                    observatory = '-81'
                elif metadata.observatorynames()[0] == 'ALMA':
                    observatory = '-7'
                metadata.close()
                ms.close()
            tb.open(vis)
            btime_vis = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
            etime_vis = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
            tb.close()
            if verbose:
                print("Beginning time of this scan " + btime_vis.iso)
                print("End time of this scan " + etime_vis.iso)

            # extend the start and end time for jpl horizons by 0.5 hr on each end
            btime = Time(btime_vis.mjd - 0.5 / 24., format='mjd')
            dur = etime_vis.mjd - btime_vis.mjd + 1.0 / 24.
        except:
            print('error in reading ms file: ' + vis + ' to obtain the ephemeris!')
            return -1

    # default the observatory to geocentric, if none provided
    if not observatory:
        observatory = '500'

    etime = Time(btime.mjd + dur, format='mjd')

    try:
        cmdstr = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&TABLE_TYPE='OBSERVER'&QUANTITIES='1,17,20'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'&CAL_FORMAT='BOTH'&SOLAR_ELONG='0,180'&CENTER='{}@399'&COMMAND='sun'&START_TIME='".format(
            observatory) + btime.iso.replace(' ', ',') + "'&STOP_TIME='" + etime.iso[:-4].replace(' ',
                                                                                                  ',') + "'&STEP_SIZE='1m'&SKIP_DAYLT='NO'&EXTRA_PREC='YES'&APPARENT='REFRACTED'"
        # cmdstr = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE='OBSERVER'&QUANTITIES='1,17,20'&CSV_FORMAT='YES'&ANG_FORMAT='DEG'&CAL_FORMAT='BOTH'&SOLAR_ELONG='0,180'&CENTER='{}@399'&COMMAND='10'&START_TIME='".format(
        #     observatory) + btime.iso.replace(' ', ',') + "'&STOP_TIME='" + etime.iso[:-4].replace(' ',
        #                                                                                           ',') + "'&STEP_SIZE='1m'&SKIP_DAYLT='NO'&EXTRA_PREC='YES'&APPARENT='REFRACTED'"
        cmdstr = cmdstr.replace("'", "%27")
        # print('1################')
        # print(cmdstr)
        try:
            context = ssl._create_unverified_context()
            f = urlopen(cmdstr, context=context)
        except:
            f = urlopen(cmdstr)
        lines = f.readlines()
        f.close()
    except:
        # todo use geocentric coordinate for the new VLA data
        import requests, collections
        params = collections.OrderedDict()
        # params['batch'] = '1'
        # params['TABLE_TYPE'] = "'OBSERVER'"
        # params['QUANTITIES'] = "'1,17,20'"
        # params['CSV_FORMAT'] = "'YES'"
        # params['ANG_FORMAT'] = "'DEG'"
        # params['CAL_FORMAT'] = "'BOTH'"
        # params['SOLAR_ELONG'] = "'0,180'"
        # if observatory == '500':
        #     params['CENTER'] = "'500'"
        # else:
        #     params['CENTER'] = "'{}@399'".format(observatory)
        # params['COMMAND'] = "'10'"
        # params['START_TIME'] = "'{}'".format(btime.iso[:-4].replace(' ', ','))
        # params['STOP_TIME'] = "'{}'".format(etime.iso[:-4].replace(' ', ','))
        # params['STEP_SIZE'] = "'1m'"
        # params['SKIP_DAYLT'] = "'NO'"
        # params['EXTRA_PREC'] = "'YES'"
        # params['APPAENT'] = "'REFRACTED'"
        # results = requests.get("https://ssd.jpl.nasa.gov/horizons_batch.cgi", params=params)

        params['EPHEM_TYPE'] = "'OBSERVER'"
        params['QUANTITIES'] = "'1,17,20'"
        params['CSV_FORMAT'] = "'YES'"
        params['ANG_FORMAT'] = "'DEG'"
        params['CAL_FORMAT'] = "'BOTH'"
        params['SOLAR_ELONG'] = "'0,180'"
        if observatory == '500':
            params['CENTER'] = "'500'"
        else:
            params['CENTER'] = "'{}@399'".format(observatory)
        params['COMMAND'] = "'sun'"
        params['START_TIME'] = "'{}'".format(btime.iso[:-4].replace(' ', ','))
        params['STOP_TIME'] = "'{}'".format(etime.iso[:-4].replace(' ', ','))
        params['STEP_SIZE'] = "'1m'"
        params['SKIP_DAYLT'] = "'NO'"
        params['EXTRA_PREC'] = "'YES'"
        params['APPARENT'] = "'REFRACTED'"
        results = requests.get("https://ssd.jpl.nasa.gov/api/horizons.api?format=text", params=params)
        # print('2################')
        # print(results)
        lines = [ll for ll in results.iter_lines()]

    # add a check for python 3
    if py3:
        lines = [l.decode('utf-8', 'backslashreplace') for l in lines]

    # print(lines)
    nline = len(lines)
    istart = 0
    for i in range(nline):
        if lines[i][0:5] == '$$SOE':  # start recording
            istart = i + 1
        if lines[i][0:5] == '$$EOE':  # end recording
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


def read_msinfo(vis=None, msinfofile=None, interp_to_scan=False, verbose=False):
    """
    Module ot read CASA measurement set and return RA, DEC, and time stamps of the phase center.
    Options for returning those from the emphemeris table (if available) or use those from the FIELD table.
    Parameters
    ----------
    vis: (required) path to the input CASA measurement set
    msinfofile: (optional) path/name of the saved numpy .npz file
    #use_ephem: (optional) if True (default), use the enclosed emphemeris table, otherwise use RA/DEC of the FIELD table
    interp_to_scan: (optional) if True, the entries of the emphemeris table
        are interpolated to the beginning of each scan. This is mainly for backward compatibility. Default is False.
    Returns
    -------
    msinfo: A dictionary contains necessary information for ephem_to_helio()
        vis: CASA measurement set
        observatory: Name of the observatory. May be used if they default to phase at the solar center.
        scans: summary info of the scans
        fieldids: a list of all FIELD ids
        scan_start_times: list of the start times of all scans, in mjd
        scan_end_times: list of the end times of all scans, in mjd
        btimes: time stamps used by ephem_to_helio() to determine the image/phase center shifts, in mjd
        ras: corresponding RA coordinates used by ephem_to_helio() to determine the image/phase center shifts, in rad
        decs: corresponding DEC coordinates used by ephem_to_helio() to determine the image/phase center shifts, in rad
    """
    import glob
    # read MS information #
    msinfo = dict.fromkeys(['vis', 'scans', 'fieldids', 'btimes', 'btimestr', 'inttimes', 'ras', 'decs', 'observatory'])
    ms.open(vis)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    if verbose:
        print('This measurement set has these fields: ', metadata.fieldnames())
    scans = ms.getscansummary()
    scanids = sorted(scans.keys(), key=lambda x: int(x))
    btimes_scan = []
    etimes_scan = []
    fieldids = []
    inttimes = []
    dirs = []
    ras_scan = []
    decs_scan = []

    for idx, scanid in enumerate(scanids):
        btimes_scan.append(scans[scanid]['0']['BeginTime'])
        etimes_scan.append(scans[scanid]['0']['EndTime'])
        fieldid = scans[scanid]['0']['FieldId']
        fieldids.append(fieldid)
        inttimes.append(scans[scanid]['0']['IntegrationTime'])
        dir = ms.getfielddirmeas('PHASE_DIR', fieldid)
        dirs.append(dir)
        ras_scan.append(dir['m0'])
        decs_scan.append(dir['m1'])

    # TODO: a caveat is that here we assume that there is only one ephemeris object in the ms
    if metadata.nfields() > 1:
        print('Warning!!! Multiple fields found in this measurement set. Will only select the first one!')

    metadata.close()
    ms.close()
    ephem_table = glob.glob(vis + '/FIELD/EPHEM*.tab')
    if len(ephem_table) == 1:
        if verbose:
            print('Found ephemeris table {}. If VLA, this is probably post-2014 data.'.format(ephem_table[0]))
            print('The visibility phase center is probably synchronized with the ephemeris. Loading data...')
        tb.open(ephem_table[0])
        col_ra = tb.getcol('RA')
        col_dec = tb.getcol('DEC')
        col_mjd = tb.getcol('MJD')
        tb.close()
        # check if the ephemeris has time that are 1 hour earlier than the measurement set
        idx, = np.where((col_mjd < (btimes_scan[0] - 1. / 24.)) | ((col_mjd > (etimes_scan[-1] + 1. / 24.))))
        if len(idx) > 0:
            print('Warning!!! I found that the ephemeris table has times that are outside that of the ms data.')
            print('Select those within +-1 hr of the ms data. Please proceed with caution...')
            col_ra = col_ra[idx]
            col_dec = col_dec[idx]
            col_mjd = col_mjd[idx]

        # Down-sample the ephemeris data to 60 s cadence (corresponding to <0.2" of shift due to solar rotation)
        dt_sec = np.median(np.diff(col_mjd)) * 24. * 3600.
        if dt_sec < 60.:
            if verbose:
                print('Ephemeris has a (median) cadence of {0:.3f}s. Down-sample to 60s.'.format(dt_sec))
            nt_res = int((col_mjd[-1] - col_mjd[0]) * 24. * 3600. / 60.)
            if nt_res > 2:
                btimes_res = np.linspace(col_mjd[0], col_mjd[-1], nt_res)
                f_ra = interp1d(col_mjd, col_ra, kind='linear')
                f_dec = interp1d(col_mjd, col_dec, kind='linear')
                ras_ = f_ra(np.array(btimes_res))
                decs_ = f_dec(np.array(btimes_res))
                ras = qa.convert(qa.quantity(ras_, 'deg'), 'rad')
                decs = qa.convert(qa.quantity(decs_, 'deg'), 'rad')
                btimes = btimes_res
        else:
            ras = qa.convert(qa.quantity(col_ra, 'deg'), 'rad')
            decs = qa.convert(qa.quantity(col_dec, 'deg'), 'rad')
            btimes = col_mjd

        ##TODO: I don't know why I was doing this. It seems it is never used.
        if interp_to_scan:
            f_ra = interp1d(col_mjd, col_ra, kind='linear')
            f_dec = interp1d(col_mjd, col_dec, kind='linear')
            ras_ = f_ra(np.array(btimes_scan))
            decs_ = f_dec(np.array(btimes_scan))
            ras = qa.convert(qa.quantity(ras_, 'deg'), 'rad')
            decs = qa.convert(qa.quantity(decs_, 'deg'), 'rad')
            btimes = btimes_scan
        msinfo['has_ephem_table'] = True
    else:
        if verbose:
            print('This measurement set does not have an ephemeris table attached. '
                  'The visibility phascenter uses RA and DEC of the FIELD attached to each scan')
        ras = ras_scan
        decs = decs_scan
        btimes = btimes_scan
        msinfo['has_ephem_table'] = False

    # put the relevent information into the dictionary
    btimestr = [qa.time(qa.quantity(btime, 'd'), form='fits', prec=10)[0] for btime in btimes]
    msinfo['vis'] = vis
    msinfo['observatory'] = observatory
    msinfo['scans'] = scans
    msinfo['fieldids'] = fieldids
    msinfo['inttimes'] = inttimes
    msinfo['scan_start_times'] = btimes_scan
    msinfo['scan_end_times'] = etimes_scan
    msinfo['btimes'] = btimes
    msinfo['btimestr'] = btimestr
    msinfo['ras'] = ras
    msinfo['decs'] = decs
    if msinfofile:
        np.savez(msinfofile, vis=vis, observatory=observatory, scans=scans, fieldids=fieldids, inttimes=inttimes,
                 btimes=btimes, btimestr=btimestr, ras=ras, decs=decs, has_ephem_table=has_ephem_table)
    return msinfo


def ephem_to_helio(vis=None, ephem=None, msinfo=None, reftime=None, dopolyfit=True,
                   usephacenter=True, geocentric=False, verbose=False):
    '''
    Module for calculating offsets of the phase center to the solar disk center using the following steps
    1. Take a solar ms database, read the scan and field information, find out the phase centers (in RA and DEC).
        This step is done with read_msinfo()
    2. Compare with the ephemeris of the solar disk center (in RA and DEC)
    3. Record RA/DEC of the phase center and offsets in RA/DEC and Helioprojective Cartesian coordinates (solar X/Y)
       inputs:
           msinfo: CASA MS information, output from read_msinfo
           ephem: solar ephem, output from read_horizons
           reftime: list of reference times (e.g., used for imaging)
                    CASA standard time format, either a single time (e.g., '2012/03/03/12:00:00'
                    or a time range (e.g., '2012/03/03/12:00:00~2012/03/03/13:00:00'. If the latter,
                    take the midpoint of the timerange for reference. If no date specified, take
                    the date of the first scan
           dopolyfit: Bool. Default: True. Works for MS database with only one source with continously tracking.
                    Disabled if usephacenter=False.
           usephacenter: Bool -- if True, correct for the RA and DEC in the ms file based on solar empheris.
                                 Otherwise assume the phasecenter is pointed to the solar disk center
                                 (EOVSA case)
           geocentric: Bool -- if True, use geocentric RA & DEC.
                        If False, use topocentric RA & DEC based on observatory location
                        Default: False

       return values:
           helio: a list of VLA pointing information
                   reftimestr: reference time, in FITS format string
                   reftime: reference time, in mjd format
                   ra: actual RA of phasecenter in the ms file at the reference time (interpolated)
                   dec: actual DEC of phasecenter in the ms file at the reference time (interpolated)
                   # CASA uses only RA and DEC of the closest field (e.g. in clean) #
                   ra_fld: RA of the CASA reference pointing direction, in radian
                   dec_fld: DEC of the CASA reference pointing direction, in radian
                   ra0: RA of the solar disk center, in radian
                   dec0: DEC of the solar disk center, in radian
                   raoff: RA offset of the phasecenter in the ms file to solar disk center, in arcsec
                   decoff: DEC offset of the phasecenter in the ms file to solar disk center, in arcsec
                   refx: heliocentric X offset of the phasecenter in the ms file to solar disk center, in arcsec
                   refy: heliocentric Y offset of the phasecenter in the ms file to solar disk center, in arcsec
                   has_ephem_table: flag to indicate if there is an ephemeris table attached.
    ######## Example #########
        from suncasa.utils import helioimage2fits as hf
        vis = vis = '22B-174_20221031_sun.1s.cal.ms'
        ephem = hf.ephem_to_helio(vis=vis, reftime='2022/10/31/20:37:10~2022/10/31/20:37:20', dopolyfit=True,
                        usephacenter=True, verbose=True)
        # Read out the ms information takes some time. To save time, one can read out the ms information first
            and supply the record here for registering multiple images. It will skip the read_msinfo() step.
        msinfo = hf.read_msinfo(vis=vis, verbose=True)
        ephem = hf.ephem_to_helio(vis=vis, msinfo=msinfo, reftime='2022/10/31/20:37:10~2022/10/31/20:37:20',
                                 dopolyfit=True, usephacenter=True, verbose=True)
    '''
    if not vis or not os.path.exists(vis):
        raise ValueError('Please provide information of the MS database!')

    if not msinfo:
        msinfo0 = read_msinfo(vis, verbose=verbose)
    else:
        if isinstance(msinfo, str):
            try:
                msinfo0 = np.load(msinfo)
            except:
                raise ValueError('The specified input msinfo file does not exist!')
        elif isinstance(msinfo, dict):
            msinfo0 = msinfo
        else:
            raise ValueError('msinfo should be either a numpy npz or a dictionary')
    if verbose:
        print('msinfo is derived from: {0:s}'.format(msinfo0['vis']))
    # scans = msinfo0['scans']
    fieldids = msinfo0['fieldids']
    btimes = msinfo0['btimes']
    btimestr = msinfo0['btimestr']
    inttimes = msinfo0['inttimes']
    ras = msinfo0['ras']
    decs = msinfo0['decs']
    scan_start_times = msinfo0['scan_start_times']
    if 'observatory' in msinfo0.keys():
        if msinfo0['observatory'] == 'EOVSA' or msinfo0['observatory'] == 'FASR':
            usephacenter = False

    if (not ephem) and usephacenter:
        if geocentric:
            # use geocentric ephemeris
            if verbose:
                print('Use geocentric RA & DEC')
            ephem = read_horizons(vis=vis, observatory='500')
        else:
            # use topocentric ephemeris
            if verbose:
                print('Use topocentric RA & DEC at observatory.')
            ephem = read_horizons(vis=vis, observatory=msinfo0['observatory'])

    if isinstance(ras, list):
        ra_rads = [ra['value'] for ra in ras]
    elif isinstance(ras, dict):
        ra_rads = ras['value']
    else:
        print('Type of msinfo0["ras"] unrecognized.')
        return 0
    if isinstance(decs, list):
        dec_rads = [dec['value'] for dec in decs]
    elif isinstance(decs, dict):
        dec_rads = decs['value']
    else:
        print('Type of msinfo0["decs"] unrecognized.')
    # fit 2nd order polynomial fits to the RAs and DECs #
    if dopolyfit:
        from suncasa.utils import fit_planet_position as fp
        if verbose:
            print('Use polynomial fit to interpolate the ms ephemeris data')
            print('Ephemeris start and end time', btimestr[0], btimestr[-1])
            print('Number of ephemeris data points used for polyfit', len(btimes))
        res = fp.fit_planet_positions(btimes, ra_rads, dec_rads, allowed_error=0.05)
        if res[0] == 0:  # success
            reftime_poly = res[1]
            cra = np.flip(np.array(res[2]))
            cdec = np.flip(np.array(res[3]))
            ra_diff = (np.polyval(cra, btimes - reftime_poly) - ra_rads) * 180. / np.pi * 3600.
            dec_diff = (np.polyval(cdec, btimes - reftime_poly) - dec_rads) * 180. / np.pi * 3600.
            if verbose:
                print('Maximum deviation of the polynomial fit in RA, DEC is {0:.2f}", {1:.2f}"'.format(
                    np.nanmax(ra_diff), np.max(dec_diff)))
        else:
            print('Polyfit cannot achieve the required accuracy. Use linear interpolation instead.')

    # find out phase center infomation in ms according to the input time or timerange #
    if not reftime:
        raise ValueError('Please specify a reference time of the image')
    if isinstance(reftime, str):
        reftime = [reftime]
    if (not isinstance(reftime, list)):
        print('input "reftime" is not a valid list. Abort...')

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
                print('Error in converting the input reftime: ' + str(reftime0) + '. Aborting...')
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
                print('Error in converting the input reftime: ' + str(reftime0) + '. Aborting...')
        helio0['reftime'] = tref_d
        # helio0['date-obs'] = qa.time(qa.quantity(tbg_d, 'd'), form='fits', prec=10)[0]
        # helio0['exptime'] = tdur_s

        # find out phase center RA and DEC in the measurement set according to the reference time
        # if dopolyfit, then use the polynomial fit to interpolate
        inttime = np.nanmean(inttimes)
        ind = bisect.bisect_left(btimes, tref_d)
        if ind > 1:
            # meaning the reference time is AFTER the start time of the ms ephemeris data
            dt = tref_d - btimes[ind - 1]
            if ind < len(btimes):
                # the reference time falls into the time range of the ms ephemeris data
                scanlen = btimes[ind] - btimes[ind - 1]
                (ra_b, ra_e) = (ra_rads[ind - 1], ra_rads[ind])
                (dec_b, dec_e) = (dec_rads[ind - 1], dec_rads[ind])
            if ind >= len(btimes):
                # the reference time falls after the last ephemeris record
                print('Warning!!! The provided reference time falls AFTER the last ephemeris record.')
                print('            I have to use the last record to register the image')
                scanlen = btimes[ind - 1] - btimes[ind - 2]
                (ra_b, ra_e) = (ra_rads[ind - 2], ra_rads[ind - 1])
                (dec_b, dec_e) = (dec_rads[ind - 2], dec_rads[ind - 1])
        if ind == 1:  # only one recorded exists (e.g., data imported from AIPS)
            print('Warning!!! Only one record exists in the ms ephemeris record.')
            print('           I have to use this only record to register the image.')
            ra_b = ra_rads[ind - 1]
            ra_e = ra_b
            dec_b = dec_rads[ind - 1]
            dec_e = dec_b
            scanlen = 10.  # radom value
            dt = 0.
        if ind < 1:
            print('Warning!!! The provided reference time falls BEFORE the first ephemeris record.')
            print('           Trying if it is within the integration time of the first record.')
            if np.abs((tref_d - btimes[0]) * 24 * 3600) < inttime / 2.0:
                ind = 1
                ra_b = ra_rads[ind - 1]
                ra_e = ra_b
                dec_b = dec_rads[ind - 1]
                dec_e = dec_b
                scanlen = 10.  # radom value
                dt = 0.
            else:
                raise ValueError('Reference time does not fall into the provided ephemeris record')
        if dopolyfit:
            ra = np.polyval(cra, tref_d - reftime_poly)
            dec = np.polyval(cdec, tref_d - reftime_poly)
        else:
            ra = ra_b + (ra_e - ra_b) / scanlen * dt
            dec = dec_b + (dec_e - dec_b) / scanlen * dt
        if ra < 0:
            ra += 2. * np.pi
        if ra_b < 0:
            ra_b += 2. * np.pi

        # compare with ephemeris from JPL Horizons
        if not usephacenter:
            # Do not need to read the information from the measurement set
            if msinfo0['observatory'] == 'EVLA':
                observatory_id = '-5'
            elif msinfo0['observatory'] == 'EOVSA' or msinfo0['observatory'] == 'FASR':
                observatory_id = '-81'
            elif msinfo0['observatory'] == 'ALMA':
                observatory_id = '-7'

            if not ephem:
                ephem = read_horizons(Time(tref_d, format='mjd'), observatory=observatory_id)

        times_ephem = ephem['time']
        ras_ephem = ephem['ra']
        decs_ephem = ephem['dec']
        p0s_ephem = ephem['p0']

        if usephacenter:
            if len(times_ephem) > 1:
                f_ra = interp1d(times_ephem, ras_ephem, kind='linear')
                f_dec = interp1d(times_ephem, decs_ephem, kind='linear')
                f_p0 = interp1d(times_ephem, p0s_ephem, kind='linear')
                ra0 = f_ra(tref_d)
                dec0 = f_dec(tref_d)
                p0 = f_p0(tref_d)
            else:
                try:
                    ra0 = ras_ephem[0]
                    dec0 = decs_ephem[0]
                    p0 = p0s_ephem[0]
                except:
                    print("Error in retrieving info from JPL solar ephemeris!")
            if ra0 < 0:
                ra0 += 2. * np.pi

            # RA and DEC offset in arcseconds
            decoff = degrees((dec - dec0)) * 3600.
            raoff = degrees((ra - ra0) * cos(dec)) * 3600.
            # Convert into heliocentric offsets
            prad = -radians(p0)
            refx = (-raoff) * cos(prad) - decoff * sin(prad)
            refy = (-raoff) * sin(prad) + decoff * cos(prad)
            helio0['raoff'] = raoff
            helio0['decoff'] = decoff
            helio0['refx'] = refx
            helio0['refy'] = refy
            helio0['ra0'] = ra0
            helio0['dec0'] = dec0
            helio0['p0'] = p0
        else:
            helio0['raoff'] = 0.
            helio0['decoff'] = 0.
            helio0['refx'] = 0.
            helio0['refy'] = 0.
            helio0['ra0'] = ras_ephem[0]
            helio0['dec0'] = decs_ephem[0]
            helio0['p0'] = p0s_ephem[0]

        helio0['ra'] = ra  # ra of the actual phase center
        helio0['dec'] = dec  # dec of the actual phase center
        ind = bisect.bisect_left(scan_start_times, tref_d)
        if ind > 1:
            dt = tref_d - scan_start_times[ind - 1]
            if ind < len(scan_start_times):
                fieldid = fieldids[ind]
                ms.open(vis)
                dir = ms.getfielddirmeas('PHASE_DIR', fieldid)
                ra_b = dir['m0']['value']
                dec_b = dir['m1']['value']
                if ra_b < 0:
                    ra_b += 2. * np.pi
                ms.close()
            if ind >= len(btimes):
                (ra_b, ra_e) = (ra_rads[ind - 2], ra_rads[ind - 1])
                (dec_b, dec_e) = (dec_rads[ind - 2], dec_rads[ind - 1])
        helio0['ra_fld'] = ra_b  # ra of the field, used as the reference in e.g., clean
        helio0['dec_fld'] = dec_b  # dec of the field, used as the refenrence in e.g., clean
        # helio['r_sun']=np.degrees(R_sun.value/(au.value*delta0))*3600. #in arcsecs
        ## add the tag whether an ephemeris table is attached or not
        helio0['has_ephem_table'] = msinfo0['has_ephem_table']
        helio.append(helio0)
    return helio


def getbeam(imagefile=None, beamfile=None):
    if not imagefile:
        raise ValueError('Please specify input images')
    bmaj = []
    bmin = []
    bpa = []
    beamunit = []
    bpaunit = []
    chans = []
    nimg = len(imagefile)
    for n in range(nimg):
        img = imagefile[n]
        bmaj_ = []
        bmin_ = []
        bpa_ = []
        if not os.path.exists(img):
            warnings.warn('{} does not exist!'.format(img))
            beamunit_ = []
            bpaunit_ = []
            chans_ = []
        else:
            ia.open(img)
            sum = ia.summary()
            ia.close()
            ia.done()
            if 'perplanebeams' in sum.keys():  # beam vary with frequency
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
            if 'restoringbeam' in sum.keys():  # only one beam
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
        print('Writing beam info to ascii file...')
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


def imreg(vis=None, imagefile=None, timerange=None,
          ephem=None, msinfo=None, fitsfile=None,
          usephacenter=True, geocentric=False, dopolyfit=True, reftime=None, offsetfile=None, beamfile=None,
          toTb=False, sclfactor=1.0, p_ang=False, overwrite=True,
          deletehistory=False, subregion='', docompress=False, verbose=False):
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
               toTb: Bool. Convert the default Jy/beam to brightness temperature?
               sclfactor: scale the image values up by its value (to compensate VLA 20 dB attenuator)
               verbose: Bool. Show more diagnostic info if True.
               usephacenter: Bool -- if True, correct for the RA and DEC in the ms file based on solar empheris.
                                     Otherwise assume the phasecenter is correctly pointed to the solar disk center
                                     (EOVSA case)
               geocentric: Bool -- if True, use geocentric RA & DEC.
                        If False, use topocentric RA & DEC based on observatory location
                        Default: False
               dopolyfit: Bool -- if True, fit the ephemeris from the measurement set using a polynomial fit.
                                  if False, just use linear interpolation
               ###### The following two parameters are only meant for temporary fixes ########################
               reftime: STRING or LIST. ONLY USED IF ANOTHER TIME (OTHER THAN TIME TO MAKE THE IMAGE)
                        IS NEEDED TO FIND RA AND DEC.
                        Each element should be in CASA standard time format, e.g., '2012/03/03/12:00:00'
               offsetfile: optionally provide an offset with a series of solar x and y offsets with timestamps
               ###############################################################################################
               subregion: Region selection. See 'help par.region' for details.
    Usage:
    >>> from suncasa.utils import helioimage2fits as hf
    >>> hf.imreg(vis='mydata.ms', imagefile='myimage.image', fitsfile='myimage.fits',
                 timerange='2017/08/21/20:21:10~2017/08/21/20:21:18')
    The output fits file is 'myimage.fits'

    History:
    BC (sometime in 2014): function was first wrote, followed by a number of edits by BC and SY
    BC (2019-07-16): Added checks for stokes parameter. Verified that for converting from Jy/beam to brightness temperature,
                     the convention of 2*k_b*T should always be used. I.e., for unpolarized source, stokes I, RR, LL, XX, YY, 
                     etc. in the output CASA images from (t)clean should all have same values of radio intensity 
                     (in Jy/beam) and brightness temperature (in K).

    '''

    if deletehistory:
        ms_clearhistory(vis)

    if not imagefile:
        raise ValueError('Please specify input image')
    if not timerange:
        raise ValueError('Please specify timerange of the input image')
    if isinstance(imagefile, str):
        imagefile = [imagefile]
    if isinstance(timerange, str):
        timerange = [timerange]
    if not fitsfile:
        fitsfile = [img + '.fits' for img in imagefile]
    if isinstance(fitsfile, str):
        fitsfile = [fitsfile]
    nimg = len(imagefile)
    if len(timerange) != nimg:
        raise ValueError('Number of input images does not equal to number of timeranges!')
    if len(fitsfile) != nimg:
        raise ValueError('Number of input images does not equal to number of output fits files!')
    nimg = len(imagefile)
    if verbose:
        print(str(nimg) + ' images to process...')

    if reftime:  # use as reference time to find solar disk RA and DEC to register the image, but not the actual timerange associated with the image
        if isinstance(reftime, str):
            reftime = [reftime] * nimg
        if len(reftime) != nimg:
            raise ValueError('Number of reference times does not match that of input images!')
        helio = ephem_to_helio(vis, ephem=ephem, msinfo=msinfo, reftime=reftime,
                               usephacenter=usephacenter, geocentric=geocentric, dopolyfit=dopolyfit, verbose=verbose)
    else:
        # use the supplied timerange to register the image
        helio = ephem_to_helio(vis, ephem=ephem, msinfo=msinfo, reftime=timerange,
                               usephacenter=usephacenter, geocentric=geocentric, dopolyfit=dopolyfit, verbose=verbose)

    if toTb:
        (bmajs, bmins, bpas, beamunits, bpaunits) = getbeam(imagefile=imagefile, beamfile=beamfile)

    for n, img in enumerate(imagefile):
        if verbose:
            print('processing image #' + str(n) + ' ' + img)
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
            print('Error in converting the input timerange: ' + str(timeran) + '. Proceeding to the next image...')
            continue

        hel = helio[n]
        if not os.path.exists(img):
            warnings.warn('{} does not existed!'.format(img))
        else:
            if os.path.exists(fitsf) and not overwrite:
                raise ValueError('Specified fits file already exists and overwrite is set to False. Aborting...')
            else:
                p0 = hel['p0']
                tb.open(img + '/logtable', nomodify=False)
                nobs = tb.nrows()
                tb.removerows([i + 1 for i in range(nobs - 1)])
                tb.close()
                ia.open(img)
                imr = ia.rotate(pa=str(-p0) + 'deg')
                if subregion != '':
                    imr = imr.subimage(region=subregion)
                imr.tofits(fitsf, history=False, overwrite=overwrite)
                imr.close()
                imsum = ia.summary()
                ia.close()
                ia.done()

            # construct the standard fits header
            # RA and DEC of the reference pixel crpix1 and crpix2
            (imra, imdec) = (imsum['refval'][0], imsum['refval'][1])

            ## When (t)clean is making an image, the default center of the image is coordinates of the associated FIELD
            ## If a new "phasecenter" is supplied in (t)clean, if
            #       CASE A (<2014-ish): No ephemeris table is attached. The visibility phase center
            #                              is the same as the RA and DEC of the FIELD
            #       CASA B (>2014-ish): Ephemeris table is attached. The visibility phase center
            #                              is the same as that from the ephemeris, and is
            #                              different from the RA and DEC of the FIELD.
            #### find out the difference between the image phase center and RA and DEC of the associated FIELD
            ddec_fld = degrees(normalise(imdec - hel['dec_fld'])) * 3600.  # in arcsec
            dra_fld = degrees(normalise(imra - hel['ra_fld']) * cos(hel['dec_fld'])) * 3600.  # in arcsec

            # Convert into image heliocentric offsets
            prad = -radians(hel['p0'])
            dx_fld = (-dra_fld) * cos(prad) - ddec_fld * sin(prad)
            dy_fld = (-dra_fld) * sin(prad) + ddec_fld * cos(prad)

            #### find out the difference between the image phase center and RA and DEC of the visibility phasecenter
            ddec_vis = degrees(normalise(imdec - hel['dec'])) * 3600.  # in arcsec
            dra_vis = degrees(normalise(imra - hel['ra']) * cos(hel['dec'])) * 3600.  # in arcsec
            # Convert into image heliocentric offsets
            dx_vis = (-dra_vis) * cos(prad) - ddec_vis * sin(prad)
            dy_vis = (-dra_vis) * sin(prad) + ddec_vis * cos(prad)
            if offsetfile:
                try:
                    offset = np.load(offsetfile)
                except:
                    raise ValueError('The specified offsetfile does not exist!')
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
                print('offset of image phase center to FIELD phase center (arcsec): dx={0:.2f}, dy={1:.2f}'.format(
                    dx_fld, dy_fld))
                print('offset of image phase center to VISIBILITY phase center (arcsec): dx={0:.2f}, dy={1:.2f}'.format(
                    dx_vis, dy_vis))
                print('offset of VISIBILITY phase center to solar disk center (arcsec): dx={0:.2f}, dy={1:.2f}'.format(
                    xoff, yoff))
            if hel['has_ephem_table']:
                if verbose:
                    print('This ms has an ephemeris table attached. Use ephemeris phase center as reference')
                (crval1, crval2) = (xoff + dx_vis, yoff + dy_vis)
            else:
                if verbose:
                    print('This ms does not have an ephemeris table attached. Use FIELD phase center as reference.')
                (crval1, crval2) = (xoff + dx_fld, yoff + dy_fld)
            # update the fits header to heliocentric coordinates

            hdu = pyfits.open(fitsf, mode='update')
            hdu[0].verify('fix')
            header = hdu[0].header
            dshape = hdu[0].data.shape
            ndim = hdu[0].data.ndim
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
            if tdur_s:
                exptime = tdur_s
            else:
                exptime = 1.
            p_angle = hel['p0']
            hgln_obs = 0.
            rsun_ref = sun.constants.radius.value
            if sunpy1:
                dsun_obs = sun.earth_distance(Time(dateobs)).to(u.meter).value
                rsun_obs = sun.angular_radius(Time(dateobs)).value
                hglt_obs = sun.B0(Time(dateobs)).value
            else:
                dsun_obs = sun.sunearth_distance(Time(dateobs)).to(u.meter).value
                rsun_obs = sun.solar_semidiameter_angular_size(Time(dateobs)).value
                hglt_obs = sun.heliographic_solar_center(Time(dateobs))[1].value
            try:
                # this works for pyfits version of CASA 4.7.0 but not CASA 4.6.0
                header.set('exptime', exptime)
                header.set('p_angle', p_angle)
                header.set('dsun_obs', dsun_obs)
                header.set('rsun_obs', rsun_obs)
                header.set('rsun_ref', rsun_ref)
                header.set('hgln_obs', hgln_obs)
                header.set('hglt_obs', hglt_obs)
            except:
                # this works for astropy.io.fits
                header.append(('exptime', exptime))
                header.append(('p_angle', p_angle))
                header.append(('dsun_obs', dsun_obs))
                header.append(('rsun_obs', rsun_obs))
                header.append(('rsun_ref', rsun_ref))
                header.append(('hgln_obs', hgln_obs))
                header.append(('hglt_obs', hglt_obs))

            # check if stokes parameter exist
            exist_stokes = False
            stokes_mapper = {1: 'I', 2: 'Q', 3: 'U', 4: 'V', -1: 'RR', -2: 'LL',
                             -3: 'RL', -4: 'LR', -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX'}
            if 'CRVAL3' in header.keys():
                if header['CTYPE3'] == 'STOKES':
                    stokenum = header['CRVAL3']
                    exist_stokes = True
            if 'CRVAL4' in header.keys():
                if header['CTYPE4'] == 'STOKES':
                    stokenum = header['CRVAL4']
                    exist_stokes = True
            if exist_stokes:
                if stokenum in stokes_mapper.keys():
                    stokesstr = stokes_mapper[stokenum]
                else:
                    print('Stokes parameter {0:d} not recognized. Assuming Stokes I'.format(stokenum))
                    stokenum = 1
                    stokesstr = 'I'
                if verbose:
                    print('This image is in Stokes ' + stokesstr)
            else:
                print('STOKES Information does not seem to exist! Assuming Stokes I')
                stokenum = 1

            # intensity units to brightness temperature
            if toTb:
                # get restoring beam info
                bmaj = bmajs[n]
                bmin = bmins[n]
                beamunit = beamunits[n]
                bpa = bpas[n]
                bpaunit = bpaunits[n]
                data = hdu[0].data  # remember the data order is reversed due to the FITS convension
                keys = list(header.keys())
                values = list(header.values())
                # which axis is frequency?
                faxis = keys[values.index('FREQ')][-1]
                faxis_ind = ndim - int(faxis)
                # find out the polarization of this image
                k_b = qa.constants('k')['value']
                c_l = qa.constants('c')['value']
                # Always use 2*kb for all polarizations
                const = 2. * k_b / c_l ** 2
                bpatmp = qa.quantity(bpa, bpaunit)['value'] - qa.convert(qa.quantity(p0, 'deg'), bpaunit)['value']
                header['BPA'] = bpatmp[0]
                if header['BUNIT'].lower() == 'jy/beam':
                    header['BUNIT'] = 'K'
                    header['BTYPE'] = 'Brightness Temperature'
                    for i in range(dshape[faxis_ind]):
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
                            bmin0 = np.radians(bmintmp / 3600.)
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
                        factor = const * nu ** 2  # SI unit
                        jy_to_si = 1e-26
                        # print(nu/1e9, beam_area, factor)
                        factor2 = sclfactor
                        # if sclfactor:
                        #     factor2 = 100.
                        if faxis == '3':
                            data[:, i, :, :] *= jy_to_si / beam_area / factor * factor2
                        if faxis == '4':
                            data[i, :, :, :] *= jy_to_si / beam_area / factor * factor2

            header = fu.headerfix(header)
            hdu.flush()
            hdu.close()

            if ndim - np.count_nonzero(np.array(dshape) == 1) > 3:
                docompress = False
                '''
                    Caveat: only 1D, 2D, or 3D images are currently supported by
                    the astropy fits compression. If a n-dimensional image data array
                    does not have at least n-3 single-dimensional entries,
                    force docompress to be False
                '''

                print(
                    'warning: The fits data contains more than 3 non squeezable dimensions. Skipping fits compression..')
            if docompress:
                fitsftmp = fitsf + ".tmp.fits"
                os.system("mv {} {}".format(fitsf, fitsftmp))
                hdu = pyfits.open(fitsftmp)
                hdu[0].verify('fix')
                header = hdu[0].header
                data = hdu[0].data
                fu.write_compressed_image_fits(fitsf, data, header, compression_type='RICE_1',
                                               quantize_level=4.0)
                os.system("rm -rf {}".format(fitsftmp))
    if deletehistory:
        ms_restorehistory(vis)
    return fitsfile


def calc_phasecenter_from_solxy(vis, timerange='', xycen=None, usemsphacenter=True, observatory=None):
    '''
    return the phase center in RA and DEC of a given solar coordinates

    :param vis: input measurement sets file
    :param timerange: can be a string or astropy.time.core.Time object, or a 2-element list of string or Time object
    :param xycen:  solar x-pos and y-pos in arcsec
    :param usemsphacenter:
    :return:
    phasecenter
    midtim: mid time of the given timerange
    '''
    ms.open(vis)
    metadata = ms.metadata()
    if not observatory:
        observatory = metadata.observatorynames()[0]

    try:
        mstrange = metadata.timerangeforobs(0)
        tst = Time(mstrange['begin']['m0']['value'], format='mjd')
        ted = Time(mstrange['end']['m0']['value'], format='mjd')
    except:
        print('Something is wrong in using metadata tool. Maybe you are using an old CASA version.')
        print('Try to use the summary outputs.')
        out = ms.summary()
        tst = Time(out['BeginTime'], format='mjd')
        ted = Time(out['EndTime'], format='mjd')
    ms.close()
    metadata.close()
    datstr = tst.iso[:10]

    if isinstance(timerange, Time):
        try:
            (sttim, edtim) = timerange
        except:
            sttim = timerange
            edtim = sttim
    else:
        if timerange == '':
            sttim = tst
            edtim = ted
        else:
            try:
                (tstart, tend) = timerange.split('~')
                if tstart[2] == ':':
                    sttim = Time(datstr + 'T' + tstart)
                    edtim = Time(datstr + 'T' + tend)
                    # timerange = '{0}/{1}~{0}/{2}'.format(datstr.replace('-', '/'), tstart, tend)
                else:
                    sttim = Time(qa.quantity(tstart, 'd')['value'], format='mjd')
                    edtim = Time(qa.quantity(tend, 'd')['value'], format='mjd')
            except:
                try:
                    if timerange[2] == ':':
                        sttim = Time(datstr + 'T' + timerange)
                        edtim = sttim
                    else:
                        sttim = Time(qa.quantity(timerange, 'd')['value'], format='mjd')
                        edtim = sttim
                except ValueError:
                    print("keyword 'timerange' in wrong format")

    midtim_mjd = (sttim.mjd + edtim.mjd) / 2.
    midtim = Time(midtim_mjd, format='mjd')
    eph = read_horizons(t0=midtim, observatory=observatory)
    if observatory == 'EOVSA' or (not usemsphacenter):
        print('Phasecenter in the ms is ignored')
        # use RA and DEC from FIELD ID 0
        tb.open(vis + '/FIELD')
        phadir = tb.getcol('PHASE_DIR').flatten()
        tb.close()
        ra0 = phadir[0]
        dec0 = phadir[1]
    else:
        ra0 = eph['ra'][0]
        dec0 = eph['dec'][0]

    if not xycen:
        # use solar disk center as default
        phasecenter = 'J2000 ' + str(ra0) + 'rad ' + str(dec0) + 'rad'
    else:
        x0 = np.radians(xycen[0] / 3600.)
        y0 = np.radians(xycen[1] / 3600.)
        p0 = np.radians(eph['p0'][0])  # p angle in radians
        raoff = -((x0) * np.cos(p0) - y0 * np.sin(p0)) / np.cos(eph['dec'][0])
        decoff = (x0) * np.sin(p0) + y0 * np.cos(p0)
        newra = ra0 + raoff
        newdec = dec0 + decoff
        phasecenter = 'J2000 ' + str(newra) + 'rad ' + str(newdec) + 'rad'
    return phasecenter, midtim
