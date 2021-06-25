from astropy.time import Time
import numpy as np


def make_ephem(vis, ephemfile=None):
    import urllib2, ssl
    from taskinit import tb
    quantities = ['1', '14', '15', '17', '19', '20', '24', '32']
    quantities = ','.join(quantities)
    tb.open(vis)
    btime = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
    etime = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
    tb.close()
    print "Beginning time of this scan " + btime.iso
    print "End time of this scan " + etime.iso

    btime = Time((btime.mjd - 1.0/60./24.), format='mjd')
    etime = Time((etime.mjd + 1.0/60./24.), format='mjd')
    startdate = btime.iso.replace(' ', ',')[:-7]
    enddate = etime.iso.replace(' ', ',')[:-7]
    cmd = ["COMMAND= '10'", "CENTER= '-5@399'", "MAKE_EPHEM= 'YES'", "TABLE_TYPE= 'OBSERVER'",
           "START_TIME= '%s'" % startdate, "STOP_TIME= '%s'" % enddate, "STEP_SIZE= '1m'", "CAL_FORMAT= 'CAL'",
           "TIME_DIGITS= 'MINUTES'", "ANG_FORMAT= 'DEG'", "OUT_UNITS= 'KM-S'", "RANGE_UNITS= 'AU'",
           "APPARENT= 'AIRLESS'", "SOLAR_ELONG= '0,180'", "SUPPRESS_RANGE_RATE= 'NO'", "SKIP_DAYLT= 'NO'",
           "EXTRA_PREC= 'NO'", "R_T_S_ONLY= 'NO'", "REF_SYSTEM= 'J2000'", "CSV_FORMAT= 'YES'", "OBJ_DATA= 'YES'",
           "TIME_DIGITS ='MIN'", "QUANTITIES= '{}'".format(quantities)]
    cmdstr = "http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&" + '&'.join(cmd)
    try:
        context = ssl._create_unverified_context()
        f = urllib2.urlopen(cmdstr, context=context)
    except:
        f = urllib2.urlopen(cmdstr)
    lines = f.readlines()
    f.close()
    istart = 0
    for i, l in enumerate(lines):
        if l[0:5] == '$$SOE':  # start recording
            istart = i + 1
        if l[0:5] == '$$EOE':  # end recording
            iend = i

    if not ephemfile:
        ephemfile = 'sun-ephem-geo.txt'
    with open(ephemfile, 'w') as fb:
        for i, l in enumerate(lines):
            if i == istart - 3:
                fb.write(
                    ' Date__(UT)__HR:MN     R.A.___(J2000.0)___DEC. Ob-lon Ob-lat Sl-lon Sl-lat   NP.ang   NP.dist               r        rdot            delta      deldot    S-T-O')
            if i >= istart and i < iend:
                l_s = l.split(',')
                l_s.pop(1)
                l_s.pop(1)
                fb.write(' '.join(l_s))
            else:
                fb.write(l)
                # with open(ephemfile,'w') as fb:
                #     for i,l in enumerate(lines):
                #         fb.write(l.replace('*m','').replace('*t',''))


def make_ephem_tb(vis, output_casa_table_name=None):
    import recipes.ephemerides.JPLephem_reader2 as jplreader
    import os.path
    ephemfile = 'sun-ephem.txt'
    make_ephem(vis, ephemfile=ephemfile)
    if os.path.isfile(ephemfile):
        output_dictionary = jplreader.readJPLephem(ephemfile)
        time_range_string = str(int(output_dictionary['earliest']['m0']['value'])) + '-' + str(
            int(output_dictionary['latest']['m0']['value'])) + 'dUTC'
        if not output_casa_table_name:
            output_casa_table_name = '_'.join([output_dictionary['NAME'], time_range_string, 'J2000.tab'])
        jplreader.ephem_dict_to_table(output_dictionary, output_casa_table_name)
        print('Ephemeris table written to {}'.format(output_casa_table_name))
    else:
        print "sun-ephem-geo.txt not found!"
