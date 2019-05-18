from astropy.time import Time
from suncasa.utils import helioimage2fits as hf
from taskinit import tb, qa, ms
import numpy as np


def msclearhistory(msfile):
    '''Clears history in the a measurement sets file

    :param msfile: string
            The name of a measurement sets file

    :return:
    '''

    tb.open(msfile + '/HISTORY', nomodify=False)
    nrows = tb.nrows()
    if nrows > 0:
        tb.removerows(range(nrows))
    tb.close()


def calc_phasecenter_from_solxy(vis, timerange='', xycen=None, usemsphacenter=True):
    '''

    :param vis:
    :param timerange:
    :param xycen:
    :param usemsphacenter:
    :return:
    phasecenter
    midtim
    '''
    tb.open(vis + '/POINTING')
    tst = Time(tb.getcell('TIME_ORIGIN', 0) / 24. / 3600., format='mjd')
    ted = Time(tb.getcell('TIME_ORIGIN', tb.nrows() - 1) / 24. / 3600., format='mjd')
    tb.close()
    datstr = tst.iso[:10]

    if timerange == '':
        sttim = tst
        edtim = ted
        # timerange = '{0}~{1}'.format(tst.iso.replace('-', '/').replace(' ', '/'),
        #                              ted.iso.replace('-', '/').replace(' ', '/'))
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
        except ValueError:
            print("keyword 'timerange' in wrong format")

    ms.open(vis)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    ms.close()

    midtim_mjd = (sttim.mjd + edtim.mjd) / 2.
    midtim = Time(midtim_mjd, format='mjd')
    eph = hf.read_horizons(t0=midtim)
    if observatory == 'EOVSA' or (not usemsphacenter):
        print('This is EOVSA data')
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
    return phasecenter,midtim
