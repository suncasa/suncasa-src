import sunpy.map as smap
from dateutil import parser
from astropy.time import Time
import numpy as np


def map2smap(sswmap):
    def mp2smap(sswmp):
        header = {}
        hd_keys = sswmp.dtype.names
        for idx, k in enumerate(hd_keys):
            if k == 'DATA':
                data = sswmp[idx]
            elif k == 'DUR':
                header['EXPTIME'] = sswmp[idx]
            elif k == 'TIME':
                header['DATE-OBS'] = Time(parser.parse(sswmp[idx])).isot
                header['DATE'] = header['DATE-OBS']
            elif k == 'XC':
                header['CRVAL1'] = sswmp[idx]
            elif k == 'DX':
                header['CDELT1'] = sswmp[idx]
            elif k == 'XUNITS':
                cunit1 = sswmp[idx]
                if isinstance(cunit1, bytes):
                    cunit1 = cunit1.decode('utf-8')
                if cunit1 == 'arcsecs':
                    header['CUNIT1'] = 'arcsec'
                else:
                    header['CUNIT1'] = cunit1
            elif k == 'L0':
                header['HGLN_OBS'] = sswmp[idx]
            elif k == 'YC':
                header['CRVAL2'] = sswmp[idx]
            elif k == 'DY':
                header['CDELT2'] = sswmp[idx]
            elif k == 'RSUN_OBS':
                header['CDELT2'] = sswmp[idx]
            elif k == 'YUNITS':
                cunit1 = sswmp[idx]
                if isinstance(cunit1, bytes):
                    cunit1 = cunit1.decode('utf-8')
                if cunit1 == 'arcsecs':
                    header['CUNIT2'] = 'arcsec'
                else:
                    header['CUNIT2'] = cunit1
            elif k == 'B0':
                header['HGLT_OBS'] = sswmp[idx]
            elif k == 'ID':
                telescop = sswmp[idx]
                if isinstance(telescop, bytes):
                    telescop = telescop.decode('utf-8')
                if 'SDO AIA' in telescop:
                    header['TELESCOP'] = 'SDO AIA'
                elif 'SDO HMI' in telescop:
                    header['TELESCOP'] = 'SDO HMI'
                elif 'EOVSA' in telescop:
                    header['TELESCOP'] = 'EOVSA'
                elif 'VLA' in telescop:
                    header['TELESCOP'] = 'VLA'
                else:
                    header['TELESCOP'] = telescop
            else:
                header[k] = sswmp[idx]
        header['CTYPE1'] = 'HPLN-TAN'
        header['CTYPE2'] = 'HPLT-TAN'
        header['CRPIX2'], header['CRPIX1'] = (np.array(data.shape) + 1.0) / 2
        # header['TELESCOP'] = 'EOVSA'
        header['RSUN_REF'] = 695508000.0
        header['P_ANGLE'] = 0.0
        return smap.Map(data, header)

    if isinstance(sswmap, np.recarray):
        res = []
        for mp in sswmap:
            res.append(mp2smap(mp))
        return res
    elif isinstance(sswmap, np.record):
        return mp2smap(sswmap)
