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
                if sswmp[idx] == 'arcsecs':
                    header['CUNIT1'] = 'arcsec'
                else:
                    header['CUNIT1'] = sswmp[idx]
            elif k == 'L0':
                header['HGLN_OBS'] = sswmp[idx]
            elif k == 'YC':
                header['CRVAL2'] = sswmp[idx]
            elif k == 'DY':
                header['CDELT2'] = sswmp[idx]
            elif k == 'RSUN_OBS':
                header['CDELT2'] = sswmp[idx]
            elif k == 'YUNITS':
                if sswmp[idx] == 'arcsecs':
                    header['CUNIT2'] = 'arcsec'
                else:
                    header['CUNIT2'] = sswmp[idx]
            elif k == 'B0':
                header['HGLT_OBS'] = sswmp[idx]
            elif k == 'ID':
                if 'SDO AIA' in sswmp[idx]:
                    header['TELESCOP'] = 'SDO AIA'
                elif 'SDO HMI' in sswmp[idx]:
                    header['TELESCOP'] = 'SDO HMI'
                elif 'EOVSA' in sswmp[idx]:
                    header['TELESCOP'] = 'EOVSA'
                elif 'VLA' in sswmp[idx]:
                    header['TELESCOP'] = 'VLA'
                else:
                    header['TELESCOP'] = sswmp[idx]
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
    elif  isinstance(sswmap, np.record):
        return mp2smap(sswmap)

