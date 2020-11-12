from dateutil import parser
from astropy.time import Time
from sunpy import map as smap
from scipy.io import readsav
import numpy as np

def idlsav2sunmap(idlsavfile):

    smaps_struc = readsav(idlsavfile)
    smaps_dict = {}
    for smap_key in smaps_struc.keys():
        smap_struc = smaps_struc[smap_key]
        smaps_dict[smap_key] = []
        smap_struc = np.squeeze(smap_struc)
        for tidx, mp in enumerate(smap_struc):
            header = {}
            hd_keys = mp.dtype.names
            for idx, k in enumerate(hd_keys):
                if k == 'DATA':
                    data = mp[k]
                elif k == 'DUR':
                    header['EXPTIME'] = mp[k]
                elif k == 'TIME':
                    header['DATE-OBS'] = Time(parser.parse(mp[k])).isot
                    header['DATE'] = header['DATE-OBS']
                elif k == 'XC':
                    header['CRVAL1'] = mp[k]
                elif k == 'DX':
                    header['CDELT1'] = mp[k]
                elif k == 'XUNITS':
                    header['CUNIT1'] = mp[k]
                elif k == 'L0':
                    header['HGLN_OBS'] = mp[k]
                elif k == 'YC':
                    header['CRVAL2'] = mp[k]
                elif k == 'DY':
                    header['CDELT2'] = mp[k]
                elif k == 'RSUN_OBS':
                    header['CDELT2'] = mp[k]
                elif k == 'YUNITS':
                    header['CUNIT2'] = mp[k]
                elif k == 'B0':
                    header['HGLT_OBS'] = mp[k]
                elif k == 'FREQ':
                    header['WAVELNTH'] = mp[k]
                elif k == 'FREQUNIT':
                    header['WAVEUNIT'] = mp[k]
                else:
                    header[k] = mp[k]
            header['CTYPE1'] = 'HPLN-TAN'
            header['CTYPE2'] = 'HPLT-TAN'
            header['CRPIX2'], header['CRPIX1'] = (np.array(data.shape) + 1.0) / 2 + 0.5
            header['TELESCOP'] = 'EOVSA'
            header['RSUN_REF'] = 695508000.0
            header['P_ANGLE'] = 0.0
            header['NAXIS'] = 2
            header['NAXIS2'], header['NAXIS1'] = data.shape
            sunmap = smap.Map(data, header)
            smaps_dict[smap_key].append(sunmap)
    return smaps_dict