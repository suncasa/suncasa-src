from astropy.time import Time
from sunpy import map as smap
import numpy as np
import astropy.units as u

def readfits(eofile):
    '''
    read eovsa image fits, adjust the date-obs to the mid time.
    :param eofile:
    :return:
    '''
    eomap_ = smap.Map(eofile)
    tbeg = Time(eomap_.date)
    newmeta = eomap_.meta
    tmid = Time(tbeg.mjd + newmeta['exptime'] / 24. / 3600 / 2, format='mjd')
    newmeta['date-obs'] = tmid.isot
    eomap = smap.Map(eomap_.data, newmeta)
    return eomap

def get_all_coordinate_from_map(sunmap):
    x, y = np.meshgrid(*[np.arange(v.value) for v in sunmap.dimensions]) * u.pix
    return sunmap.pixel_to_world(x, y)