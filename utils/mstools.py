# from astropy.time import Time
# from suncasa.utils import helioimage2fits as hf
from taskinit import tb, qa, ms
# import numpy as np


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
