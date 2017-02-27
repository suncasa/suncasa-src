import numpy as np


def polsfromfitsheader(header):
    '''
    get polarisation information from fits header
    :param header: fits header
    :return pols: polarisation stokes
    '''
    try:
        stokeslist = ['{}'.format(int(ll)) for ll in
                      (header["CRVAL4"] + np.arange(header["NAXIS4"]) * header["CDELT4"])]
        stokesdict = {'1': 'I', '2': 'Q', '3': 'U', '4': 'V', '-1': 'RR', '-2': 'LL', '-3': 'RL', '-4': 'LR',
                      '-5': 'XX', '-6': 'YY', '-7': 'XY', '-8': 'YX'}
        pols = map(lambda x: stokesdict[x], stokeslist)
    except:
        print("error in fits header!")
    return pols


def fitcltodf(datain, gauss=True):
    '''
    convert the results from pimfit or pmaxfit tasks to pandas DataFrame structure.
    :param datain: The component list from pimfit or pmaxfit tasks
    :param gauss: True if the results is from pimfit, otherwise False.
    :return: the pandas DataFrame structure.
    '''
    import pandas as pd

    ra2arcsec = 180. * 3600. / np.pi

    for ll in datain['timestamps']:
        tidx = datain['timestamps'].index(ll)
        if datain['succeeded'][tidx]:
            pols = datain['outputs'][tidx].keys()
            freqstrs = []
            fits_local = []
            for comp in datain['outputs'][tidx][pols[0]]['results'].keys():
                if comp.startswith('component'):
                    freqstrs.append('{:.3f}'.format(
                        datain['outputs'][tidx][pols[0]]['results'][comp]['spectrum']['frequency']['m0']['value']))
                    fits_local.append(datain['imagenames'][tidx].split('/')[-1])
            dspecDF = pd.DataFrame({'freqstr': freqstrs, 'fits_local': fits_local})
            for ppit in pols:
                shape_latitude = []
                shape_longitude = []
                shape_latitude_err = []
                shape_longitude_err = []
                shape_majoraxis = []
                shape_minoraxis = []
                shape_positionangle = []
                peak = []
                beam_major = []
                beam_minor = []
                beam_positionangle = []
                for comp in datain['outputs'][tidx][pols[0]]['results'].keys():
                    if comp.startswith('component'):
                        if gauss:
                            majoraxis = datain['outputs'][tidx][ppit]['results'][comp]['shape']['majoraxis']['value']
                            minoraxis = datain['outputs'][tidx][ppit]['results'][comp]['shape']['minoraxis']['value']
                            positionangle = datain['outputs'][tidx][ppit]['results'][comp]['shape']['positionangle'][
                                'value']
                            bmajor = datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['major'][
                                'value']
                            bminor = datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['minor'][
                                'value']
                            bpositionangle = \
                                datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['positionangle'][
                                    'value']
                            shape_majoraxis.append(majoraxis)
                            shape_minoraxis.append(minoraxis)
                            shape_positionangle.append(positionangle)
                            beam_major.append(bmajor)
                            beam_minor.append(bminor)
                            beam_positionangle.append(bpositionangle)
                            fluxpeak = datain['outputs'][tidx][ppit]['results'][comp]['peak']['value']
                        else:
                            fluxpeak = datain['outputs'][tidx][ppit]['results'][comp]['flux']['value'][0]
                        longitude = datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['m0'][
                                        'value'] * ra2arcsec
                        latitude = datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['m1'][
                                       'value'] * ra2arcsec
                        longitude_err = \
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['longitude'][
                                'value']
                        latitude_err = \
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['latitude'][
                                'value']
                        shape_longitude.append(longitude)
                        shape_latitude.append(latitude)
                        shape_longitude_err.append(longitude_err)
                        shape_latitude_err.append(latitude_err)
                        peak.append(fluxpeak)
                if gauss:
                    dspecDF['shape_latitude{}'.format(ppit)] = shape_latitude
                    dspecDF['shape_longitude{}'.format(ppit)] = shape_longitude
                    dspecDF['shape_latitude_err{}'.format(ppit)] = shape_latitude_err
                    dspecDF['shape_longitude_err{}'.format(ppit)] = shape_longitude_err
                    dspecDF['peak{}'.format(ppit)] = peak
                    dspecDF['shape_majoraxis{}'.format(ppit)] = shape_majoraxis
                    dspecDF['shape_minoraxis{}'.format(ppit)] = shape_minoraxis
                    dspecDF['shape_positionangle{}'.format(ppit)] = shape_positionangle
                    dspecDF['beam_major{}'.format(ppit)] = beam_major
                    dspecDF['beam_minor{}'.format(ppit)] = beam_minor
                    dspecDF['beam_positionangle{}'.format(ppit)] = beam_positionangle
                else:
                    dspecDF['shape_latitude{}'.format(ppit)] = shape_latitude
                    dspecDF['shape_longitude{}'.format(ppit)] = shape_longitude
                    dspecDF['shape_latitude_err{}'.format(ppit)] = shape_latitude_err
                    dspecDF['shape_longitude_err{}'.format(ppit)] = shape_longitude_err
                    dspecDF['peak{}'.format(ppit)] = peak

    return dspecDF


def getcolctinDF(dspecDF, col):
    '''
    return the count of how many times of the element starts with col occurs in columns of dspecDF
    :param dspecDF:
    :param col: the start string
    :return: the count and items
    '''
    itemset1 = col
    itemset2 = dspecDF.columns.tolist()
    items = []
    for ll in itemset2:
        if ll.startswith(itemset1):
            items.append(ll)
    itemct = [ll.startswith(itemset1) for ll in itemset2].count(True)
    columns = items
    return [itemct, columns]


def dspecDFfilter(dspecDF, pol):
    '''
    filter the unselect polarisation from dspecDF
    :param dspecDF: the original dspecDF
    :param pol: selected polarisation, dtype = string
    :return: the output dspecDF
    '''
    colnlistcom = ['shape_latitude', 'shape_longitude', 'peak','shape_latitude_err','shape_longitude_err']
    colnlistgaus = ['shape_majoraxis', 'shape_minoraxis', 'shape_positionangle', 'beam_major', 'beam_minor',
                    'beam_positionangle']
    ## above are the columns to filter
    colnlistess = dspecDF.columns.tolist()
    if getcolctinDF(dspecDF, 'peak')[0] > 1:
        for ll in colnlistcom + colnlistgaus:
            colinfo = getcolctinDF(dspecDF, ll)
            if colinfo[0] > 0:
                for cc in colinfo[1]:
                    if cc in colnlistess:
                        colnlistess.remove(cc)
        dspecDF1 = dspecDF.copy()[colnlistess]
        for ll in colnlistcom:
            dspecDF1[ll] = dspecDF.copy()[ll + pol]
        if getcolctinDF(dspecDF, 'shape_majoraxis')[0] > 0:
            for ll in colnlistgaus:
                dspecDF1[ll] = dspecDF.copy()[ll + pol]
        return dspecDF1
    else:
        return dspecDF
