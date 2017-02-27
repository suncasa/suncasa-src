import numpy as np
__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

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
                for comp in datain['outputs'][tidx][ppit]['results'].keys():
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







def get_contour_data(X, Y, Z):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from bokeh.models import (ColumnDataSource)
    cs = plt.contour(X, Y, Z, levels=(np.arange(5, 10, 2) / 10.0 * np.nanmax(Z)).tolist(), cmap=cm.Greys_r)
    xs = []
    ys = []
    xt = []
    yt = []
    col = []
    text = []
    isolevelid = 0
    for isolevel in cs.collections:
        # isocol = isolevel.get_color()[0]
        # thecol = 3 * [None]
        theiso = '{:.0f}%'.format(cs.get_array()[isolevelid] / Z.max() * 100)
        isolevelid += 1
        # for i in xrange(3):
        # thecol[i] = int(255 * isocol[i])
        thecol = '#%02x%02x%02x' % (220, 220, 220)
        # thecol = '#03FFF9'

        for path in isolevel.get_paths():
            v = path.vertices
            x = v[:, 0]
            y = v[:, 1]
            xs.append(x.tolist())
            ys.append(y.tolist())
            xt.append(x[len(x) / 2])
            yt.append(y[len(y) / 2])
            text.append(theiso)
            col.append(thecol)

    source = ColumnDataSource(data={'xs': xs, 'ys': ys, 'line_color': col, 'xt': xt, 'yt': yt, 'text': text})
    return source


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
    return g.ravel()


# def maxfit(image):
#     # NAME:
#     #    maxfit
#     # PURPOSE:
#     #    find maximum position of an image using 2D Gaussian fit
#     # INPUTS:
#     #    image: sunpy
#
#     data = image.data
#     nx, ny = data.shape
#     dx, dy = image.scale.x.value, image.scale.y.value
#     xc, yc = image.center.x.value, image.center.y.value
#     mapx, mapy = (np.linspace(0, nx - 1, nx) - image.reference_pixel.x.value + 1 + 0.5) * dx + xc, (
#         np.linspace(0, ny - 1, ny) - image.reference_pixel.y.value + 1 + 0.5) * dy + yc
#     mapx, mapy = np.meshgrid(mapx, mapy)
#     idxmax = np.where(data == np.amax(data))
#     idxhm = np.where(data >= np.amax(data) / 2)
#     hmfw = np.sqrt(len(idxhm[0]) / np.pi) * dx
#     xmax, ymax = mapx[idxmax[0][0], idxmax[1][0]], mapy[idxmax[0][0], idxmax[1][0]]
#     theta = np.arctan(
#         abs(idxhm[1] - image.reference_pixel.y.value).sum() / abs(idxhm[0] - image.reference_pixel.x.value).sum())
#     initial_guess = (np.amax(data), xmax, ymax, hmfw, hmfw, theta, data.std())
#
#     try:
#         popt, pcov = opt.curve_fit(twoD_Gaussian, (mapx, mapy), data.ravel(), p0=initial_guess)
#     except:
#         popt = np.empty((7))
#         popt[:] = np.nan
#
#     return popt

