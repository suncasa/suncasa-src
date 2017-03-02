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


def freqsfromfitsheader(header):
    '''
    get frequency in GHz from fits header
    :param header: fits header
    :return pols: polarisation stokes
    '''
    try:
        freqs = ['{:.3f}'.format(ll / 1e9) for ll in
                 (header["CRVAL3"] + np.arange(header["NAXIS3"]) * header["CDELT3"])]
        return freqs
    except:
        print("error in fits header!")
        raise ValueError


def transfitdict2DF(datain, gaussfit=True):
    '''
    convert the results from pimfit or pmaxfit tasks to pandas DataFrame structure.
    :param datain: The component list from pimfit or pmaxfit tasks
    :param gaussfit: True if the results is from pimfit, otherwise False.
    :return: the pandas DataFrame structure.
    '''
    import pandas as pd

    ra2arcsec = 180. * 3600. / np.pi
    dspecDF0 = pd.DataFrame()
    for ll in datain['timestamps']:
        tidx = datain['timestamps'].index(ll)
        if datain['succeeded'][tidx]:
            pols = datain['outputs'][tidx].keys()
            dspecDFlist = []
            dspecDF = pd.DataFrame()
            for ppit in pols:
                dspecDFtmp = pd.DataFrame()
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
                freqstrs = []
                fits_local = []
                for comp in datain['outputs'][tidx][ppit]['results'].keys():
                    if comp.startswith('component'):
                        if gaussfit:
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
                        freqstrs.append('{:.3f}'.format(
                            datain['outputs'][tidx][ppit]['results'][comp]['spectrum']['frequency']['m0']['value']))
                        fits_local.append(datain['imagenames'][tidx].split('/')[-1])
                if gaussfit:
                    dspecDFtmp['shape_latitude{}'.format(ppit)] = pd.Series(shape_latitude)
                    dspecDFtmp['shape_longitude{}'.format(ppit)] = pd.Series(shape_longitude)
                    dspecDFtmp['shape_latitude_err{}'.format(ppit)] = pd.Series(shape_latitude_err)
                    dspecDFtmp['shape_longitude_err{}'.format(ppit)] = pd.Series(shape_longitude_err)
                    dspecDFtmp['peak{}'.format(ppit)] = pd.Series(peak)
                    dspecDFtmp['shape_majoraxis{}'.format(ppit)] = pd.Series(shape_majoraxis)
                    dspecDFtmp['shape_minoraxis{}'.format(ppit)] = pd.Series(shape_minoraxis)
                    dspecDFtmp['shape_positionangle{}'.format(ppit)] = pd.Series(shape_positionangle)
                    dspecDFtmp['beam_major{}'.format(ppit)] = pd.Series(beam_major)
                    dspecDFtmp['beam_minor{}'.format(ppit)] = pd.Series(beam_minor)
                    dspecDFtmp['beam_positionangle{}'.format(ppit)] = pd.Series(beam_positionangle)
                else:
                    dspecDFtmp['shape_latitude{}'.format(ppit)] = pd.Series(shape_latitude)
                    dspecDFtmp['shape_longitude{}'.format(ppit)] = pd.Series(shape_longitude)
                    dspecDFtmp['shape_latitude_err{}'.format(ppit)] = pd.Series(shape_latitude_err)
                    dspecDFtmp['shape_longitude_err{}'.format(ppit)] = pd.Series(shape_longitude_err)
                    dspecDFtmp['peak{}'.format(ppit)] = pd.Series(peak)
                dspecDFtmp['freqstr'.format(ppit)] = pd.Series(freqstrs)
                dspecDFtmp['fits_local'.format(ppit)] = pd.Series(fits_local)
                dspecDFlist.append(dspecDFtmp)
            for DFidx, DFit in enumerate(dspecDFlist):
                if DFidx == 0:
                    dspecDF = dspecDFlist[0]
                else:
                    dspecDF = pd.merge(dspecDF.copy(), DFit, how='outer', on=['freqstr', 'fits_local'])
            dspecDF0 = dspecDF0.append(dspecDF, ignore_index=True)

    return dspecDF0


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
    colnlistcom = ['shape_latitude', 'shape_longitude', 'peak', 'shape_latitude_err', 'shape_longitude_err']
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


def c_correlate(a, v):
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)
    return np.correlate(a, v, mode='same')


def XcorrMap(z, x, y, doxscale = True):
    '''
    get the cross correlation map along y axis
    :param z: data
    :param x: x axis
    :param y: y axis
    :return:
    '''
    from scipy.interpolate import splev, splrep
    if doxscale:
        xfit = np.linspace(x[0], x[-1], 10 * len(x) + 1)
        zfit = np.zeros((len(y), len(xfit)))
        for yidx1, yq in enumerate(y):
            xx = x
            yy = z[yidx1, :]
            s = len(yy)  # (len(yy) - np.sqrt(2 * len(yy)))*2
            tck = splrep(xx, yy, s=s)
            ys = splev(xfit, tck)
            zfit[yidx1, :] = ys
    else:
        xfit=x
        zfit=z
    ny, nxfit = zfit.shape
    ccpeak = np.empty((ny - 1, ny - 1))
    ccpeak[:] = np.nan
    ccmax = ccpeak.copy()
    ya = ccpeak.copy()
    yv = ccpeak.copy()
    yidxa = ccpeak.copy()
    yidxv = ccpeak.copy()
    for idx1 in xrange(1, ny):
        for idx2 in xrange(0, idx1):
            lightcurve1 = zfit[idx1, :]
            lightcurve2 = zfit[idx2, :]
            ccval = c_correlate(lightcurve1, lightcurve2)
            if sum(lightcurve1) != 0 and sum(lightcurve2) != 0:
                cmax = np.amax(ccval)
                cpeak = np.argmax(ccval) - (nxfit - 1) / 2
            else:
                cmax = 0
                cpeak = 0
            ccmax[idx2, idx1 - 1] = cmax
            ccpeak[idx2, idx1 - 1] = cpeak
            ya[idx2, idx1 - 1] = y[idx1 - 1]
            yv[idx2, idx1 - 1] = y[idx2]
            yidxa[idx2, idx1 - 1] = idx1 - 1
            yidxv[idx2, idx1 - 1] = idx2
            if idx1 - 1 != idx2:
                ccmax[idx1 - 1, idx2] = cmax
                ccpeak[idx1 - 1, idx2] = cpeak
                ya[idx1 - 1, idx2] = y[idx2]
                yv[idx1 - 1, idx2] = y[idx1 - 1]
                yidxa[idx1 - 1, idx2] = idx2
                yidxv[idx1 - 1, idx2] = idx1 - 1

    return {'zfit': zfit, 'ccmax': ccmax, 'ccpeak': ccpeak, 'x': x, 'nx': len(x), 'xfit': xfit, 'nxfit': nxfit, 'y': y,
            'ny': ny, 'yv': yv, 'ya': ya, 'yidxv': yidxv, 'yidxa': yidxa}

