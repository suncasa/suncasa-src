def polsfromfitsheader(header):
    import numpy as np
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
    ra2arcsec = 180. * 3600. / np.pi

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
    # freq = []
    freqstrs = []
    fits_local = []

    for ll in datain['timestamps']:
        tidx = datain['timestamps'].index(ll)
        if datain['succeeded'][tidx]:
            pols = datain['outputs'][tidx].keys()
            for ff in xrange(datain['outputs'][tidx]['results']['nelements']):
                comp = 'component{}'.format(ff)
                if datain['outputs'][tidx]['results'][comp]['peak']['value'] > 0.0:
                    if gauss:
                        longitude = [datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['m0'][
                                         'value'] * ra2arcsec for ppit in pols]
                        latitude = [datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['m1'][
                                        'value'] * ra2arcsec for ppit in pols]
                        longitude_err = [
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['longitude'][
                                'value'] for ppit in pols]
                        latitude_err = [
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['latitude'][
                                'value'] for ppit in pols]
                        majoraxis = [datain['outputs'][tidx][ppit]['results'][comp]['shape']['majoraxis']['value'] for
                                     ppit in pols]
                        minoraxis = [datain['outputs'][tidx][ppit]['results'][comp]['shape']['minoraxis']['value'] for
                                     ppit in pols]
                        positionangle = [
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['positionangle']['value'] for ppit
                            in pols]
                        fluxpeak = [datain['outputs'][tidx][ppit]['results'][comp]['peak']['value'] for ppit in pols]
                        bmajor = [datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['major']['value']
                                  for ppit in pols]
                        bminor = [datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['minor']['value']
                                  for ppit in pols]
                        bpositionangle = [
                            datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['positionangle'][
                                'value'] for ppit in pols]
                        shape_majoraxis.append(majoraxis)
                        shape_minoraxis.append(minoraxis)
                        shape_positionangle.append(positionangle)
                        beam_major.append(bmajor)
                        beam_minor.append(bminor)
                        beam_positionangle.append(bpositionangle)
                    else:
                        longitude = [datain['outputs'][tidx][ppit][comp]['shape']['direction']['m0'][
                                         'value'] * ra2arcsec for ppit in pols]
                        latitude = [datain['outputs'][tidx][ppit][comp]['shape']['direction']['m1'][
                                        'value'] * ra2arcsec for ppit in pols]
                        longitude_err = [
                            datain['outputs'][tidx][ppit][comp]['shape']['direction']['error']['longitude']['value'] for
                            ppit in pols]
                        latitude_err = [
                            datain['outputs'][tidx][ppit][comp]['shape']['direction']['error']['latitude']['value'] for
                            ppit in pols]
                        fluxpeak = [datain['outputs'][tidx][ppit][comp]['flux']['value'] for ppit in pols]
                    shape_longitude.append(longitude)
                    shape_latitude.append(latitude)
                    shape_longitude_err.append(longitude_err)
                    shape_latitude_err.append(latitude_err)
                    peak.append(fluxpeak)
                    freqstrs.append(
                        '{:.3f}'.format(
                            datain['outputs'][tidx]['results'][comp]['spectrum']['frequency']['m0']['value']))
                    fits_local.append(datain['imagenames'][tidx].split('/')[-1])
