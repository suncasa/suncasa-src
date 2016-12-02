import json
import os
from  casac import *
import pickle
import numpy as np
import pandas as pd

database_dir = "${SUNCASADB}"
database_dir = os.path.expandvars(database_dir) + '/'
if os.path.exists('CASA_imfit_args.json'):
    with open('CASA_imfit_args.json', 'r') as fp:
        CASA_imfit_args = json.load(fp)
    for key, val in CASA_imfit_args.items():
        exec (key + '= {}'.format(val))

    if 'struct_id' in locals():
        structure_id = struct_id
    else:
        raise ValueError('define a struct_id!!!')
    print 'Script for imfit --- {} in {}'.format(structure_id, event_id)
    print ''
    if not 'ncpu' in locals():
        ncpu = 10

    default('pimfit')
    with open('CASA_imfit_args.json', 'r') as fp:
        CASA_imfit_args = json.load(fp)
    for key, val in CASA_imfit_args.items():
        exec (key + '= {}'.format(val))
    # inp(pimfit)
    out = pimfit()

    imgdir = database_dir + event_id + '/' + struct_id + '/Synthesis_Image/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    with open(imgdir + 'CASA_imfit_out', 'w') as fp:
        pickle.dump(out, fp)

    # todo add RR LL check. if stokes='RRLL', do pimfit separately with RR and LL.
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
    ra2arcsec = 180. * 3600. / np.pi
    for ll in out['timestamps']:
        tidx = out['timestamps'].index(ll)
        if out['succeeded'][tidx]:
            for ff in xrange(out['outputs'][tidx]['results']['nelements']):
                comp = 'component{}'.format(ff)
                if out['outputs'][tidx]['results'][comp]['peak']['value'] != 0.0:
                    shape_longitude.append(
                        out['outputs'][tidx]['results'][comp]['shape']['direction']['m0']['value'] * ra2arcsec)
                    shape_latitude.append(
                        out['outputs'][tidx]['results'][comp]['shape']['direction']['m1']['value'] * ra2arcsec)
                    shape_latitude_err.append(
                        out['outputs'][tidx]['results'][comp]['shape']['direction']['error']['latitude']['value'])
                    shape_longitude_err.append(
                        out['outputs'][tidx]['results'][comp]['shape']['direction']['error']['longitude']['value'])
                    shape_majoraxis.append(out['outputs'][tidx]['results'][comp]['shape']['majoraxis']['value'])
                    shape_minoraxis.append(out['outputs'][tidx]['results'][comp]['shape']['minoraxis']['value'])
                    shape_positionangle.append(out['outputs'][tidx]['results'][comp]['shape']['positionangle']['value'])
                    peak.append(out['outputs'][tidx]['results'][comp]['peak']['value'])
                    beam_major.append(out['outputs'][tidx]['results'][comp]['beam']['beamarcsec']['major']['value'])
                    beam_minor.append(out['outputs'][tidx]['results'][comp]['beam']['beamarcsec']['minor']['value'])
                    beam_positionangle.append(
                        out['outputs'][tidx]['results'][comp]['beam']['beamarcsec']['positionangle']['value'])
                    freqstrs.append(
                        '{:.3f}'.format(out['outputs'][tidx]['results'][comp]['spectrum']['frequency']['m0']['value']))
                    # freq.append(float(
                    #     '{:.3f}'.format(out['outputs'][tidx]['results'][comp]['spectrum']['frequency']['m0']['value'])))
                    fits_local.append(out['imagenames'][tidx].split('/')[-1])
    dspecDF2 = pd.DataFrame({
        'shape_latitude': shape_latitude,
        'shape_longitude': shape_longitude,
        'shape_latitude_err': shape_latitude_err,
        'shape_longitude_err': shape_longitude_err,
        'shape_majoraxis': shape_majoraxis,
        'shape_minoraxis': shape_minoraxis,
        'shape_positionangle': shape_positionangle,
        'peak': peak,
        'beam_major': beam_major,
        'beam_minor': beam_minor,
        'beam_positionangle': beam_positionangle,
        'freqstr': freqstrs,
        'fits_local': fits_local
    })
    with open(database_dir + event_id +'/'+ struct_id + '/dspecDF-save', 'rb') as fp:
        dspecDF1 = pickle.load(fp)
    for ll in dspecDF1.index:
        tmp = dspecDF1.loc[ll, 'freq']
        dspecDF1.loc[ll, 'freq'] = float('{:.3f}'.format(tmp))
    dspecDF = pd.merge(dspecDF1, dspecDF2, how='left', on=['freqstr', 'fits_local'])
    with open(database_dir + event_id +'/'+ struct_id + '/dspecDF-save', 'wb') as fp:
        pickle.dump(dspecDF,fp)

else:
    print 'CASA arguments config file not found!!'

    # CASA <43>: out.keys()
    #   Out[43]: ['timestamps', 'imagenames', 'succeeded', 'outputs']

    # CASA <44>: out['outputs'][0].keys()
    #   Out[44]: ['converged', 'deconvolved', 'results']

    # CASA <63>: out['outputs'][0]['deconvolved']['nelements']
    #   Out[63]: 256

    # CASA <45>: out['outputs'][0]['deconvolved'].keys()
    #   Out[45]:
    # ['component66',
    #  'component39',
    #  'component38',
    #  'component37',
    #  'component36',
    #  .............]

    # CASA <42>: out['outputs'][0]['deconvolved']['component50']
    #   Out[42]:
    # {'beam': {'beamarcsec': {'major': {'unit': 'arcsec',
    #                                    'value': 80.77617645263672},
    #                          'minor': {'unit': 'arcsec',
    #                                    'value': 54.635799407958984},
    #                          'positionangle': {'unit': 'deg',
    #                                            'value': -20.626220703125}},
    #           'beampixels': 200.02533455860916,
    #           'beamster': 1.175370395548373e-07},
    #  'flux': {'error': array([ 3.29035864,  0.        ,  0.        ,  0.        ]),
    #           'polarisation': 'Stokes',
    #           'unit': 'Jy',
    #           'value': array([ 22.37381579,   0.        ,   0.        ,   0.        ])},
    #  'ispoint': False,
    #  'label': '',
    #  'peak': {'error': 2.2380235668176343,
    #           'unit': 'Jy/beam',
    #           'value': 20.558514700697245},
    #  'shape': {'direction': {'error': {'latitude': {'unit': 'arcsec',
    #                                                 'value': 4.268624648944596},
    #                                    'longitude': {'unit': 'arcsec',
    #                                                  'value': 5.503895283292848}},
    #                          'm0': {'unit': 'rad',
    #                                 'value': -0.003066262239961642}, * 180/pi*3600
    #                          'm1': {'unit': 'rad',
    #                                 'value': -0.0006359638389362692}, * 180/pi*3600
    #                          'refer': 'J2000',
    #                          'type': 'direction'},
    #            'majoraxis': {'unit': 'arcsec', 'value': 101.22624749838064},
    #            'majoraxiserror': {'unit': 'arcsec', 'value': 16.8409753696521},
    #            'minoraxis': {'unit': 'arcsec', 'value': 47.44776658169024},
    #            'minoraxiserror': {'unit': 'arcsec', 'value': 20.731554584768958},
    #            'positionangle': {'unit': 'deg', 'value': 70.38611353707329},
    #            'positionangleerror': {'unit': 'deg', 'value': 20.07539719062812},
    #            'type': 'Gaussian'},
    #  'spectrum': {'channel': 50,
    #               'frequency': {'m0': {'unit': 'GHz', 'value': 1.096},
    #                             'refer': 'TOPO',
    #                             'type': 'frequency'},
    #               'type': 'Constant'},
    #  'sum': {'unit': 'Jy/beam', 'value': 2194.9042081832886}}

    # CASA <69>: out['outputs'][0]['results']['component50']
    #   Out[69]:
    # {'beam': {'beamarcsec': {'major': {'unit': 'arcsec',
    #                                    'value': 80.77617645263672},
    #                          'minor': {'unit': 'arcsec',
    #                                    'value': 54.635799407958984},
    #                          'positionangle': {'unit': 'deg',
    #                                            'value': -20.626220703125}},
    #           'beampixels': 200.02533455860916,
    #           'beamster': 1.175370395548373e-07},
    #  'flux': {'error': array([ 3.29035864,  0.        ,  0.        ,  0.        ]),
    #           'polarisation': 'Stokes',
    #           'unit': 'Jy',
    #           'value': array([ 22.37381579,   0.        ,   0.        ,   0.        ])},
    #  'ispoint': False,
    #  'label': '',
    #  'peak': {'error': 0.9975390755164335,
    #           'unit': 'Jy/beam',
    #           'value': 9.163407415626928},
    #  'shape': {'direction': {'error': {'latitude': {'unit': 'arcsec',
    #                                                 'value': 4.268624648944596},
    #                                    'longitude': {'unit': 'arcsec',
    #                                                  'value': 5.503895283292848}},
    #                          'm0': {'unit': 'rad',
    #                                 'value': -0.003066262239961642},
    #                          'm1': {'unit': 'rad',
    #                                 'value': -0.0006359638389362692},
    #                          'refer': 'J2000',
    #                          'type': 'direction'},
    #            'majoraxis': {'unit': 'arcsec', 'value': 115.03827860906479},
    #            'majoraxiserror': {'unit': 'arcsec', 'value': 13.29532969982013},
    #            'minoraxis': {'unit': 'arcsec', 'value': 93.67016307294118},
    #            'minoraxiserror': {'unit': 'arcsec', 'value': 9.604845436573132},
    #            'positionangle': {'unit': 'deg', 'value': 71.18955299149792},
    #            'positionangleerror': {'unit': 'deg', 'value': 20.07539719062812},
    #            'type': 'Gaussian'},
    #  'spectrum': {'channel': 50,
    #               'frequency': {'m0': {'unit': 'GHz', 'value': 1.096},
    #                             'refer': 'TOPO',
    #                             'type': 'frequency'},
    #               'type': 'Constant'},
    #  'sum': {'unit': 'Jy/beam', 'value': 2194.9042081832886}}
