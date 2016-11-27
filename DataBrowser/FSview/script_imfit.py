import json
import os
from  casac import *
import pickle

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

# dspecDF_tmp = pd.DataFrame({'time': xx - xx[0],
#                             'freq': yy,
#                             'timestr': timestrs,
#                             'dspec': tab2_spec_plt.flatten(),
#                             'fits_local': fits_local,
#                             'fits_global': fits_global})
# 'x_pos': [ll[0][1] for ll in Gauss_params],
# 'y_pos': [ll[0][2] for ll in Gauss_params],
# 'x_width': [ll[0][3] for ll in Gauss_params],
# 'y_width': [ll[0][4] for ll in Gauss_params],
# 'amp_gaus': [ll[0][0] for ll in Gauss_params],
# 'theta': [ll[0][5] for ll in Gauss_params],
# 'amp_offset': [ll[0][6] for ll in Gauss_params]})


else:
    print 'CASA arguments config file not found!!'

# CASA <43>: out.keys()
#   Out[43]: ['timestamps', 'imagenames', 'succeeded', 'outputs']

# CASA <44>: out['outputs'][0].keys()
#   Out[44]: ['converged', 'deconvolved', 'results']

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
