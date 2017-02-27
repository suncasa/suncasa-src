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

    if Dgaussfit:
        default('pimfit')
        with open('CASA_imfit_args.json', 'r') as fp:
            CASA_imfit_args = json.load(fp)
        for key, val in CASA_imfit_args.items():
            exec (key + '= {}'.format(val))
        out = pimfit()
    else:
        default('pmaxfit')
        with open('CASA_imfit_args.json', 'r') as fp:
            CASA_imfit_args = json.load(fp)
        for key, val in CASA_imfit_args.items():
            exec (key + '= {}'.format(val))
        out = pmaxfit()

    imgdir = database_dir + event_id + '/' + struct_id + '/Synthesis_Image/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    with open(imgdir + 'CASA_imfit_out', 'w') as fp:
        pickle.dump(out, fp)

    # todo add RR LL check. if stokes='RRLL', do pimfit separately with RR and LL.
    # todo add deconvolved results
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
                if out['outputs'][tidx]['results'][comp]['peak']['value'] > 0.0:
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
    print 'imfit results saved to '+database_dir + event_id +'/'+ struct_id + '/dspecDF-save'

else:
    print 'CASA arguments config file not found!!'
