import json
import os
from  casac import *
import pickle

# import numpy as np
# import pandas as pd
# from suncasa.utils import DButil

database_dir = "${SUNCASADB}"
database_dir = os.path.expandvars(database_dir) + '/'
if os.path.exists('CASA_imfit_args.json'):
    with open('CASA_imfit_args.json', 'r') as fp:
        CASA_imfit_args = json.load(fp)
    for key, val in CASA_imfit_args.items():
        exec (key + '= {}'.format(val))

    if 'struct_id' in locals():
        structure_id = struct_id
        if 'clean_id' in locals():
            if 'imfit_id' in locals():
                imfitIDdir = database_dir + event_id + '/' + struct_id + '/' + clean_id + '/' + imfit_id + '/'
            else:
                raise ValueError('define a imfit_id!!!')
        else:
            raise ValueError('define a clean_id!!!')
    else:
        raise ValueError('define a struct_id!!!')
    print 'Script for imfit --- {} in {}'.format(structure_id, event_id)
    print ''
    if not 'ncpu' in locals():
        ncpu = 10

    if gaussfit:
        default('pimfit')
        with open('CASA_imfit_args.json', 'r') as fp:
            CASA_imfit_args = json.load(fp)
        for key, val in CASA_imfit_args.items():
            exec (key + '= {}'.format(val))
        out = pimfit()
        os.system('cp pimfit.last {}/'.format(imfitIDdir))
    else:
        default('pmaxfit')
        with open('CASA_imfit_args.json', 'r') as fp:
            CASA_imfit_args = json.load(fp)
        for key, val in CASA_imfit_args.items():
            exec (key + '= {}'.format(val))
        out = pmaxfit()
        os.system('cp pmaxfit.last {}/'.format(imfitIDdir))

    if not os.path.exists(imfitIDdir):
        os.mkdir(imfitIDdir)
    with open(imfitIDdir + 'CASA_imfit_out', 'w') as fp:
        pickle.dump(out, fp)

        # todo add deconvolved results
        # dspecDF2 = DButil.transfitdict2DF(out, gaussfit=gaussfit)
        # with open(imfitIDdir + '/dspecDF-save', 'rb') as fp:
        #     dspecDF1 = pickle.load(fp)
        # for ll in dspecDF1.index:
        #     tmp = dspecDF1.loc[ll, 'freq']
        #     dspecDF1.loc[ll, 'freq'] = float('{:.3f}'.format(tmp))
        # dspecDF = pd.merge(dspecDF1, dspecDF2, how='left', on=['freqstr', 'fits_local'])
        # with open(imfitIDdir + '/dspecDF-save', 'wb') as fp:
        #     pickle.dump(dspecDF, fp)
        # print 'imfit results saved to ' + imfitIDdir + '/dspecDF-save'

else:
    print 'CASA arguments config file not found!!'
