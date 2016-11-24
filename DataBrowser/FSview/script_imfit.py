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


else:
    print 'CASA arguments config file not found!!'
