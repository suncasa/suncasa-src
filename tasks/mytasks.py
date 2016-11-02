#
# User defined tasks setup.
# Generated from buildmytask.
#

<<<<<<< HEAD
if sys.path[1] != '/local/software/suncasa/tasks':
  sys.path.insert(1, '/local/software/suncasa/tasks')
=======
if sys.path[1] != '/afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks':
  sys.path.insert(1, '/afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks')
>>>>>>> origin/master
from odict import odict
if not globals().has_key('mytasks') :
  mytasks = odict()

mytasks['pimfit'] = 'Fit one or more elliptical Gaussian components on an image region(s)'
mytasks['ptclean'] = 'Parallelized clean in consecutive time steps'
mytasks['subvs'] = 'Vector-subtraction in UV using selected time ranges and spectral channels as background'

if not globals().has_key('task_location') :
  task_location = odict()

<<<<<<< HEAD
task_location['pimfit'] = '/local/software/suncasa/tasks'
task_location['ptclean'] = '/local/software/suncasa/tasks'
task_location['subvs'] = '/local/software/suncasa/tasks'
=======
task_location['pimfit'] = '/afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks'
task_location['ptclean'] = '/afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks'
task_location['subvs'] = '/afs/cad.njit.edu/research/physics/binchen/1/bchen/Dropbox/bc_python/suncasa/tasks'
>>>>>>> origin/master
import inspect
myglobals = sys._getframe(len(inspect.stack())-1).f_globals
tasksum = myglobals['tasksum'] 
for key in mytasks.keys() :
  tasksum[key] = mytasks[key]

from pimfit_cli import  pimfit_cli as pimfit
from ptclean_cli import  ptclean_cli as ptclean
from subvs_cli import  subvs_cli as subvs
