import os
import numpy as np
import vla_prep
import shutil
import pdb
import glob
from scipy.signal import medfilt2d
import json

with open(path_database + 'StrID_list.dat', 'rb') as f:
    StructureIdList = pickle.load(f)

for ll in StructureIdList.index:
    tmp = StructureIdList.iloc[ll]['timeran']
    StructureIdList.iloc[ll]['timeran'] = [tmp]
    tmp = StructureIdList.iloc[ll]['freqran']
    StructureIdList.iloc[ll]['freqran'] = [tmp]
    tmp = StructureIdList.iloc[ll]['str_id']
    StructureIdList.iloc[ll]['str_id'] = [tmp]
    tmp = StructureIdList.iloc[ll]['date']
    StructureIdList.iloc[ll]['date'] = [tmp]

database_dir = './database/'
in_json = glob.glob(database_dir + '*.json')
# Reading data back
with open(in_json[0], 'r') as f:
    StructureId = json.load(f)
tim0, tim1 = StructureId['time'][0:2]

with open('dspecDF-save', 'rb') as f:
    dspecDF = pickle.load(f)
for ll in dspecDF.index:
    tmp = dspecDF.iloc[ll]['thumbnail']
    tmp = 'FSview/' + '/'.join(tmp.split('/')[1:])
    dspecDF.loc[ll, 'thumbnail'] = tmp
with open('dspecDF-save', 'wb') as f:
    pickle.dump(dspecDF, f)

thumbnail = [vla_local_thumbnailpath + ll[0:-4] + 'small.jpg' for ll in fits_local]
with open('dspecDF-save', 'rb') as f:
    dspecDF = pickle.load(f)
for ll in dspecDF.index:
    tmp = timestrs[ll]
    dspecDF.loc[ll, 'timestr'] = tmp
for ll in dspecDF.index:
    tmp = fits_local[ll]
    dspecDF.loc[ll, 'fits_local'] = tmp
    tmp = fits_global[ll]
    dspecDF.loc[ll, 'fits_global'] = tmp
    tmp = thumbnail[ll]
    dspecDF.loc[ll, 'thumbnail'] = tmp

for ll in dspecDF.index:
    # tmp = dspecDF.iloc[ll]['fits_local']
    # idx = tmp.find('U01_2014-11-01') + 21
    # tmp = tmp[:idx]+'.'+tmp[idx:]
    # dspecDF.loc[ll, 'fits_local'] = tmp
    # tmp = dspecDF.iloc[ll]['fits_global']
    # idx = tmp.find('U01_2014-11-01') + 21
    # tmp = tmp[:idx]+'.'+tmp[idx:]
    # dspecDF.loc[ll, 'fits_global'] = tmp
    # tmp = dspecDF.iloc[ll]['thumbnail']
    tmp = dspecDF.loc[ll, 'thumbnail'].replace('QLook', 'FSview')
    dspecDF.loc[ll, 'thumbnail'] = tmp

dspecDFtmp.loc[np.asarray([76, 77, 79]), 'x-dummy'] = dspecDF.loc[np.asarray([76, 77, 79]), 'x_pos']

nrows = len(dspecDF.index)
dspecDF['x_dummy'] = pd.Series([None] * nrows, index=dspecDF.index)
dspecDF['y_dummy'] = pd.Series([None] * nrows, index=dspecDF.index)
dspecDF['amp_dummy'] = pd.Series([None] * nrows, index=dspecDF.index)

with open('dspecDF-save', 'wb') as f:
    pickle.dump(dspecDF, f)

# # Writing JSON data
# with open('data.json', 'w') as f:
#     json.dump(data, f)

clearcal(vis=msfile)
tb.open(msfile, nomodify=False)
colnames = tb.colnames()
cols2rm = ["MODEL_DATA", 'CORRECTED_DATA']
for l in range(len(cols2rm)):
    if cols2rm[l] in colnames:
        tb.removecols(cols2rm[l])
tb.close()

#  clean :: Invert and deconvolve images with selected algorithm
vis = 'SUN01_20141101.T163940-164700.50ms.cal.ms'  # Name of input visibility file
imagename = 'U04.test1'  # Pre-name of output images
outlierfile = ''  # Text file with image names, sizes, centers for outliers
field = ''  # Field Name or id
spw = '0~3'  # Spectral windows e.g. '0~3', '' is all
selectdata = True  # Other data selection parameters
timerange = '16:46:21.425~16:46:21.475'  # Range of time to select from data
uvrange = ''  # Select data within uvrange
antenna = ''  # Select data based on antenna/baseline
scan = ''  # Scan number range
observation = ''  # Observation ID range
intent = ''  # Scan Intent(s)

mode = 'channel'  # Spectral gridding type (mfs, channel, velocity, frequency)
nchan = -1  # Number of channels (planes) in output image; -1 = all
start = ''  # Begin the output cube at the frequency of this channel in the MS
width = 1  # Width of output channel relative to MS channel (# to average)
interpolation = 'linear'  # Spectral interpolation (nearest, linear, cubic).
resmooth = False  # Re-restore the cube image to a common beam when True
chaniter = False  # Clean each channel to completion (True), or all channels each cycle (False)
outframe = ''  # default spectral frame of output image

gridmode = ''  # Gridding kernel for FFT-based transforms, default='' None
niter = 500  # Maximum number of iterations
gain = 0.1  # Loop gain for cleaning
threshold = '0.0mJy'  # Flux level to stop cleaning, must include units: '1.0mJy'
psfmode = 'clark'  # Method of PSF calculation to use during minor cycles
imagermode = 'csclean'  # Options: 'csclean' or 'mosaic', '', uses psfmode
cyclefactor = 1.5  # Controls how often major cycles are done. (e.g. 5 for frequently)
cyclespeedup = -1  # Cycle threshold doubles in this number of iterations

multiscale = []  # Deconvolution scales (pixels); [] = standard clean
interactive = True  # Use interactive clean (with GUI viewer)
npercycle = 100  # Clean iterations before interactive prompt (can be changed)

mask = []  # Cleanbox(es), mask image(s), region(s), or a level
imsize = [512, 512]  # x and y image size in pixels.  Single value: same for both
cell = ['4.0arcsec']  # x and y cell size(s). Default unit arcsec.
phasecenter = 'J2000 14h26m22.7351 -14d29m29.801'  # Image center: direction or field index
restfreq = ''  # Rest frequency to assign to image (see help)
stokes = 'LL'  # Stokes params to image (eg I,IV,IQ,IQUV)
weighting = 'natural'  # Weighting of uv (natural, uniform, briggs, ...)
uvtaper = False  # Apply additional uv tapering of visibilities
modelimage = ''  # Name of model image(s) to initialize cleaning
restoringbeam = ['']  # Output Gaussian restoring beam for CLEAN image
pbcor = False  # Output primary beam-corrected image
minpb = 0.2  # Minimum PB level to use
usescratch = True  # True if to save model visibilities in MODEL_DATA column
allowchunk = False  # Divide large image cubes into channel chunks for deconvolution

for t0 in tim2:
    timestr0 = jdutil.jd_to_datetime(t0 / 3600. / 24.)
    # print timestr0.microsecond,'{:03d}'.format(int(round(timestr0.microsecond / 1e3)))
    timestr = timestr0.strftime('%Y-%m-%dT%H%M%S') + '.{:03d}'.format(int(round(timestr0.microsecond / 1e3)))
    print timestr

tofits = True
if tofits:
    import json
    import numpy as np
    import glob
    import os
    from  casac import *
    import pickle
    from astropy.time import Time

    database_dir = "${SUNCASADB}"
    database_dir = os.path.expandvars(database_dir) + '/'
    if os.path.exists('CASA_CLN_args.json'):
        with open('CASA_CLN_args.json', 'r') as fp:
            CASA_CLN_args = json.load(fp)
        for key, val in CASA_CLN_args.items():
            exec (key + '= {}'.format(val))
    import suncasa.vla.vla_prep as vla_prep

    ephem = vla_prep.read_horizons(ephemfile=ephemfile)
    imagenames = glob.glob('slfcal/U04/*.image')
    vlafits = ['/'.join(img.split('/')[0:-1]) + '/' + img.split('/')[-1][0:-6] + '.fits' for img in imagenames]
    timeran_tmp = [img.split('/')[-1][0:-6] for img in imagenames]
    # timeran = ['{}:{}:{}.{:03d}~{}:{}:{}.{:03d}'.format(ll[9:11],ll[11:13],ll[13:15],int(ll[16:])-25,ll[9:11],ll[11:13],ll[13:15],int(ll[16:])+25) for ll in timeran_tmp]
    timeran0 = [
        '{}-{}-{}T{}:{}:{}.{:03d}'.format(ll[0:4], ll[4:6], ll[6:8], ll[9:11], ll[11:13], ll[13:15], int(ll[16:])) for
        ll in timeran_tmp]
    timeran1 = Time(timeran0, format='isot', scale='utc')
    dt = 25. / 1000 / 24 / 3600
    timeranjd0 = timeran1.jd - dt
    timeranjd1 = timeran1.jd + dt
    timeran0 = list(Time(timeranjd0, format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso)
    timeran1 = list(Time(timeranjd1, format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso)
    timeran0 = [ll.replace(' ', 'T') for ll in timeran0]
    timeran1 = [ll.replace(' ', 'T') for ll in timeran1]
    timeran=['{}~{}'.format(ll[0],ll[1]) for ll in zip(timeran0,timeran1)]
    reftime = timeran
    helio = vla_prep.ephem_to_helio(msinfo=msinfofile, ephem=ephem, reftime=reftime)
    vla_prep.imreg(imagefile=imagenames, fitsfile=vlafits, helio=helio, toTb=False, scl100=True)
