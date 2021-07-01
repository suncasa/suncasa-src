import multiprocessing
from suncasa.utils import qlookplot as ql

msfile = 'IDB20210507_1840-1920XXYY.cal.10s.slfcaled.ms'
specfile = msfile + '.dspec.npz'
## set the time interval
timerange = '2021/05/07/19:02:00~2021/05/07/19:02:30'
## Bin width for time averaging
twidth = 1
## frequency range
spw = ['{}'.format(l) for l in range(48)]
## image center for clean in solar X-Y in arcsec
xycen = [-900, 280]
## pixel scale
cell = ['2.5arcsec']
## number of pixels in X and Y. If only one value is provided, NX = NY
imsize = [128]
## field of view of the zoomed-in panels in unit of arcsec
fov = [300, 300]
## select stokes XX
stokes = 'XX'
## set True if make a movie
mkmovie = True
## set True if generate compressed fits
docompress = True
## numbers of CPU threads for computing
ncpu = multiprocessing.cpu_count()
## set False to plot EOVSA images as open contours, set True to plot as filled contours
opencontour = False
clevels = [0.9, 1.0]

#### ---- Control knobs for AIA plotting ---- ####
## set True if plot AIA images as the background
plotaia = True
## AIA passband. The options are [171,131,304,335,211,193,94,1600,1700]
aiawave = 1600

ql.qlookplot(vis=msfile, specfile=specfile, timerange=timerange, spw=spw,
             ncpu=ncpu, xycen=xycen, imsize=imsize, fov=fov, cell=cell,
             restoringbeam=['30arcsec'],
             opencontour=opencontour, clevels=clevels,
             plotaia=plotaia, aiawave=aiawave,
             mkmovie=mkmovie, twidth=twidth, docompress=docompress, stokes=stokes, overwrite=False)

