import lmfit
import numpy as np
from suncasa.pyGSFIT import gsutils
from importlib import reload

reload(gsutils)

# msfile = 'IDB20210507_1840-1920XXYY.cal.10s.slfcaled.ms'
## SDO/AIA image fits
aiafile = 'aia.lev1_euv_12s.2021-05-07T190209Z.171.image_lev1.fits'
# aiafile = 'aia.lev1_uv_24s.2021-05-07T190214Z.1600.image_lev1.fits'
## EOVSA image fits created above
eofiles = ['EOVSA_20210507T190205.000000.outim.image.allbd.fits']
# eofiles = ['EOVSA_20210507T190135.000000.outim.image.allbd.fits','EOVSA_20210507T190205.000000.outim.image.allbd.fits']
spws = ['{}'.format(l) for l in range(48)]
## image center ("xycen"), pixel scale ("cell"), and image field of view ("fov")
## for the plots
xycen = np.array([-895, 290])
fov = np.array([200, 200])
## The alpha blending value for , between 0 (transparent) and 1 (opaque).
calpha = 0.5
## The respective maximum intensity of EOVSA images
clevels = np.array([0.3, 1.0])
freqghz_bound = [2., 14.0]
# freqghz_bound = [-1, 100]

gt = gsutils.GStool(eofiles, aiafile=aiafile, xycen=xycen, fov=fov, freqghz_bound=freqghz_bound,
                    calpha=calpha, clevels=clevels)

dep = 20.  # total source depth, 1e8 cm
Bmag = 30.
Bmagerr = 0.0
delta = 4.0
deltaerr = 0.0
# em = region.em
# nth = np.sqrt(em / (dep * 725.e5)) / 1e10
nth = 1.0
ntherr = 0.0
nrl = 1.0e6
nrlerr = 0.0
E_hi = 0.1
Emin = 0.1
Emax = 5.
nrlh = nrl * (E_hi ** (1. - delta) - Emax ** (1. - delta)) / (Emin ** (1. - delta) - Emax ** (1. - delta)) + 1e-8
nrlherr = 0.0
lognrlh = np.log10(nrlh)
lognrlherr = nrlherr / nrlh / np.log(10.)
Tth = 30.  # Thermal temperature in MK
# Tth = region.wT / 1e6  # Thermal temperature in MK
theta = 70.

params = lmfit.Parameters()
params.add('lnf', value=-2., vary=False)  # ln of fractional error adjustment

# params.add('ssz2', value=ssz2, vary=False)  # pixel size in arcsec
params.add('depth', value=dep, min=5.0, max=300., vary=False)  # Column depth, 1e8 cm (1 Mm)
params.add('Bmag', value=Bmag, min=20., max=100., vary=True)  # Magnetic field strength in G
params.add('Tth', value=Tth, min=5, max=30., vary=False)  # Thermal temperature in MK
params.add('nth', value=nth, min=0.1, max=20.0, vary=False)  # Thermal density in 1e10 cm^{-3}
params.add('lognrlh', value=np.log10(nrlh), min=1., max=10., vary=True)  # Non-thermal density in cm^{-3}
params.add('delta', value=delta, min=1.5, max=10.0, vary=True)  # Power law index for nonthermal electrons
params.add('theta', value=theta, min=20., max=85, vary=False)  # theta
params.add('Emin', value=Emin, min=0.001, max=0.05, vary=False)  # Emin
params.add('Emax', value=Emax, min=0.1, max=300, vary=False)  # Emax

gt.set_params(params=params)

gt.fit()
