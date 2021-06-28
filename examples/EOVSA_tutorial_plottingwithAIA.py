import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from suncasa.utils import fitsutils as fu
from sunpy import map as smap
import copy

## SDO/AIA image fits
aiafile = 'aia.lev1_euv_12s.2021-05-07T190209Z.171.image_lev1.fits'
# aiafile = 'aia.lev1_uv_24s.2021-05-07T190214Z.1600.image_lev1.fits'
## EOVSA image fits created above
eofile = 'EOVSA_20210507T190205.000000.outim.image.allbd.fits'

## image center ("xycen"), pixel scale ("cell"), and image field of view ("fov")
## for the plots
xycen = np.array([-895, 290])
fov = np.array([200, 200])
xlim = xycen[0] + np.array([-1, 1]) * 0.5 * fov[0]
ylim = xycen[1] + np.array([-1, 1]) * 0.5 * fov[1]
## Select the lower 48 spectral windows for plotting
fidxs = np.arange(0, 48)
## The alpha blending value for , between 0 (transparent) and 1 (opaque).
alpha = 0.5
## The respective maximum intensity of EOVSA images
vmin = 0.3


## Get solar coordinates of every pixels in the radio image.
rmap, ndim, npol_fits, stokaxis, rcfreqs, rdata, rheader = fu.read_compressed_image_fits(eofile)
nspw = len(rcfreqs)
eodate = Time(rmap.date.mjd + rmap.exposure_time.value / 2. / 24 / 3600, format='mjd')
ny, nx = rmap.data.shape
x0, x1 = (np.array([1, rmap.meta['NAXIS1']]) - rmap.meta['CRPIX1']) * rmap.meta['CDELT1'] + \
         rmap.meta['CRVAL1']
y0, y1 = (np.array([1, rmap.meta['NAXIS2']]) - rmap.meta['CRPIX2']) * rmap.meta['CDELT2'] + \
         rmap.meta['CRVAL2']
dx = rmap.meta['CDELT1']
dy = rmap.meta['CDELT2']
mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)

## Plot the AIA maps on the background
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.8), sharex=True, sharey=True)
aiacmap = plt.get_cmap('gray_r')
aiacmap.set_under('white')
aiacmap.set_over('k')
aiamap = smap.Map(aiafile)
ax = axs[0]
aiamap.plot(axes=ax, cmap=aiacmap)
ax = axs[1]
aiamap.plot(axes=ax, cmap=aiacmap)

## Plot EOVSA images as filled contour on top of the AIA image
vmins = np.array([vmin] * len(rcfreqs))
spwcmap = plt.get_cmap('RdYlBu')
colors_spws = spwcmap(np.linspace(0, 1, nspw))
cts = []

for idx, fidx in enumerate(fidxs):
    data = rdata[0, fidx, ...]
    clevels1 = np.linspace(vmins[idx], 1.0, 2) * np.nanmax(data)
    cts.append(ax.contourf(mapx, mapy, data, levels=clevels1,
                           # colors=[colors_spws[fidx]] * len(clevels1), ## contour
                           colors=[colors_spws[fidx]],  ## contourf
                           alpha=alpha))
ax.set_xlim(xlim)
ax.set_ylim(ylim)
for ax in axs:
    ax.set_xlabel('Solar-X [arcsec]')
    ax.set_ylabel('Solar-y [arcsec]')
    ax.set_title('')
ax.text(0.02, 0.98, 'SOL2020-11-29 M4.4 Event', ha='left', va='top', color='k', transform=ax.transAxes, fontsize=12,
        fontweight='bold')
ax.text(0.02, 0.01, ' '.join(['AIA {:.0f} Å'.format(aiamap.wavelength.value), aiamap.date.datetime.strftime('%Y-%m-%dT%H:%M:%S')]), ha='left',
        va='bottom',
        color='k', transform=ax.transAxes)
ax.text(0.02, 0.05, ' '.join(['EOVSA     ', eodate.datetime.strftime('%Y-%m-%dT%H:%M:%S')]), ha='left', va='bottom',
        color='k', transform=ax.transAxes)

## Add the color bar for the EOVSA multi-frequency images
ax1_pos = ax.get_position().extents
caxwidth = (ax1_pos[2] - ax1_pos[0]) * 0.1
caxcenter = ax1_pos[2] + caxwidth * 0.75
cayheight = ax1_pos[3] - 0.06 - ax1_pos[1]
rcfreqsplt = rcfreqs
bdwidth = np.nanmean(np.diff(rcfreqsplt))
bounds = np.linspace(rcfreqsplt[0] - bdwidth / 2.0, rcfreqsplt[-1] + bdwidth / 2.0, len(rcfreqsplt) + 1)
ticks = rcfreqsplt[1:][::3]
cax = plt.axes((caxcenter - caxwidth / 2.0, ax1_pos[1], caxwidth, cayheight))

cb = colorbar.ColorbarBase(cax, norm=colors.Normalize(vmax=rcfreqsplt[-1], vmin=rcfreqsplt[0]), cmap=plt.cm.RdYlBu,
                           orientation='vertical', boundaries=bounds, spacing='uniform', ticks=ticks,
                           format='%4.1f  ', alpha=1)
plt.text(0.5, 1.05, 'MW', ha='center', va='bottom', transform=cax.transAxes, color='k', fontweight='normal')
plt.text(0.5, 1.01, '[GHz]', ha='center', va='bottom', transform=cax.transAxes, color='k', fontweight='normal')

cax.xaxis.set_visible(False)
cax.tick_params(axis="y", direction="in", pad=-21.5, length=0, colors='k', labelsize=8)

fig.tight_layout()