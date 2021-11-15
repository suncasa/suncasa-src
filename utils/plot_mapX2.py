# import os
# import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.colorbar as colorbar
# import sunpy.cm.cm as cm  ## to bootstrap sdoaia color map
# import matplotlib.cm as cm
# import matplotlib.colors as colors
import astropy.units as u
# from astropy.io import fits
# import matplotlib.dates as mdates
# from astropy.time import Time
from sunpy import map as smap
# import matplotlib.gridspec as gridspec
# import numpy.ma as ma
# import matplotlib.patches as patches
# from suncasa.utils import stackplot as stp
# from IPython import embed
# from astropy.coordinates import SkyCoord
import numpy as np
from suncasa.utils import DButil
import warnings
import matplotlib.patches as patches
import numpy.ma as ma


class Sunmap():

    def __init__(self, sunmap, aia=False):
        if aia:
            try:
                sunmap = DButil.normalize_aiamap(sunmap)
            except:
                pass
            data = sunmap.data
            data[data < 1.0] = 1.0
            self.sunmap = smap.Map(data, sunmap.meta)
        else:
            self.sunmap = sunmap

    def map2wcsgrids(self, sunpymap=None, cell=False):
        '''

        :param sunpymap:
        :param cell: if True, return the coordinates of the pixel centers. if False, return the coordinates of the pixel boundaries
        :return:
        '''
        # embed()
        import astropy.units as u
        from suncasa.utils import stputils as stpu
        if sunpymap is None:
            sunpymap = self.sunmap
        ny, nx = sunpymap.data.shape
        x0, x1, y0, y1 = stpu.get_map_corner_coord(sunpymap).value

        dx = sunpymap.scale.axis1.to(u.arcsec / u.pix).value
        dy = sunpymap.scale.axis2.to(u.arcsec / u.pix).value

        if cell:
            mapx, mapy = np.linspace(x0, x1, nx) - dx / 2.0, np.linspace(y0, y1, ny) - dy / 2.0
            mapx = np.tile(mapx, ny).reshape(ny, nx)
            mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
        else:
            nx += 1
            ny += 1
            mapx, mapy = np.linspace(x0 - dx, x1 + dx, nx), np.linspace(y0 - dy, y1 + dy, ny)
            mapx = np.tile(mapx, ny).reshape(ny, nx)
            mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
        return mapx, mapy

    def get_map_extent(self, sunpymap=None, rot=0, rangereverse=False):
        from suncasa.utils import stputils as stpu
        if sunpymap is None:
            sunpymap = self.sunmap
        rot = rot % 360
        x0, x1, y0, y1 = stpu.get_map_corner_coord(sunpymap).value
        if rot == 90:
            extent = np.array([y1, y0, x0, x1])
            extent = extent - np.array([sunpymap.scale.axis2.value] * 2 + [sunpymap.scale.axis1.value] * 2) / 2.0
        elif rot == 180:
            extent = np.array([x1, x0, y1, y0])
            extent = extent - np.array([sunpymap.scale.axis1.value] * 2 + [sunpymap.scale.axis2.value] * 2) / 2.0
        elif rot == 270:
            extent = np.array([y0, y1, x1, x0])
            extent = extent - np.array([sunpymap.scale.axis1.value] * 2 + [sunpymap.scale.axis2.value] * 2) / 2.0
        else:
            extent = np.array([x0, x1, y0, y1])
            extent = extent - np.array([sunpymap.scale.axis1.value] * 2 + [sunpymap.scale.axis2.value] * 2) / 2.0
        if rangereverse:
            x0, x1, y0, y1 = extent
            extent = -np.array([x1, x0, y1, y0])
        return extent

    def imshow(self, axes=None, rot=0, rangereverse=False, maskon=False, image_enhance=False, **kwargs):
        '''
        :param sunpymap:
        :param axes:
        :param rot: rotation angle in degrees counter-clockwise. Must be an integer multiple of 90.
        :param kwargs:
        :return:
        '''
        sunpymap = self.sunmap
        if axes is None:
            axes = plt.subplot()
        rot = rot % 360
        if rot == 0:
            imdata = sunpymap.data
        elif rot == 90:
            imdata = sunpymap.data.transpose()[:, ::-1]
        elif rot == 180:
            imdata = sunpymap.data[::-1, ::-1]
        elif rot == 270:
            imdata = sunpymap.data.transpose()[::-1, :]
        else:
            warnings.warn('rot must be an integer multiple of 90. rot not implemented!')
            imdata = sunpymap.data
            rot = 0
        extent = self.get_map_extent(rot=rot, rangereverse=rangereverse)

        if maskon:
            if isinstance(maskon, bool):
                imdataplt = ma.masked_invalid(imdata)
                immask = imdataplt.mask
            elif isinstance(maskon, dict):
                if 'masked_equal' in maskon.keys():
                    imdataplt = ma.masked_equal(imdata, maskon['masked_equal'])
                elif 'masked_greater' in maskon.keys():
                    imdataplt = ma.masked_greater(imdata, maskon['masked_greater'])
                elif 'masked_less' in maskon.keys():
                    imdataplt = ma.masked_less(imdata, maskon['masked_less'])
                elif 'masked_greater_equal' in maskon.keys():
                    imdataplt = ma.masked_greater_equal(imdata, maskon['masked_greater_equal'])
                elif 'masked_less_equal' in maskon.keys():
                    imdataplt = ma.masked_less_equal(imdata, maskon['masked_less_equal'])
                elif 'masked_outside' in maskon.keys():
                    v1, v2 = maskon['masked_outside']
                    imdataplt = ma.masked_outside(imdata, v1, v2)
                elif 'masked_inside' in maskon.keys():
                    v1, v2 = maskon['masked_inside']
                    imdataplt = ma.masked_inside(imdata, v1, v2)
                elif 'masked_invalid' in maskon.keys():
                    imdataplt = ma.masked_invalid(imdata)
                else:
                    raise ValueError('maskon key wrong.')
                immask = imdataplt.mask
            else:
                raise TypeError('maskon must be bool or dict type.')
        else:
            imdataplt = imdata.copy()

        if image_enhance:
            dmax = np.nanmax(imdataplt)
            dmin = np.nanmin(imdataplt)
            from skimage.exposure import equalize_adapthist
            if isinstance(image_enhance, dict):
                imdataplt = equalize_adapthist(imdataplt, **image_enhance) * (dmax - dmin) + dmin
            else:
                imdataplt = equalize_adapthist(imdataplt) * (dmax - dmin) + dmin

        if maskon:
            imdataplt = ma.masked_array(imdataplt, immask)

        if isinstance(axes, list):
            ims = []
            for ax in axes:
                im = ax.imshow(imdataplt, extent=extent, origin='lower', **kwargs)
                ims.append(im)
                if rot == 0:
                    ax.set_xlabel('Solar X [arcsec]')
                    ax.set_ylabel('Solar Y [arcsec]')
                elif rot == 90:
                    ax.set_xlabel('-Solar Y [arcsec]')
                    ax.set_ylabel('Solar X [arcsec]')
                elif rot == 180:
                    ax.set_xlabel('-Solar X [arcsec]')
                    ax.set_ylabel('-Solar Y [arcsec]')
                elif rot == 270:
                    ax.set_xlabel('Solar Y [arcsec]')
                    ax.set_ylabel('-Solar X [arcsec]')
            return ims
        else:
            im = axes.imshow(imdataplt, extent=extent, origin='lower', **kwargs)

            if rot == 0:
                axes.set_xlabel('Solar X [arcsec]')
                axes.set_ylabel('Solar Y [arcsec]')
            elif rot == 90:
                axes.set_xlabel('-Solar Y [arcsec]')
                axes.set_ylabel('Solar X [arcsec]')
            elif rot == 180:
                axes.set_xlabel('-Solar X [arcsec]')
                axes.set_ylabel('-Solar Y [arcsec]')
            elif rot == 270:
                axes.set_xlabel('Solar Y [arcsec]')
                axes.set_ylabel('-Solar X [arcsec]')
            return im

    def contour(self, axes=None, rot=0, mapx=None, mapy=None, rangereverse=False, **kwargs):
        sunpymap = self.sunmap
        if axes is None:
            axes = plt.subplot()
        rot = rot % 360
        if (mapx is None) or (mapy is None):
            if rot == 0:
                mapx, mapy = self.map2wcsgrids(cell=True)
            elif rot == 90:
                mapy, mapx = self.map2wcsgrids(cell=True)
            elif rot == 180:
                mapx, mapy = self.map2wcsgrids(cell=True)
            elif rot == 270:
                mapy, mapx = self.map2wcsgrids(cell=True)

        if isinstance(axes, list):
            ims = []
            for ax in axes:
                im = ax.contour(mapx, mapy, sunpymap.data, **kwargs)
                ims.append(im)
                extent = self.get_map_extent(rot=rot, rangereverse=rangereverse)
                ax.set_xlim(extent[:2])
                ax.set_ylim(extent[2:])
            return ims
        else:
            im = axes.contour(mapx, mapy, sunpymap.data, **kwargs)
            extent = self.get_map_extent(rot=rot, rangereverse=rangereverse)
            axes.set_xlim(extent[:2])
            axes.set_ylim(extent[2:])
            return im

    def contourf(self, axes=None, rot=0, mapx=None, mapy=None, rangereverse=False, **kwargs):
        sunpymap = self.sunmap
        if axes is None:
            axes = plt.subplot()
        rot = rot % 360
        if (mapx is None) or (mapy is None):
            if rot == 0:
                mapx, mapy = self.map2wcsgrids(cell=True)
            elif rot == 90:
                mapy, mapx = self.map2wcsgrids(cell=True)
            elif rot == 180:
                mapx, mapy = self.map2wcsgrids(cell=True)
            elif rot == 270:
                mapy, mapx = self.map2wcsgrids(cell=True)

        if isinstance(axes, list):
            ims = []
            for ax in axes:
                im = ax.contourf(mapx, mapy, sunpymap.data, **kwargs)
                ims.append(im)
                extent = self.get_map_extent(rot=rot, rangereverse=rangereverse)
                ax.set_xlim(extent[:2])
                ax.set_ylim(extent[2:])
            return ims
        else:
            im = axes.contourf(mapx, mapy, sunpymap.data, **kwargs)
            extent = self.get_map_extent(rot=rot, rangereverse=rangereverse)
            axes.set_xlim(extent[:2])
            axes.set_ylim(extent[2:])
            return im

    def draw_limb(self, axes=None, rangereverse=False, **kwargs):
        if 'c' not in kwargs and 'color' not in kwargs:
            kwargs['c'] = 'w'
        if 'ls' not in kwargs and 'linestyle' not in kwargs:
            kwargs['ls'] = 'solid'
        sunpymap = self.sunmap
        if axes is None:
            axes = plt.gca()

        rsun = sunpymap.rsun_obs
        phi = np.linspace(-180, 180, num=181) * u.deg
        x = np.cos(phi) * rsun
        y = np.sin(phi) * rsun
        if isinstance(axes, list):
            ims = []
            for ax in axes:
                ax.set_autoscale_on(False)
                im = ax.plot(x, y, **kwargs)
                ims.append(im)
            return ims
        else:
            axes.set_autoscale_on(False)
            im = axes.plot(x, y, **kwargs)
            return im

    def draw_grid(self, axes=None, rot=0, grid_spacing=None, **kwargs):
        sunpymap = self.sunmap
        if grid_spacing is None:
            grid_spacing = 15 * u.deg

        def hgs2hcc(rsun, lon, lat, B0, L0):
            lon_L0 = lon - L0
            x = rsun * np.cos(lat) * np.sin(lon)
            y = rsun * (np.sin(lat) * np.cos(B0) - np.cos(lat) * np.cos(lon_L0) * np.sin(B0))
            z = rsun * (np.sin(lat) * np.sin(B0) + np.cos(lat) * np.cos(lon_L0) * np.cos(B0))
            return x, y, z

        def hcc2hpc(x, y, z, dsun):
            d = np.sqrt(x ** 2 + y ** 2 + (dsun - z) ** 2)
            Tx = np.arctan2(x, dsun - z)
            Ty = np.arcsin(y / d)
            return Tx, Ty

        if 'c' not in kwargs and 'color' not in kwargs:
            kwargs['c'] = 'w'
        if 'ls' not in kwargs and 'linestyle' not in kwargs:
            kwargs['ls'] = 'dotted'
        dsun = sunpymap.dsun
        rsun = sunpymap.rsun_meters
        if axes is None:
            axes = plt.gca()
        im = []
        b0 = sunpymap.heliographic_latitude.to(u.deg)
        l0 = sunpymap.heliographic_longitude.to(u.deg)
        hg_longitude_deg = np.linspace(-90, 90, num=91) * u.deg
        hg_latitude_deg = np.arange(0, 90, grid_spacing.to(u.deg).value)
        hg_latitude_deg = np.hstack([-hg_latitude_deg[1:][::-1], hg_latitude_deg]) * u.deg
        for lat in hg_latitude_deg:
            c = hgs2hcc(rsun, hg_longitude_deg, lat * np.ones(91), b0, l0)
            coords = hcc2hpc(c[0], c[1], c[2], dsun)
            if rot in [90, 270]:
                coords_ = [coords[1], coords[0]]
            else:
                coords_ = coords
            if isinstance(axes, list):
                for ax in axes:
                    im += ax.plot(coords_[0].to(u.arcsec), coords_[1].to(u.arcsec), **kwargs)
            else:
                im += axes.plot(coords_[0].to(u.arcsec), coords_[1].to(u.arcsec), **kwargs)

        hg_longitude_deg = np.arange(0, 90, grid_spacing.to(u.deg).value)
        hg_longitude_deg = np.hstack([-hg_longitude_deg[1:][::-1], hg_longitude_deg]) * u.deg
        hg_latitude_deg = np.linspace(-90, 90, num=91) * u.deg

        for lon in hg_longitude_deg:
            c = hgs2hcc(rsun, lon * np.ones(91), hg_latitude_deg, b0, l0)
            coords = hcc2hpc(c[0], c[1], c[2], dsun)
            if rot in [90, 270]:
                coords_ = [coords[1], coords[0]]
            else:
                coords_ = coords
            if isinstance(axes, list):
                for ax in axes:
                    im += ax.plot(coords_[0].to(u.arcsec), coords_[1].to(u.arcsec), **kwargs)
            else:
                im += axes.plot(coords_[0].to(u.arcsec), coords_[1].to(u.arcsec), **kwargs)
        return im

    def draw_rectangle(self, bottom_left, width, height, axes=None, **kwargs):
        if 'ec' not in kwargs and 'edgecolor' not in kwargs:
            kwargs['ec'] = 'w'
        if 'ls' not in kwargs and 'linestyle' not in kwargs:
            kwargs['ls'] = 'solid'
        if 'fill' not in kwargs:
            kwargs['fill'] = False
        if axes is None:
            axes = plt.gca()

        if isinstance(axes, list):
            ims = []
            for ax in axes:
                ax.set_autoscale_on(False)
                im = ax.add_patch(patches.Rectangle(bottom_left, width, height, **kwargs))
                ims.append(im)
            return ims
        else:
            axes.set_autoscale_on(False)
            im = axes.add_patch(patches.Rectangle(bottom_left, width, height, **kwargs))
            return im

    def imshow_RGB(self, maps, axes=None, returndataonly=False, rangereverse=False):
        from scipy import ndimage
        from astropy.coordinates import SkyCoord
        mapR = maps[0]
        znewR = mapR.data
        aiamapx, aiamapy = self.map2wcsgrids(sunpymap=mapR, cell=False)
        mapG = maps[1]
        XX, YY = mapG.data_to_pixel(SkyCoord(aiamapx * u.arcsec, aiamapy * u.arcsec, frame=mapG.coordinate_frame))
        znewG = ndimage.map_coordinates(mapG.data, [YY, XX], order=1)
        mapB = maps[2]
        XX, YY = mapB.data_to_pixel(SkyCoord(aiamapx * u.arcsec, aiamapy * u.arcsec, frame=mapB.coordinate_frame))
        znewB = ndimage.map_coordinates(mapB.data, [YY, XX], order=1)

        znewR = np.sqrt(znewR)
        znewG = np.sqrt(znewG)
        znewB = np.sqrt(znewB)

        vmax, vmin = np.sqrt(5000), np.sqrt(10)
        # clrange=DButil.sdo_aia_scale_dict(304)
        znewR[znewR > vmax] = vmax
        znewR[znewR < vmin] = vmin
        # clrange=DButil.sdo_aia_scale_dict(94)
        vmax, vmin = np.sqrt(20000), np.sqrt(200)
        znewG[znewG > vmax] = vmax
        znewG[znewG < vmin] = vmin
        # clrange=DButil.sdo_aia_scale_dict(211)
        vmax, vmin = np.sqrt(5000), np.sqrt(100)
        znewB[znewB > vmax] = vmax
        znewB[znewB < vmin] = vmin
        znewR = (znewR - np.nanmin(znewR)) / (np.nanmax(znewR) - np.nanmin(znewR))
        znewG = (znewG - np.nanmin(znewG)) / (np.nanmax(znewG) - np.nanmin(znewG))
        znewB = (znewB - np.nanmin(znewB)) / (np.nanmax(znewB) - np.nanmin(znewB))
        # znew1 = np.sqrt(znew1)
        # znew2 = np.sqrt(znew2)
        # imshow(np.sqrt(np.stack([znew0, znew1, znew2], axis=-1)), extent=list(aiamap.xrange.value) + list(aiamap.yrange.value),origin='lower')
        if returndataonly:
            return np.stack([znewR, znewG, znewB], axis=-1)
        else:
            if axes:
                pass
            else:
                axes = plt.subplot()
            extent = self.get_map_extent(sunpymap=mapR, rangereverse=rangereverse)
            if isinstance(axes, list):
                ims = []
                for ax in axes:
                    im = ax.imshow(np.stack([znewR, znewG, znewB], axis=-1),
                                   extent=extent, origin='lower')
                    ims.append(im)
                return ims
            else:
                return axes.imshow(np.stack([znewR, znewG, znewB], axis=-1),
                                   extent=extent, origin='lower')
