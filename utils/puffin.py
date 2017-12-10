import astropy.units as u
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import sunpy.coordinates
import sunpy.map
import sunpy.cm
import sunpy.wcs as wcs
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from bokeh.models.mappers import LogColorMapper, LinearColorMapper
from skimage.io import imread
from suncasa.utils import DButil

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

# define the colormaps
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]


def getAIApalette(wavelength):
    clmap = cm.get_cmap("sdoaia" + wavelength)  # choose any matplotlib colormap here
    palette = [colors.rgb2hex(m).encode("ascii").upper() for m in clmap(np.arange(clmap.N))]
    return palette


class PuffinMap:
    """
    Produce A Generic Column data source from map for Bokeh plot

    Parameters
    ----------
    map : Sunpy map
    """
    __slots__ = ['smap', 'plot_height', 'plot_width', 'x', 'y', 'dw', 'dh', 'x_range', 'y_range']

    def __init__(self, data=None, header=None, smap=None, plot_height=None, plot_width=None, *args,
                 **kwargs):
        if not smap:
            smap = sunpy.map.Map((data, header))
        self.smap = smap
        if not plot_height:
            plot_height = 400.0
        if not plot_width:
            plot_width = 400.0
        self.plot_height = plot_height
        self.plot_width = plot_width

        """
        define the plot region based the plot_width and plot_height
        ------------------
        x (DataSpecProperty) - The x-coordinates to locate the image anchors.
        y (DataSpecProperty) - The y-coordinates to locate the image anchors.
        dw (UnitsSpecProperty) - The widths of the plot regions that the images will occupy.
        dh (UnitsSpecProperty) - The height of the plot region that the image will occupy.
        """
        try:
            dx, dy = self.smap.scale.x.value, self.smap.scale.y.value
        except:
            dx, dy = self.smap.scale[0].value, self.smap.scale[1].value

        self.dx, self.dy = dx, dy
        x0, x1 = self.smap.xrange.value[0] - dx / 2.0, self.smap.xrange.value[1] + dx / 2.0
        y0, y1 = self.smap.yrange.value[0] - dy / 2.0, self.smap.yrange.value[1] + dy / 2.0
        self.x, self.y = [x0], [y0]
        self.dw, self.dh = [x1 - x0], [y1 - y0]
        self.x_range = [x0, x1]
        # self.y_range = [y0, y0 + (y1 - y0) * float(self.plot_height) / float(self.plot_width)]
        self.y_range = [y0, y1]

    def meshgrid(self, rescale=1.0, *args, **kwargs):
        XX, YY = np.meshgrid(np.arange(self.smap.data.shape[1] * rescale), np.arange(self.smap.data.shape[0] * rescale))
        mesh =  self.smap.pixel_to_world(XX / rescale * u.pix, YY / rescale * u.pix)
        x,y = mesh.Tx,mesh.Tx
        return x, y

    def meshgridpix(self, rescale=1.0, *args, **kwargs):
        XX, YY = np.meshgrid(np.arange(self.smap.data.shape[1] * rescale) / rescale,
                             np.arange(self.smap.data.shape[0] * rescale) / rescale)
        x, y = XX * u.pix, YY * u.pix
        return x, y

    def ImageSource(self, *args, **kwargs):
        """maps the Sunpy map to Bokeh DataSource
        """
        # XX, YY = (np.arange(self.smap.data.shape[0]), np.zeros(self.smap.data.shape[0]))
        # x = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)[0]
        # XX, YY = (np.zeros(self.smap.data.shape[1]), np.arange(self.smap.data.shape[1]))
        # y = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)[1]
        data = self.smap.data.copy()
        data[~np.isnan(data)] = data[~np.isnan(data)] / self.smap.exposure_time.value
        # return {'data': [data], 'xx': [x], 'yy': [y]}
        return {'data': [data]}

    def DrawGridSource(self, grid_spacing=15 * u.deg, *args, **kwargs):
        """maps the Longitude and Latitude grids to Bokeh DataSource
        """
        # XX, YY = np.meshgrid(np.arange(self.smap.data.shape[0]), np.arange(self.smap.data.shape[1]))
        # x, y = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)
        dsun = self.smap.dsun

        b0 = self.smap.heliographic_latitude.to(u.deg).value
        l0 = self.smap.heliographic_longitude.to(u.deg).value
        units = self.smap.spatial_units

        xs = []
        ys = []

        hg_longitude_deg = np.linspace(-180, 180, num=361) + l0
        hg_latitude_deg = np.arange(-90, 90, grid_spacing.to(u.deg).value)

        # draw the latitude lines
        for lat in hg_latitude_deg:
            x, y = wcs.convert_hg_hpc(hg_longitude_deg, lat * np.ones(361), b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units.axis1, occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            xs.append(x.tolist())
            ys.append(y.tolist())

        hg_longitude_deg = np.arange(-180, 180, grid_spacing.to(u.deg).value) + l0
        hg_latitude_deg = np.linspace(-90, 90, num=181)

        # draw the longitude lines
        for lon in hg_longitude_deg:
            x, y = wcs.convert_hg_hpc(lon * np.ones(181), hg_latitude_deg, b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units.axis1, occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            xs.append(x.tolist())
            ys.append(y.tolist())

        return ColumnDataSource(data={'xs': xs, 'ys': ys})

    def DrawLimbSource(self, *args, **kwargs):
        """maps the solar limb to Bokeh DataSource
        """
        radius = self.smap.rsun_obs.value

        phi = np.linspace(-180, 180, num=361) * np.pi / 180.
        x = np.cos(phi) * radius
        y = np.sin(phi) * radius

        return ColumnDataSource(data={'x': x, 'y': y})

    def PlotMap(self, DrawLimb=True, DrawGrid=True, grid_spacing=15 * u.deg, ignore_coord=False, title=None, tools=None,
                x_range=None, y_range=None,
                palette=None, imagetype='image', *args,
                **kwargs):
        """Plot the map using the bokeh.plotting interface
        """

        source_image = self.ImageSource()
        if not title:
            title = self.smap.name

        if not tools:
            tools = 'pan,wheel_zoom,save,reset'
        if not x_range:
            x_range = self.x_range
        if not y_range:
            y_range = self.y_range

        if not palette:
            if self.smap.observatory == 'SDO' and self.smap.instrument[0:3] == 'AIA':
                wavelngth = '{:.0f}'.format(self.smap.wavelength.value)
                # clmap = cm.get_cmap("sdoaia" + wavelngth)  # choose any matplotlib colormap here
                # palette = [colors.rgb2hex(m).encode("ascii").upper() for m in clmap(np.arange(clmap.N))]
                palette = getAIApalette(wavelngth)
                clrange = DButil.sdo_aia_scale_dict(wavelength=wavelngth, imagetype=imagetype)
                print clrange
                if clrange['log']:
                    colormapper = LogColorMapper(palette=palette, low=clrange['low'], high=clrange['high'])
                else:
                    colormapper = LinearColorMapper(palette=palette, low=clrange['low'], high=clrange['high'])
            else:
                palette = bokehpalette_jet
                colormapper = LinearColorMapper(palette=palette)
        else:
            colormapper = LinearColorMapper(palette=palette)

        if ignore_coord:
            x_range = [0 - 0.5, self.smap.data.shape[0] + 0.5]
            y_range = [0 - 0.5, self.smap.data.shape[1] + 0.5]
            p_xaxis_visible = False
            p_yaxis_visible = False
            p_xaxis_axislabel = ''
            p_yaxis_axislabel = ''
            x0 = 0 - 0.5
            y0 = 0 - 0.5
            dw = self.smap.data.shape[0]
            dh = self.smap.data.shape[1]
        else:
            p_xaxis_visible = True
            p_yaxis_visible = True
            p_xaxis_axislabel = 'X-position [arcsec]'
            p_yaxis_axislabel = 'Y-position [arcsec]'
            x0 = self.x
            y0 = self.y
            dw = self.dw
            dh = self.dh

        p_image = figure(tools=tools, x_range=x_range,
                         y_range=y_range, title=title, plot_height=self.plot_height,
                         plot_width=self.plot_width, *args, **kwargs)
        p_image.xaxis.axis_label = p_xaxis_axislabel
        p_image.yaxis.axis_label = p_yaxis_axislabel
        p_image.xaxis.visible = p_xaxis_visible
        p_image.yaxis.visible = p_yaxis_visible
        r_img = p_image.image(image=source_image['data'], x=x0, y=y0, dw=dw, dh=dh, color_mapper=colormapper)
        if DrawLimb:
            p_image.line(x='x', y='y', line_color='white', line_dash='solid', source=self.DrawLimbSource())
            if DrawGrid:
                p_image.multi_line(xs='xs', ys='ys', line_color='white', line_dash='dotted',
                                   source=self.DrawGridSource(grid_spacing=grid_spacing))
        return p_image, r_img


class PuffinImage:
    """
    Produce a A Generic Column data source from image for Bokeh plot

    Parameters
    ----------
    image : image data
    imfile: image file
    """
    __slots__ = ['plot_height', 'plot_width', 'img', 'x', 'y', 'dw', 'dh', 'x_range', 'y_range']

    def __init__(self, image=None, imfile=None, plot_height=None, plot_width=None, *args, **kwargs):
        if image is None and imfile is None:
            raise IOError("Missing input image or image file")
        else:
            if imfile is None:
                im = image
            else:
                im = imread(imfile)
            if not plot_height:
                plot_height = im.shape[0]
            if not plot_width:
                plot_width = im.shape[1]

            self.plot_height = plot_height
            self.plot_width = plot_width
            img = np.empty((plot_width, plot_height), dtype=np.uint32)
            view = img.view(dtype=np.uint8).reshape((plot_width, plot_height, 4))
            view[:, :, 0:3] = im[:, :, 0:3]
            self.img = np.flipud(img)

        """
        defind the plot region based the plot_width and plot_height
        ------------------
        x (DataSpecProperty) - The x-coordinates to locate the image anchors.
        y (DataSpecProperty) - The y-coordinates to locate the image anchors.
        dw (UnitsSpecProperty) - The widths of the plot regions that the images will occupy.
        dh (UnitsSpecProperty) - The height of the plot region that the image will occupy.
        """

        self.x, self.y = 0, 0
        self.dw, self.dh = plot_width, plot_height
        self.x_range = [0, plot_width]
        self.y_range = [0, plot_height]

    def ImageSource(self, *args, **kwargs):
        """maps the jpg to Bokeh DataSource
        """
        x = np.arange(self.img.shape[0])
        y = np.arange(self.img.shape[1])
        return ColumnDataSource(data={'data': [self.img], 'xx': [x], 'yy': [y]})

    def PlotImage(self, title=None, x_range=None, y_range=None, *args, **kwargs):
        """Plot the image using the bokeh.plotting interface
        """

        source_image = self.ImageSource()
        if not title:
            title = ''

        if not x_range:
            x_range = self.x_range
        if not y_range:
            y_range = self.y_range

        # plot the global vla image
        p_image = figure(tools='', x_range=x_range, y_range=y_range, title=title, toolbar_location=None,
                         plot_height=self.plot_height, plot_width=self.plot_width)
        p_image.xaxis.visible = False
        p_image.yaxis.visible = False
        r_img = p_image.image_rgba(image="data", x=self.x, y=self.y, dw=self.dw, dh=self.dh, source=source_image)

        return p_image, r_img
