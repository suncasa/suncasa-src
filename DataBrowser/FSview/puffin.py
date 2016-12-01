import astropy.units as u
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import sunpy.coordinates
import sunpy.map
import sunpy.wcs as wcs
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from skimage.io import imread

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

# define the colormaps
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]


class PuffinMap:
    """
    Produce A Generic Column data source from map for Bokeh plot

    Parameters
    ----------
    map : Sunpy map
    """
    __slots__ = ['smap', 'plot_height', 'plot_width', 'x', 'y', 'dw', 'dh', 'x_range', 'y_range']

    def __init__(self, data=None, header=None, smap=None, plot_height=None, plot_width=None, webgl=False, *args,
                 **kwargs):
        if not smap:
            smap = sunpy.map.Map((data, header))
        self.smap = smap
        if not plot_height:
            plot_height = 400.0
        if not plot_width:
            plot_width = 400.0
        self.webgl = webgl
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
        x0, x1 = self.smap.xrange.value[0], self.smap.xrange.value[1]
        y0, y1 = self.smap.yrange.value[0], self.smap.yrange.value[1]
        self.x, self.y = [x0], [y0]
        self.dw, self.dh = [x1 - x0], [y1 - y0]
        self.x_range = [x0, x1]
        self.y_range = [y0, y0 + (y1 - y0) * float(self.plot_height) / float(self.plot_width)]

    def meshgrid(self, *args, **kwargs):
        XX, YY = np.meshgrid(np.arange(self.smap.data.shape[0]), np.arange(self.smap.data.shape[1]))
        x, y = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)
        return x, y

    def meshgridpix(self, *args, **kwargs):
        XX, YY = np.meshgrid(np.arange(self.smap.data.shape[0]), np.arange(self.smap.data.shape[1]))
        x, y = XX * u.pix, YY * u.pix
        return x, y

    def ImageSource(self, *args, **kwargs):
        """maps the Sunpy map to Bokeh DataSource
        """
        XX, YY = (np.arange(self.smap.data.shape[0]), np.zeros(self.smap.data.shape[0]))
        x = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)[0]
        XX, YY = (np.zeros(self.smap.data.shape[1]), np.arange(self.smap.data.shape[1]))
        y = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)[1]
        return ColumnDataSource(data={'data': [self.smap.data], 'xx': [x], 'yy': [y]})

    def DrawGridSource(self, grid_spacing=15 * u.deg, *args, **kwargs):
        """maps the Longitude and Latitude grids to Bokeh DataSource
        """
        XX, YY = np.meshgrid(np.arange(self.smap.data.shape[0]), np.arange(self.smap.data.shape[1]))
        x, y = self.smap.pixel_to_data(XX * u.pix, YY * u.pix)
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
                                      angle_units=units.x, occultation=True)
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
                                      angle_units=units[0], occultation=True)
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

    def PlotMap(self, DrawLimb=True, DrawGrid=True, grid_spacing=15 * u.deg, title=None, x_range=None, y_range=None,
                palette=None, *args,
                **kwargs):
        """Plot the map using the bokeh.plotting interface
        """

        source_image = self.ImageSource()
        if not title:
            title = self.smap.name

        if not x_range:
            x_range = self.x_range
        if not y_range:
            y_range = self.y_range
        if not palette:
            palette = bokehpalette_jet
        # plot the global vla image
        p_image = figure(tools='pan,wheel_zoom,save,reset', webgl=self.webgl, x_range=x_range, y_range=y_range,
                         title=title,
                         plot_height=self.plot_height,
                         plot_width=self.plot_width, *args, **kwargs)
        p_image.xaxis.axis_label = 'X-position [arcsec]'
        p_image.yaxis.axis_label = 'Y-position [arcsec]'
        r_img = p_image.image(image="data", x=self.x, y=self.y, dw=self.dw, dh=self.dh, source=source_image,
                              palette=palette)
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
