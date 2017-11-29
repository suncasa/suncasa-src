import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.time import Time
import matplotlib.dates as mdates
from astropy.io import fits
import signalsmooth as ss
import numpy.ma as ma
import matplotlib.cm as cm
import sunpy.map
import astropy.units as u
import pickle
import glob
# from astropy import units as u
# import sunpy.map as smap
from scipy.interpolate import griddata
from scipy import signal
from suncasa.utils import DButil
from scipy.interpolate import splev, splrep
import scipy.ndimage
from IPython import embed
import time
from tqdm import *
from copy import deepcopy
from functools import partial
import multiprocessing as mp
import gc
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate


def resettable(f):
    import copy

    def __init_and_copy__(self, *args, **kwargs):
        f(self, *args)
        self.__original_dict__ = copy.deepcopy(self.__dict__)

        def reset(o=self):
            o.__dict__ = o.__original_dict__

        self.reset = reset

    return __init_and_copy__


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.lfilter(b, a, data)
    return y


def b_filter(data, lowcut, highcut, fs, ix):
    x = data[ix]
    y = butter_bandpass_filter(x, lowcut, highcut, fs, order=3)
    return {'idx': ix, 'y': y}


def FitSlit(xx, yy, cutwidth, cutang, cutlength, s=None, method='Polyfit', ascending=True):
    if len(xx) <= 3 or method == 'Polyfit':
        '''polynomial fit'''
        out = DButil.polyfit(xx, yy, cutlength, len(xx) - 1 if len(xx) <= 3 else 2)
        xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    else:
        if method == 'Param_Spline':
            '''parametic spline fit'''
            out = DButil.paramspline(xx, yy, cutlength, s=s)
            xs, ys, posangs = out['xs'], out['ys'], out['posangs']
        else:
            '''spline fit'''
            out = DButil.spline(xx, yy, cutlength, s=s)
            xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    if not ascending and (method != 'Param_Spline' or len(xx) <= 3):
        xs, ys = xs[::-1], ys[::-1]
        posangs = posangs[::-1]
    dist = DButil.findDist(xs, ys)
    dists = np.cumsum(dist)
    posangs2 = posangs + np.pi / 2
    cutwidths = dists * np.tan(cutang) + cutwidth
    xs0 = xs - cutwidths / 2. * np.cos(posangs2)
    ys0 = ys - cutwidths / 2. * np.sin(posangs2)
    xs1 = xs + cutwidths / 2. * np.cos(posangs2)
    ys1 = ys + cutwidths / 2. * np.sin(posangs2)
    return {'xcen': xs, 'ycen': ys, 'xs0': xs0, 'ys0': ys0, 'xs1': xs1, 'ys1': ys1, 'cutwidth': cutwidths,
            'posangs': posangs, 'posangs2': posangs2, 'dist': dists}


def MakeSlit(pointDF):
    pointDFtmp = pointDF
    xx = pointDFtmp.loc[:, 'xx'].values
    yy = pointDFtmp.loc[:, 'yy'].values
    if len(pointDFtmp.index) <= 1:
        cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [], 'posangs': [],
                      'posangs2': [], 'dist': []}
    else:
        # if len(pointDFtmp.index) <= 3:
        cutslitplt = FitSlit(xx, yy, 10, 0.0, 200, method='Polyfit')
    return cutslitplt


def getimprofile(data, cutslit, xrange=None, yrange=None):
    num = len(cutslit['xcen'])
    if num > 1:
        intens = np.zeros(num)
        ndy, ndx = data.shape
        if xrange is not None and yrange is not None:
            xs0 = (cutslit['xs0'] - xrange[0]) / (xrange[1] - xrange[0]) * ndx
            xs1 = (cutslit['xs1'] - xrange[0]) / (xrange[1] - xrange[0]) * ndx
            ys0 = (cutslit['ys0'] - yrange[0]) / (yrange[1] - yrange[0]) * ndy
            ys1 = (cutslit['ys1'] - yrange[0]) / (yrange[1] - yrange[0]) * ndy
        else:
            xs0 = cutslit['xs0']
            xs1 = cutslit['xs1']
            ys0 = cutslit['ys0']
            ys1 = cutslit['ys1']
        for ll in range(num):
            inten = DButil.improfile(data, [xs0[ll], xs1[ll]], [ys0[ll], ys1[ll]], interp='nearest')
            intens[ll] = np.mean(inten)
        intensdist = {'x': cutslit['dist'], 'y': intens}
        return intensdist


def plot_map(smap, dspec=None, diff=False, SymLogNorm=False, linthresh=0.5, returnImAx=False, uni_cm=False, *args,
             **kwargs):
    import sunpy.cm.cm as cm  ## to bootstrap sdoaia color map
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    from suncasa.utils import DButil
    try:
        clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
    except:
        clrange = {'high': None, 'log': False, 'low': None}
    plt.clf()
    if dspec:
        ax = plt.subplot(121)
    else:
        ax = plt.subplot()
    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = clrange['low']
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = clrange['high']
    if diff:
        if SymLogNorm:
            norm = colors.SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)
        else:
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    try:
        imshow_args = {'cmap': cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth'])), 'norm': norm,
                       'interpolation': 'nearest', 'origin': 'lower'}
    except:
        imshow_args = {'cmap': 'gray', 'norm': norm, 'interpolation': 'nearest', 'origin': 'lower'}
    try:
        if smap.coordinate_system.x == 'HG':
            xlabel = 'Longitude [{lon}]'.format(lon=smap.spatial_units.x)
        else:
            xlabel = 'X-position [{xpos}]'.format(xpos=smap.spatial_units.x)
        if smap.coordinate_system.y == 'HG':
            ylabel = 'Latitude [{lat}]'.format(lat=smap.spatial_units.y)
        else:
            ylabel = 'Y-position [{ypos}]'.format(ypos=smap.spatial_units.y)
    except:
        if smap.coordinate_system.axis1 == 'HG':
            xlabel = 'Longitude [{lon}]'.format(lon=smap.spatial_units.axis1)
        else:
            xlabel = 'X-position [{xpos}]'.format(xpos=smap.spatial_units.axis1)
        if smap.coordinate_system.axis2 == 'HG':
            ylabel = 'Latitude [{lat}]'.format(lat=smap.spatial_units.axis2)
        else:
            ylabel = 'Y-position [{ypos}]'.format(ypos=smap.spatial_units.axis2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    imshow_args.update({'extent': list(smap.xrange.value) + list(smap.yrange.value)})
    if smap.detector == 'HMI':
        im1 = ax.imshow(np.rot90(smap.data, 2), **imshow_args)
    else:
        im1 = ax.imshow(smap.data, **imshow_args)
    plt.title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['date-obs']))
    plt.colorbar(im1, ax=ax, label='DN counts per second')
    ax.set_autoscale_on(False)
    if dspec:
        fig = plt.gcf()
        # if figsize:
        #     fig.set_size_inches(figsize)
        # else:
        #     fig.set_size_inches(13, 5)
        if uni_cm:
            dspec['args']['norm'] = norm
        ax2 = plt.subplot(122)
        im2 = plt.pcolormesh(dspec['x'], dspec['y'], dspec['dspec'], **dspec['args'])
        date_format = mdates.DateFormatter('%H:%M:%S')
        ax2.xaxis_date()
        ax2.xaxis.set_major_formatter(date_format)
        fig.autofmt_xdate(rotation=30)
        ax2.yaxis.set_label_text(dspec['ytitle'])
        plt.colorbar(im2, ax=ax2, label=dspec['ctitle'])
        ax2.set_autoscale_on(False)
        if 'axvspan' in dspec.keys():
            vspan = ax2.axvspan(dspec['axvspan'][0], dspec['axvspan'][1], alpha=0.5, color='white')
        if 'xs' in dspec.keys() and 'ys' in dspec.keys():
            ax2.plot(dspec['xs'], dspec['ys'], '--', lw=2.0, alpha=0.7, c='black')
        if 'xlim' in dspec.keys():
            ax2.set_xlim(dspec['xlim'])
        if 'ylim' in dspec.keys():
            ax2.set_ylim(dspec['ylim'])
        if returnImAx:
            return ax, im1, ax2, im2, vspan
        else:
            return ax, ax2
    else:
        if returnImAx:
            return ax, im1
        else:
            return ax
            # ax.autoscale(True, 'both', True)
            # ax.autoscale_view(True, True, True)
            # ax.relim(visible_only=True)


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')
    y = y[window_len - 1:-(window_len - 1)]
    return y


def grid(x, y, z, resX=20, resY=40):
    "Convert 3 column data to matplotlib grid"
    xi = np.linspace(np.nanmin(x), np.nanmax(x), resX)
    yi = np.linspace(np.nanmin(y), np.nanmax(y), resY)
    Z = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z


def polyfit(x, y, length, deg):
    xs = np.linspace(x.min(), x.max(), length)
    z = np.polyfit(x=x, y=y, deg=deg)
    p = np.poly1d(z)
    ys = p(xs)
    rms = np.sqrt(np.sum((np.polyval(z, x) - y) ** 2) / len(x))
    return {'xs': xs, 'ys': ys, 'rms': rms}


class CutslitBuilder:
    def __init__(self, axes, cutwidth=5, cutang=0, cutlength=80, scale=1.0):
        self.axes = axes
        self.clickedpoints, = self.axes.plot([], [], 'o', color='white')
        self.slitline, = self.axes.plot([], [], color='white', ls='solid')
        self.slitline0, = self.axes.plot([], [], color='white', ls='dotted')
        self.slitline1, = self.axes.plot([], [], color='white', ls='dotted')
        self.cutlength = cutlength
        self.cutwidth = cutwidth
        self.cutang = cutang
        self.scale = scale
        self.xx = list(self.clickedpoints.get_xdata())
        self.yy = list(self.clickedpoints.get_ydata())
        self.cid = self.clickedpoints.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        tmode = '{}'.format(self.clickedpoints.figure.canvas.toolbar.mode)
        if tmode == '':
            if event.inaxes != self.axes:
                return
            if event.button == 1:
                self.xx.append(event.xdata)
                self.yy.append(event.ydata)
            elif event.button == 3:
                if len(self.xx) > 0:
                    self.xx.pop()
                    self.yy.pop()
            self.clickedpoints.set_data(self.xx, self.yy)
            self.clickedpoints.figure.canvas.draw()
            self.update()
        else:
            if event.inaxes != self.axes:
                return
            if event.button == 1 or event.button == 3:
                self.clickedpoints.figure.canvas.toolbar.set_message('Uncheck toolbar button {} first!'.format(tmode))

    def update(self):
        xx = np.array(self.xx, dtype=np.float64)
        yy = np.array(self.yy, dtype=np.float64)

        if len(self.xx) <= 1:
            cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                          'posangs': [], 'posangs2': [], 'dist': []}
        else:
            if len(self.xx) <= 3:
                cutslitplt = FitSlit(xx, yy, self.cutwidth * self.scale, self.cutang, self.cutlength, method='Polyfit')
            else:
                cutslitplt = FitSlit(xx, yy, self.cutwidth * self.scale, self.cutang, self.cutlength, s=len(xx),
                                     method='Param_Spline')
        self.cutslitplt = cutslitplt
        self.slitline.set_data(cutslitplt['xcen'], cutslitplt['ycen'])
        self.slitline0.set_data(cutslitplt['xs0'], cutslitplt['ys0'])
        self.slitline1.set_data(cutslitplt['xs1'], cutslitplt['ys1'])
        self.slitline.figure.canvas.draw()
        self.slitline0.figure.canvas.draw()
        self.slitline1.figure.canvas.draw()


class Stackplot:
    instrum_meta = {'SDO/AIA': {'scale': 0.6 * u.arcsec / u.pix}}
    suncasadb = os.getenv('SUNCASADB')
    if not suncasadb:
        print('Environmental variable for SUNCASA database path not defined')
        print('Use default path')
        suncasadb = './'
    fitsdir = suncasadb + '/aiaBrowserData/Download/'
    mapcube = None
    mapcube_diff = None
    mapcube_plot = None
    cutslitbd = None
    stackplt = None
    trange = None
    wavelength = None
    fitsfile = None
    fov = None
    binpix = None
    dt_data = None

    @resettable
    def __init__(self, infile=None):
        if infile:
            if isinstance(infile, sunpy.map.mapcube.MapCube):
                self.mapcube = infile
                self.mapcube_info()
            else:
                self.mapcube_fromfile(infile)

    def make_mapcube(self, trange, outfile=None, fov=None, wavelength='171', binpix=1, dt_data=1, derotate=False,
                     tosave=True):
        if isinstance(trange, list):
            if isinstance(trange[0], Time):
                trange = Time([trange[0], trange[-1]])
                fitsfile = DButil.readsdofile(datadir=self.fitsdir, wavelength=wavelength, jdtime=trange.jd)
            else:
                fitsfile = trange
        elif isinstance(trange, Time):
            fitsfile = DButil.readsdofile(datadir=self.fitsdir, wavelength=wavelength, jdtime=trange.jd)
        else:
            print(
                'Input trange format not recognized. trange can either be a file list or a timerange of astropy Time object')

        maplist = []
        print 'Loading fits files....'
        for ll in tqdm(fitsfile[::dt_data]):
            maptmp = sunpy.map.Map(ll)
            if fov:
                x0, x1, y0, y1 = fov
                submaptmp = maptmp.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                          u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
            else:
                submaptmp = maptmp
            submaptmp = submaptmp.resample(u.Quantity(submaptmp.dimensions) / binpix)
            if submaptmp.detector == 'HMI':
                pass
            else:
                submaptmp = DButil.normalize_aiamap(submaptmp)
            maplist.append(submaptmp)
        if derotate:
            mapcube = mapcube_solar_derotate(sunpy.map.Map(maplist, cube=True))
        else:
            mapcube = sunpy.map.Map(maplist, cube=True)
        trange = Time([mapcube[0].date, mapcube[-1].date])
        self.fitsfile = fitsfile
        self.dt_data = dt_data
        self.mapcube = mapcube
        self.mapcube_info()

        if tosave:
            if not outfile:
                outfile = 'mapcube_{0}_bin{3}_dtdata{4}_{1}_{2}'.format(mapcube[0].meta['wavelnth'],
                                                                        trange[0].isot[:-4].replace(':', ''),
                                                                        trange[1].isot[:-4].replace(':', ''), binpix,
                                                                        dt_data)
            self.mapcube_tofile(outfile)
        gc.collect()

    def mapcube_fromfile(self, infile):
        with open(infile, 'rb') as sf:
            print('Loading mapcube....')
            tmp = pickle.load(sf)
            if isinstance(tmp, dict):
                if not isinstance(tmp['mp'], sunpy.map.mapcube.MapCube):
                    print('Load failed. mapcube must be a instance of sunpy.map.mapcube.MapCube')
                    return
                self.mapcube = tmp['mp']
                self.dt_data = tmp['dt_data']
                self.fitsfile = tmp['fitsfile']
            else:
                if not isinstance(tmp, sunpy.map.mapcube.MapCube):
                    print('Load failed. mapcube must be a instance of sunpy.map.mapcube.MapCube')
                    return
                self.mapcube = tmp
            self.mapcube_info()

    def mapcube_tofile(self, outfile=None, mapcube=None):
        if not mapcube:
            mapcube = self.mapcube
        mp_info = self.mapcube_info(mapcube)
        if not outfile:
            outfile = 'mapcube_{0}_{1}_{2}'.format(mapcube[0].meta['wavelnth'],
                                                   self.trange[0].isot[:-4].replace(':', ''),
                                                   self.trange[1].isot[:-4].replace(':', ''))
        with open(outfile, 'wb') as sf:
            print('Saving mapcube to {}'.format(outfile))
            pickle.dump({'mp': mapcube, 'trange': mp_info['trange'], 'fov': mp_info['fov'], 'binpix': mp_info['binpix'],
                         'dt_data': self.dt_data, 'fitsfile': self.fitsfile}, sf)

    def mapcube_drot(self):
        self.mapcube = mapcube_solar_derotate(self.mapcube)
        return self.mapcube

    def mapcube_resample(self, binpix=1):
        print 'resampling mapcube.....'
        maplist = []
        for idx, ll in enumerate(tqdm(self.mapcube)):
            maplist.append(deepcopy(ll.resample(u.Quantity(ll.dimensions) / binpix)))
        self.mapcube = sunpy.map.Map(maplist, cube=True)
        self.binpix *= binpix

    def mapcube_mkdiff(self, mode=0, dt_frm=3, medfilt=None, bfilter=False, lowcut=0.1, highcut=50, outfile=None,
                       tosave=False):
        modes = {0: 'rdiff', 1: 'rratio', 2: 'bdiff', 3: 'bratio'}
        maplist = []
        datacube = self.mapcube.as_array()
        if medfilt:
            print 'median filtering map.....'
            for idx, ll in enumerate(tqdm(self.mapcube)):
                datacube[:, :, idx] = signal.medfilt(datacube[:, :, idx], medfilt)
        print 'making the running diff mapcube.....'
        for idx, ll in enumerate(tqdm(self.mapcube)):
            maplist.append(deepcopy(ll))
            if idx - dt_frm < 0:
                sidx = 0
            else:
                sidx = idx - dt_frm
            if modes[mode] == 'rdiff':
                maplist[idx].data = datacube[:, :, idx] - datacube[:, :, sidx]
            elif modes[mode] == 'rratio':
                maplist[idx].data = datacube[:, :, idx] / datacube[:, :, sidx]
            elif modes[mode] == 'bdiff':
                maplist[idx].data = datacube[:, :, idx] - datacube[:, :, 0]
            elif modes[mode] == 'bratio':
                maplist[idx].data = datacube[:, :, idx] / datacube[:, :, 0]
        mapcube_diff = sunpy.map.Map(maplist, cube=True)

        if bfilter:
            datacube = mapcube_diff.as_array()
            datacube_ft = np.zeros_like(datacube)
            ny, nx, nt = datacube_ft.shape
            fs = len(mapcube_diff)
            ncpu = mp.cpu_count() - 1
            print 'filtering the mapcube in time domain.....'
            for ly in tqdm(xrange(ny)):
                b_filter_partial = partial(b_filter, datacube[ly], lowcut, highcut, fs)
                pool = mp.Pool(ncpu)
                res = pool.map(b_filter_partial, xrange(nx))
                pool.close()
                pool.join()
                for lx in xrange(nx):
                    datacube_ft[ly, lx] = res[lx]['y']

            for idx, ll in enumerate(tqdm(mapcube_diff)):
                mapcube_diff[idx].data = datacube_ft[:, :, idx]

        if tosave:
            if not outfile:
                outfile = 'mapcube_{5}_{0}_bin{3}_dtdata{4}_{1}_{2}'.format(self.mapcube[0].meta['wavelnth'],
                                                                            self.trange[0].isot[:-4].replace(':', ''),
                                                                            self.trange[1].isot[:-4].replace(':', ''),
                                                                            self.binpix, self.dt_data, modes[mode])
            self.mapcube_tofile(outfile=outfile, mapcube=mapcube_diff)
        self.mapcube_diff = mapcube_diff
        return mapcube_diff

    def plot_mapcube(self, mapcube=None, hdr=False, vmax=None, vmin=None, diff=False, sav_img=False, out_dir=None,
                     dpi=100, anim=False):
        '''

        :param mapcube:
        :param hdr:
        :param vmax:
        :param vmin:
        :param diff:
        :param sav_img:
        :param out_dir:
        :param dpi:
        :param anim:
        :return:
        '''
        if mapcube:
            mapcube_plot = deepcopy(mapcube)
        else:
            if diff:
                mapcube_plot = deepcopy(self.mapcube_diff)
            else:
                mapcube_plot = deepcopy(self.mapcube)
        if mapcube_plot is None:
            print('No mapcube found. Load a mapcube first!')
            return
        if not isinstance(mapcube_plot, sunpy.map.mapcube.MapCube):
            print('mapcube must be a instance of sunpy.map.mapcube.MapCube')
            return
        if hdr:
            for idx, smap in enumerate(tqdm(mapcube_plot)):
                smap = DButil.sdo_aia_scale_hdr(smap)
                mapcube_plot[idx].data = smap.data
        if not diff:
            for idx, smap in enumerate(tqdm(mapcube_plot)):
                mapcube_plot[idx].data[np.where(smap.data < 1)] = 1
        self.mapcube_plot = mapcube_plot
        # sp = stackplot(parent_obj = self, mapcube = mapcube_plot)
        fig_mapcube = plt.figure(figsize=(8, 7))
        try:
            clrange = DButil.sdo_aia_scale_dict(mapcube_plot[0].meta['wavelnth'])
        except:
            clrange = {'high': None, 'log': False, 'low': None}
        if not vmax:
            vmax = clrange['high']
        if not vmin:
            vmin = clrange['low']
        if sav_img:
            if out_dir is None:
                out_dir = './'

            ax, im1 = plot_map(mapcube_plot[0], vmax=vmax, vmin=vmin, diff=diff, returnImAx=True)
            if anim:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.animation as animation
                nframe = len(mapcube_plot)

                def update_frame(num):
                    smap = mapcube_plot[int(num)]
                    # smap.data[smap.data<1]=1
                    im1.set_data(smap.data)
                    # im1.set_extent(list(smap.xrange.value) + list(smap.yrange.value))
                    ax.set_title(
                        '{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['t_obs']))
                    fig_mapcube.canvas.draw()
                    return

                ani = animation.FuncAnimation(fig_mapcube, update_frame, nframe, interval=50, blit=False)

            prompt = ''
            while not (prompt.lower() in ['y', 'n']):
                prompt = raw_input('Satisfied with current FOV? [y/n]')
            if prompt.lower() == 'n':
                return
            if anim:
                print 'Saving movie to {}'.format(out_dir)
                Writer = animation.writers['ffmpeg']
                writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
                ani.save('{0}/{2}{1}.mp4'.format(out_dir, mapcube_plot[0].meta['wavelnth'], mapcube_plot[0].detector),
                         writer=writer)
            else:
                plt.ioff()
                print 'Saving images to {}'.format(out_dir)
                for smap in tqdm(mapcube_plot):
                    im1.set_data(smap.data)
                    # im1.set_extent(list(smap.xrange.value) + list(smap.yrange.value))
                    # ax.set_xlim(list(smap.xrange.value))
                    # ax.set_ylim(list(smap.yrange.value))
                    if smap.meta.has_key('t_obs'):
                        tstr = smap.meta['t_obs']
                    else:
                        tstr = smap.meta['date-obs']
                    ax.set_title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, tstr))
                    t_map = Time(tstr)
                    fig_mapcube.canvas.draw()
                    fig_mapcube.savefig('{0}/{3}{1}-{2}.png'.format(out_dir, smap.meta['wavelnth'],
                                                                    t_map.iso.replace(' ', 'T').replace(':',
                                                                                                        '').replace('-',
                                                                                                                    '')[
                                                                    :-4], smap.detector), format='png', dpi=dpi)
                plt.ion()
        else:
            ax, im1 = plot_map(mapcube_plot[0], vmax=vmax, vmin=vmin, diff=diff, returnImAx=True)
            plt.subplots_adjust(bottom=0.10)
            dims = mapcube_plot[0].dimensions
            diagpix = int(np.sqrt(dims[0] ** 2 + dims[1] ** 2).value)
            axcolor = 'lightgoldenrodyellow'
            # axStackplt = plt.axes([0.8, 0.02, 0.10, 0.05], facecolor=axcolor)
            # bStackplt = Button(axStackplt, 'StackPlt')
            pixscale = ((self.fov[1] - self.fov[0]) / dims[0].value + (self.fov[3] - self.fov[2]) / dims[1].value) / 2.0
            axFrame = plt.axes([0.10, 0.03, 0.40, 0.02], facecolor=axcolor)
            sFrame = Slider(axFrame, 'frame', 0, len(mapcube_plot) - 1, valinit=0, valfmt='%0.0f')
            axCutwdth = plt.axes([0.65, 0.02, 0.10, 0.01], facecolor=axcolor)
            sCutwdth = Slider(axCutwdth, 'Width[pix]', 1, int(diagpix / 4.0), valinit=5, valfmt='%0.0f')
            axCutang = plt.axes([0.65, 0.04, 0.10, 0.01], facecolor=axcolor)
            sCutang = Slider(axCutang, 'Angle[deg]', 0.0, 90.0, valinit=3.0, valfmt='%.1f')
            axCutlngth = plt.axes([0.65, 0.06, 0.10, 0.01], facecolor=axcolor)
            sCutlngth = Slider(axCutlngth, 'Length[pix]', 20, int(diagpix * 4), valinit=150, valfmt='%0.0f')
            self.cutslitbd = CutslitBuilder(ax, cutwidth=sCutwdth.val, cutang=sCutang.val / 180. * np.pi,
                                            cutlength=sCutlngth.val, scale=pixscale)

            # def bStackplt_update(event):
            #     # print bStackplt.val
            #     print('button clicked')
            #
            # bStackplt.on_clicked(bStackplt_update)

            def sFrame_update(val):
                frm = sFrame.val
                smap = mapcube_plot[int(frm)]
                # smap.data[smap.data<1]=1
                im1.set_data(smap.data)
                if smap.meta.has_key('t_obs'):
                    tstr = smap.meta['t_obs']
                else:
                    tstr = smap.meta['date-obs']
                ax.set_title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, tstr))
                fig_mapcube.canvas.draw()

            sFrame.on_changed(sFrame_update)

            def sCutwdth_update(val):
                wdth = sCutwdth.val
                self.cutslitbd.cutwidth = wdth
                self.cutslitbd.update()

            sCutwdth.on_changed(sCutwdth_update)

            def sCutang_update(val):
                ang = sCutang.val / 180. * np.pi
                self.cutslitbd.cutang = ang
                self.cutslitbd.update()

            sCutang.on_changed(sCutang_update)

            def sCutlngth_update(val):
                lngth = sCutlngth.val
                self.cutslitbd.cutlength = lngth
                self.cutslitbd.update()

            sCutlngth.on_changed(sCutlngth_update)
        return

    def cutslit_fromfile(self, infile):
        if self.cutslitbd:
            with open('{}'.format(infile), 'rb') as sf:
                cutslit = pickle.load(sf)
            self.cutslitbd.clickedpoints.set_data(cutslit['x'], cutslit['y'])
            self.cutslitbd.clickedpoints.figure.canvas.draw()
            self.cutslitbd.cutslitplt = cutslit['cutslit']
            self.cutslitbd.slitline.set_data(cutslit['cutslit']['xcen'], cutslit['cutslit']['ycen'])
            self.cutslitbd.slitline0.set_data(cutslit['cutslit']['xs0'], cutslit['cutslit']['ys0'])
            self.cutslitbd.slitline1.set_data(cutslit['cutslit']['xs1'], cutslit['cutslit']['ys1'])
            self.cutslitbd.slitline.figure.canvas.draw()
            self.cutslitbd.slitline0.figure.canvas.draw()
            self.cutslitbd.slitline1.figure.canvas.draw()
        else:
            print('plot_mapcube first before loading cutslit from file!')

    def cutslit_tofile(self, outfile=None, cutslit=None):
        if not cutslit:
            cutslit = {'x': self.cutslitbd.clickedpoints.get_xdata(), 'y': self.cutslitbd.clickedpoints.get_ydata(),
                       'cutslit': self.cutslitbd.cutslitplt}
        with open('{}'.format(outfile), 'wb') as sf:
            pickle.dump(cutslit, sf)

    def make_stackplot(self, mapcube):
        stackplt = []
        print 'making the stack plot...'
        for idx, smap in enumerate(tqdm(mapcube)):
            intens = getimprofile(smap.data, self.cutslitbd.cutslitplt, xrange=smap.xrange.value,
                                  yrange=smap.yrange.value)
            stackplt.append(intens['y'])
        if len(stackplt) > 1:
            stackplt = np.vstack(stackplt)
            self.stackplt = stackplt.transpose()
        else:
            print('Too few timestamps. Failed to make a stack plot map.')

    def plot_stackplot(self, mapcube=None, hdr=False, vmax=None, vmin=None, diff=False, uni_cm=False, cutslitplt=None):
        if mapcube:
            mapcube_plot = deepcopy(mapcube)
        else:
            mapcube_plot = self.mapcube_plot
        if mapcube_plot is None:
            print('No mapcube found. Load a mapcube first!')
            return
        if not isinstance(mapcube_plot, sunpy.map.mapcube.MapCube):
            print('mapcube must be a instance of sunpy.map.mapcube.MapCube')
            return
        self.make_stackplot(mapcube_plot)
        fig_mapcube = plt.figure(figsize=(14, 7))
        try:
            clrange = DButil.sdo_aia_scale_dict(mapcube_plot[0].meta['wavelnth'])
        except:
            clrange = {'high': None, 'log': False, 'low': None}
        if not vmax:
            vmax = clrange['high']
        if not vmin:
            vmin = clrange['low']
        norm = colors.Normalize(vmin=np.min(self.stackplt), vmax=np.max(self.stackplt))
        cutslitplt = self.cutslitbd.cutslitplt
        try:
            dspec = {'dspec': self.stackplt, 'x': self.tplt.plot_date, 'y': cutslitplt['dist'],
                     'ytitle': 'Distance [arcsec]', 'ctitle': 'DN counts per second',
                     'args': {'norm': norm, 'cmap': cm.get_cmap('sdoaia{}'.format(mapcube_plot[0].meta['wavelnth']))}}
        except:
            dspec = {'dspec': self.stackplt, 'x': self.tplt.plot_date, 'y': cutslitplt['dist'],
                     'ytitle': 'Distance [arcsec]', 'ctitle': '',
                     'args': {'norm': norm, 'cmap': 'gray_r'}}

        dtplot = np.mean(np.diff(self.tplt.plot_date))
        dspec['axvspan'] = [self.tplt[0].plot_date, self.tplt[0].plot_date + dtplot]
        ax, im1, ax2, im2, vspan = plot_map(mapcube_plot[0], dspec, vmax=vmax, vmin=vmin, diff=diff, returnImAx=True,
                                            uni_cm=uni_cm)
        plt.subplots_adjust(bottom=0.10)
        ax.plot(cutslitplt['xcen'], cutslitplt['ycen'], color='white', ls='solid')
        ax.plot(cutslitplt['xs0'], cutslitplt['ys0'], color='white', ls='dotted')
        ax.plot(cutslitplt['xs1'], cutslitplt['ys1'], color='white', ls='dotted')
        axcolor = 'lightgoldenrodyellow'
        axframe2 = plt.axes([0.1, 0.03, 0.40, 0.02], facecolor=axcolor)
        sframe2 = Slider(axframe2, 'frame', 0, len(mapcube_plot) - 1, valinit=0, valfmt='%0.0f')

        def update2(val):
            frm = int(sframe2.val)
            smap = mapcube_plot[frm]
            im1.set_data(smap.data)
            if smap.meta.has_key('t_obs'):
                tstr = smap.meta['t_obs']
            else:
                tstr = smap.meta['date-obs']
            ax.set_title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, tstr))
            vspan_xy = vspan.get_xy()
            vspan_xy[np.array([0, 1, 4]), 0] = self.tplt[frm].plot_date
            if frm < len(self.tplt) - 1:
                vspan_xy[np.array([2, 3]), 0] = self.tplt[frm + 1].plot_date
            else:
                vspan_xy[np.array([2, 3]), 0] = self.tplt[frm].plot_date
            vspan.set_xy(vspan_xy)
            fig_mapcube.canvas.draw()

        sframe2.on_changed(update2)

    def mapcube_info(self, mapcube=None):
        if mapcube:
            trange = Time([mapcube[0].date, mapcube[-1].date])
            fov = np.hstack([mapcube[0].xrange.value, mapcube[0].yrange.value])
            binpix = int(
                np.round(np.mean([ll.value for ll in mapcube[0].scale]) / self.instrum_meta['SDO/AIA']['scale'].value))
            return {'trange': trange, 'fov': fov, 'binpix': binpix}
        else:
            self.trange = Time([self.mapcube[0].date, self.mapcube[-1].date])
            self.fov = np.hstack([self.mapcube[0].xrange.value, self.mapcube[0].yrange.value])
            self.binpix = int(np.round(
                np.mean([ll.value for ll in self.mapcube[0].scale]) / self.instrum_meta['SDO/AIA']['scale'].value))
            return {'trange': self.trange, 'fov': self.fov, 'binpix': self.binpix}

    # def mapcube2image(self,mapcube=None,figsize=(7,5)):
    #     if mapcube:
    #         pass
    #     else:
    #         mapcube = self.mapcube_plot

    @property
    def tplt(self, mapcube=None):
        if not mapcube:
            mapcube = self.mapcube
        t = []
        for idx, smap in enumerate(mapcube):
            if smap.meta.has_key('t_obs'):
                tstr = smap.meta['t_obs']
            else:
                tstr = smap.meta['date-obs']
            t.append(tstr)
        return Time(t)

    @classmethod
    def set_fits_dir(cls, fitsdir):
        cls.fitsdir = fitsdir
        # def __repr__(self):
        #     if self.mpcube:
        #         print('')

        # @property
        # def mapcubefile(self):
        #     return self.mapcubefile
        #
        # @mapcubefile.setter
        # def mapcubefile(self, infile):
        #     with open(infile, 'rb') as sf:
        #         print('Loading mapcube....')
        #         self.mpcube = pickle.load(sf)
        #
        # @mapcubefile.deleter
        # def mapcubefile(self):
        #     print('Delete mapcube')
        #     self.mpcube = None
