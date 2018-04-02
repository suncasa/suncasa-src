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
from scipy.signal import butter, lfilter
from scipy.signal import medfilt
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


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


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
        cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                      'posangs': [], 'posangs2': [], 'dist': []}
    else:
        # if len(pointDFtmp.index) <= 3:
        cutslitplt = FitSlit(xx, yy, 10, 0.0,
                             200, method='Polyfit')
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
            inten = DButil.improfile(data, [xs0[ll], xs1[ll]],
                                     [ys0[ll], ys1[ll]], interp='nearest')
            intens[ll] = np.mean(inten)
        intensdist = {'x': cutslit['dist'], 'y': intens}
        return intensdist


def plot_map(smap, dspec=None, diff=False, SymLogNorm=False, linthresh=0.5, returnImAx=False, *args, **kwargs):
    import sunpy.cm.cm as cm  ## to bootstrap sdoaia color map
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    from suncasa.utils import DButil
    clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
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
        imshow_args = {'cmap': 'gray', 'norm': norm,
                       'interpolation': 'nearest', 'origin': 'lower'}
    if smap.coordinate_system.x == 'HG':
        xlabel = 'Longitude [{lon}]'.format(lon=smap.spatial_units.x)
    else:
        xlabel = 'X-position [{xpos}]'.format(xpos=smap.spatial_units.x)
    if smap.coordinate_system.y == 'HG':
        ylabel = 'Latitude [{lat}]'.format(lat=smap.spatial_units.y)
    else:
        ylabel = 'Y-position [{ypos}]'.format(ypos=smap.spatial_units.y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    imshow_args.update({'extent': list(smap.xrange.value) + list(smap.yrange.value)})
    if smap.detector == 'HMI':
        im1 = ax.imshow(np.rot90(smap.data, 2), **imshow_args)
    else:
        im1 = ax.imshow(smap.data, **imshow_args)
    plt.title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['t_obs']))
    plt.colorbar(im1, ax=ax, label='DN counts per second')
    ax.set_autoscale_on(False)
    if dspec:
        fig = plt.gcf()
        # if figsize:
        #     fig.set_size_inches(figsize)
        # else:
        #     fig.set_size_inches(13, 5)
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


datadir = '/Users/fisher/Desktop/work/2017/LLL/aia/'
mkmc = 1
if mkmc:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate

    t1, t2 = Time('2014-11-09T10:00:00'), Time('2014-11-09T13:00:00')
    binwdth = 2
    dt_data = 1
    trange = [t1.jd, t2.jd]
    SDOdir = '/Volumes/NAOC-001/work/database/aiaBrowserData/Download/'
    sdofile = DButil.readsdofile(datadir=SDOdir, wavelength='171', jdtime=trange)
    # SDOdir = '/Volumes/NAOC-001/work/2015/20141111/data/prep/0171/'
    # sdofile = glob.glob(SDOdir + 'AIA20141109_*_0171.fits')
    # sdofile = sorted(sdofile)
    # x0, x1, y0, y1 = -660., -575., -200., -115.
    x0, x1, y0, y1 = -850, -500., 730., 1200.  # large
    maplist = []
    print 'Loading fits files....'
    for ll in tqdm(sdofile[::dt_data]):
        maptmp = sunpy.map.Map(ll)
        submaptmp = maptmp.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                  u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        submaptmp = submaptmp.resample(u.Quantity(submaptmp.dimensions) / binwdth)
        if submaptmp.detector == 'HMI':
            pass
        else:
            submaptmp = DButil.normalize_aiamap(submaptmp)
        maplist.append(submaptmp)
        # print '{}/{} loaded.'.format(idx + 1, len(sdofile))
    doderotate = False
    if doderotate:
        mapcube = mapcube_solar_derotate(sunpy.map.Map(maplist, cube=True))
    else:
        mapcube = sunpy.map.Map(maplist, cube=True)
    outfile = '{0}/mapcube_{1}_bin{4}_dtdata{5}_{2}_{3}'.format(datadir, mapcube[0].meta['wavelnth'],
                                                                t1.isot[:-4].replace(':', ''),
                                                                t2.isot[:-4].replace(':', ''), binwdth, dt_data)
    with open(outfile, 'wb') as sf:
        print 'Saving mapcube....'
        pickle.dump(mapcube, sf)
    gc.collect()
else:
    with open('{}/mapcube_{}'.format(datadir, 304), 'rb') as sf:
        print 'Loading mapcube....'
        mapcube = pickle.load(sf)
        # submaplist = []
        # for sidx, smap in enumerate(mapcube0):
        #     submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
        #                             u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        #     submaplist.append(submaptmp)
        # mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))

rdiff_ratio = 1
mapcubename = 'mapcube_171_bin2_dtdata1_2014-11-09T140000_2014-11-09T200000'
mapcubename = 'mapcube_304_bin2_dtdata1_2014-11-09T150000_2014-11-09T200000'
mapcubename = 'mapcube_171_bin2_dtdata1_2014-11-09T150000_2014-11-09T200000'
mapcubename = 'mapcube_171_bin2_dtdata1_2014-11-09T100000_2014-11-09T130000'
loadmapcube = True
if loadmapcube:
    with open('{}/{}'.format(datadir, mapcubename), 'rb') as sf:
        print 'Loading mapcube ....'
        mapcube = pickle.load(sf)
if mapcube[0].meta['wavelnth'] == 304:
    for idx, smap in enumerate(mapcube):
        mapcube[idx].data[mapcube[idx].data < 1.0] = 1.0

if rdiff_ratio:
    dt_frm = 6
    domedfilt = 0
    # figure()
    # plt.subplot(121)
    # plt.imshow(mapcube[0].data)
    # plt.subplot(122)
    # plt.imshow(medfilt(mapcube[0].data, 3))
    maps_rdiff = []
    datacube = mapcube.as_array()
    if domedfilt:
        for idx, ll in enumerate(tqdm(mapcube)):
            datacube[:, :, idx] = medfilt(datacube[:, :, idx], 3)
    print 'making the rdiff mapcube.....'
    for idx, ll in enumerate(tqdm(mapcube)):
        maps_rdiff.append(deepcopy(ll))
        if idx - dt_frm < 0:
            sidx = 0
        else:
            sidx = idx - dt_frm
        maps_rdiff[idx].data = datacube[:, :, idx] - datacube[:, :, sidx]

    mapcube_rdiff = sunpy.map.Map(maps_rdiff, cube=True)
    outdir = '{0}img_{1}_diff'.format(datadir, mapcubename)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        os.system('rm -rf {}/*.png'.format(outdir))
    datacube = mapcube_rdiff.as_array()
    # fig_cutslit = plt.figure(5, figsize=(7, 5))
    # ax = plot_map(mapcube_rdiff[0], vmax=500, vmin=-500, diff=True)
    datacube_ft = np.zeros_like(datacube)
    ny, nx, nt = datacube_ft.shape


    def b_filter(data, lowcut, highcut, fs, ix):
        x = data[ix]
        y = butter_bandpass_filter(x, lowcut, highcut, fs, order=3)
        return {'idx': ix, 'y': y}


    dofilter = True
    if dofilter:
        fs = len(mapcube_rdiff)
        lowcut = 0.1
        highcut = 50
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

        for idx, ll in enumerate(tqdm(mapcube_rdiff)):
            mapcube_rdiff[idx].data = datacube_ft[:, :, idx]
        outfile = '{0}/mapcube_{1}_bin{4}_dtdata{5}_{2}_{3}_diff_filter'.format(datadir, mapcube[0].meta['wavelnth'],
                                                                                t1.isot[:-4].replace(':', ''),
                                                                                t2.isot[:-4].replace(':', ''), binwdth,
                                                                                dt_data)
        with open(outfile, 'wb') as sf:
            print 'Saving mapcube....'
            pickle.dump(mapcube_rdiff, sf)
    else:
        mapcube_rdiff_name = 'mapcube_171_bin2_dtdata1_2014-11-09T100000_2014-11-09T130000_diff_filter'
        with open('{}/{}'.format(datadir, mapcube_rdiff_name), 'rb') as sf:
            print 'Loading mapcube ....'
            mapcube_rdiff = pickle.load(sf)
    # mapcube_rdiff = mapcube_rdiff
    mapcube_plot = mapcube_rdiff[30:301]

    # x = datacube[300, 300]
    # y = butter_bandpass_filter(x, lowcut, highcut, fs, order=6)
    # clf()
    # plot(x)
    # plot(y)
    doplot = False
    if doplot:
        fig = plt.figure(1, figsize=(7, 5))
        plt.ioff()
        if mapcube_plot[0].meta['wavelnth'] == 304:
            # vmax, vmin, diff = 2, -2, True
            vmax, vmin, diff = 3, 0.8, False
        elif mapcube_plot[0].meta['wavelnth'] == 193:
            vmax, vmin = 1000, -1000
        elif mapcube_plot[0].meta['wavelnth'] == 171:
            vmax, vmin, diff = 10, -10, True
        elif mapcube_plot[0].meta['wavelnth'] == 335:
            vmax, vmin = 50, -50
        else:
            vmax, vmin = 1000, -1000
        print 'making images.....'
        for smap in tqdm(mapcube_plot):
            ax = plot_map(smap, vmax=vmax, vmin=vmin, diff=diff)
            t_map = Time(smap.meta['t_obs'])
            fig.savefig('{0}/{3}{1}-{2}.png'.format(outdir, smap.meta['wavelnth'],
                                                    t_map.iso.replace(' ', 'T').replace(':',
                                                                                        '').replace('-',
                                                                                                    '')[
                                                    :-4], smap.detector), format='png', dpi=100)
        plt.ion()
else:
    mapcube_plot = deepcopy(mapcube[30:301])
    for idx, smap in enumerate(tqdm(mapcube_plot)):
        smap = DButil.sdo_aia_scale_hdr(smap)
        mapcube_plot[idx].data = smap.data

docutslit = 1
if docutslit:
    class CutslitBuilder:
        def __init__(self, axes, cutwidth=5, cutang=0, cutlength=80):
            self.axes = axes
            self.clickedpoints, = self.axes.plot([], [], 'o', color='white')
            self.slitline, = self.axes.plot([], [], color='white', ls='solid')
            self.slitline0, = self.axes.plot([], [], color='white', ls='dotted')
            self.slitline1, = self.axes.plot([], [], color='white', ls='dotted')
            self.cutlength = cutlength
            self.cutwidth = cutwidth
            self.cutang = cutang
            self.xx = list(self.clickedpoints.get_xdata())
            self.yy = list(self.clickedpoints.get_ydata())
            self.cid = self.clickedpoints.figure.canvas.mpl_connect('button_press_event', self)

        def __call__(self, event):
            # print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata))
            if event.inaxes != self.axes:
                return
            if event.button == 1:
                self.xx.append(event.xdata)
                self.yy.append(event.ydata)
            elif event.button == 3:
                self.xx.pop()
                self.yy.pop()
            xx = np.array(self.xx, dtype=np.float64)
            yy = np.array(self.yy, dtype=np.float64)
            self.clickedpoints.set_data(self.xx, self.yy)
            self.clickedpoints.figure.canvas.draw()

            if len(self.xx) <= 1:
                cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                              'posangs': [], 'posangs2': [], 'dist': []}
            else:
                if len(self.xx) <= 3:
                    cutslitplt = FitSlit(xx, yy, self.cutwidth, self.cutang, self.cutlength, method='Polyfit')
                else:
                    cutslitplt = FitSlit(xx, yy, self.cutwidth, self.cutang, self.cutlength, s=len(xx),
                                         method='Param_Spline')
            self.cutslitplt = cutslitplt
            self.slitline.set_data(cutslitplt['xcen'], cutslitplt['ycen'])
            self.slitline0.set_data(cutslitplt['xs0'], cutslitplt['ys0'])
            self.slitline1.set_data(cutslitplt['xs1'], cutslitplt['ys1'])
            self.slitline.figure.canvas.draw()
            self.slitline0.figure.canvas.draw()
            self.slitline1.figure.canvas.draw()


    dt_frm = 3
    domedfilt = 0
    # mapcubename = 'mapcube_171_bin2_dtdata5_2014-11-09T150000_2014-11-09T170000'

    tplot = []
    for idx, smap in enumerate(tqdm(mapcube_plot)):
        # maps_rdiff.append(deepcopy(smap))
        tplot.append(smap.meta['t_obs'])
    tplot = Time(tplot)
    fig_cutslit = plt.figure(5, figsize=(7, 5))
    clrange = DButil.sdo_aia_scale_dict(mapcube_plot[0].meta['wavelnth'])
    if mapcube_plot[0].meta['wavelnth'] == 304:
        # vmax, vmin, diff = 2, -2, True
        vmax, vmin, diff = 10, 1.0, False
    elif mapcube_plot[0].meta['wavelnth'] == 193:
        vmax, vmin = 1000, -1000
    elif mapcube_plot[0].meta['wavelnth'] == 171:
        vmax, vmin, diff = 10, -10, True
        vmax, vmin, diff = clrange['high'], clrange['low'], False
    elif mapcube_plot[0].meta['wavelnth'] == 335:
        vmax, vmin = 50, -50
    else:
        vmax, vmin = 1000, -1000
    ax, im1 = plot_map(mapcube_plot[0], vmax=vmax, vmin=vmin, diff=diff, returnImAx=True)
    # ax, im1 = plot_map(mapcube_plot[0], vmax=3, vmin=-3, diff=True, returnImAx=True)
    cutslitbuilder = CutslitBuilder(ax, cutwidth=20, cutang=np.pi / 60, cutlength=150)
    axcolor = 'lightgoldenrodyellow'
    axframe = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
    sframe = Slider(axframe, 'frame', 0, len(mapcube_plot) - 1, valinit=0, valfmt='%0.0f')


    def update(val):
        frm = sframe.val
        smap = mapcube_plot[int(frm)]
        # smap.data[smap.data<1]=1
        im1.set_data(smap.data)
        ax.set_title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['t_obs']))
        fig_cutslit.canvas.draw()


    sframe.on_changed(update)

    to_save = 0
    to_restore = 1
    if to_save:
        cutslit = {'x': cutslitbuilder.clickedpoints.get_xdata(), 'y': cutslitbuilder.clickedpoints.get_ydata(),
                   'cutslit': cutslitbuilder.cutslitplt}
        with open('{}cutslit2'.format(datadir), 'wb') as sf:
            pickle.dump(cutslit, sf)
    if to_restore:
        with open('{}cutslit2'.format(datadir), 'rb') as sf:
            cutslit = pickle.load(sf)
        cutslitbuilder.clickedpoints.set_data(cutslit['x'], cutslit['y'])
        cutslitbuilder.clickedpoints.figure.canvas.draw()
        cutslitbuilder.cutslitplt = cutslit['cutslit']
        cutslitbuilder.slitline.set_data(cutslit['cutslit']['xcen'], cutslit['cutslit']['ycen'])
        cutslitbuilder.slitline0.set_data(cutslit['cutslit']['xs0'], cutslit['cutslit']['ys0'])
        cutslitbuilder.slitline1.set_data(cutslit['cutslit']['xs1'], cutslit['cutslit']['ys1'])
        cutslitbuilder.slitline.figure.canvas.draw()
        cutslitbuilder.slitline0.figure.canvas.draw()
        cutslitbuilder.slitline1.figure.canvas.draw()
    dostackplot = 1
    if dostackplot:
        stackplt = []
        print 'making the stack plot...'
        for idx, smap in enumerate(tqdm(mapcube_plot)):
            intens = getimprofile(smap.data, cutslitbuilder.cutslitplt,
                                  xrange=smap.xrange.value,
                                  yrange=smap.yrange.value)
            stackplt.append(intens['y'])
        if len(stackplt) > 1:
            stackplt = np.vstack(stackplt)
            stackplt = stackplt.transpose()
        # if mapcube_plot[0].meta['wavelnth'] == 304:
        #     # vmax, vmin, diff = 2, -2, True
        #     vmax, vmin, diff = 3, 0.8, False
        # elif mapcube_plot[0].meta['wavelnth'] == 193:
        #     vmax, vmin = 1000, -1000
        # elif mapcube_plot[0].meta['wavelnth'] == 171:
        #     vmax, vmin, diff = 5, -5, True
        # elif mapcube_plot[0].meta['wavelnth'] == 335:
        #     vmax, vmin = 50, -50
        # else:
        #     vmax, vmin = 1000, -1000
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        # norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        cutslitplt = cutslitbuilder.cutslitplt
        dspec = {'dspec': stackplt, 'x': tplot.plot_date, 'y': cutslitplt['dist'],
                 'ytitle': 'Distance [arcsec]',
                 'ctitle': 'DN counts per second',
                 'args': {'norm': norm, 'cmap': cm.get_cmap('sdoaia{}'.format(mapcube_plot[0].meta['wavelnth']))}}
        with open(outfile + '_dspec', 'wb') as sf:
            print 'Saving dspec....'
            pickle.dump(dspec, sf)

        with open(outfile + '_dspec', 'rb') as sf:
            print 'Loading dspec....'
            dspec = pickle.load(sf)

        fig_stackplt = plt.figure(6, figsize=(12, 5))
        dtplot = np.mean(np.diff(tplot.plot_date))
        dspec['axvspan'] = [tplot[0].plot_date, tplot[0].plot_date + dtplot]
        ax, im1, ax2, im2, vspan = plot_map(mapcube_plot[0], dspec, vmax=vmax, vmin=vmin, diff=diff,
                                            returnImAx=True)
        ax.plot(cutslitplt['xcen'], cutslitplt['ycen'], color='white', ls='solid')
        ax.plot(cutslitplt['xs0'], cutslitplt['ys0'], color='white', ls='dotted')
        ax.plot(cutslitplt['xs1'], cutslitplt['ys1'], color='white', ls='dotted')
        axcolor = 'lightgoldenrodyellow'
        axframe2 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
        sframe2 = Slider(axframe2, 'frame', 0, len(mapcube_plot) - 1, valinit=0, valfmt='%0.0f')


        def update2(val):
            frm = int(sframe2.val)
            smap = mapcube_plot[frm]
            im1.set_data(smap.data)
            ax.set_title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['t_obs']))
            vspan_xy = vspan.get_xy()
            vspan_xy[np.array([0, 1, 4]), 0] = tplot[frm].plot_date
            vspan_xy[np.array([2, 3]), 0] = tplot[frm].plot_date + dtplot
            vspan.set_xy(vspan_xy)
            fig_stackplt.canvas.draw()


        sframe2.on_changed(update2)

        doplot = False
        if doplot:
            outdir = '{0}img_{1}_diff_stackplt'.format(datadir, mapcubename)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            else:
                os.system('rm -rf {}/*.png'.format(outdir))
            print 'producing images...'
            for xidx, smap in enumerate(tqdm(mapcube_plot)):
                dspec['axvspan'] = [tplot[xidx].plot_date, tplot[xidx].plot_date + dtplot]
                ax, ax2 = plot_map(smap, dspec, vmax=vmax, vmin=vmin, diff=True)
                ax.plot(cutslitplt['xcen'], cutslitplt['ycen'], color='white', ls='solid')
                ax.plot(cutslitplt['xs0'], cutslitplt['ys0'], color='white', ls='dotted')
                ax.plot(cutslitplt['xs1'], cutslitplt['ys1'], color='white', ls='dotted')
                fig_stackplt.savefig('{0}/{3}{1}-{2}.png'.format(outdir, smap.meta['wavelnth'],
                                                                 tplot[xidx].iso.replace(' ', 'T').replace(':',
                                                                                                           '').replace(
                                                                     '-',
                                                                     '')[
                                                                 :-4], smap.detector), format='png', dpi=100)
