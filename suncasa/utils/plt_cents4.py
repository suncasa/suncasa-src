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
import sunpy.cm.cm as cm  ## to bootstrap sdoaia color map
import matplotlib.cm as cm
import sunpy.map
import astropy.units as u
import pickle
import glob
from suncasa.utils import ctplot as ct
# from astropy import units as u
# import sunpy.map as smap
from scipy.interpolate import griddata
from suncasa.utils import DButil
from scipy.interpolate import splev, splrep
import scipy.ndimage
# import pickle
from IPython import embed
import time
from tqdm import *


# def smooth(y, box_pts):
#     box = np.ones(box_pts) / box_pts
#     y_smooth = np.convolve(y, box, mode='same')
#     return y_smooth

# def plot_map(smap, position=None):
#     import sunpy.cm.cm as cm  ## to bootstrap sdoaia color map
#     import matplotlib.cm as cm
#     import matplotlib.colors as colors
#     from suncasa.utils import DButil
#     clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
#     plt.clf()
#     ax = plt.subplot()
#     if position:
#         ax.set_position(position)
#     norm = colors.LogNorm(vmin=clrange['low'], vmax=clrange['high'])
#     imshow_args = {'cmap': cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth'])), 'norm': norm,
#                    'interpolation': 'nearest', 'origin': 'lower'}
#     if smap.coordinate_system.x == 'HG':
#         xlabel = 'Longitude [{lon}]'.format(lon=smap.spatial_units.x)
#     else:
#         xlabel = 'X-position [{xpos}]'.format(xpos=smap.spatial_units.x)
#     if smap.coordinate_system.y == 'HG':
#         ylabel = 'Latitude [{lat}]'.format(lat=smap.spatial_units.y)
#     else:
#         ylabel = 'Y-position [{ypos}]'.format(ypos=smap.spatial_units.y)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     imshow_args.update({'extent': list(smap.xrange.value) + list(smap.yrange.value)})
#     cax = ax.imshow(smap.data, **imshow_args)
#     plt.title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, str(smap.date)[:-7]))
#     plt.colorbar(cax, ax=ax)
#     ax.set_autoscale_on(False)
#     # ax.autoscale(True, 'both', True)
#     # ax.autoscale_view(True, True, True)
#     # ax.relim(visible_only=True)
#     return ax

def FitSlit(xx, yy, cutwidth, cutang, cutlength, s=None, method='Polyfit'):
    # if len(xx) <= 3 or method == 'Polyfit':
    #     '''polynomial fit'''
    out = DButil.polyfit(xx, yy, cutlength, len(xx) - 1 if len(xx) <= 3 else 2)
    xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    # else:
    #     if method == 'Param_Spline':
    #         '''parametic spline fit'''
    #         out = DButil.paramspline(xx, yy, cutlength, s=s)
    #         xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    #     else:
    #         '''spline fit'''
    #         out = DButil.spline(xx, yy, cutlength, s=s)
    #         xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    # if not ascending and (fitmethod != 'Param_Spline' or len(xx) <= 3):
    #     xs, ys = xs[::-1], ys[::-1]
    #     posangs = posangs[::-1]
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


def plot_map(smap, dspec=None, diff=False, SymLogNorm=False, linthresh=0.5, smap_colorbar=True, *args, **kwargs):
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
        if 'cmap' in kwargs.keys():
            cmap = kwargs['cmap']
        else:
            cmap = cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth']))
        imshow_args = {'cmap': cmap, 'norm': norm, 'interpolation': 'nearest', 'origin': 'lower'}
    except:
        imshow_args = {'cmap': 'gray', 'norm': norm, 'interpolation': 'nearest', 'origin': 'lower'}
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
        cax = ax.imshow(np.rot90(smap.data, 2), **imshow_args)
    else:
        cax = ax.imshow(smap.data, **imshow_args)
    plt.title('{} {} {} {}'.format(smap.observatory, smap.detector, smap.wavelength, smap.meta['t_obs']))
    if smap_colorbar:
        plt.colorbar(cax, ax=ax, label='DN counts per second')
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
            ax2.axvspan(dspec['axvspan'][0], dspec['axvspan'][1], alpha=0.5, color='white')
        if 'xs' in dspec.keys() and 'ys' in dspec.keys():
            ax2.plot(dspec['xs'], dspec['ys'], '--', lw=2.0, alpha=0.7, c='black')
        if 'xlim' in dspec.keys():
            ax2.set_xlim(dspec['xlim'])
        if 'ylim' in dspec.keys():
            ax2.set_ylim(dspec['ylim'])
        return ax, ax2
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


class LineBuilder:
    def __init__(self, scatters, lines, cutlength=50):
        self.line = []
        self.scatter = []
        self.cutlength = cutlength
        self.xs = list(scatters[0].get_xdata())
        self.ys = list(scatters[0].get_ydata())
        self.xsline, self.ysline = [], []
        for idx, s in enumerate(scatters):
            self.line.append(lines[idx])
            self.scatter.append(s)
            self.cid = s.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        # print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata))
        axes = [scatter.axes for scatter in self.scatter]
        if scatter.figure.canvas.toolbar.mode == '':
            if event.inaxes not in axes:
                return
            if event.button == 1:
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
            elif event.button == 3:
                if len(self.xs) > 0:
                    self.xs.pop()
                    self.ys.pop()
            xs = np.array(self.xs, dtype=np.float64)
            # ys = np.array(self.ys, dtype=np.float64)
            for scatter in self.scatter:
                scatter.set_data(self.xs, self.ys)
                scatter.figure.canvas.draw()
            xs0 = np.float64(np.nanmin(xs))
            nxs = len(self.xs)

            if nxs >= 2:
                self.cutlength = np.ceil((self.xs[-1] - self.xs[0]) * 24 * 3600 / 0.1)
                if nxs <= 3:
                    out = DButil.polyfit((xs - xs0) * 24. * 3600., self.ys, self.cutlength,
                                         len(self.xs) - 1 if len(self.xs) <= 3 else 2)
                else:
                    out = DButil.paramspline((xs - xs0) * 24. * 3600., self.ys, self.cutlength, s=0.00001)
                out['xs'] = out['xs'] / 24. / 3600. + xs0
                for line in self.line:
                    line.set_data(out['xs'], out['ys'])
                    line.figure.canvas.draw()
                self.xsline, self.ysline = out['xs'], out['ys']
            else:
                for line in self.line:
                    line.set_data([], [])
                    line.figure.canvas.draw()


def pltspec(x, y, specs, vrange, pltline=True, figsize=(10, 8)):
    if len(specs) == 3:
        date_format = mdates.DateFormatter('%H:%M:%S')
        fig, axs = plt.subplots(3, 1, figsize=figsize, sharey=True, sharex=True)
        im1 = axs[0].pcolormesh(x, y, specs[0], vmin=vrange['spec'][0], vmax=vrange['spec'][1], cmap=cm.jet)
        axs[0].xaxis_date()
        axs[0].xaxis.set_major_formatter(date_format)
        axs[0].yaxis.set_label_text('Frequency [GHz]')
        fig.colorbar(im1, orientation='vertical', ax=axs[0], label='flux [sfu]')
        im2 = axs[1].pcolormesh(x, y, specs[1], vmin=vrange['x'][0], vmax=vrange['x'][1], cmap=cm.jet)
        axs[1].xaxis_date()
        axs[1].xaxis.set_major_formatter(date_format)
        axs[1].yaxis.set_label_text('Frequency [GHz]')
        fig.colorbar(im2, orientation='vertical', ax=axs[1], label='x [arcsec]')
        im3 = axs[2].pcolormesh(x, y, specs[2], vmin=vrange['y'][0], vmax=vrange['y'][1], cmap=cm.jet)
        axs[2].xaxis_date()
        axs[2].xaxis.set_major_formatter(date_format)
        axs[2].yaxis.set_label_text('Frequency [GHz]')
        fig.colorbar(im3, orientation='vertical', ax=axs[2], label='y [arcsec]')
        fig.autofmt_xdate(rotation=45)
        fig.tight_layout()
        if pltline:
            line1, = axs[0].plot([], [], '-', c='gray', alpha=0.7)  # empty line
            line2, = axs[1].plot([], [], '-', c='gray', alpha=0.7)  # empty line
            line3, = axs[2].plot([], [], '-', c='gray', alpha=0.7)  # empty line
            lines = [line1, line2, line3]
            scatter1, = axs[0].plot([], [], 'o', c='magenta', alpha=0.7)  # empty line
            scatter2, = axs[1].plot([], [], 'o', c='magenta', alpha=0.7)  # empty line
            scatter3, = axs[2].plot([], [], 'o', c='magenta', alpha=0.7)  # empty line
            scatters = [scatter1, scatter2, scatter3]
            linebuilder = LineBuilder(scatters, lines, 40)
            return [fig, axs, lines, scatters, linebuilder]
        else:
            return [fig, axs]


vrange = {'U01': {'spec': [3, 20], 'x': [-645, -630], 'y': [-140, -125], 'f': [1.1, 1.4], 'xlim': None, 'ylim': None},
          'U02': {'spec': [10, 25], 'x': [-640, -635], 'y': [-135, -115], 'f': [1.1, 1.35],
                  'xlim': [735538.69837392925, 735538.69847811805], 'ylim': [0.994, 1.426]},
          'U04': {'spec': [10, 50], 'x': [-645, -630], 'y': [-135, -120], 'f': [1.0, 1.35],
                  'xlim': [735538.69878972694, 735538.69896720233], 'ylim': [0.994, 1.470]}}
dspec_cmaps = ['jet'] * 3
branch_colors = ['cyan', 'green', 'royalblue']
doshift = 0
domovie = 0
domovie2 = 0
dopltvdspec = 0
mkmc = 0
dosmooth = 0
docorrect = 0
mknewfits = 0
rdiff_ratio = 0
docutslit = 1
mkmovie = 0
mkmovie2 = 0
mkmovie3 = 0
mkmovie4 = 0
dopltmass = 0
dopltmass_mlab = 0
# --------------
timelineplt = 0
timcontour = 0
tbin = 8
mask = ''
thrshd = 10
dmax = 55
box_pts = 11
# --------------
burstname = 'U04'
os.chdir('/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/{}-centroid'.format(burstname))
datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'

tstep = tbin * 2

with open('rms-{}'.format(burstname), 'rb') as sf:
    rms = pickle.load(sf)
cents = np.load('centroidsLL.npy')
# cents = np.load('centroidsLL.npy')
x = cents.item()['x']
y = cents.item()['y']
# f = cents.item()['freq']
# t = cents.item()['time']
p = cents.item()['peak']
f = cents.item()['freq']
fitsarr = cents.item()['fitsarr']
ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
t = Time(cents.item()['time'] / 3600. / 24., format='jd')
# tdatetime = mdates.date2num(t.datetime)
tdatetime = t.plot_date
r = np.empty_like(p)
r[:] = np.NAN
f0 = rms['freq'].index('{:.3f}'.format(f[0]))
f1 = len(f) + f0

for idx, ll in enumerate(t):
    dt = np.min(np.abs(ll.jd - rms['time'].jd)) * 24 * 3600 * 1000
    if dt < 30:
        try:
            r[:, idx] = rms['rms'][f0:f1, np.argmin(np.abs(ll.jd - rms['time'].jd))]
        except:
            r[:-(len(f) - len(rms['freq'])), idx] = rms['rms'][:, np.argmin(np.abs(ll.jd - rms['time'].jd))]
# trange = [Time('2014-11-01T16:45:38.725').plot_date,
#           Time('2014-11-01T16:45:52.525').plot_date]
trange = []
if trange:
    tidx = np.where((tdatetime >= trange[0]) & (tdatetime <= trange[1]))[0]
    tdatetime = tdatetime[tidx[0]:tidx[-1]]
    p = p[:, tidx[0]:tidx[-1]]
    x = x[:, tidx[0]:tidx[-1]]
    y = y[:, tidx[0]:tidx[-1]]
    ff = ff[:, tidx[0]:tidx[-1]]
    fitsarr = fitsarr[:, tidx[0]:tidx[-1]]
# checked that time and frequency axes are all the same
nfreq, ntime = p.shape

# logp=np.log10(p)
# mlogp=ma.masked_outside(logp,dmin,dmax)
# mx=ma.array(x,mask=mlogp.mask)
# my=ma.array(y,mask=mlogp.mask)
# tplt=tdatetime
# logp = np.log10(p)

if burstname == 'U04':
    docorrect = 0
if docorrect:
    bandedge = [63, 127, 191]
    nf = bandedge[2] - bandedge[1] + 1
    xcorr = np.tile(np.linspace(0, -4, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
    # xcorr = np.tile(-(np.exp(np.log(4+1) * np.linspace(0, 1, nf)) - 1), ntime).reshape(ntime, nf).swapaxes(0, 1)
    x[bandedge[1]:bandedge[2] + 1, :] += xcorr
    ycorr = np.tile(np.linspace(0, -4.5, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
    y[bandedge[1]:bandedge[2] + 1, :] += ycorr

if dosmooth:
    # first smooth along the frequency axis
    for i in range(len(tdatetime)):
        # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
        x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
        y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
if mask:
    if mask == 'threshold':
        mp = ma.masked_outside(p, thrshd, dmax)
    elif mask == 'native':
        mp = p
    mx = ma.array(x, mask=mp.mask)
    my = ma.array(y, mask=mp.mask)
    # mp = ma.array(p, mask=mp.mask)
    mf = ma.array(ff, mask=mp.mask)
else:
    # mlogp = logp
    mx = x
    my = y
    mp = p
    mf = f

if dopltvdspec:
    # embed()


    fig1, axs1, lines, scatters, linebuilder = pltspec(tdatetime, f, [mp, mx, my], vrange[burstname])

    to_save = 0
    to_restore = 1
    if to_save:
        sct = [{'x': ll.get_xdata(), 'y': ll.get_ydata()} for ll in scatters]
        lns = [{'x': ll.get_xdata(), 'y': ll.get_ydata()} for ll in lines]
        with open('branch1_status', 'wb') as sf:
            pickle.dump([sct, lns], sf)

    SDOdir = os.getenv('SUNCASADB') + 'aiaBrowserData/Download/'
    aiamap = DButil.readsdofile(datadir=SDOdir, wavelength='304', jdtime=t[0].jd, timtol=20. / 3600. / 24.)
    # x0, x1, y0, y1 = -660., -620., -150., -110.
    x0, x1, y0, y1 = -660., -600., -170., -110.
    submap = aiamap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]), u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
    fig2 = plt.figure(2, figsize=(10, 8))
    ax = plot_map(submap)

    fig3, axs3 = pltspec(tdatetime, f, [mp, mx, my], vrange[burstname], pltline=False, figsize=(8, 8))
    # ax.set_xlim(-660., -620.)
    # ax.set_ylim(-150., -110.)
    # statusfiles = glob.glob('./branch?_status')
    # statusfiles = ['./branch2_status']*2
    statusfiles = ['', 'branch2_status']
    statusfiles = ['branch1_status']
    if to_restore:
        for sidx, statusfile in enumerate(statusfiles):
            if statusfile:
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)

                lines = []
                for idx, line in enumerate(linebuilder.line):
                    line.set_data(lns[idx]['x'], lns[idx]['y'])
                    line.figure.canvas.draw()
                    lines.append(line)
                scatters = []
                for idx, scatter in enumerate(linebuilder.scatter):
                    scatter.set_data(sct[idx]['x'], sct[idx]['y'])
                    scatter.figure.canvas.draw()
                    scatters.append(scatter)

                linebuilder.scatter = scatters
                linebuilder.line = lines
                linebuilder.mapx = list(scatters[0].get_xdata())
                linebuilder.mapy = list(scatters[0].get_ydata())
                linebuilder.xsline, linebuilder.ysline = lns[0]['x'], lns[0]['y']

                xs, ys = linebuilder.xsline, linebuilder.ysline
                for ll in xrange(3):
                    axs3[ll].plot(xs, ys, '-', lw=2.0, c=branch_colors[sidx], alpha=0.7)
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                s0 = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                r0 = r[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                s0 = ma.array(s0, mask=xinterp.mask)
                r0 = ma.array(r0, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                s = ma.compressed(s)
                s0 = ma.compressed(s0)
                r0 = ma.compressed(r0)
                c = ma.compressed(c)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                xplot = xplot  # [::3]
                yplot = yplot  # [::3]
                s = s  # [::3]
                s0 = s0  # [::3]
                r0 = r0  # [::3]
                c = c  # [::3]
                err = (30. * r0 / s0)
                ax.errorbar(xplot, yplot, xerr=err, yerr=err, fmt='none', c='white', alpha=0.3, errorevery=3)
                # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                #           yplot[1:] - yplot[:-1], angles='uv', scale_units=None, scale=None, alpha=0.6,
                #           color=branch_colors[sidx])
                ax.quiver(xplot[:-1], yplot[:-1], (xplot[1:] - xplot[:-1]), (yplot[1:] - yplot[:-1]), angles='uv',
                          scale_units=None, scale=10, alpha=0.6, color=branch_colors[sidx])
                scats = ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                                   vmax=vrange[burstname]['f'][1], alpha=1.0)
                plt.colorbar(scats, ax=ax, label='Frequency')
    else:
        xs, ys = linebuilder.xsline, linebuilder.ysline
        xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
        yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
        xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
        yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
        s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
        s0 = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
        r0 = r[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
        c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
        xinterp = ma.masked_invalid(xinterp)
        yinterp = ma.array(yinterp, mask=xinterp.mask)
        s = ma.array(s, mask=xinterp.mask)
        s0 = ma.array(s0, mask=xinterp.mask)
        r0 = ma.array(r0, mask=xinterp.mask)
        c = ma.array(c, mask=xinterp.mask)
        xinterp = ma.compressed(xinterp)
        yinterp = ma.compressed(yinterp)
        s = ma.compressed(s)
        s0 = ma.compressed(s0)
        r0 = ma.compressed(r0)
        c = ma.compressed(c)
        xplot = smooth(xinterp, box_pts)
        yplot = smooth(yinterp, box_pts)
        # xplot = xinterp
        # yplot = yinterp
        xplot = xplot  # [::3]
        yplot = yplot  # [::3]
        s = s  # [::3]
        s0 = s0  # [::3]
        r0 = r0  # [::3]
        c = c  # [::3]

        err = (30. * r0 / s0)
        ax.errorbar(xplot, yplot, xerr=err, yerr=err, fmt='none', c='white', alpha=0.3, errorevery=3)
        ax.quiver(xplot[:-1], yplot[:-1], (xplot[1:] - xplot[:-1]) / 2, (yplot[1:] - yplot[:-1]) / 2, angles='uv',
                  scale_units=None, scale=None, alpha=0.6, color=branch_colors[0])
        scats = ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                           vmax=vrange[burstname]['f'][1], alpha=1.0)
        plt.colorbar(scats, ax=ax, label='Frequency')

    fig2.savefig('AIA.png', format='png', dpi=300)
    fig3.savefig('dspec.png', format='png', dpi=200)
    # fig1.close()
    # fig2.close()
    # fig3.close()


if mkmc:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate

    trange = [Time('2014-11-01T16:38:00').jd, Time('2014-11-01T16:52:00').jd]
    SDOdir = os.getenv('SUNCASADB') + 'aiaBrowserData/Download/'
    sdofile = DButil.readsdofile(datadir=SDOdir, wavelength='211', jdtime=trange)
    # sdofile = glob.glob('/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/database/aiaBrowserData/Download/2014/11/01/hmi.M_45s.20141101_*_TAI.2.magnetogram.fits')
    # sdofile = sorted(sdofile)
    # x0, x1, y0, y1 = -660., -575., -200., -115.
    x0, x1, y0, y1 = -690., -510., -270., -90.  # large
    maplist = []
    for ll in tqdm(sdofile):
        maptmp = sunpy.map.Map(ll)
        if maptmp.detector == 'HMI':
            pass
        else:
            maptmp = DButil.normalize_aiamap(maptmp)
        submaptmp = maptmp.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                  u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        maplist.append(submaptmp)
        # print '{}/{} loaded.'.format(idx + 1, len(sdofile))
    mapcube = mapcube_solar_derotate(sunpy.map.Map(maplist, cube=True))
    with open('{}/img_centroid_new/mapcube_{}_large'.format(datadir, mapcube[0].meta['wavelnth']), 'wb') as sf:
        pickle.dump(mapcube, sf)
else:
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube = pickle.load(sf)
        # submaplist = []
        # for sidx, smap in enumerate(mapcube0):
        #     submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
        #                             u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        #     submaplist.append(submaptmp)
        # mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))

# make an AIA movie
if mkmovie:
    burstnames = ['', 'U02', 'U04']
    t_bursts = Time(['2014-11-01T16:43:53', '2014-11-01T16:45:43', '2014-11-01T16:46:19'])
    mapcubename = 'mapcube_hmi_45s'
    from copy import deepcopy

    with open('{}/img_centroid_new/{}'.format(datadir, mapcubename), 'rb') as sf:
        mapcube = pickle.load(sf)

    outdir = '{0}img_centroid_new/img_{1}'.format(datadir, mapcubename)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        os.system('rm -rf {}/*.png'.format(outdir))
    fig = plt.figure(1, figsize=(7, 5))
    plt.ioff()
    for smap in tqdm(mapcube):
        time_tol = 6.
        if mapcube[0].meta['wavelnth'] in [1600, 1700]:
            time_tol = 12.
        if mapcubename == 'mapcube_hmi_45s':
            ax = plot_map(smap, vmax=500, vmin=-500, diff=True)
            t_map = Time(smap.meta['t_obs'])
            time_tol = 22.5
        else:
            ax = plot_map(smap)
            t_map = Time(smap.meta['t_obs'])
        tdiff = np.abs(t_map.jd - t_bursts.jd)
        if np.min(tdiff) < time_tol / 24 / 3600:
            bidx = np.argmin(tdiff)
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstnames[bidx]))
            # cents = np.load('centroidsLL.npy')
            x = cents.item()['x']
            y = cents.item()['y']
            # f = cents.item()['freq']
            # t = cents.item()['time']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstnames[bidx]))
            if burstnames[bidx] == 'U04':
                docorrect = 0
            if docorrect:
                bandedge = [63, 127, 191]
                nf = bandedge[2] - bandedge[1] + 1
                xcorr = np.tile(np.linspace(0, -4, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                # xcorr = np.tile(-(np.exp(np.log(4+1) * np.linspace(0, 1, nf)) - 1), ntime).reshape(ntime, nf).swapaxes(0, 1)
                x[bandedge[1]:bandedge[2] + 1, :] += xcorr
                ycorr = np.tile(np.linspace(0, -4.5, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                y[bandedge[1]:bandedge[2] + 1, :] += ycorr
            if mask:
                if mask == 'threshold':
                    mp = ma.masked_outside(p, thrshd, dmax)
                elif mask == 'native':
                    mp = p
                mx = ma.array(x, mask=mp.mask)
                my = ma.array(y, mask=mp.mask)
                # mp = ma.array(p, mask=mp.mask)
                mf = ma.array(ff, mask=mp.mask)
            else:
                # mlogp = logp
                mx = x
                my = y
                mp = p
                mf = f
            if dosmooth:
                # first smooth along the frequency axis
                for i in range(len(tdatetime)):
                    # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                    x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                    y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
            for sidx, statusfile in enumerate(statusfiles):
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)
                xs, ys = lns[sidx]['x'], lns[sidx]['y']
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                s = ma.compressed(s)
                c = ma.compressed(c)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                xplot = xplot  # [::3]
                yplot = yplot  # [::3]
                s = s  # [::3]
                c = c  # [::3]

                # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                #           yplot[1:] - yplot[:-1], angles='uv', scale=None, alpha=0.6,
                #           color=color[sidx])
                ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                           vmax=vrange[burstname]['f'][1], alpha=1.0)
        fig.savefig('{0}/{3}{1}-{2}.png'.format(outdir, smap.meta['wavelnth'],
                                                t_map.iso.replace(' ', 'T').replace(':', '').replace('-', '')[:-4],
                                                smap.detector), format='png', dpi=100)
    plt.ion()
    # fig.close()

if rdiff_ratio:
    d_frm = 3
    domedfilt = 0
    burstnames = ['', 'U02', 'U04']
    t_bursts = Time(['2014-11-01T16:43:53', '2014-11-01T16:45:43', '2014-11-01T16:46:19'])
    mapcubename = 'mapcube_211_large'
    from copy import deepcopy

    with open('{}/img_centroid_new/{}'.format(datadir, mapcubename), 'rb') as sf:
        mapcube = pickle.load(sf)
    from scipy.signal import medfilt

    # figure()
    # plt.subplot(121)
    # plt.imshow(mapcube[0].data)
    # plt.subplot(122)
    # plt.imshow(medfilt(mapcube[0].data, 3))
    maps_rdiff = []
    for idx, ll in enumerate(mapcube):
        maps_rdiff.append(deepcopy(ll))
        if idx - d_frm < 0:
            sidx = 0
        else:
            sidx = idx - d_frm
        if domedfilt:
            maps_rdiff[idx].data = medfilt(maps_rdiff[idx].data, 5) - medfilt(mapcube[sidx].data, 5)
        else:
            maps_rdiff[idx].data = maps_rdiff[idx].data - mapcube[sidx].data
    mapcube_rdiff = sunpy.map.Map(maps_rdiff, cube=True)
    outdir = '{0}img_centroid_new/img_{1}_diff'.format(datadir, mapcubename)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        os.system('rm -rf {}/*.png'.format(outdir))
    # fig_cutslit = plt.figure(5, figsize=(7, 5))
    # ax = plot_map(mapcube_rdiff[0], vmax=500, vmin=-500, diff=True)

    fig = plt.figure(1, figsize=(7, 5))
    plt.ioff()
    for smap in tqdm(mapcube_rdiff):
        if smap.meta['wavelnth'] == 304:
            ax = plot_map(smap, vmax=500, vmin=-500, diff=True)
        elif smap.meta['wavelnth'] == 193:
            ax = plot_map(smap, vmax=1000, vmin=-1000, diff=True)
        elif smap.meta['wavelnth'] == 171:
            ax = plot_map(smap, vmax=1000, vmin=-1000, diff=True)
        elif smap.meta['wavelnth'] == 335:
            ax = plot_map(smap, vmax=50, vmin=-50, diff=True)
        elif smap.meta['wavelnth'] == 211:
            ax = plot_map(smap, vmax=500, vmin=-500, diff=True, SymLogNorm=True)
        else:
            ax = plot_map(smap, vmax=1000, vmin=-1000, diff=True)
        t_map = Time(smap.meta['t_obs'])
        tdiff = np.abs(t_map.jd - t_bursts.jd)
        time_tol = 6.5
        if mapcube[0].meta['wavelnth'] in [1600, 1700]:
            time_tol = 12.
        if np.min(tdiff) < time_tol / 24 / 3600:
            bidx = np.argmin(tdiff)
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstnames[bidx]))
            # cents = np.load('centroidsLL.npy')
            x = cents.item()['x']
            y = cents.item()['y']
            # f = cents.item()['freq']
            # t = cents.item()['time']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstnames[bidx]))
            if burstnames[bidx] == 'U04':
                docorrect = 0
            if docorrect:
                bandedge = [63, 127, 191]
                nf = bandedge[2] - bandedge[1] + 1
                xcorr = np.tile(np.linspace(0, -4, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                # xcorr = np.tile(-(np.exp(np.log(4+1) * np.linspace(0, 1, nf)) - 1), ntime).reshape(ntime, nf).swapaxes(0, 1)
                x[bandedge[1]:bandedge[2] + 1, :] += xcorr
                ycorr = np.tile(np.linspace(0, -4.5, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                y[bandedge[1]:bandedge[2] + 1, :] += ycorr
            if mask:
                if mask == 'threshold':
                    mp = ma.masked_outside(p, thrshd, dmax)
                elif mask == 'native':
                    mp = p
                mx = ma.array(x, mask=mp.mask)
                my = ma.array(y, mask=mp.mask)
                # mp = ma.array(p, mask=mp.mask)
                mf = ma.array(ff, mask=mp.mask)
            else:
                # mlogp = logp
                mx = x
                my = y
                mp = p
                mf = f
            if dosmooth:
                # first smooth along the frequency axis
                for i in range(len(tdatetime)):
                    # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                    x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                    y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
            for sidx, statusfile in enumerate(statusfiles):
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)
                xs, ys = lns[sidx]['x'], lns[sidx]['y']
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                s = ma.compressed(s)
                c = ma.compressed(c)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                xplot = xplot  # [::3]
                yplot = yplot  # [::3]
                s = s  # [::3]
                c = c  # [::3]

                # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                #           yplot[1:] - yplot[:-1], angles='uv', scale=None, alpha=0.6,
                #           color=color[sidx])
                ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                           vmax=vrange[burstname]['f'][1], alpha=1.0)
        fig.savefig('{0}/{3}{1}-{2}.png'.format(outdir, smap.meta['wavelnth'],
                                                t_map.iso.replace(' ', 'T').replace(':', '').replace('-', '')[:-4],
                                                smap.detector), format='png', dpi=100)
        # time.sleep(0.01)
    plt.ion()


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
                # if len(pointDFtmp.index) <= 3:
                cutslitplt = FitSlit(xx, yy, self.cutwidth, self.cutang, self.cutlength, method='Polyfit')
            self.cutslitplt = cutslitplt
            self.slitline.set_data(cutslitplt['xcen'], cutslitplt['ycen'])
            self.slitline0.set_data(cutslitplt['xs0'], cutslitplt['ys0'])
            self.slitline1.set_data(cutslitplt['xs1'], cutslitplt['ys1'])
            self.slitline.figure.canvas.draw()
            self.slitline0.figure.canvas.draw()
            self.slitline1.figure.canvas.draw()


    d_frm = 3
    domedfilt = 0
    burstnames = ['U01', 'U02', 'U04']
    t_bursts = Time(['2014-11-01T16:43:53', '2014-11-01T16:45:43', '2014-11-01T16:46:19'])
    mapcubename = 'mapcube_171'
    from copy import deepcopy

    with open('{}/img_centroid_new/{}'.format(datadir, mapcubename), 'rb') as sf:
        mapcube = pickle.load(sf)
    from scipy.signal import medfilt

    maps_rdiff = []
    tplot = []
    for idx, smap in enumerate(tqdm(mapcube)):
        maps_rdiff.append(deepcopy(smap))
        tplot.append(smap.meta['t_obs'])
        if idx - d_frm < 0:
            sidx = 0
        else:
            sidx = idx - d_frm
        if domedfilt:
            maps_rdiff[idx].data = medfilt(maps_rdiff[idx].data, 5) - medfilt(mapcube[sidx].data, 5)
        else:
            maps_rdiff[idx].data = maps_rdiff[idx].data - mapcube[sidx].data
    tplot = Time(tplot)
    mapcube_rdiff = sunpy.map.Map(maps_rdiff, cube=True)
    fig_cutslit = plt.figure(5, figsize=(7, 5))
    ax = plot_map(mapcube_rdiff[39], vmax=500, vmin=-500, diff=True)
    cutslitbuilder = CutslitBuilder(ax, cutwidth=1.8, cutang=0, cutlength=80)

    to_save = 0
    to_restore = 0
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
        outdir = '{0}img_centroid_new/img_{1}_diff_stackplt'.format(datadir, mapcubename)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        else:
            os.system('rm -rf {}/*.png'.format(outdir))
        stackplt = []
        for idx, smap in enumerate(tqdm(mapcube_rdiff)):
            intens = getimprofile(smap.data, cutslitbuilder.cutslitplt, xrange=smap.xrange.value,
                                  yrange=smap.yrange.value)
            stackplt.append(intens['y'])
        if len(stackplt) > 1:
            stackplt = np.vstack(stackplt)
            stackplt = stackplt.transpose()
        if mapcube_rdiff[0].meta['wavelnth'] == 304:
            vmax, vmin = 500, -500
        elif mapcube_rdiff[0].meta['wavelnth'] == 193:
            vmax, vmin = 1000, -1000
        elif mapcube_rdiff[0].meta['wavelnth'] == 171:
            vmax, vmin = 1000, -1000
        elif mapcube_rdiff[0].meta['wavelnth'] == 335:
            vmax, vmin = 50, -50
        else:
            vmax, vmin = 1000, -1000
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cutslitplt = cutslitbuilder.cutslitplt
        dspec = {'dspec': stackplt, 'x': tplot.plot_date, 'y': cutslitplt['dist'], 'ytitle': 'Distance [arcsec]',
                 'ctitle': 'DN counts per second',
                 'args': {'norm': norm, 'cmap': cm.get_cmap('sdoaia{}'.format(mapcube_rdiff[0].meta['wavelnth']))}}
        fig_stackplt = plt.figure(6, figsize=(11, 5))
        dtplot = np.mean(np.diff(tplot.plot_date))
        for xidx, smap in enumerate(tqdm(mapcube_rdiff)):
            dspec['axvspan'] = [tplot[xidx].plot_date, tplot[xidx].plot_date + dtplot]
            ax, ax2 = plot_map(smap, dspec, vmax=vmax, vmin=vmin, diff=True)
            ax.plot(cutslitplt['xcen'], cutslitplt['ycen'], color='white', ls='solid')
            ax.plot(cutslitplt['xs0'], cutslitplt['ys0'], color='white', ls='dotted')
            ax.plot(cutslitplt['xs1'], cutslitplt['ys1'], color='white', ls='dotted')
            time_tol = 6.
            t_map = Time(smap.meta['t_obs'])
            if mapcube[0].meta['wavelnth'] in [1600, 1700]:
                time_tol = 12.
            if mapcubename == 'mapcube_hmi_45s':
                time_tol = 22.5
            tdiff = np.abs(t_map.jd - t_bursts.jd)
            if np.min(tdiff) < time_tol / 24 / 3600:
                bidx = np.argmin(tdiff)
                cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstnames[bidx]))
                # cents = np.load('centroidsLL.npy')
                x = cents.item()['x']
                y = cents.item()['y']
                # f = cents.item()['freq']
                # t = cents.item()['time']
                p = cents.item()['peak']
                f = cents.item()['freq']
                fitsarr = cents.item()['fitsarr']
                ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
                t = Time(cents.item()['time'] / 3600. / 24., format='jd')
                # tdatetime = mdates.date2num(t.datetime)
                tdatetime = t.plot_date
                nfreq, ntime = p.shape
                statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstnames[bidx]))
                if burstnames[bidx] == 'U04':
                    docorrect = 0
                if docorrect:
                    bandedge = [63, 127, 191]
                    nf = bandedge[2] - bandedge[1] + 1
                    xcorr = np.tile(np.linspace(0, -4, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                    # xcorr = np.tile(-(np.exp(np.log(4+1) * np.linspace(0, 1, nf)) - 1), ntime).reshape(ntime, nf).swapaxes(0, 1)
                    x[bandedge[1]:bandedge[2] + 1, :] += xcorr
                    ycorr = np.tile(np.linspace(0, -4.5, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                    y[bandedge[1]:bandedge[2] + 1, :] += ycorr
                if mask:
                    if mask == 'threshold':
                        mp = ma.masked_outside(p, thrshd, dmax)
                    elif mask == 'native':
                        mp = p
                    mx = ma.array(x, mask=mp.mask)
                    my = ma.array(y, mask=mp.mask)
                    # mp = ma.array(p, mask=mp.mask)
                    mf = ma.array(ff, mask=mp.mask)
                else:
                    # mlogp = logp
                    mx = x
                    my = y
                    mp = p
                    mf = f
                if dosmooth:
                    # first smooth along the frequency axis
                    for i in range(len(tdatetime)):
                        # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                        x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                        y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
                for sidx, statusfile in enumerate(statusfiles):
                    with open(statusfile, 'rb') as sf:
                        sct, lns = pickle.load(sf)
                    xs, ys = lns[sidx]['x'], lns[sidx]['y']
                    xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                    yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                    xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                    yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                    s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                    c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                    xinterp = ma.masked_invalid(xinterp)
                    yinterp = ma.array(yinterp, mask=xinterp.mask)
                    s = ma.array(s, mask=xinterp.mask)
                    c = ma.array(c, mask=xinterp.mask)
                    xinterp = ma.compressed(xinterp)
                    yinterp = ma.compressed(yinterp)
                    s = ma.compressed(s)
                    c = ma.compressed(c)
                    xplot = smooth(xinterp, box_pts)
                    yplot = smooth(yinterp, box_pts)
                    xplot = xplot  # [::3]
                    yplot = yplot  # [::3]
                    s = s  # [::3]
                    c = c  # [::3]

                    # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                    #           yplot[1:] - yplot[:-1], angles='uv', scale=None, alpha=0.6,
                    #           color=color[sidx])
                    ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                               vmax=vrange[burstname]['f'][1], alpha=1.0)
            fig_stackplt.savefig('{0}/{3}{1}-{2}.png'.format(outdir, smap.meta['wavelnth'],
                                                             tplot[xidx].iso.replace(' ', 'T').replace(':', '').replace(
                                                                 '-', '')[:-4], smap.detector), format='png', dpi=100)

# do a movie of burst motion
if mkmovie2:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    import matplotlib.animation as animation

    fig = plt.figure(1, figsize=(12, 5))
    # fig2 = plt.figure(2)
    plt.ioff()
    burstnames = ['U01', 'U02', 'U04']
    burstnames = ['', 'U02', 'U04']
    datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube0 = pickle.load(sf)
    submaplist = []

    t_maps = [ll['t_obs'] for ll in mapcube0.all_meta()]
    t_maps = Time(t_maps)
    for bidx, burstname in enumerate(burstnames):
        if burstname:
            if burstname in ['U02', 'U04']:
                x0, x1, y0, y1 = -647.5, -630, -135, -117.5
            else:
                x0, x1, y0, y1 = -652.5, -622.5, -147.5, -117.5
            submaplist = []
            for sidx, smap in enumerate(mapcube0):
                submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                        u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
                submaplist.append(submaptmp)
            mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstname))
            x = cents.item()['x']
            y = cents.item()['y']
            # f = cents.item()['freq']
            # t = cents.item()['time']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstname))
            if burstname == 'U04':
                docorrect = 0
            if docorrect:
                bandedge = [63, 127, 191]
                nf = bandedge[2] - bandedge[1] + 1
                xcorr = np.tile(np.linspace(0, -4, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                # xcorr = np.tile(-(np.exp(np.log(4+1) * np.linspace(0, 1, nf)) - 1), ntime).reshape(ntime, nf).swapaxes(0, 1)
                x[bandedge[1]:bandedge[2] + 1, :] += xcorr
                ycorr = np.tile(np.linspace(0, -4.5, nf), ntime).reshape(ntime, nf).swapaxes(0, 1)
                y[bandedge[1]:bandedge[2] + 1, :] += ycorr
            if mask:
                if mask == 'threshold':
                    mp = ma.masked_outside(p, thrshd, dmax)
                elif mask == 'native':
                    mp = p
                mx = ma.array(x, mask=mp.mask)
                my = ma.array(y, mask=mp.mask)
                # mp = ma.array(p, mask=mp.mask)
                mf = ma.array(ff, mask=mp.mask)
            else:
                # mlogp = logp
                mx = x
                my = y
                mp = p
                mf = f
            if dosmooth:
                # first smooth along the frequency axis
                for i in range(len(tdatetime)):
                    # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                    x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                    y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
            dspec = {'dspec': p, 'x': tdatetime, 'y': f, 'xlim': vrange[burstname]['xlim'],
                     'ylim': vrange[burstname]['ylim'], 'ytitle': 'Frequency [GHz]', 'ctitle': 'Jy/Beam',
                     'args': {'vmin': vrange[burstname]['spec'][0], 'vmax': vrange[burstname]['spec'][1],
                              'cmap': 'jet'}}
            for sidx, statusfile in enumerate(statusfiles):
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)

                xs, ys = lns[sidx]['x'], lns[sidx]['y']
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                finterp = f[np.floor(yspix, out=None).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                s0 = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                r0 = r[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                finterp = ma.array(finterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                s0 = ma.array(s0, mask=xinterp.mask)
                r0 = ma.array(r0, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xs = ma.array(xs, mask=xinterp.mask)
                ys = ma.array(ys, mask=xinterp.mask)
                xspix = ma.array(xspix, mask=xinterp.mask)
                yspix = ma.array(yspix, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                finterp = ma.compressed(finterp)
                s = ma.compressed(s)
                s0 = ma.compressed(s0)
                r0 = ma.compressed(r0)
                c = ma.compressed(c)
                xs = ma.compressed(xs)
                ys = ma.compressed(ys)
                xspix = ma.compressed(xspix)
                yspix = ma.compressed(yspix)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                fplot = smooth(finterp, box_pts)
                xplot_dif = np.diff(xplot)
                yplot_dif = np.diff(yplot)

                # xplot = xplot#[::3]
                # yplot = yplot#[::3]
                # s = s#[::3]
                # c = c#[::3]
                dspec['xs'] = xs
                dspec['ys'] = ys
                err = (30. * r0 / s0)
                tplts = Time(xs, format='plot_date')
                tplts_dif = np.diff(tplts.jd * 24 * 3600)
                vxplot, vyplot = xplot_dif / tplts_dif, yplot_dif / tplts_dif
                vxplot = np.append(vxplot, 0)
                vyplot = np.append(vyplot, 0)
                outdir = '{0}{1}-centroid/{2}/'.format(datadir, burstname, os.path.basename(statusfile)[:7])
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                # os.system('rm -rf {}*.png'.format(outdir))
                # for xidx, tplt in enumerate(tqdm(tplts)):
                #     plt.figure(1)
                #     fig.clf()
                #     smap = mapcube[np.argmin(np.abs(tplt.jd - t_maps.jd))]
                #     if xidx == 0:
                #         dspec['axvspan'] = [tplt.plot_date, (tplt.plot_date + tplts[xidx + 1].plot_date) / 2]
                #     elif xidx == len(tplts.plot_date) - 1:
                #         dspec['axvspan'] = [(tplt.plot_date + tplts[xidx - 1].plot_date) / 2, tplt.plot_date]
                #     else:
                #         dspec['axvspan'] = [(tplt.plot_date + tplts[xidx - 1].plot_date) / 2,
                #                             (tplt.plot_date + tplts[xidx + 1].plot_date) / 2]
                #     ax, ax2 = plot_map(smap, dspec, cmap='gray_r')
                #
                #     ax.errorbar(xplot[xidx], yplot[xidx], xerr=err[xidx], yerr=err[xidx], fmt='none', c='white',
                #                 alpha=0.3)
                #     if xidx < len(tplts) - 1:
                #         ax.quiver(xplot[:-1][xidx], yplot[:-1][xidx], xplot[1:][xidx] - xplot[:-1][xidx],
                #                   yplot[1:][xidx] - yplot[:-1][xidx], angles='uv', scale_units=None, scale=None,
                #                   alpha=0.6, color=branch_colors[sidx])
                #     ax.scatter(xplot[xidx], yplot[xidx], c=c[xidx], cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                #                vmax=vrange[burstname]['f'][1], alpha=1.0)
                #     scalebartext = '3,000 km'
                #     scalelength = 3.e6 * u.m * smap.rsun_obs / smap.rsun_meters
                #     ct.add_scalebar(scalelength.value, 0.95, 0.10, ax, label=scalebartext, yoff=-0.008, color='white',
                #                     alpha=1.0)
                #     tstring = tplt.iso.replace(' ', 'T').replace(':', '').replace('-', '')
                #     fig.savefig('{0}{1}-{2}.png'.format(outdir, os.path.basename(statusfile)[:7], tstring),
                #                 format='png', dpi=150)
                # DButil.img2movie('{0}{1}'.format(outdir, os.path.basename(statusfile)[:7]),
                #                  outname='{0}{1}-movie'.format(outdir, os.path.basename(statusfile)[:7]),
                #                  overwrite=True)
                fig3 = plt.figure(3, figsize=(6, 6))
                fig3.clf()
                circle1 = plt.Circle((0, 0), 1.5, fill=False)
                plt.quiver(vxplot[:-1], vyplot[:-1], vxplot[1:] - vxplot[:-1], vyplot[1:] - vyplot[:-1], angles='uv',
                           scale_units=None, scale=None, alpha=0.3, color='black')
                ax = plt.gca()
                ax.add_artist(circle1)
                ax.scatter(vxplot, vyplot, alpha=1.0, c=xs, cmap=dspec_cmaps[0])
                ax.set_xlim(-15, 15)
                ax.set_ylim(-15, 15)
                ax.set_xlabel('Vx [arcsec/s]')
                ax.set_ylabel('Vy [arcsec/s]')
                fig3.savefig('{0}velocity.png'.format(outdir), format='png', dpi=300)

                fig2 = plt.figure(2, figsize=(4, 6))
                fig2.clf()
                tplt = tplts[len(tplts) / 2]
                smap = mapcube[np.argmin(np.abs(tplt.jd - t_maps.jd))]
                ax = plot_map(smap, cmap='gray_r', smap_colorbar=False)
                ax.errorbar(xplot, yplot, xerr=err, yerr=err, fmt='none', c='white', alpha=0.3, errorevery=3)
                ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1], yplot[1:] - yplot[:-1], angles='uv',
                          scale_units=None, scale=None, alpha=0.6, color=branch_colors[sidx])
                ax.scatter(xplot, yplot, c=c, cmap=dspec_cmaps[0], vmin=vrange[burstname]['f'][0],
                           vmax=vrange[burstname]['f'][1], alpha=1.0)
                ax.set_xlim(-641.78044977746868, -635.17863477746869)
                ax.set_ylim(-131., -119.)
                ax.set_axis_off()
                ax.set_title = ''
                scalebartext = '1,500 km'
                scalelength = 1.5e6 * u.m * smap.rsun_obs / smap.rsun_meters
                ct.add_scalebar(scalelength.value, 0.95, 0.10, ax, label=scalebartext, yoff=-0.008, color='white',
                                alpha=1.0)
                fig2.savefig('{0}overview.png'.format(outdir), format='png', dpi=150)
    plt.ion()
    # fig.close()
    # fig2.close()

# do a movie of 3D trajectories of the burst
if mkmovie3:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    fig2 = plt.figure(2, figsize=(9, 6))
    fig2.clf()
    ax = fig2.add_subplot(111, projection='3d')
    burstnames = ['U01', 'U02', 'U04']
    burstnames = ['', 'U02', '']
    datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube0 = pickle.load(sf)
    submaplist = []

    t_maps = [ll['t_obs'] for ll in mapcube0.all_meta()]
    t_maps = Time(t_maps)
    for bidx, burstname in enumerate(burstnames):
        if burstname:
            if burstname in ['U02', 'U04']:
                if burstname == 'U04':
                    tmean = Time('2014-11-01 16:46:18.900')
                    tspan = 2.5 / 24 / 3600
                else:
                    tmean = Time('2014-11-01 16:45:44')
                    tspan = 5. / 24 / 3600
                x0, x1, y0, y1 = -647.5, -630, -135, -117.5
            else:
                x0, x1, y0, y1 = -652.5, -622.5, -147.5, -117.5
            tmin = Time(tmean.jd - tspan, format='jd')
            tmax = Time(tmean.jd + tspan, format='jd')
            submaplist = []
            for sidx, smap in enumerate(mapcube0):
                submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                        u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
                submaplist.append(submaptmp)
            mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstname))
            x = cents.item()['x']
            y = cents.item()['y']
            # f = cents.item()['freq']
            # t = cents.item()['time']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstname))

            if mask:
                if mask == 'threshold':
                    mp = ma.masked_outside(p, thrshd, dmax)
                elif mask == 'native':
                    mp = p
                mx = ma.array(x, mask=mp.mask)
                my = ma.array(y, mask=mp.mask)
                # mp = ma.array(p, mask=mp.mask)
                mf = ma.array(ff, mask=mp.mask)
            else:
                # mlogp = logp
                mx = x
                my = y
                mp = p
                mf = f
            if dosmooth:
                # first smooth along the frequency axis
                for i in range(len(tdatetime)):
                    # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                    x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                    y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
            for sidx, statusfile in enumerate(statusfiles):
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)

                xs, ys = lns[sidx]['x'], lns[sidx]['y']
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                finterp = f[np.floor(yspix, out=None).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                s0 = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                r0 = r[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                # c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                c = xs
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                finterp = ma.array(finterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                s0 = ma.array(s0, mask=xinterp.mask)
                r0 = ma.array(r0, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xs = ma.array(xs, mask=xinterp.mask)
                ys = ma.array(ys, mask=xinterp.mask)
                xspix = ma.array(xspix, mask=xinterp.mask)
                yspix = ma.array(yspix, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                finterp = ma.compressed(finterp)
                s = ma.compressed(s)
                s0 = ma.compressed(s0)
                r0 = ma.compressed(r0)
                c = ma.compressed(c)
                xs = ma.compressed(xs)
                ys = ma.compressed(ys)
                xspix = ma.compressed(xspix)
                yspix = ma.compressed(yspix)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                fplot = smooth(finterp, box_pts)
                neplot = (fplot * 1e9 / 8.98e3) ** 2 / 1e10
                xplot_dif = np.diff(xplot)
                yplot_dif = np.diff(yplot)

                err = (30. * r0 / s0)
                tplts = Time(xs, format='plot_date')
                tplts_dif = np.diff(tplts.jd * 24 * 3600)
                vxplot, vyplot = xplot_dif / tplts_dif, yplot_dif / tplts_dif
                vxplot = np.append(vxplot, 0)
                vyplot = np.append(vyplot, 0)
                outdir = '{0}{1}-centroid/{2}/'.format(datadir, burstname, os.path.basename(statusfile)[:7])
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

                tplt = tplts[len(tplts) / 2]
                # ax = plot_map(smap)
                # ax.errorbar(xplot, yplot, xerr=err, yerr=err, fmt='none', c='white', alpha=0.3, errorevery=3)
                # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                #           yplot[1:] - yplot[:-1], angles='uv', scale_units=None, scale=None, alpha=0.6,
                #           color=branch_colors[sidx])
                cax = ax.scatter(xplot, yplot, neplot, c=c, cmap='jet', vmin=tmin.plot_date, vmax=tmax.plot_date,
                                 alpha=1.0, s=80)
            smap = mapcube[np.argmin(np.abs(tplt.jd - t_maps.jd))]
            nx, ny = smap.dimensions
            x_smap = np.linspace(smap.xrange[0].value, smap.xrange[1].value, nx.value)
            y_smap = np.linspace(smap.yrange[0].value, smap.yrange[1].value, ny.value)
            xx, yy = np.meshgrid(x_smap, y_smap)
            clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
            vmax = clrange['high']
            vmin = clrange['low']
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth']))
            imshow_args = {'cmap': cmap}
            ax.view_init(elev=195., azim=60)
            if ax.zaxis_inverted():
                ax.invert_zaxis()
            # ax.plot_surface(xx, yy, np.ones(smap.data.shape) * 2.4, rstride=1, cstride=1,
            #                 facecolors=cmap(smap.data.astype(float) / smap.data.max()), norm=norm, shade=False)
            ax.contourf(xx, yy, smap.data, 256, zdir='z', offset=2.4, zorder=1, **imshow_args)
            cbar = plt.colorbar(cax, ax=ax, label='Time')
            cbar.set_ticks(
                [cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin) for t in cbar.ax.get_yticks()])
            cbar_labels = map(lambda x: x.strftime('%H:%M:%S'), Time(
                [(cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin)) for t in cbar.ax.get_yticks()],
                format='plot_date').to_datetime())
            cbar.set_ticklabels(cbar_labels)
            ax.set_zlim(1.2, 2.4)
            ax.set_xlabel('Solar X [{xpos}]'.format(xpos=smap.spatial_units.x))
            ax.set_ylabel('Solar Y [{ypos}]'.format(ypos=smap.spatial_units.y))
            ax.set_zlabel(r'Density [$\times 10^{10} cm^{-3}$]')
            # fig2.savefig('{0}overview.png'.format(outdir), format='png', dpi=150)
            # fig.close()
            # fig2.close()

            nframe = 361
            az, el = np.linspace(-300, 60, nframe), [195] * nframe


            def update_frame(num, az, el):
                ax.view_init(elev=el[num], azim=az[num])
                ax.set_zlim(1.2, 2.4)
                return


            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
            burst_ani = animation.FuncAnimation(fig2, update_frame, nframe, fargs=(az, el), interval=50, blit=False)
            outdir = '{0}{1}-centroid/'.format(datadir, burstname)
            burst_ani.save(outdir + 'burst.mp4', writer=writer)

if dopltmass:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    fig2 = plt.figure(2, figsize=(9, 6))
    fig2.clf()
    ax = fig2.add_subplot(111, projection='3d')
    burstnames = ['', '', 'U04']
    datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube0 = pickle.load(sf)
    submaplist = []

    t_maps = [ll['t_obs'] for ll in mapcube0.all_meta()]
    t_maps = Time(t_maps)
    for bidx, burstname in enumerate(burstnames):
        if burstname:
            if burstname in ['U02', 'U04']:
                if burstname == 'U04':
                    tmean = Time('2014-11-01 16:46:18.900')
                    tspan = 2.5 / 24 / 3600
                else:
                    tmean = Time('2014-11-01 16:45:44')
                    tspan = 5. / 24 / 3600
                x0, x1, y0, y1 = -647.5, -630, -135, -117.5
            else:
                x0, x1, y0, y1 = -652.5, -622.5, -147.5, -117.5
            tmin = Time(tmean.jd - tspan, format='jd')
            tmax = Time(tmean.jd + tspan, format='jd')
            submaplist = []
            for sidx, smap in enumerate(mapcube0):
                submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                        u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
                submaplist.append(submaptmp)
            mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstname))
            x = cents.item()['x']
            y = cents.item()['y']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstname))

            thrshd = 40
            mp = ma.masked_less(p, thrshd)
            mx = ma.array(x, mask=mp.mask)
            my = ma.array(y, mask=mp.mask)
            mf = ma.array(ff, mask=mp.mask)
            mt = ma.array(np.vstack([t.jd] * nfreq), mask=mp.mask)
            neplot = (mf * 1e9 / 8.98e3) ** 2 / 1e10
            c = Time(mt, format='jd').plot_date
            xplot = mx
            yplot = my
            tplt = np.mean(mt)
            cax = ax.scatter(xplot, yplot, neplot, c=c, cmap='jet', vmin=tmin.plot_date, vmax=tmax.plot_date, alpha=1.0,
                             s=80)
            smap = mapcube[np.argmin(np.abs(tplt - t_maps.jd))]
            nx, ny = smap.dimensions
            x_smap = np.linspace(smap.xrange[0].value, smap.xrange[1].value, nx.value)
            y_smap = np.linspace(smap.yrange[0].value, smap.yrange[1].value, ny.value)
            xx, yy = np.meshgrid(x_smap, y_smap)
            clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
            vmax = clrange['high']
            vmin = clrange['low']
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth']))
            imshow_args = {'cmap': cmap}
            ax.view_init(elev=195., azim=60)
            if ax.zaxis_inverted():
                ax.invert_zaxis()
            ax.contourf(xx, yy, smap.data, 256, zdir='z', offset=2.4, zorder=1, **imshow_args)
            cbar = plt.colorbar(cax, ax=ax, label='Time')
            cbar.set_ticks(
                [cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin) for t in cbar.ax.get_yticks()])
            cbar_labels = map(lambda x: x.strftime('%H:%M:%S'), Time(
                [(cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin)) for t in cbar.ax.get_yticks()],
                format='plot_date').to_datetime())
            cbar.set_ticklabels(cbar_labels)
            ax.set_zlim(1.2, 2.4)
            ax.set_xlabel('Solar X [{xpos}]'.format(xpos=smap.spatial_units.x))
            ax.set_ylabel('Solar Y [{ypos}]'.format(ypos=smap.spatial_units.y))
            ax.set_zlabel(r'Density [$\times 10^{10} cm^{-3}$]')
            # fig2.savefig('{0}overview.png'.format(outdir), format='png', dpi=150)
            # fig.close()
            # fig2.close()

            nframe = 361
            az, el = np.linspace(-300, 60, nframe), [195] * nframe


            def update_frame(num, az, el):
                ax.view_init(elev=el[num], azim=az[num])
                ax.set_zlim(1.2, 2.4)
                return


            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
            burst_ani = animation.FuncAnimation(fig2, update_frame, nframe, fargs=(az, el), interval=50, blit=False)
            outdir = '{0}{1}-centroid/'.format(datadir, burstname)
            burst_ani.save(outdir + 'burst.mp4', writer=writer)


# same as mkmovie4, using mayavi to plot. to be completed
if mkmovie4:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    from mayavi import mlab

    burstnames = ['U01', 'U02', 'U04']
    burstnames = ['', '', 'U04']
    datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube0 = pickle.load(sf)
    submaplist = []

    t_maps = [ll['t_obs'] for ll in mapcube0.all_meta()]
    t_maps = Time(t_maps)
    for bidx, burstname in enumerate(burstnames):
        if burstname:
            if burstname in ['U02', 'U04']:
                if burstname == 'U04':
                    tmean = Time('2014-11-01 16:46:18.900')
                    tspan = 2.5 / 24 / 3600
                else:
                    tmean = Time('2014-11-01 16:45:44')
                    tspan = 5. / 24 / 3600
                x0, x1, y0, y1 = -647.5, -630, -135, -117.5
            else:
                x0, x1, y0, y1 = -652.5, -622.5, -147.5, -117.5
            tmin = Time(tmean.jd - tspan, format='jd')
            tmax = Time(tmean.jd + tspan, format='jd')
            submaplist = []
            for sidx, smap in enumerate(mapcube0):
                submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                        u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
                submaplist.append(submaptmp)
            mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstname))
            x = cents.item()['x']
            y = cents.item()['y']
            # f = cents.item()['freq']
            # t = cents.item()['time']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstname))

            if mask:
                if mask == 'threshold':
                    mp = ma.masked_outside(p, thrshd, dmax)
                elif mask == 'native':
                    mp = p
                mx = ma.array(x, mask=mp.mask)
                my = ma.array(y, mask=mp.mask)
                # mp = ma.array(p, mask=mp.mask)
                mf = ma.array(ff, mask=mp.mask)
            else:
                # mlogp = logp
                mx = x
                my = y
                mp = p
                mf = f
            if dosmooth:
                # first smooth along the frequency axis
                for i in range(len(tdatetime)):
                    # p[:, i] = ss.smooth(p[:, i], window_len=4, window='flat')
                    x[:, i] = ss.smooth(x[:, i], window_len=4, window='flat')
                    y[:, i] = ss.smooth(y[:, i], window_len=4, window='flat')
            for sidx, statusfile in enumerate(statusfiles):
                with open(statusfile, 'rb') as sf:
                    sct, lns = pickle.load(sf)

                xs, ys = lns[sidx]['x'], lns[sidx]['y']
                xspix = [(ll - tdatetime[0]) / (tdatetime[-1] - tdatetime[0]) * (ntime - 1) for ll in xs]
                yspix = [(ll - f[0]) / (f[-1] - f[0]) * (nfreq - 1) for ll in ys]
                xinterp = x[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                yinterp = y[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                finterp = f[np.floor(yspix, out=None).astype(np.int)]
                s = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)] ** 2 / 10.
                s0 = p[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                r0 = r[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                # c = ff[np.floor(yspix, out=None).astype(np.int), np.floor(xspix).astype(np.int)]
                c = xs
                xinterp = ma.masked_invalid(xinterp)
                yinterp = ma.array(yinterp, mask=xinterp.mask)
                finterp = ma.array(finterp, mask=xinterp.mask)
                s = ma.array(s, mask=xinterp.mask)
                s0 = ma.array(s0, mask=xinterp.mask)
                r0 = ma.array(r0, mask=xinterp.mask)
                c = ma.array(c, mask=xinterp.mask)
                xs = ma.array(xs, mask=xinterp.mask)
                ys = ma.array(ys, mask=xinterp.mask)
                xspix = ma.array(xspix, mask=xinterp.mask)
                yspix = ma.array(yspix, mask=xinterp.mask)
                xinterp = ma.compressed(xinterp)
                yinterp = ma.compressed(yinterp)
                finterp = ma.compressed(finterp)
                s = ma.compressed(s)
                s0 = ma.compressed(s0)
                r0 = ma.compressed(r0)
                c = ma.compressed(c)
                xs = ma.compressed(xs)
                ys = ma.compressed(ys)
                xspix = ma.compressed(xspix)
                yspix = ma.compressed(yspix)
                xplot = smooth(xinterp, box_pts)
                yplot = smooth(yinterp, box_pts)
                fplot = smooth(finterp, box_pts)
                neplot = (fplot * 1e9 / 8.98e3) ** 2 / 1e10
                xplot_dif = np.diff(xplot)
                yplot_dif = np.diff(yplot)

                err = (30. * r0 / s0)
                tplts = Time(xs, format='plot_date')
                tplts_dif = np.diff(tplts.jd * 24 * 3600)
                vxplot, vyplot = xplot_dif / tplts_dif, yplot_dif / tplts_dif
                vxplot = np.append(vxplot, 0)
                vyplot = np.append(vyplot, 0)
                outdir = '{0}{1}-centroid/{2}/'.format(datadir, burstname, os.path.basename(statusfile)[:7])
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

                tplt = tplts[len(tplts) / 2]
                # ax = plot_map(smap)
                # ax.errorbar(xplot, yplot, xerr=err, yerr=err, fmt='none', c='white', alpha=0.3, errorevery=3)
                # ax.quiver(xplot[:-1], yplot[:-1], xplot[1:] - xplot[:-1],
                #           yplot[1:] - yplot[:-1], angles='uv', scale_units=None, scale=None, alpha=0.6,
                #           color=branch_colors[sidx])
                pts = mlab.points3d(xplot, yplot, neplot * 10, colormap='jet', vmin=tmin.plot_date, vmax=tmax.plot_date,
                                    scale_factor=0.2)
                # cax = ax.scatter(xplot, yplot, neplot, c=c, cmap='jet',
                #                  vmin=tmin.plot_date,
                #                  vmax=tmax.plot_date, alpha=1.0, s=80)
            smap = mapcube[np.argmin(np.abs(tplt.jd - t_maps.jd))]
            nx, ny = smap.dimensions
            x_smap = np.linspace(smap.xrange[0].value, smap.xrange[1].value, nx.value)
            y_smap = np.linspace(smap.yrange[0].value, smap.yrange[1].value, ny.value)
            xx, yy = np.meshgrid(x_smap, y_smap)
            clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
            vmax = clrange['high']
            vmin = clrange['low']
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth']))
            imshow_args = {'cmap': cmap}
            # ax.view_init(elev=195., azim=60)
            # if ax.zaxis_inverted():
            #     ax.invert_zaxis()
            # ax.plot_surface(xx, yy, np.ones(smap.data.shape) * 2.4, rstride=1, cstride=1,
            #                 facecolors=cmap(smap.data.astype(float) / smap.data.max()), norm=norm, shade=False)
            mlab.imshow(smap.data,
                        extent=[smap.xrange[0].value, smap.xrange[1].value, smap.yrange[0].value, smap.yrange[1].value,
                                2.4 * 10, 2.4 * 10], colormap='gray', vmax=vmax, vmin=vmin)
            # ax.contourf(xx, yy, smap.data, 256, zdir='z', offset=2.4, zorder=1, **imshow_args)
            # cbar = plt.colorbar(cax, ax=ax, label='Time')
            # cbar.set_ticks(
            #     [cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin) for t in cbar.ax.get_yticks()])
            # cbar_labels = map(lambda x: x.strftime('%H:%M:%S'),
            #                   Time([(cax.colorbar.vmin + t * (cax.colorbar.vmax - cax.colorbar.vmin)) for t in
            #                         cbar.ax.get_yticks()], format='plot_date').to_datetime())
            # cbar.set_ticklabels(cbar_labels)
            ax.set_zlim(1.2, 2.4)
            ax.set_xlabel('Solar X [{xpos}]'.format(xpos=smap.spatial_units.x))
            ax.set_ylabel('Solar Y [{ypos}]'.format(ypos=smap.spatial_units.y))
            ax.set_zlabel(r'Density [$\times 10^{10} cm^{-3}$]')
            # fig2.savefig('{0}overview.png'.format(outdir), format='png', dpi=150)
            # fig.close()
            # fig2.close()

            nframe = 361
            az, el = np.linspace(-300, 60, nframe), [195] * nframe


            def update_frame(num, az, el):
                ax.view_init(elev=el[num], azim=az[num])
                ax.set_zlim(1.2, 2.4)
                return


            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
            burst_ani = animation.FuncAnimation(fig2, update_frame, nframe, fargs=(az, el), interval=50, blit=False)
            outdir = '{0}{1}-centroid/'.format(datadir, burstname)
            burst_ani.save(outdir + 'burst.mp4', writer=writer)

if dopltmass_mlab:
    from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    from mayavi import mlab

    burstnames = ['', '', 'U04']
    datadir = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/Uburst/'
    with open('{}/img_centroid_new/mapcube_{}'.format(datadir, '304'), 'rb') as sf:
        mapcube0 = pickle.load(sf)
    submaplist = []

    t_maps = [ll['t_obs'] for ll in mapcube0.all_meta()]
    t_maps = Time(t_maps)
    for bidx, burstname in enumerate(burstnames):
        if burstname:
            if burstname in ['U02', 'U04']:
                if burstname == 'U04':
                    tmean = Time('2014-11-01 16:46:18.900')
                    tspan = 2.5 / 24 / 3600
                    # tmean = Time('2014-11-01 16:46:22')
                    # tspan = 8.0 / 24 / 3600
                else:
                    tmean = Time('2014-11-01 16:45:44')
                    tspan = 5. / 24 / 3600
                x0, x1, y0, y1 = -647.5, -630, -135, -117.5
            else:
                x0, x1, y0, y1 = -652.5, -622.5, -147.5, -117.5
            tmin = Time(tmean.jd - tspan, format='jd')
            tmax = Time(tmean.jd + tspan, format='jd')
            submaplist = []
            for sidx, smap in enumerate(mapcube0):
                submaptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                        u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
                submaplist.append(submaptmp)
            mapcube = mapcube_solar_derotate(sunpy.map.Map(submaplist, cube=True))
            cents = np.load('{}{}-centroid/centroidsLL.npy'.format(datadir, burstname))
            x = cents.item()['x']
            y = cents.item()['y']
            p = cents.item()['peak']
            f = cents.item()['freq']
            fitsarr = cents.item()['fitsarr']
            ff = np.tile(f, p.shape[1]).reshape(p.shape[1], p.shape[0]).swapaxes(0, 1)
            t = Time(cents.item()['time'] / 3600. / 24., format='jd')
            # tdatetime = mdates.date2num(t.datetime)
            tdatetime = t.plot_date
            nfreq, ntime = p.shape
            statusfiles = glob.glob('{}{}-centroid/branch?_status'.format(datadir, burstname))

            thrshd = 32
            # thrshd = 27
            # thrshd = 10
            mp = ma.masked_less_equal(p, thrshd)
            mx = ma.array(x, mask=mp.mask)
            my = ma.array(y, mask=mp.mask)
            mf = ma.array(ff, mask=mp.mask)
            fig4, axs4, lines, scatters, linebuilder = pltspec(tdatetime, f, [mp, mx, my], vrange[burstname])
            mt = ma.array(np.vstack([t.jd] * nfreq), mask=mp.mask)
            neplot = (ma.compressed(mf) * 1e9 / 8.98e3) ** 2 / 1e10
            hgtplot = 30. * np.log(2.5 / neplot) - 5.0
            pplot = ma.compressed(mp) ** 3
            pplot = pplot / np.max(pplot)  # / 2.0
            c = Time(ma.compressed(mt), format='jd').plot_date
            xplot = ma.compressed(mx)
            yplot = ma.compressed(my)
            tplt = np.mean(mt)
            mfig = mlab.figure(size=(900, 700), bgcolor=(0.99, 0.99, 0.99), fgcolor=(0, 0, 0))
            pts = mlab.points3d(xplot, yplot, hgtplot, colormap='jet', vmin=tmin.plot_date, vmax=tmax.plot_date,
                                scale_factor=0.5)
            # extent=[np.min(xplot), np.max(xplot), np.min(yplot), np.max(yplot), 1,20])
            pts.glyph.color_mode = "color_by_scalar"
            pts.mlab_source.dataset.point_data.scalars = c
            pts.glyph.scale_mode = 'scale_by_vector'
            pts.mlab_source.dataset.point_data.vectors = np.transpose(np.tile(pplot, (3, 1)))
            # mlab_clb = mlab.colorbar(pts, title='Time', orientation='vertical')
            smap = mapcube[np.argmin(np.abs(tplt - t_maps.jd))]
            nx, ny = smap.dimensions
            x_smap = np.linspace(smap.xrange[0].value, smap.xrange[1].value, nx.value)
            y_smap = np.linspace(smap.yrange[0].value, smap.yrange[1].value, ny.value)
            xx, yy = np.meshgrid(x_smap, y_smap)
            clrange = DButil.sdo_aia_scale_dict(wavelength=smap.meta['wavelnth'])
            vmax = clrange['high']
            vmin = clrange['low']
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap('sdoaia{}'.format(smap.meta['wavelnth']))
            imshow_args = {'cmap': cmap}
            mlab_im = mlab.imshow(np.transpose(smap.data),
                                  extent=[smap.xrange[0].value, smap.xrange[1].value, smap.yrange[0].value,
                                          smap.yrange[1].value, 0, 0], colormap='gray', vmax=vmax, vmin=vmin)
            mlab_outline = mlab.outline(pts, opacity=0.2,
                                        extent=[smap.xrange[0].value, smap.xrange[1].value, smap.yrange[0].value,
                                                smap.yrange[1].value, 0, np.nanmax(hgtplot)])
            mlab_axes = mlab.axes(mlab_outline)
            mlab_axes.axes.font_factor = 1.
            mlab_xlabel = mlab.xlabel('Solar X [{xpos}]'.format(xpos=smap.spatial_units.x), object=mlab_axes)
            mlab_ylabel = mlab.ylabel('Solar Y [{ypos}]'.format(ypos=smap.spatial_units.y), object=mlab_axes)
            # mlab_zlabel = mlab.zlabel(r'Density [$\times 10^{10} cm^{-3}$]', object=mlab_axes)
            mlab_zlabel = mlab.zlabel('Pseudo-height', object=mlab_axes)
            mlab.orientation_axes()
            v = (
                -22.922482347017496, 82.773599226676239, 42.199999999998845, [-639.82400909, -127.15131564, 5.98696341])
            # mlab.view(azimuth=-22.733470835484969, elevation=78.510163166936607, distance=42.2,
            #           focalpoint=[-639.80083297, -127.26821974, 5.43944316])
            mlab.view(*v)
            mlab.roll(-86.8584361987099)
            v = mlab.view()
            outdir = '{0}{1}-centroid/3Dplot/'.format(datadir, burstname)
            ext = '.png'
            filename = os.path.join(outdir, 'overview{}'.format(ext))
            mlab.savefig(filename)

            doanimation = False

            if doanimation:
                mspts = pts.mlab_source
                if burstname == 'U02':
                    tidxs, = np.where((t >= Time('2014-11-01 16:45:39')) & (t <= Time('2014-11-01 16:45:48')))
                elif burstname == 'U04':
                    tidxs, = np.where((t >= Time('2014-11-01 16:46:16')) & (t <= Time('2014-11-01 16:46:22')))
                prefix = 'ani'


                @mlab.animate(delay=10)
                def anim(tidxs):
                    mfig = mlab.gcf()
                    for i, tt in enumerate(tidxs):
                        # animation updates here
                        mp_ = mp[:, tt]
                        mx_ = mx[:, tt]
                        my_ = my[:, tt]
                        mf_ = mf[:, tt]
                        mt_ = mt[:, tt]
                        if mp_.count():
                            neplot = (ma.compressed(mf_) * 1e9 / 8.98e3) ** 2 / 1e10
                            hgtplot = 30. * np.log(2.5 / neplot) - 5.0
                            pplot = ma.compressed(mp_) ** 3
                            pplot = pplot / np.max(pplot)  # / 2.0
                            c = Time(ma.compressed(mt_), format='jd').plot_date
                            xplot = ma.compressed(mx_)
                            yplot = ma.compressed(my_)
                            # pts.mlab_source.dataset.point_data.vectors = np.transpose(np.tile(pplot, (3, 1)))
                            pts.mlab_source.reset(x=xplot, y=yplot, z=hgtplot, u=pplot, v=pplot, w=pplot, scalars=c)
                            mfig.scene.reset_zoom()
                            mlab.view(*v)
                            tstring = t[tt].iso.replace(' ', 'T').replace(':', '').replace('-', '')
                            filename = os.path.join(outdir, '{}_{}{}'.format(prefix, tstring, ext))
                            mlab.savefig(filename)
                            yield


                anim(tidxs)
                # mlab.show()

                DButil.img2movie(imgprefix=outdir + prefix + '_', outname=outdir + 'burst_3Dplot', size='900x664',
                                 overwrite=True)
