import json
import os, sys
import pickle
import time
from collections import OrderedDict
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
from sys import platform
import scipy.ndimage as sn
from math import radians, cos, sin
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Slider, Button, TextInput, RadioButtonGroup, CheckboxGroup,
                          BoxSelectTool, LassoSelectTool, HoverTool, Spacer, TapTool, LabelSet, Div)
from bokeh.models.mappers import LinearColorMapper
from bokeh.models.widgets import Select, Toggle
from bokeh.palettes import Spectral11
from bokeh.plotting import figure, curdoc
import glob
from astropy.time import Time
from suncasa.utils.puffin import PuffinMap
from suncasa.utils import DButil

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

if platform == "linux" or platform == "linux2":
    print 'Runing QLook in Linux platform'
    for ll in xrange(5100, 5100 + 10):
        os.system('fuser -n tcp -k {}'.format(ll))
elif platform == "darwin":
    print 'Runing QLook in OS X platform'
    for ll in xrange(5100, 5100 + 10):
        os.system(
            'port=($(lsof -i tcp:{}|grep python2.7 |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
        os.system('port=($(lsof -i tcp:{}|grep Google |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
elif platform == "win32":
    print 'Runing QLook in Windows platform'

'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
DButil.initconfig(suncasa_dir)
'''load config file'''
config_main = DButil.loadjsonfile(suncasa_dir + 'DataBrowser/config.json')
database_dir = config_main['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/'
config_EvtID = DButil.loadjsonfile('{}config_EvtID_curr.json'.format(database_dir))
SDOdir = DButil.getSDOdir(config_main, database_dir + '/aiaBrowserData/', suncasa_dir)
spec_square_rs_tmax = config_main['plot_config']['tab_FSview_base']['spec_square_rs_tmax']
spec_square_rs_fmax = config_main['plot_config']['tab_FSview_base']['spec_square_rs_fmax']
spec_image_rs_ratio = config_main['plot_config']['tab_FSview_base']['spec_image_rs_ratio']

# do_spec_regrid = False

'''define the colormaps'''
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]
colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
colormap_viridis = cm.get_cmap("viridis")  # choose any matplotlib colormap here
bokehpalette_viridis = [colors.rgb2hex(m) for m in colormap_viridis(np.arange(colormap_viridis.N))]


def downsample_dspecDF(spec_square_rs_tmax=None, spec_square_rs_fmax=None):
    global dspecDF0_rs, dspecDF0, spec_rs_step
    spec_sz = len(dspecDF0.index)
    spec_sz_max = spec_square_rs_tmax * spec_square_rs_fmax
    if spec_sz > spec_sz_max:
        spec_rs_step = next(i for i in xrange(1, 1000) if spec_sz / i < spec_sz_max)
    else:
        spec_rs_step = 1
    dspecDF0_rs = dspecDF0.loc[::spec_rs_step, :]


def make_spec_plt(spec_plt_R, spec_plt_L):
    spec_I = (spec_plt_R + spec_plt_L) / 2
    spec_V = (spec_plt_R - spec_plt_L) / 2
    spec_max_IRL = int(
        max(spec_plt_R.max(), spec_plt_L.max(), spec_I.max())) * 1.2
    spec_min_IRL = (int(min(spec_plt_R.min(), spec_plt_L.min(), spec_I.min())) / 10) * 10
    spec_max_V = max(abs(int(spec_V.max())), abs(int(spec_V.min()))) * 1.2
    spec_min_V = -spec_max_V
    spec_max_pol = {'RR': spec_max_IRL, 'LL': spec_max_IRL, 'I': spec_max_IRL, 'V': spec_max_V}
    spec_min_pol = {'RR': spec_min_IRL, 'LL': spec_min_IRL, 'I': spec_min_IRL, 'V': spec_min_V}
    spec_pol = {'RR': spec_plt_R, 'LL': spec_plt_L}
    spec_pol['I'] = (spec_pol['RR'] + spec_pol['LL']) / 2
    spec_pol['V'] = (spec_pol['RR'] - spec_pol['LL']) / 2
    return {'spec': spec_pol, 'max': spec_max_pol, 'min': spec_min_pol}


ports = []


def tab2_panel_XCorr_update():
    global clickmode
    if tab2_BUT_XCorr.label == 'XCorr':
        But_ClickMode.label = 'ClickMode: Double'
        tab2_BUT_XCorr.label = 'GoXCorr'
        tab2_panel_Div_exit.text = '<p><b>click two points in dynamic spectrum to select time and frequency range.</b></p>'
    elif 'GoXCorr':
        But_ClickMode.label = 'ClickMode: Single'
        tab2_BUT_XCorr.label = 'XCorr'
        global dspecDF_select
        time0, time1 = Time((dspecDF_select['time'].min() + timestart) / 3600. / 24., format='jd'), Time(
            (dspecDF_select['time'].max() + timestart) / 3600. / 24., format='jd')
        freq0, freq1 = dspecDF_select['freq'].min(), dspecDF_select['freq'].max()
        timeidx0 = next(i for i in xrange(tab2_ntim) if tab2_tim[i] >= time0.mjd * 24. * 3600.)
        timeidx1 = next(i for i in xrange(tab2_ntim - 1, -1, -1) if tab2_tim[i] <= time1.mjd * 24. * 3600.) + 1
        freqidx0 = next(i for i in xrange(tab2_nfreq) if tab2_freq[i] >= freq0)
        freqidx1 = next(i for i in xrange(tab2_nfreq - 1, -1, -1) if tab2_freq[i] <= freq1) + 1
        dspecSel = tab2_r_dspec.data_source.data['image'][0][freqidx0:(freqidx1 + 1), timeidx0:(timeidx1 + 1)]
        freqSel = tab2_freq[freqidx0:(freqidx1 + 1)]
        timSel = tab2_tim[timeidx0:(timeidx1 + 1)]
        CC_dict = DButil.XCorrMap(dspecSel, timSel, freqSel)
        CC_save = struct_dir + 'CC_save.npz'
        np.savez(CC_save, spec=dspecSel, specfit=CC_dict['zfit'], ccmax=CC_dict['ccmax'], ccpeak=CC_dict['ccpeak'],
                 tim=CC_dict['x'], ntim=CC_dict['nx'], timfit=CC_dict['xfit'], ntimfit=CC_dict['nxfit'],
                 freq=CC_dict['y'], nfreq=CC_dict['ny'], freqv=CC_dict['yv'], freqa=CC_dict['ya'],
                 fidxv=CC_dict['yidxv'], fidxa=CC_dict['yidxa'])
        try:
            tab2_panel_Div_exit.text = '<p><b>{}</b> saved.</p>'.format(CC_save)
        except:
            pass
        port = DButil.getfreeport()
        print 'bokeh serve {}DataBrowser/XCorr --show --port {} &'.format(suncasa_dir, port)
        os.system('bokeh serve {}DataBrowser/XCorr --show --port {} &'.format(suncasa_dir, port))
        ports.append(port)
        tab2_r_dspec_patch.data_source.data = {'xx': [], 'yy': []}


def tab2_ClickMode_handler():
    global clickmode
    if But_ClickMode.label == "ClickMode: Single":
        But_ClickMode.label = "ClickMode: Double"
        clickmode = 'singleclick'
    elif But_ClickMode.label == "ClickMode: Double":
        But_ClickMode.label = "ClickMode: Single"
        clickmode = 'doubleclick'


def tab2_dspec_selection_change(attrname, old, new):
    global tapPointDF_dspec, dspecDF_select
    global SRC_dspec_quadselround
    if clickmode == 'singleclick':
        tab2_r_dspec_patch.data_source.data = {'xx': [], 'yy': []}
        tab2_SRC_dspec_quadx_selected = tab2_SRC_dspec_quadx.selected['1d']['indices']
        tab2_SRC_dspec_quady_selected = tab2_SRC_dspec_quady.selected['1d']['indices']
        if len(tab2_SRC_dspec_quadx_selected) > 0 and len(tab2_SRC_dspec_quady_selected) > 0:
            quadx_selected = tab2_SRC_dspec_quadx_selected[0]
            quady_selected = tab2_SRC_dspec_quady_selected[0]
            idx = quady_selected * tab2_ntim + quadx_selected
            tmp = tab2_r_dspec.data_source.data['image'][0]
            timstr = dspecDF0['timestr'].iloc[idx]
            freqstr = '{:.3f} [GHz]'.format(dspecDF0['freq'].iloc[idx])
            tooltips = ['{} {}'.format(timstr, freqstr)]
            tab2_dspec_xPro_data = {'x': tab2_dtim, 'y': tmp[quady_selected, :]}
            r_dspec_xPro.data_source.data = tab2_dspec_xPro_data
            tab2_dspec_xPro_hover_data = {'x': [tab2_dtim[quadx_selected]], 'y':
                [tmp[quady_selected, quadx_selected]],
                                          'tooltips': tooltips}
            r_dspec_xPro_hover.data_source.data = tab2_dspec_xPro_hover_data
            tab2_dspec_yPro_data = {'x': tmp[:, quadx_selected], 'y': tab2_freq}
            r_dspec_yPro.data_source.data = tab2_dspec_yPro_data
            tab2_dspec_yPro_hover_data = {'x': [tmp[quady_selected, quadx_selected]],
                                          'y': [tab2_freq[quady_selected]],
                                          'tooltips': tooltips}
            r_dspec_yPro_hover.data_source.data = tab2_dspec_yPro_hover_data
            tab2_r_dspec_line_x.data_source.data = {'time': [tab2_dtim[0], tab2_dtim[-1]],
                                                    'freq': [tab2_freq[quady_selected]] * 2}
            tab2_r_dspec_line_y.data_source.data = {'time': [tab2_dtim[quadx_selected]] * 2,
                                                    'freq': [tab2_freq[0], tab2_freq[-1]]}
        elif len(tab2_SRC_dspec_quadx_selected) == 0 and len(tab2_SRC_dspec_quady_selected) == 0:
            r_dspec_xPro.data_source.data = {'x': [], 'y': []}
            r_dspec_xPro_hover.data_source.data = {'x': [], 'y': [], 'tooltips': []}
            r_dspec_yPro.data_source.data = {'x': [], 'y': []}
            r_dspec_yPro_hover.data_source.data = {'x': [], 'y': [], 'tooltips': []}
            tab2_r_dspec_line_x.data_source.data = {'time': [], 'freq': []}
            tab2_r_dspec_line_y.data_source.data = {'time': [], 'freq': []}
    elif clickmode == 'doubleclick':
        SRC_dspec_quadselround += 1
        if SRC_dspec_quadselround == 2:
            SRC_dspec_quadselround = 0
            tab2_SRC_dspec_quadx_selected = tab2_SRC_dspec_quadx.selected['1d']['indices']
            tab2_SRC_dspec_quady_selected = tab2_SRC_dspec_quady.selected['1d']['indices']
            tab2_SRC_dspec_quadx.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
            tab2_SRC_dspec_quady.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
            if len(tab2_SRC_dspec_quadx_selected) > 0 and len(tab2_SRC_dspec_quady_selected) > 0:
                quadx_selected = tab2_SRC_dspec_quadx_selected[0]
                quady_selected = tab2_SRC_dspec_quady_selected[0]
                tapPointDF_dspec = tapPointDF_dspec.append(pd.Series(
                    {'xx': tab2_dtim[quadx_selected], 'yy': tab2_freq[quady_selected]}), ignore_index=True)
                if len(tapPointDF_dspec.index) == 1:
                    tab2_r_dspec_lines.data_source.data = {
                        'xs': [[tab2_dtim[quadx_selected]] * 2, [tab2_dtim[0], tab2_dtim[-1]]],
                        'ys': [[tab2_freq[0], tab2_freq[-1]], [tab2_freq[quady_selected]] * 2]}
                    tab2_r_dspec_patch.data_source.data = {'xx': [], 'yy': []}
                elif len(tapPointDF_dspec.index) == 2:
                    x0, x1 = tapPointDF_dspec['xx'].min(), tapPointDF_dspec['xx'].max()
                    y0, y1 = tapPointDF_dspec['yy'].min(), tapPointDF_dspec['yy'].max()
                    tab2_r_dspec_patch.data_source.data = {'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]}
                    tab2_r_dspec_lines.data_source.data = {'xs': [], 'ys': []}
                    dspecDF_select = \
                        dspecDF0_rs[dspecDF0_rs['time'] >= x0][dspecDF0_rs['time'] <= x1][
                            dspecDF0_rs['freq'] >= y0][
                            dspecDF0_rs['freq'] <= y1]
                    tapPointDF_dspec = pd.DataFrame({'xx': [], 'yy': []})
                    time0, time1 = dspecDF_select['time'].min() + timestart, dspecDF_select['time'].max() + timestart
                    t0_char = Time(time0 / 3600. / 24., format='jd', scale='utc', precision=3,
                                   out_subfmt='date_hms').iso
                    date0_char = t0_char.split(' ')[0].replace('-', '/')
                    time0_char = t0_char.split(' ')[1]
                    t1_char = Time(time1 / 3600. / 24., format='jd', scale='utc', precision=3,
                                   out_subfmt='date_hms').iso
                    date1_char = t1_char.split(' ')[0].replace('-', '/')
                    time1_char = t1_char.split(' ')[1]
                    tab2_tCLN_Param_dict['timerange'] = "'{}/{}~{}/{}'".format(date0_char, time0_char, date1_char,
                                                                               time1_char)
                    freq0, freq1 = dspecDF_select['freq'].min(), dspecDF_select['freq'].max()
                    freqrange = "'{:.3f}~{:.3f} GHz'".format(freq0, freq1)
                    tab2_tCLN_Param_dict['freqrange'] = freqrange
                    tab2_Div_tCLN_text = ' '.join(
                        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
                    tab2_Div_tCLN.text = tab2_Div_tCLN_text


def tab2_update_dspec_image(attrname, old, new):
    global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
    select_pol = tab2_Select_pol.value
    select_bl = tab2_Select_bl.value
    bl_index = tab2_bl.index(select_bl)
    spec_plt_R = tab2_spec[0, bl_index, :, :]
    spec_plt_L = tab2_spec[1, bl_index, :, :]
    spec_plt_I = (tab2_spec[0, bl_index, :, :] + tab2_spec[1, bl_index, :, :]) / 2.
    spec_plt_V = (tab2_spec[0, bl_index, :, :] - tab2_spec[1, bl_index, :, :]) / 2.
    spec_plt_max_IRL = int(
        max(spec_plt_R.max(), spec_plt_L.max(), spec_plt_I.max())) * 1.2
    spec_plt_min_IRL = (int(min(spec_plt_R.min(), spec_plt_L.min(), spec_plt_I.min())) / 10) * 10
    spec_plt_max_V = max(abs(int(spec_plt_V.max())), abs(int(spec_plt_V.min()))) * 1.2
    spec_plt_min_V = -spec_plt_max_V
    spec_plt_max_pol = {'RR': spec_plt_max_IRL, 'LL': spec_plt_max_IRL, 'I': spec_plt_max_IRL, 'V': spec_plt_max_V}
    spec_plt_min_pol = {'RR': spec_plt_min_IRL, 'LL': spec_plt_min_IRL, 'I': spec_plt_min_IRL, 'V': spec_plt_min_V}
    spec_plt_pol = {'RR': spec_plt_R, 'LL': spec_plt_L}
    spec_plt_pol['I'] = (spec_plt_pol['RR'] + spec_plt_pol['LL']) / 2
    spec_plt_pol['V'] = (spec_plt_pol['RR'] - spec_plt_pol['LL']) / 2
    spec_plt_max = spec_plt_max_pol[select_pol]
    spec_plt_min = spec_plt_min_pol[select_pol]
    spec_plt = spec_plt_pol[select_pol]
    if select_pol == 'V':
        tab2_Select_colorspace.value = 'linear'
    if tab2_Select_colorspace.value == 'log' and select_pol != 'V':
        tab2_r_dspec.data_source.data['image'] = [np.log(spec_plt)]
    else:
        tab2_r_dspec.data_source.data['image'] = [spec_plt]
    tab2_p_dspec_xPro.y_range.start = spec_plt_min
    tab2_p_dspec_xPro.y_range.end = spec_plt_max
    tab2_p_dspec_yPro.x_range.start = spec_plt_min
    tab2_p_dspec_yPro.x_range.end = spec_plt_max


def tab2_BUT_tCLN_param_add():
    global tab2_tCLN_Param_dict
    tab2_Div_tCLN2.text = ' '
    txts = tab2_input_tCLN.value.strip()
    txts = txts.split(';')
    for txt in txts:
        txt = txt.strip()
        txt = txt.split('=')
        if len(txt) == 2:
            key, val = txt
            tab2_tCLN_Param_dict[key.strip()] = val.strip()
        else:
            tab2_Div_tCLN2.text = '<p>Input syntax: <b>uvtaper</b>=True; <b>niter</b>=200; ' \
                                  '<b>cell</b>=["5.0arcsec", "5.0arcsec"]. Any spaces will be ignored.</p>'
    tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
    tab2_Div_tCLN.text = tab2_Div_tCLN_text


def tab2_BUT_tCLN_param_delete():
    global tab2_tCLN_Param_dict
    tab2_Div_tCLN2.text = ' '
    txts = tab2_input_tCLN.value.strip()
    txts = txts.split(';')
    for key in txts:
        try:
            tab2_tCLN_Param_dict.pop(key)
        except:
            tab2_Div_tCLN2.text = '<p>Input syntax: <b>uvtaper</b>; <b>niter</b>; ' \
                                  '<b>cell</b>. Any spaces will be ignored.</p>'
    tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
    tab2_Div_tCLN.text = tab2_Div_tCLN_text


def tab2_BUT_tCLN_param_default():
    global tab2_tCLN_Param_dict, dspecDF0
    tab2_tCLN_Param_dict = OrderedDict()
    time0, time1 = dspecDF0['time'].min() + timestart, dspecDF0['time'].max() + timestart
    t0_char = Time(time0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
    date0_char = t0_char.split(' ')[0].replace('-', '/')
    time0_char = t0_char.split(' ')[1]
    t1_char = Time(time1 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
    date1_char = t1_char.split(' ')[0].replace('-', '/')
    time1_char = t1_char.split(' ')[1]
    freq0, freq1 = dspecDF0['freq'].min(), dspecDF0['freq'].max()
    freqrange = "'{:.3f}~{:.3f} GHz'".format(freq0, freq1)
    tab2_tCLN_Param_dict['freqrange'] = freqrange
    tab2_tCLN_Param_dict['workdir'] = "'./'"
    tab2_tCLN_Param_dict['vis'] = "''"
    tab2_tCLN_Param_dict['imageprefix'] = "'slfcal/{}{}'".format(struct_id, CleanID)
    tab2_tCLN_Param_dict['ncpu'] = "10"
    tab2_tCLN_Param_dict['twidth'] = "1"
    tab2_tCLN_Param_dict['doreg'] = "False"
    tab2_tCLN_Param_dict['timerange'] = "''"
    tab2_tCLN_Param_dict['uvrange'] = "''"
    tab2_tCLN_Param_dict['antenna'] = "''"
    tab2_tCLN_Param_dict['ephemfile'] = "'horizons_sun_20141101.radecp'"
    tab2_tCLN_Param_dict['msinfofile'] = "'SUN01_20141101.T163940-164700.50ms.cal.msinfo.npz'"
    tab2_tCLN_Param_dict['event_id'] = "'{}'".format(event_id.replace("/", ""))
    tab2_tCLN_Param_dict['struct_id'] = "'{}'".format(struct_id.replace("/", ""))
    tab2_tCLN_Param_dict['clean_id'] = "'{}'".format(CleanID.replace("/", ""))
    tab2_tCLN_Param_dict['mode'] = "'channel'"
    tab2_tCLN_Param_dict['imagermode'] = "'csclean'"
    tab2_tCLN_Param_dict['weighting'] = "'natural'"
    tab2_tCLN_Param_dict['gain'] = '0.1'
    tab2_tCLN_Param_dict['psfmode'] = "'clark'"
    tab2_tCLN_Param_dict['imsize'] = ""'[256, 256]'""
    tab2_tCLN_Param_dict['cell'] = "['3.0arcsec', '3.0arcsec']"
    tab2_tCLN_Param_dict['phasecenter'] = "'J2000 14h26m22.7351 -14d29m29.801'"
    tab2_tCLN_Param_dict['mask'] = "''"
    tab2_tCLN_Param_dict['stokes'] = "'RRLL'"
    tab2_tCLN_Param_dict['uvrange'] = "''"
    tab2_tCLN_Param_dict['niter'] = "500"
    tab2_tCLN_Param_dict['usescratch'] = "False"
    tab2_tCLN_Param_dict['interactive'] = "False"
    tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
    tab2_Div_tCLN.text = tab2_Div_tCLN_text
    tab2_Div_tCLN2.text = '<p><b>Default parameter Restored.</b></p>'


def tab2_BUT_tCLN_param_reload():
    global tab2_tCLN_Param_dict
    with open(CleanID_dir + 'CASA_CLN_args.json', 'r') as fp:
        tab2_tCLN_Param_dict = json.load(fp)
    tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
    tab2_Div_tCLN.text = tab2_Div_tCLN_text
    tab2_Div_tCLN2.text = '<p>CASA arguments reload from config file in <b>{}</b>.</p>'.format(
        CleanID_dir)


def tab2_BUT_tCLN_param_save():
    with open(CleanID_dir + 'CASA_CLN_args.json', 'w') as fp:
        json.dump(tab2_tCLN_Param_dict, fp)
    tab2_Div_tCLN2.text = '<p>CASA script and arguments config file saved to <b>{}</b>.</p>'.format(
        CleanID_dir)
    timestrs = []
    fits_local = []
    fits_global = []
    if 'twidth' in tab2_tCLN_Param_dict.keys():
        val = tab2_tCLN_Param_dict['twidth']
        exec ('twidth = int({})'.format(val))
    else:
        twidth = 1
    if 'workdir' in tab2_tCLN_Param_dict.keys():
        val = tab2_tCLN_Param_dict['workdir']
        exec ('workdir = {}'.format(val))
    else:
        workdir = './'
    # os.system('cp {} {}'.format(CleanID_dir + 'CASA_CLN_args.json', workdir))
    os.system('cp {}/DataBrowser/ToClean/script_clean.py {}'.format(suncasa_dir, CleanID_dir))

    for ii in range(tab2_ntim):
        iit = int(ii) / twidth * twidth
        t0 = xx[iit] - tab2_dt / 2
        datestr = Time(t0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date').iso
        timestr0 = Time(t0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        timestr0 = timestr0.split(' ')[1]
        timestr = datestr.replace("-", "") + 'T' + timestr0.replace(":", "")
        timestrs.append(timestr0)
        fits_local.append(timestr + '.fits')
        fits_global.append(timestr + '.fits')
    timestrs = timestrs * int(tab2_nfreq)
    fits_local = fits_local * int(tab2_nfreq)
    fits_global = fits_global * int(tab2_nfreq)
    freqstrs = ['{:.3f}'.format(ll) for ll in yy]
    dspecDFout = pd.DataFrame({'time': xx - xx[0],
                               'freq': yy,
                               'timestr': timestrs,
                               'freqstr': freqstrs,
                               'dspec': tab2_spec_plt.flatten(),
                               'fits_local': fits_local,
                               'fits_global': fits_global})
    with open(FS_dspecDF, 'wb') as fp:
        pickle.dump(dspecDFout, fp)
    tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-base saved to <b>{}</b>. '.format(
        CleanID_dir) + 'Click the <b>clean</b> button to clean. When finished, \
        go back to <b>QLook</b> window, select StrID <b>{}</b> and \
        click <b>FSview</b> button again.</p>'.format(
        CleanID_dir, struct_id[0:-1])


def tab2_BUT_tCLN_clean():
    cwd = os.getcwd()
    try:
        tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-base saved to <b>{}.</b></p>\
        <p>CASA clean is in processing.</p>'.format(CleanID_dir)
        os.chdir(CleanID_dir)
        suncasapy46 = config_main['core']['casapy46']
        suncasapy46 = os.path.expandvars(suncasapy46)
        os.system('{} -c script_clean.py'.format(suncasapy46))
        tab2_Div_tCLN2.text = '<p>Clean finished, go back to <b>QLook</b> window, select StrID <b>{}</b> and \
            click <b>FSview</b> button again.</p>'.format(
            struct_id[0:-1])
    except:
        tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-base saved to <b>{}.</b></p>\
        <p>Do image clean with CASA manually.</p>'.format(CleanID_dir) + '<p>When finished,\
         go back to <b>QLook</b> window, select StrID <b>{}</b>\
          and click <b>FSview</b> button again.</p>'.format(struct_id[0:-1])
    os.chdir(cwd)


def tab2_panel2_exit():
    tab2_panel_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


event_id = config_EvtID['datadir']['event_id']
try:
    with open(database_dir + event_id + 'CurrFS.json', 'r') as fp:
        FS_config = json.load(fp)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']
struct_dir = database_dir + event_id + struct_id
CleanID = FS_config['datadir']['clean_id']
CleanID_dir = struct_dir + CleanID
FS_specfile = FS_config['datadir']['FS_specfile']
FS_dspecDF = CleanID_dir + 'dspecDF-base'
tab2_specdata = np.load(FS_specfile)
tab2_spec = tab2_specdata['spec']
tab2_npol = tab2_specdata['npol']
tab2_nbl = tab2_specdata['nbl']
tab2_tim = tab2_specdata['tim']
tab2_dt = np.median(np.diff(tab2_tim))
tab2_freq = tab2_specdata['freq'] / 1e9
tab2_freq = [float('{:.03f}'.format(ll)) for ll in tab2_freq]
tab2_df = np.median(np.diff(tab2_freq))
tab2_ntim = len(tab2_tim)
tab2_nfreq = len(tab2_freq)

if isinstance(tab2_specdata['bl'].tolist(), str):
    tab2_bl = tab2_specdata['bl'].item().split(';')
elif isinstance(tab2_specdata['bl'].tolist(), list):
    tab2_bl = ['&'.join(ll) for ll in tab2_specdata['bl'].tolist()]
else:
    raise ValueError('Please check the data of {}'.format(FS_specfile))
SRC_dspec_quadselround = 0
clickmode = 'singleclick'
bl_index = 0
tab2_pol = 'I'
sz_spec = tab2_spec.shape
tab2_spec_plt_R = tab2_spec[0, bl_index, :, :]
tab2_spec_plt_L = tab2_spec[1, bl_index, :, :]
spec_pol_dict = make_spec_plt(tab2_spec_plt_R, tab2_spec_plt_L)
tab2_spec_plt_pol = spec_pol_dict['spec']
spec_plt_max_pol = spec_pol_dict['max']
spec_plt_min_pol = spec_pol_dict['min']
tab2_spec_plt = tab2_spec_plt_pol[tab2_pol]
spec_plt_max = spec_plt_max_pol[tab2_pol]
spec_plt_min = spec_plt_min_pol[tab2_pol]

tab2_dtim = tab2_tim - tab2_tim[0]
tab2_dur = tab2_dtim[-1] - tab2_dtim[0]
tim_map = ((np.tile(tab2_tim, tab2_nfreq).reshape(tab2_nfreq, tab2_ntim) / 3600. / 24. + 2400000.5)) * 86400.
freq_map = np.tile(tab2_freq, tab2_ntim).reshape(tab2_ntim, tab2_nfreq).swapaxes(0, 1)
xx = tim_map.flatten()
yy = freq_map.flatten()
timestart = xx[0]

tab2_panel_Div_exit = Div(
    text="""<p><b>Warning</b>: Click the <b>Exit ToClean</b>
            first before closing the tab</p></b>""",
    width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'])
timestrs = []
for ll in tab2_tim:
    timestr = Time(ll / 3600. / 24., format='mjd', scale='utc', precision=3, out_subfmt='date_hms').iso
    timestrs.append(timestr.split(' ')[1])
timestrs = timestrs * tab2_nfreq
dspecDF0 = pd.DataFrame({'time': xx - xx[0],
                         'freq': yy,
                         'dspec': tab2_spec_plt.flatten(),
                         'timestr': timestrs})

rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()
# timestart = xx[0]
TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
# downsample_dspecDF(spec_square_rs_tmax=spec_square_rs_tmax, spec_square_rs_fmax=spec_square_rs_fmax)
# tab2_SRC_dspec_square = ColumnDataSource(dspecDF0_rs)
dspecDF0_rs = dspecDF0.copy()

'''create the dynamic spectrum plot'''
tab2_p_dspec = figure(tools=TOOLS, webgl=config_main['plot_config']['WebGL'],
                      plot_width=config_main['plot_config']['tab_ToClean']['dspec_wdth'],
                      plot_height=config_main['plot_config']['tab_ToClean']['dspec_hght'],
                      x_range=(tab2_dtim[0] - tab2_dt / 2.0, tab2_dtim[-1] + tab2_dt / 2.0),
                      y_range=(tab2_freq[0] - tab2_df / 2.0, tab2_freq[-1] + tab2_df / 2.0),
                      toolbar_location="above")
tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
tab2_p_dspec.axis.visible = True
tab2_p_dspec.title.text = "Dynamic spectrum"
tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
tab2_r_dspec = tab2_p_dspec.image(image=[tab2_spec_plt], x=tab2_dtim[0] - tab2_dt / 2.0, y=tab2_freq[0] - tab2_df / 2.0,
                                  dw=tab2_dtim[-1] - tab2_dtim[0] + tab2_dt,
                                  dh=tab2_freq[-1] - tab2_freq[0] + tab2_df, palette=bokehpalette_jet)

tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.0, line_alpha=0.8,
                                        line_color='white', source=tab2_source_idx_line_x)
tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.0, line_alpha=0.8,
                                        line_color='white', source=tab2_source_idx_line_y)

xLFE = tab2_dtim - tab2_dt / 2.0
xRTE = np.append(xLFE[1:], xLFE[-1] + tab2_dt)
yBTE = tab2_freq - tab2_df / 2.0
yTPE = np.append(yBTE[1:], yBTE[-1] + tab2_df)
tab2_SRC_dspec_quadx = ColumnDataSource(
    {'left': xLFE.ravel(), 'right': xRTE.ravel(), 'bottom': np.repeat(yBTE[0], tab2_ntim),
     'top': np.repeat(yBTE[-1] + tab2_df, tab2_ntim)})
tab2_SRC_dspec_quady = ColumnDataSource(
    {'left': np.repeat(xLFE[0], tab2_nfreq), 'right': np.repeat(xLFE[-1] + tab2_dt, tab2_nfreq),
     'bottom': yBTE.ravel(),
     'top': yTPE.ravel()})
tab2_r_dspec_quadx = tab2_p_dspec.quad('left', 'right', 'top', 'bottom', source=tab2_SRC_dspec_quadx,
                                       fill_alpha=0.0, fill_color=None,
                                       line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                       selection_fill_color=None,
                                       nonselection_fill_alpha=0.0,
                                       selection_line_alpha=0.0, selection_line_color=None,
                                       nonselection_line_alpha=0.0)
tab2_r_dspec_quady = tab2_p_dspec.quad('left', 'right', 'top', 'bottom', source=tab2_SRC_dspec_quady,
                                       fill_alpha=0.0, fill_color=None,
                                       line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                       selection_fill_color=None,
                                       nonselection_fill_alpha=0.0,
                                       selection_line_alpha=0.0, selection_line_color=None,
                                       nonselection_line_alpha=0.0)
tab2_p_dspec.add_tools(TapTool(renderers=[tab2_r_dspec_quadx, tab2_r_dspec_quady]))
tab2_SRC_dspec_quadx.on_change('selected', tab2_dspec_selection_change)
tab2_SRC_dspec_quady.on_change('selected', tab2_dspec_selection_change)

tapPointDF_dspec = pd.DataFrame({'xx': [], 'yy': []})
tab2_SRC_dspec_lines = ColumnDataSource({'xs': [], 'ys': []})
tab2_r_dspec_lines = tab2_p_dspec.multi_line('xs', 'ys', source=tab2_SRC_dspec_lines, line_color='Magenta',
                                             line_width=1, alpha=0.8)
tab2_SRC_dspec_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
tab2_r_dspec_patch = tab2_p_dspec.patch('xx', 'yy', source=tab2_SRC_dspec_Patch,
                                        fill_color=None, fill_alpha=0.5, line_color="Magenta",
                                        line_alpha=1.0, line_width=1)
tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                        line_color='white', source=tab2_source_idx_line_x)
tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                        line_color='white', source=tab2_source_idx_line_y)
# tab2_p_dspec.add_tools(BoxSelectTool())
# tab2_p_dspec.select(BoxSelectTool).select_every_mousemove = False
# tab2_p_dspec.border_fill_color = "silver"
tab2_p_dspec.border_fill_alpha = 0.4
tab2_p_dspec.axis.major_tick_out = 0
tab2_p_dspec.axis.major_tick_in = 5
tab2_p_dspec.axis.minor_tick_out = 0
tab2_p_dspec.axis.minor_tick_in = 3
tab2_p_dspec.axis.major_tick_line_color = "white"
tab2_p_dspec.axis.minor_tick_line_color = "white"

# tab2_dspec_selected = None



# tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)

tab2_Select_pol = Select(title="Polarization:", value='I', options=['RR', 'LL', 'I', 'V'],
                         width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'])
tab2_Select_bl = Select(title="Baseline:", value=tab2_bl[0], options=tab2_bl,
                        width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'])
tab2_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                                width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'])

tab2_p_dspec_xPro = figure(tools='', webgl=config_main['plot_config']['WebGL'],
                           plot_width=config_main['plot_config']['tab_ToClean']['dspec_xPro_wdth'],
                           plot_height=config_main['plot_config']['tab_ToClean']['dspec_xPro_hght'],
                           x_range=tab2_p_dspec.x_range, y_range=(spec_plt_min, spec_plt_max),
                           title="Time profile", toolbar_location=None)
tab2_SRC_dspec_xPro = ColumnDataSource({'x': [], 'y': []})
tab2_SRC_dspec_xPro_hover = ColumnDataSource({'x': [], 'y': [], 'tooltips': []})
r_dspec_xPro = tab2_p_dspec_xPro.line(x='x', y='y', alpha=1.0, line_width=1, source=tab2_SRC_dspec_xPro)
r_dspec_xPro_c = tab2_p_dspec_xPro.circle(x='x', y='y', size=5, fill_alpha=0.2, fill_color='grey',
                                          line_color=None,
                                          source=tab2_SRC_dspec_xPro)
r_dspec_xPro_hover = tab2_p_dspec_xPro.circle(x='x', y='y', size=5, fill_alpha=0.5, fill_color='firebrick',
                                              line_color='firebrick', source=tab2_SRC_dspec_xPro_hover)
l_dspec_xPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                              source=tab2_SRC_dspec_xPro_hover,
                              render_mode='canvas')
l_dspec_xPro_hover.text_font_size = '10pt'
tab2_p_dspec_xPro.add_layout(l_dspec_xPro_hover)
tab2_p_dspec_xPro.title.text_font_size = '6pt'
tab2_p_dspec_xPro.background_fill_color = "beige"
tab2_p_dspec_xPro.background_fill_alpha = 0.4
tab2_p_dspec_xPro.xaxis.axis_label = 'Seconds since ' + tim0_char
tab2_p_dspec_xPro.yaxis.axis_label_text_font_size = '5px'
tab2_p_dspec_xPro.yaxis.axis_label = 'Intensity [sfu]'
# tab2_p_dspec_xPro.border_fill_color = "silver"
tab2_p_dspec_xPro.border_fill_alpha = 0.4
tab2_p_dspec_xPro.axis.major_tick_out = 0
tab2_p_dspec_xPro.axis.major_tick_in = 5
tab2_p_dspec_xPro.axis.minor_tick_out = 0
tab2_p_dspec_xPro.axis.minor_tick_in = 3
tab2_p_dspec_xPro.axis.major_tick_line_color = "black"
tab2_p_dspec_xPro.axis.minor_tick_line_color = "black"

tab2_p_dspec_yPro = figure(tools='', webgl=config_main['plot_config']['WebGL'],
                           plot_width=config_main['plot_config']['tab_ToClean']['dspec_yPro_wdth'],
                           plot_height=config_main['plot_config']['tab_ToClean']['dspec_yPro_hght'],
                           x_range=(spec_plt_min, spec_plt_max), y_range=tab2_p_dspec.y_range,
                           title="Frequency profile", toolbar_location=None)
tab2_SRC_dspec_yPro = ColumnDataSource({'x': [], 'y': []})
tab2_SRC_dspec_yPro_hover = ColumnDataSource({'x': [], 'y': [], 'tooltips': []})
r_dspec_yPro = tab2_p_dspec_yPro.line(x='x', y='y', alpha=1.0, line_width=1, source=tab2_SRC_dspec_yPro)
r_dspec_yPro_c = tab2_p_dspec_yPro.circle(x='x', y='y', size=5, fill_alpha=0.2, fill_color='grey',
                                          line_color=None,
                                          source=tab2_SRC_dspec_yPro)
r_dspec_yPro_hover = tab2_p_dspec_yPro.circle(x='x', y='y', size=5, fill_alpha=0.5, fill_color='firebrick',
                                              line_color='firebrick', source=tab2_SRC_dspec_yPro_hover)
l_dspec_yPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph', source=tab2_SRC_dspec_yPro_hover,
                              render_mode='canvas')
l_dspec_yPro_hover.text_font_size = '10pt'
tab2_p_dspec_yPro.add_layout(l_dspec_yPro_hover)
tab2_p_dspec_yPro.title.text_font_size = '6pt'
tab2_p_dspec_yPro.yaxis.visible = False
tab2_p_dspec_yPro.background_fill_color = "beige"
tab2_p_dspec_yPro.background_fill_alpha = 0.4
tab2_p_dspec_yPro.xaxis.axis_label = 'Intensity [sfu]'
tab2_p_dspec_yPro.xaxis.axis_label_text_font_size = '5px'
tab2_p_dspec_yPro.min_border_bottom = 0
tab2_p_dspec_yPro.min_border_left = 0
# tab2_p_dspec_yPro.border_fill_color = "silver"
tab2_p_dspec_yPro.border_fill_alpha = 0.4
tab2_p_dspec_yPro.axis.major_tick_out = 0
tab2_p_dspec_yPro.axis.major_tick_in = 5
tab2_p_dspec_yPro.axis.minor_tick_out = 0
tab2_p_dspec_yPro.axis.minor_tick_in = 3
tab2_p_dspec_yPro.axis.major_tick_line_color = "black"
tab2_p_dspec_yPro.axis.minor_tick_line_color = "black"

tab2_ctrls = [tab2_Select_bl, tab2_Select_pol, tab2_Select_colorspace]
for ctrl in tab2_ctrls:
    ctrl.on_change('value', tab2_update_dspec_image)

tab2_input_tCLN = TextInput(value="Input the param here", title="Clean task parameters:",
                            width=config_main['plot_config']['tab_ToClean']['input_tCLN_wdth'])
tab2_Div_tCLN = Div(text='', width=config_main['plot_config']['tab_ToClean']['tab2_Div_tCLN_wdth'])
tab2_Div_tCLN2 = Div(text='', width=config_main['plot_config']['tab_ToClean']['input_tCLN_wdth'])

tab2_BUT_tCLN_param_default()

tab2_BUT_tCLN_param_ADD = Button(label='Add to Param',
                                 width=config_main['plot_config']['tab_ToClean']['button_wdth'],
                                 button_type='primary')
tab2_BUT_tCLN_param_ADD.on_click(tab2_BUT_tCLN_param_add)
tab2_BUT_tCLN_param_DEL = Button(label='Delete Param',
                                 width=config_main['plot_config']['tab_ToClean']['button_wdth'],
                                 button_type='warning')
tab2_BUT_tCLN_param_DEL.on_click(tab2_BUT_tCLN_param_delete)
tab2_BUT_tCLN_param_Default = Button(label='Default Param',
                                     width=config_main['plot_config']['tab_ToClean']['button_wdth'])
tab2_BUT_tCLN_param_Default.on_click(tab2_BUT_tCLN_param_default)
tab2_SPCR_LFT_BUT_tCLN_param_DEL = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_SPCR_LFT_BUT_tCLN_param_RELOAD = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_SPCR_LFT_BUT_tCLN_param_SAVE = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_SPCR_LFT_BUT_tCLN_param_DEFAULT = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_SPCR_LFT_BUT_tCLN_param_FSVIEW = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_SPCR_LFT_Div_tCLN2 = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth50'],
                                 height=config_main['plot_config']['tab_ToClean']['space_hght10'])

tab2_BUT_tCLN_param_RELOAD = Button(label='reload Param',
                                    width=config_main['plot_config']['tab_ToClean']['button_wdth'])
tab2_BUT_tCLN_param_RELOAD.on_click(tab2_BUT_tCLN_param_reload)
tab2_BUT_tCLN_param_SAVE = Button(label='save Param',
                                  width=config_main['plot_config']['tab_ToClean']['button_wdth'],
                                  button_type='success')
tab2_BUT_tCLN_param_SAVE.on_click(tab2_BUT_tCLN_param_save)
tab2_BUT_tCLN_CLEAN = Button(label='clean',
                             width=config_main['plot_config']['tab_ToClean']['button_wdth'],
                             button_type='success')
tab2_SPCR_ABV_BUT_tCLN = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'],
                                height=config_main['plot_config']['tab_ToClean']['space_hght18'])
tab2_BUT_tCLN_CLEAN.on_click(tab2_BUT_tCLN_clean)
tab2_SPCR_LFT_BUT_CLEAN = Spacer(width=config_main['plot_config']['tab_ToClean']['space_wdth10'])
tab2_BUT_FSview = Button(label='FSview', width=config_main['plot_config']['tab_ToClean']['button_wdth'],
                         button_type='primary')

tab2_BUT_XCorr = Button(label='XCorr',
                        width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'],
                        button_type='warning')
tab2_BUT_XCorr.on_click(tab2_panel_XCorr_update)

But_ClickMode = Button(label="ClickMode: Double",
                       width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'], button_type="primary")
But_ClickMode.on_click(tab2_ClickMode_handler)

tab2_panel2_BUT_exit = Button(label='Exit ToClean',
                              width=config_main['plot_config']['tab_ToClean']['widgetbox_wdth1'],
                              button_type='danger')
tab2_panel2_BUT_exit.on_click(tab2_panel2_exit)
panel2 = row(column(row(tab2_p_dspec, tab2_p_dspec_yPro),
                    row(column(row(column(tab2_p_dspec_xPro,
                                          row(tab2_input_tCLN),
                                          row(tab2_BUT_tCLN_param_ADD, tab2_SPCR_LFT_BUT_tCLN_param_DEL,
                                              tab2_BUT_tCLN_param_DEL, tab2_SPCR_LFT_BUT_tCLN_param_DEFAULT,
                                              tab2_BUT_tCLN_param_Default, tab2_SPCR_LFT_BUT_tCLN_param_RELOAD,
                                              tab2_BUT_tCLN_param_RELOAD, tab2_SPCR_LFT_BUT_tCLN_param_SAVE,
                                              tab2_BUT_tCLN_param_SAVE, tab2_SPCR_LFT_BUT_CLEAN,
                                              tab2_BUT_tCLN_CLEAN)),
                                   widgetbox(tab2_Select_pol, tab2_Select_bl, tab2_Select_colorspace,
                                             But_ClickMode, tab2_BUT_XCorr, tab2_panel2_BUT_exit, tab2_panel_Div_exit,
                                             width=config_main['plot_config']['tab_ToClean'][
                                                 'widgetbox_wdth2'])), tab2_Div_tCLN2))), tab2_Div_tCLN)

curdoc().add_root(panel2)
curdoc().title = "ToClean"
