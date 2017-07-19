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
                          BoxSelectTool, LassoSelectTool, HoverTool, TapTool, Spacer, LabelSet, Div)
from bokeh.models.mappers import LinearColorMapper
from bokeh.models.widgets import Select, RangeSlider
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
tidx_prev = None

# do_spec_regrid = False

'''define the colormaps'''
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]
colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
colormap_viridis = cm.get_cmap("viridis")  # choose any matplotlib colormap here
bokehpalette_viridis = [colors.rgb2hex(m) for m in colormap_viridis(np.arange(colormap_viridis.N))]
'''
-------------------------- panel 2,3   --------------------------
'''


def read_fits(fname):
    hdulist = fits.open(fname)
    hdu = hdulist[0]
    # hdulist.close()
    return hdu


def goodchan(hdu):
    ndx = hdu.header["NAXIS1"]
    ndy = hdu.header["NAXIS2"]
    xc = ndx / 2
    yc = ndy / 2
    hdu_goodchan = \
        np.where(np.nanmean(hdu.data[0, :, yc - ndy / 16:yc + ndy / 16, xc - ndx / 16:xc + ndx / 16], axis=(-1, -2)))[0]
    return hdu_goodchan


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


def tab2_vdspec_update():
    global tab2_Select_pol_opt, spec_pol_dict
    select_pol = tab2_Select_pol.value
    if tab2_BUT_vdspec.label == "VEC Dyn Spec":
        if len(tab2_r_vla_ImgRgn_patch.data_source.data['xx']) > 0:
            tab2_Select_pol_opt = tab2_Select_pol.options
            tab2_Select_pol.options = pols
            tab2_BUT_vdspec.label = "Dyn Spec"
            spec_plt_R = np.zeros((tab2_nfreq, tab2_ntim))
            spec_plt_L = np.zeros((tab2_nfreq, tab2_ntim))
            if len(pols) > 1:
                for ll in xrange(tab2_ntim):
                    hdufile = fits_LOCL_dir + dspecDF0.loc[ll, :]['fits_local']
                    if os.path.exists(hdufile):
                        hdu = read_fits(hdufile)
                        hdu_goodchan = goodchan(hdu)
                        nfreq_hdu = hdu_goodchan[-1] - hdu_goodchan[0] + 1
                        freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
                        freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
                        try:
                            idxfreq = freq.index(freq_ref)
                        except:
                            idxfreq = (float(freq_ref) - float(freq[0])) / float('{:.3f}'.format(tab2_df))
                            idxfreq = int(np.round(idxfreq))
                        vla_l = hdu.data[0, :, y0pix:y1pix + 1, x0pix:x1pix + 1]
                        vla_r = hdu.data[1, :, y0pix:y1pix + 1, x0pix:x1pix + 1]
                        if idxfreq >= 0:
                            spec_plt_R[idxfreq:idxfreq + nfreq_hdu, ll] = \
                                np.nanmean(vla_l, axis=(-1, -2))[hdu_goodchan[0]:hdu_goodchan[-1] + 1]
                            spec_plt_L[idxfreq:idxfreq + nfreq_hdu, ll] = \
                                np.nanmean(vla_r, axis=(-1, -2))[hdu_goodchan[0]:hdu_goodchan[-1] + 1]
                        else:
                            spec_plt_R[0:idxfreq + nfreq_hdu, ll] = \
                                np.nanmean(vla_l, axis=(-1, -2))[hdu_goodchan[0] - idxfreq:hdu_goodchan[-1] + 1]
                            spec_plt_L[0:idxfreq + nfreq_hdu, ll] = \
                                np.nanmean(vla_r, axis=(-1, -2))[hdu_goodchan[0] - idxfreq:hdu_goodchan[-1] + 1]
                    tab2_Div_LinkImg_plot.text = """<p><b>Vec Dspec in calculating...</b></p><p>{}</p>""".format(
                        DButil.ProgressBar(ll + 1, tab2_ntim + 1, decimals=0, length=16, empfill='=', fill='#'))
                spec_plt_R[spec_plt_R < 0] = 0
                spec_plt_L[spec_plt_L < 0] = 0
                tab2_Div_LinkImg_plot.text = """<p><b>Vec Dspec in calculating...</b></p><p>{}</p>""".format(
                    DButil.ProgressBar(tab2_ntim + 1, tab2_ntim + 1, decimals=0, length=16, empfill='=',
                                       fill='#'))
            elif len(pols) == 1:
                for ll in xrange(tab2_ntim):
                    hdufile = fits_LOCL_dir + dspecDF0.loc[ll, :]['fits_local']
                    if os.path.exists(hdufile):
                        hdu = read_fits(hdufile)
                        hdu_goodchan = goodchan(hdu)
                        nfreq_hdu = hdu_goodchan[-1] - hdu_goodchan[0] + 1
                        freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
                        freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
                        idxfreq = freq.index(freq_ref)
                        vladata = hdu.data[0, :, y0pix:y1pix + 1, x0pix:x1pix + 1]
                        vlaflux = np.nanmean(vladata, axis=(-1, -2))[hdu_goodchan[0]:hdu_goodchan[-1] + 1]
                        spec_plt_R[idxfreq:idxfreq + nfreq_hdu, ll] = vlaflux
                    tab2_Div_LinkImg_plot.text = """<p><b>Vec Dspec in calculating...</b></p><p>{}</p>""".format(
                        DButil.ProgressBar(ll + 1, tab2_ntim + 1, decimals=0, length=16, empfill='=', fill='#'))
                spec_plt_R[spec_plt_R < 0] = 0
                spec_plt_L = spec_plt_R
                tab2_Div_LinkImg_plot.text = """<p><b>Vec Dspec in calculating...</b></p><p>{}</p>""".format(
                    DButil.ProgressBar(tab2_ntim + 1, tab2_ntim + 1, decimals=0, length=16, empfill='=',
                                       fill='#'))
            tab2_Div_LinkImg_plot.text = '<p><b>Vec Dspec calculated.</b></p>'
            spec_pol_dict = make_spec_plt(spec_plt_R, spec_plt_L)
            tab2_bl_pol_cls_change(None, select_pol)
            tab2_p_dspec.title.text = "Vector Dynamic spectrum"
        else:
            tab2_Div_LinkImg_plot.text = '<p><b>Warning:</b> select a region first.</p>'
    else:
        tab2_Select_pol.options = tab2_Select_pol_opt
        select_bl = tab2_Select_bl.value
        tab2_BUT_vdspec.label = "VEC Dyn Spec"
        bl_index = tab2_bl.index(select_bl)
        spec_plt_R = tab2_spec[0, bl_index, :, :]
        spec_plt_L = tab2_spec[1, bl_index, :, :]
        spec_pol_dict = make_spec_plt(spec_plt_R, spec_plt_L)
        # spec_plt_pol = spec_pol_dict['spec']
        tab2_bl_pol_cls_change(bl_index, select_pol)
        tab2_p_dspec.title.text = "Dynamic spectrum"
        tab2_Div_LinkImg_plot.text = ''


ports = []


def tab2_panel_XCorr_update():
    global clickmode
    if tab2_BUT_XCorr.label == 'XCorr':
        clickmode = 'doubleclick'
        But_ClickMode.label = 'ClickMode: Double'
        tab2_BUT_XCorr.label = 'GoXCorr'
        tab2_Div_LinkImg_plot.text = '<p><b>click two points in dynamic spectrum to select time and frequency range.</b></p>'
    elif 'GoXCorr':
        clickmode = 'singleclick'
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
            tab2_Div_LinkImg_plot.text = '<p><b>{}</b> saved.</p>'.format(CC_save)
        except:
            pass
        port = DButil.getfreeport()
        print 'bokeh serve {}DataBrowser/XCorr --show --port {} &'.format(suncasa_dir, port)
        os.system('bokeh serve {}DataBrowser/XCorr --show --port {} &'.format(suncasa_dir, port))
        ports.append(port)
        tab2_r_dspec_patch.data_source.data = {'xx': [], 'yy': []}


def tab2_update_dspec_image(attrname, old, new):
    global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
    select_pol = tab2_Select_pol.value
    select_bl = tab2_Select_bl.value
    bl_index = tab2_bl.index(select_bl)
    if tab2_BUT_vdspec.label == "VEC Dyn Spec":
        tab2_bl_pol_cls_change(bl_index, select_pol)
    else:
        tab2_bl_pol_cls_change(None, select_pol)


def tab2_bl_pol_cls_change(bl_index, select_pol):
    global spec_pol_dict, dspecDF0_rs
    global tab2_SRC_dspec_square, tab2_p_dspec
    global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
    if tab2_BUT_vdspec.label == "VEC Dyn Spec":
        spec_plt_R = tab2_spec[0, bl_index, :, :]
        spec_plt_L = tab2_spec[1, bl_index, :, :]
        spec_pol_dict = make_spec_plt(spec_plt_R, spec_plt_L)
    # else:
    #     spec_plt_R = spec_pol_dict['spec']['RR']
    #     spec_plt_L = spec_pol_dict['spec']['LL']

    if select_pol == 'V':
        tab2_Select_colorspace.value = 'linear'
    if tab2_Select_colorspace.value == 'log' and select_pol != 'V':
        tmp = spec_pol_dict['spec'][select_pol]
        tmp[tmp < 1.0] = 1.0
        tab2_r_dspec.data_source.data['image'] = [np.log(tmp)]
    else:
        tab2_r_dspec.data_source.data['image'] = [spec_pol_dict['spec'][select_pol]]
    # tab2_SRC_dspec_square.data['dspec'] = spec_pol_dict['spec'][select_pol].flatten()
    tab2_p_dspec_xPro.y_range.start = spec_pol_dict['min'][select_pol]
    tab2_p_dspec_xPro.y_range.end = spec_pol_dict['max'][select_pol]
    tab2_p_dspec_yPro.x_range.start = spec_pol_dict['min'][select_pol]
    tab2_p_dspec_yPro.x_range.end = spec_pol_dict['max'][select_pol]


def slider_LinkImg_update(polonly=False):
    global hdu, select_vla_pol, dspecDF0, tidx_prev
    select_vla_pol = tab2_Select_vla_pol.value
    tab2_Slider_time_LinkImg.start = next(
        i for i in xrange(tab2_ntim) if tab2_dtim[i] >= tab2_p_dspec.x_range.start)
    tab2_Slider_time_LinkImg.end = next(
        i for i in xrange(tab2_ntim - 1, -1, -1) if tab2_dtim[i] <= tab2_p_dspec.x_range.end) + 1
    tab2_Slider_freq_LinkImg.start = next(
        i for i in xrange(tab2_nfreq) if tab2_freq[i] >= tab2_p_dspec.y_range.start)
    tab2_Slider_freq_LinkImg.end = next(
        i for i in xrange(tab2_nfreq - 1, -1, -1) if tab2_freq[i] <= tab2_p_dspec.y_range.end) + 1
    tidx = int(tab2_Slider_time_LinkImg.value)
    fidx = int(tab2_Slider_freq_LinkImg.value)
    tab2_r_dspec_line_x.data_source.data = ColumnDataSource(
        pd.DataFrame({'time': [tab2_dtim[tidx], tab2_dtim[tidx]],
                      'freq': [tab2_freq[0], tab2_freq[-1]]})).data
    tab2_r_dspec_line_y.data_source.data = ColumnDataSource(
        pd.DataFrame({'time': [tab2_dtim[0], tab2_dtim[-1]],
                      'freq': [tab2_freq[fidx], tab2_freq[fidx]]})).data
    timeidxs = np.unique(dspecDF0['time'])
    fitsfile = dspecDF0[dspecDF0.time == timeidxs[tidx]].iloc[0]['fits_local']
    hdufile = fits_LOCL_dir + fitsfile
    print hdufile
    if os.path.exists(hdufile):
        if tidx != tidx_prev:
            if not polonly or not 'hdu' in globals():
                hdu = read_fits(hdufile)
        hdu_goodchan = goodchan(hdu)
        freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
        freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
        try:
            idxfreq = freq.index(freq_ref)
        except:
            idxfreq = (float(freq_ref) - float(freq[0])) / float('{:.3f}'.format(tab2_df))
            idxfreq = int(np.round(idxfreq))
        fidx_hdu = fidx - idxfreq
        print 'tidx,tidx_prev,fidx:', tidx, tidx_prev, fidx_hdu
        if hdu_goodchan[0] <= fidx_hdu <= hdu_goodchan[-1]:
            if select_vla_pol == 'RR':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :]
            elif select_vla_pol == 'LL':
                vladata = hdu.data[pols.index('LL'), fidx_hdu, :, :]
            elif select_vla_pol == 'I':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :] + hdu.data[pols.index('1'), fidx_hdu, :, :]
            elif select_vla_pol == 'V':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :] - hdu.data[pols.index('1'), fidx_hdu, :, :]
            pfmap = PuffinMap(vladata, hdu.header, plot_height=tab2_LinkImg_HGHT,
                              plot_width=tab2_LinkImg_WDTH, webgl=config_main['plot_config']['WebGL'])
            tab2_r_vla.data_source.data['image'] = pfmap.ImageSource()['data']
            mapx, mapy = pfmap.meshgrid()
            mapx, mapy = mapx.value, mapy.value
            SRC_contour = DButil.get_contour_data(mapx, mapy, pfmap.smap.data)
            tab2_r_vla_multi_line.data_source.data = SRC_contour.data
            tab2_Div_LinkImg_plot.text = '<p><b>{}</b> loaded.</p>'.format(fitsfile)
        else:
            tab2_Div_LinkImg_plot.text = '<p><b>freq idx</b> out of range.</p>'
    else:
        tab2_Div_LinkImg_plot.text = '<p><b>{}</b> not found.</p>'.format(fitsfile)
    tidx_prev = tidx


def tab3_slider_LinkImg_update(attrname, old, new):
    slider_LinkImg_update()


def tab2_Select_vla_pol_update(attrname, old, new):
    global hdu, select_vla_pol, tidx_prev
    select_vla_pol = tab2_Select_vla_pol.value
    slider_LinkImg_update()


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
            tab2_Slider_time_LinkImg.value = quadx_selected
            tab2_Slider_freq_LinkImg.value = quady_selected
            idx = quady_selected * tab2_ntim + quadx_selected
            tmp = tab2_r_dspec.data_source.data['image'][0]
            timstr = dspecDF0['timestr'].iloc[idx]
            freqstr = dspecDF0['freqstr'].iloc[idx] + ' [GHz]'
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


def tab2_update_MapRES(attrname, old, new):
    start_timestamp = time.time()
    select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
    dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
    aia_resampled_map = aiamap.resample(dimensions)
    aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                    plot_height=config_main['plot_config']['tab_FSview_base']['aia_hght'],
                                    plot_width=config_main['plot_config']['tab_FSview_base']['aia_wdth'],
                                    webgl=config_main['plot_config']['WebGL'])
    tab2_r_aia.data_source.data['image'] = aia_resampled_pfmap.ImageSource()['data']
    hmi_resampled_map = hmimap.resample(dimensions)
    hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                    plot_height=config_main['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_main['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_main['plot_config']['WebGL'])
    # SRC_HMI = hmi_resampled_pfmap.ImageSource()
    tab2_r_hmi.data_source.data['image'] = hmi_resampled_pfmap.ImageSource()['data']
    print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


def tab2_save_region():
    if len(tab2_r_vla_ImgRgn_patch.data_source.data['xx']) > 0:
        pangle = hdu.header['p_angle']
        x0Deg, x1Deg, y0Deg, y1Deg = (x0 - hdu.header['CRVAL1']) / 3600., (
            x1 - hdu.header['CRVAL1']) / 3600., (
                                         y0 - hdu.header['CRVAL2']) / 3600., (
                                         y1 - hdu.header['CRVAL2']) / 3600.
        p0 = -pangle
        prad = radians(p0)
        dx0 = (x0Deg) * cos(prad) - y0Deg * sin(prad)
        dy0 = (x0Deg) * sin(prad) + y0Deg * cos(prad)
        dx1 = (x1Deg) * cos(prad) - y1Deg * sin(prad)
        dy1 = (x1Deg) * sin(prad) + y1Deg * cos(prad)
        x0Deg, x1Deg, y0Deg, y1Deg = (dx0 + hdu.header['CRVAL1'] / 3600.), (
            dx1 + hdu.header['CRVAL1'] / 3600.), (dy0 + hdu.header['CRVAL2'] / 3600.), (
                                         dy1 + hdu.header['CRVAL2'] / 3600.)
        c0fits = SkyCoord(ra=(x0Deg) * u.degree, dec=(y0Deg) * u.degree)
        c1fits = SkyCoord(ra=(x1Deg) * u.degree, dec=(y1Deg) * u.degree)
        rgnfits = '#CRTFv0 CASA Region Text Format version 0\n\
        box [[{}], [{}]] coord=J2000, linewidth=1, \
        linestyle=-, symsize=1, symthick=1, color=magenta, \
        font="DejaVu Sans", fontsize=11, fontstyle=normal, \
        usetex=false'.format(', '.join(c0fits.to_string('hmsdms').split(' ')),
                             ', '.join(c1fits.to_string('hmsdms').split(' ')))
        with open(rgnfitsfile, "w") as fp:
            fp.write(rgnfits)
        tab2_Div_LinkImg_plot.text = '<p>region saved to <b>{}</b>.</p>'.format(rgnfitsfile)
    else:
        tab2_Div_LinkImg_plot.text = '<p><b>Warning:</b> select a region first.</p>'


def tab2_panel_exit():
    tab2_panel2_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    for ll in ports:
        if platform == "linux" or platform == "linux2":
            os.system('fuser -n tcp -k {}'.format(ll))
        elif platform == "darwin":
            os.system(
                'port=($(lsof -i tcp:{}|grep python2.7 |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
            os.system(
                'port=($(lsof -i tcp:{}|grep Google |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
        print 'port {} killed'.format(ll)
    raise SystemExit


def tab2_vla_region_select(attrname, old, new):
    global tapPointDF_vla, SRC_vla_quadselround
    global x0, x1, y0, y1
    global x0pix, x1pix, y0pix, y1pix
    SRC_vla_quadselround += 1
    if SRC_vla_quadselround == 2:
        vla_quadx_selected = SRC_vla_quadx.selected['1d']['indices']
        vla_quady_selected = SRC_vla_quady.selected['1d']['indices']
        SRC_vla_quadx.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        SRC_vla_quady.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        if len(vla_quadx_selected) > 0 and len(vla_quady_selected) > 0:
            tapPointDF_vla = tapPointDF_vla.append(
                pd.Series({'xx': mapx[vla_quady_selected[0], vla_quadx_selected[0]],
                           'yy': mapy[vla_quady_selected[0], vla_quadx_selected[0]],
                           'idx': vla_quadx_selected[0], 'idy': vla_quady_selected[0]}),
                ignore_index=True)
            if len(tapPointDF_vla.index) == 1:
                r_vla_line.data_source.data = {
                    'xs': [[mapx[vla_quady_selected[0], vla_quadx_selected[0]]] * 2,
                           [mapx[vla_quady_selected[0], 0],
                            mapx[vla_quady_selected[0], -1]]],
                    'ys': [[mapy[0, vla_quadx_selected[0]],
                            mapy[-1, vla_quadx_selected[0]]],
                           [mapy[vla_quady_selected[0], vla_quadx_selected[0]]] * 2]}
            elif len(tapPointDF_vla.index) == 2:
                x0, x1 = tapPointDF_vla['xx'].min(), tapPointDF_vla['xx'].max()
                y0, y1 = tapPointDF_vla['yy'].min(), tapPointDF_vla['yy'].max()
                if x1 > x0 + mapx[0, 1] - mapx[0, 0] and y1 > y0 + mapy[0, 1] - \
                        mapy[0, 0]:
                    ## select at least 4 points
                    tab2_r_vla_ImgRgn_patch.data_source.data = {'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]}
                    x0pix, x1pix = int(tapPointDF_vla['idx'].min()), int(tapPointDF_vla['idx'].max())
                    y0pix, y1pix = int(tapPointDF_vla['idy'].min()), int(tapPointDF_vla['idy'].max())
                    tab2_tImfit_Param_dict['box'] = "'{},{},{},{}'".format(x0pix, y0pix, x1pix, y1pix)
                    if tab2_tImfit_Param_dict['gaussfit'] == 'True':
                        tab2_BUT_tImfit.label = 'pimfit'
                        tab2_maxfit_checkbox.active = []
                        tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
                        on an image region(s)</b></p>' + ' '.join(
                            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
                    else:
                        tab2_tImfit_Param_dict['gaussfit'] = 'False'
                        tab2_BUT_tImfit.label = 'pmaxfit'
                        tab2_maxfit_checkbox.active = [0]
                        tab2_Div_tImfit_text = '<p><b>#  maxfit :: do one parabolic fit components \
                        on an image region(s)</b></p>' + ' '.join(
                            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
                    tab2_Div_tImfit.text = tab2_Div_tImfit_text
                tapPointDF_vla = pd.DataFrame({'xx': [], 'yy': [], 'idx': [], 'idy': []})
                r_vla_line.data_source.data = {'xs': [], 'ys': []}
        else:
            tab2_r_vla_ImgRgn_patch.data_source.data = ColumnDataSource(
                pd.DataFrame({'xx': [], 'yy': []})).data
        SRC_vla_quadselround = 0


def Domaxfit(new):
    global tab2_tImfit_Param_dict
    if len(tab2_maxfit_checkbox.active) == 0:
        tab2_tImfit_Param_dict['gaussfit'] = 'True'
        tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
        tab2_BUT_tImfit.label = 'pimfit'
    else:
        tab2_tImfit_Param_dict['gaussfit'] = 'False'
        tab2_Div_tImfit_text = '<p><b>#  maxfit :: do one parabolic fit components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
        tab2_BUT_tImfit.label = 'pmaxfit'
    tab2_Div_tImfit.text = tab2_Div_tImfit_text


def tab2_BUT_tImfit_param_add():
    tab2_Div_tImfit2.text = ' '
    txts = tab2_input_tImfit.value.strip()
    txts = txts.split(';')
    for txt in txts:
        txt = txt.strip()
        txt = txt.split('=')
        if len(txt) == 2:
            key, val = txt
            tab2_tImfit_Param_dict[key.strip()] = val.strip()
        else:
            tab2_Div_tImfit2.text = '<p>Input syntax: <b>stokes</b>="LL"; \
            <b>ncpu</b>=10; Any spaces will be ignored.</p>'
    if tab2_tImfit_Param_dict['gaussfit'] == 'True':
        tab2_maxfit_checkbox.active = []
        tab2_BUT_tImfit.label = 'pimfit'
        tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    else:
        tab2_tImfit_Param_dict['gaussfit'] = 'False'
        tab2_maxfit_checkbox.active = [0]
        tab2_BUT_tImfit.label = 'pmaxfit'
        tab2_Div_tImfit_text = '<p><b>#  maxfit :: do one parabolic fit components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    tab2_Div_tImfit.text = tab2_Div_tImfit_text


def tab2_BUT_tImfit_param_delete():
    global tab2_tImfit_Param_dict
    tab2_Div_tImfit2.text = ' '
    txts = tab2_input_tImfit.value.strip()
    txts = txts.split(';')
    for key in txts:
        try:
            if key in ['gaussfit', 'event_id', 'struct_id', 'clean_id', 'imfit_id', 'imagefiles']:
                pass
            else:
                tab2_tImfit_Param_dict.pop(key)
        except:
            tab2_Div_tImfit2.text = '<p>Input syntax: <b>stokes</b>; <b>ncpu</b>; ' \
                                    '<b>region</b>. Any spaces will be ignored.</p>'
    tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
    on an image region(s)</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    tab2_Div_tImfit.text = tab2_Div_tImfit_text


def tab2_BUT_tImfit_param_default():
    global tab2_tImfit_Param_dict
    tab2_tImfit_Param_dict = OrderedDict()
    vlafileliststr = "'" + "','".join(vlafile) + "'"
    tab2_tImfit_Param_dict['gaussfit'] = "False"
    tab2_tImfit_Param_dict['event_id'] = "'{}'".format(event_id.replace("/", ""))
    tab2_tImfit_Param_dict['struct_id'] = "'{}'".format(struct_id.replace("/", ""))
    tab2_tImfit_Param_dict['clean_id'] = "'{}'".format(CleanID.replace("/", ""))
    tab2_tImfit_Param_dict['imfit_id'] = "'{}'".format(DButil.getcurtimstr(prefix='ImfitID_', suffix=''))
    tab2_tImfit_Param_dict['ncpu'] = "10"
    tab2_tImfit_Param_dict['box'] = "''"
    tab2_tImfit_Param_dict['width'] = "5"
    tab2_tImfit_Param_dict['getcentroid'] = "False"
    # tab2_tImfit_Param_dict['stokes'] = "'{}'".format(tab2_Select_vla_pol.value)
    tab2_tImfit_Param_dict['mask'] = "''"
    tab2_tImfit_Param_dict['imagefiles'] = "[{}]".format(vlafileliststr)
    tab2_maxfit_checkbox.active = [0]
    tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
    on an image region(s)</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    tab2_Div_tImfit.text = tab2_Div_tImfit_text
    tab2_Div_tImfit2.text = '<p><b>Default parameter Restored.</b></p>'


def tab2_BUT_timfit_param_reload():
    global tab2_tImfit_Param_dict
    infile = CleanID_dir + 'CASA_imfit_args.json'
    try:
        tab2_tImfit_Param_dict = DButil.loadjsonfile(infile)
        if tab2_tImfit_Param_dict['gaussfit'] == 'True':
            tab2_BUT_tImfit.label = 'pimfit'
            tab2_maxfit_checkbox.active = []
            tab2_Div_tImfit_text = '<p><b>#  imfit :: Fit one or more elliptical Gaussian components \
            on an image region(s)</b></p>' + ' '.join(
                "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
        else:
            tab2_BUT_tImfit.label = 'pmaxfit'
            tab2_maxfit_checkbox.active = [0]
            tab2_Div_tImfit_text = '<p><b>#  maxfit :: do one parabolic fit components \
            on an image region(s)</b></p>' + ' '.join(
                "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
        tab2_Div_tImfit.text = tab2_Div_tImfit_text
        tab2_Div_tImfit2.text = '<p>CASA pimfit arguments reload from config file in <b>{}</b>.</p>'.format(CleanID_dir)
    except:
        tab2_Div_tImfit2.text = '<p>{} not found!!!</p>'.format(infile)


def tab2_BUT_tImfit_update():
    global ImfitID_dir
    exec ('imfit_id= {}'.format(tab2_tImfit_Param_dict['imfit_id']))
    ImfitID_dir = CleanID_dir + imfit_id + '/'
    if not os.path.exists(ImfitID_dir):
        os.makedirs(ImfitID_dir)
    outfile = ImfitID_dir + 'CASA_imfit_args.json'
    DButil.updatejsonfile(outfile, tab2_tImfit_Param_dict)
    os.system('cp {}/DataBrowser/FSview/script_imfit.py {}'.format(suncasa_dir,
                                                                   ImfitID_dir))
    tab2_Div_tImfit2.text = '<p>CASA pimfit script and arguments config\
     file saved to <b>{}</b>.</p>'.format(ImfitID_dir)
    cwd = os.getcwd()
    # try:
    tab2_Div_tImfit2.text = '<p>CASA imfit script and arguments config file saved to <b>{}.</b></p>\
    <p>CASA imfit is <b>in processing</b>.</p>'.format(ImfitID_dir)
    os.chdir(ImfitID_dir)
    exec ('gaussfit = {}'.format(tab2_tImfit_Param_dict['gaussfit']))
    exec ('getcentroid = {}'.format(tab2_tImfit_Param_dict['getcentroid']))
    if gaussfit:
        suncasapy = config_main['core']['casapy47']
    else:
        suncasapy = config_main['core']['casapy46']
    suncasapy = os.path.expandvars(suncasapy)
    print suncasapy
    os.system('{} -c script_imfit.py'.format(suncasapy))
    with open(ImfitID_dir + '/CASA_imfit_out', 'rb') as f:
        out = pickle.load(f)

    dspecDF2 = DButil.transfitdict2DF(out, gaussfit=gaussfit, getcentroid=getcentroid)
    with open(CleanID_dir + '/dspecDF-base', 'rb') as fp:
        dspecDF1 = pickle.load(fp)
    for ll in dspecDF1.index:
        tmp = dspecDF1.loc[ll, 'freq']
        dspecDF1.loc[ll, 'freq'] = float('{:.3f}'.format(tmp))
    dspecDF = pd.merge(dspecDF1, dspecDF2, how='left', on=['freqstr', 'fits_local'])
    with open(ImfitID_dir + '/dspecDF-save', 'wb') as fp:
        pickle.dump(dspecDF, fp)
    print 'imfit results saved to ' + ImfitID_dir + '/dspecDF-save'
    tab2_Div_tImfit2.text = '<p>imfit finished, go back to <b>QLook</b> \
    window, select StrID <b>{}</b> and click <b>FSview</b> button again.</p>'.format(struct_id[0:-1])
    # except:
    #     tab2_Div_tImfit2.text = '<p>CASA imfit script and arguments config file \
    #     saved to <b>{}.</b></p><p>Do imfit with CASA manually.</p>'.format(
    #         ImfitID_dir) + '<p>When finished, go back to <b>QLook</b> window, \
    #         select StrID <b>{}</b> and click <b>FSview</b> button again.</p>'.format(struct_id[0:-1])
    os.chdir(cwd)


event_id = config_EvtID['datadir']['event_id']
event_dir = database_dir + event_id
try:
    infile = event_dir + 'CurrFS.json'
    FS_config = DButil.loadjsonfile(infile)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']
struct_dir = database_dir + event_id + struct_id
CleanID = FS_config['datadir']['clean_id']
CleanID_dir = struct_dir + CleanID
FS_dspecDF = CleanID_dir + 'dspecDF-base'
FS_specfile = FS_config['datadir']['FS_specfile']
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

clickmode = 'singleclick'
SRC_dspec_quadselround = 0
SRC_vla_quadselround = 0
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
fits_LOCL = config_EvtID['datadir']['fits_LOCL']
fits_GLOB = config_EvtID['datadir']['fits_GLOB']
fits_LOCL_dir = CleanID_dir + fits_LOCL
fits_GLOB_dir = CleanID_dir + fits_GLOB

if os.path.exists(FS_dspecDF):
    with open(FS_dspecDF, 'rb') as f:
        dspecDF0 = pickle.load(f)
    vlafile = glob.glob(fits_LOCL_dir + '*.fits')
    if len(vlafile) > 0:
        tab2_panel2_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>\
                                first before closing the tab</p></b>""",
                                   width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()

        dspecDF_select = dspecDF0
        '''create the regridded dynamic spectrum plot'''
        TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
        # downsample_dspecDF(spec_square_rs_tmax=spec_square_rs_tmax, spec_square_rs_fmax=spec_square_rs_fmax)
        dspecDF0_rs = dspecDF0.copy()
        '''create the dynamic spectrum plot'''
        TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
        tab2_SRC_dspec_square = ColumnDataSource(dspecDF0_rs)
        tab2_p_dspec = figure(tools=TOOLS, webgl=config_main['plot_config']['WebGL'],
                              plot_width=config_main['plot_config']['tab_FSview_base']['dspec_wdth'],
                              plot_height=config_main['plot_config']['tab_FSview_base']['dspec_hght'],
                              x_range=(tab2_dtim[0] - tab2_dt / 2.0, tab2_dtim[-1] + tab2_dt / 2.0),
                              y_range=(tab2_freq[0] - tab2_df / 2.0, tab2_freq[-1] + tab2_df / 2.0),
                              toolbar_location="above")
        tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        tab2_p_dspec.axis.visible = True
        tab2_p_dspec.title.text = "Dynamic spectrum"
        tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
        tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
        tab2_r_dspec = tab2_p_dspec.image(image=[tab2_spec_plt], x=tab2_dtim[0] - tab2_dt / 2.0,
                                          y=tab2_freq[0] - tab2_df / 2.0,
                                          dw=tab2_dur + tab2_dt,
                                          dh=tab2_freq[-1] - tab2_freq[0] + tab2_df, palette=bokehpalette_jet)

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
                                                fill_color=None, fill_alpha=0.0, line_color="Magenta",
                                                line_alpha=0.8, line_width=1)

        tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                line_color='white', source=tab2_source_idx_line_x)
        tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                line_color='white', source=tab2_source_idx_line_y)
        # tab2_p_dspec.border_fill_color = "silver"
        tab2_p_dspec.border_fill_alpha = 0.4
        tab2_p_dspec.axis.major_tick_out = 0
        tab2_p_dspec.axis.major_tick_in = 5
        tab2_p_dspec.axis.minor_tick_out = 0
        tab2_p_dspec.axis.minor_tick_in = 3
        tab2_p_dspec.axis.major_tick_line_color = "white"
        tab2_p_dspec.axis.minor_tick_line_color = "white"

        tab2_Select_pol = Select(title="Polarization:", value='I', options=['RR', 'LL', 'I', 'V'],
                                 width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_Select_bl = Select(title="Baseline:", value=tab2_bl[0], options=tab2_bl,
                                width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                                        width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_p_dspec_xPro = figure(tools='', webgl=config_main['plot_config']['WebGL'],
                                   plot_width=config_main['plot_config']['tab_FSview_base']['dspec_xPro_wdth'],
                                   plot_height=config_main['plot_config']['tab_FSview_base']['dspec_xPro_hght'],
                                   x_range=tab2_p_dspec.x_range, y_range=(spec_plt_min, spec_plt_max),
                                   title="Time profile", toolbar_location=None)
        tab2_SRC_dspec_xPro = ColumnDataSource({'x': [], 'y': []})
        tab2_SRC_dspec_xPro_hover = ColumnDataSource({'x': [], 'y': [], 'tooltips': []})
        r_dspec_xPro = tab2_p_dspec_xPro.line(x='x', y='y', alpha=1.0, line_width=1, line_color='black',
                                              source=tab2_SRC_dspec_xPro)
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
        tab2_p_dspec_xPro.xaxis.axis_label = 'Seconds since ' + tim0_char
        tab2_p_dspec_xPro.yaxis.axis_label = 'Intensity [sfu]'
        # tab2_p_dspec_xPro.border_fill_color = "silver"
        tab2_p_dspec_xPro.border_fill_alpha = 0.4
        tab2_p_dspec_xPro.axis.major_tick_out = 0
        tab2_p_dspec_xPro.axis.major_tick_in = 5
        tab2_p_dspec_xPro.axis.minor_tick_out = 0
        tab2_p_dspec_xPro.axis.minor_tick_in = 3

        tab2_p_dspec_yPro = figure(tools='', webgl=config_main['plot_config']['WebGL'],
                                   plot_width=config_main['plot_config']['tab_FSview_base']['dspec_yPro_wdth'],
                                   plot_height=config_main['plot_config']['tab_FSview_base']['dspec_yPro_hght'],
                                   x_range=(spec_plt_min, spec_plt_max), y_range=tab2_p_dspec.y_range,
                                   title="Frequency profile", toolbar_location=None)
        tab2_SRC_dspec_yPro = ColumnDataSource({'x': [], 'y': []})
        tab2_SRC_dspec_yPro_hover = ColumnDataSource({'x': [], 'y': [], 'tooltips': []})
        r_dspec_yPro = tab2_p_dspec_yPro.line(x='x', y='y', alpha=1.0, line_width=1, line_color='black',
                                              source=tab2_SRC_dspec_yPro)
        r_dspec_yPro_c = tab2_p_dspec_yPro.circle(x='x', y='y', size=5, fill_alpha=0.2, fill_color='grey',
                                                  line_color=None,
                                                  source=tab2_SRC_dspec_yPro)
        r_dspec_yPro_hover = tab2_p_dspec_yPro.circle(x='x', y='y', size=5, fill_alpha=0.5, fill_color='firebrick',
                                                      line_color='firebrick', source=tab2_SRC_dspec_yPro_hover)
        l_dspec_yPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                      source=tab2_SRC_dspec_yPro_hover,
                                      render_mode='canvas')
        l_dspec_yPro_hover.text_font_size = '10pt'
        tab2_p_dspec_yPro.add_layout(l_dspec_yPro_hover)
        tab2_p_dspec_yPro.yaxis.visible = False
        tab2_p_dspec_yPro.xaxis.axis_label = 'Intensity [sfu]'
        # tab2_p_dspec_yPro.border_fill_color = "silver"
        tab2_p_dspec_yPro.border_fill_alpha = 0.4
        tab2_p_dspec_yPro.min_border_bottom = 0
        tab2_p_dspec_yPro.min_border_left = 0
        tab2_p_dspec_yPro.axis.major_tick_out = 0
        tab2_p_dspec_yPro.axis.major_tick_in = 5
        tab2_p_dspec_yPro.axis.minor_tick_out = 0
        tab2_p_dspec_yPro.axis.minor_tick_in = 3

        tab2_BUT_vdspec = Button(label="VEC Dyn Spec",
                                 width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                 button_type="success")
        tab2_Select_pol_opt = ['RR', 'LL', 'I', 'V']

        tab2_BUT_vdspec.on_click(tab2_vdspec_update)
        tab2_ctrls = [tab2_Select_bl, tab2_Select_pol, tab2_Select_colorspace]
        for ctrl in tab2_ctrls:
            ctrl.on_change('value', tab2_update_dspec_image)

        # initial the VLA map contour source
        tab2_SRC_vlamap_contour = ColumnDataSource(
            data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})
        tab2_SRC_vlamap_peak = ColumnDataSource(
            data={'dspec': [], 'shape_longitude': [], 'shape_latitude': [], 'peak': []})

        # import the vla image
        hdu = read_fits(vlafile[0])
        hdu_goodchan = goodchan(hdu)
        vla_local_pfmap = PuffinMap(hdu.data[0, hdu_goodchan[0], :, :], hdu.header,
                                    plot_height=config_main['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_main['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_main['plot_config']['WebGL'])
        # plot the contour of vla image
        mapx, mapy = vla_local_pfmap.meshgrid(rescale=1.0)
        mapx, mapy = mapx.value, mapy.value
        mapvlasize = mapy.shape
        tab2_SRC_vlamap_contour = DButil.get_contour_data(mapx, mapy, vla_local_pfmap.smap.data)
        # ImgDF0 = pd.DataFrame({'xx': mapx.ravel(), 'yy': mapy.ravel()})
        # tab2_SRC_vla_square = ColumnDataSource(ImgDF0)
        colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
        bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        tab2_SRC_ImgRgn_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))

        # aiamap = sunpy.map.Map(filepath)
        print xx[0] / 3600. / 24., SDOdir
        aiamap = DButil.readsdofile(datadir=SDOdir, wavelength='171', jdtime=xx[0] / 3600. / 24.,
                                    timtol=tab2_dur / 3600. / 24.)
        MapRES = 256
        dimensions = u.Quantity([MapRES, MapRES], u.pixel)
        aia_resampled_map = aiamap.resample(dimensions)

        # plot the global AIA image
        aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                        plot_height=config_main['plot_config']['tab_FSview_base']['aia_hght'],
                                        plot_width=config_main['plot_config']['tab_FSview_base']['aia_wdth'],
                                        webgl=config_main['plot_config']['WebGL'])

        tab2_p_aia, tab2_r_aia = aia_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg)
        tab2_p_aia.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_aia.title.text_font_size = '6pt'
        # tab2_p_aia.border_fill_color = "silver"
        tab2_p_aia.border_fill_alpha = 0.4
        tab2_p_aia.axis.major_tick_out = 0
        tab2_p_aia.axis.major_tick_in = 5
        tab2_p_aia.axis.minor_tick_out = 0
        tab2_p_aia.axis.minor_tick_in = 3
        tab2_p_aia.axis.major_tick_line_color = "white"
        tab2_p_aia.axis.minor_tick_line_color = "white"
        tab2_p_aia.grid.grid_line_color = None
        tab2_p_aia.background_fill_color = "black"
        tab2_r_aia_ImgRgn_patch = tab2_p_aia.patch('xx', 'yy', source=tab2_SRC_ImgRgn_Patch,
                                                   fill_color=None, fill_alpha=0.5, line_color="white",
                                                   line_alpha=1, line_width=1)

        colormap = cm.get_cmap("gray")  # choose any matplotlib colormap here
        bokehpalette_sdohmimag = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        # hmimap = sunpy.map.Map(filepath)
        # todo fix this bug
        hmimap = aiamap
        # plot the global HMI image
        hmi_resampled_map = hmimap.resample(dimensions)
        hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                        plot_height=config_main['plot_config']['tab_FSview_base']['vla_hght'],
                                        plot_width=config_main['plot_config']['tab_FSview_base']['vla_wdth'],
                                        webgl=config_main['plot_config']['WebGL'])

        tab2_p_hmi, tab2_r_hmi = hmi_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                             x_range=tab2_p_aia.x_range,
                                                             y_range=tab2_p_aia.y_range,
                                                             palette=bokehpalette_sdohmimag)
        tab2_p_hmi.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_hmi.yaxis.visible = False
        # tab2_p_hmi.border_fill_color = "silver"
        tab2_p_hmi.border_fill_alpha = 0.4
        tab2_p_hmi.axis.major_tick_out = 0
        tab2_p_hmi.axis.major_tick_in = 5
        tab2_p_hmi.axis.minor_tick_out = 0
        tab2_p_hmi.axis.minor_tick_in = 3
        tab2_p_hmi.axis.major_tick_line_color = "white"
        tab2_p_hmi.axis.minor_tick_line_color = "white"
        tab2_p_hmi.grid.grid_line_color = None
        tab2_p_hmi.background_fill_color = "black"
        tab2_r_hmi_ImgRgn_patch = tab2_p_hmi.patch('xx', 'yy', source=tab2_SRC_ImgRgn_Patch,
                                                   fill_color=None, fill_alpha=0.5, line_color="white",
                                                   line_alpha=1, line_width=1)
        # plot the vla image
        tab2_p_vla, tab2_r_vla = vla_local_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                         palette=bokehpalette_SynthesisImg,
                                                         x_range=tab2_p_aia.x_range,
                                                         y_range=tab2_p_aia.y_range)
        tab2_p_vla.title.text_font_size = '6pt'
        tab2_p_vla.yaxis.visible = False
        # tab2_p_vla.border_fill_color = "silver"
        tab2_p_vla.border_fill_alpha = 0.4
        tab2_p_vla.axis.major_tick_out = 0
        tab2_p_vla.axis.major_tick_in = 5
        tab2_p_vla.axis.minor_tick_out = 0
        tab2_p_vla.axis.minor_tick_in = 3
        tab2_p_vla.axis.major_tick_line_color = "white"
        tab2_p_vla.axis.minor_tick_line_color = "white"
        tab2_p_vla.grid.grid_line_color = None
        tab2_p_vla.background_fill_color = "black"

        tapPointDF_vla = pd.DataFrame(
            {'xx': [], 'yy': [], 'idx': [], 'idy': []})  ## the selected point to fit in aia submap
        dx = np.mean(np.diff(mapx[0, :]))
        dy = np.mean(np.diff(mapy[:, 0]))
        xLFE = mapx[0, :] - dx / 2.0
        xRTE = np.append(xLFE[1:], xLFE[-1] + dx)
        yBTE = mapy[:, 0] - dy / 2.0
        yTPE = np.append(yBTE[1:], yBTE[-1] + dy)
        SRC_vla_quadx = ColumnDataSource(
            {'left': xLFE.ravel(), 'right': xRTE.ravel(), 'bottom': np.repeat(yBTE[0], mapvlasize[1]),
             'top': np.repeat(yBTE[-1] + dy, mapvlasize[1])})
        SRC_vla_quady = ColumnDataSource(
            {'left': np.repeat(xLFE[0], mapvlasize[0]), 'right': np.repeat(xLFE[-1] + dx, mapvlasize[0]),
             'bottom': yBTE.ravel(),
             'top': yTPE.ravel()})
        r_vla_quadx = tab2_p_vla.quad('left', 'right', 'top', 'bottom', source=SRC_vla_quadx,
                                      fill_alpha=0.0, fill_color=None,
                                      line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                      selection_fill_color=None,
                                      nonselection_fill_alpha=0.0, nonselection_fill_color=None,
                                      selection_line_alpha=0.0, selection_line_color=None,
                                      nonselection_line_alpha=0.0)
        r_vla_quady = tab2_p_vla.quad('left', 'right', 'top', 'bottom', source=SRC_vla_quady,
                                      fill_alpha=0.0, fill_color=None,
                                      line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                      selection_fill_color=None,
                                      nonselection_fill_alpha=0.0, nonselection_fill_color=None,
                                      selection_line_alpha=0.0, selection_line_color=None,
                                      nonselection_line_alpha=0.0)
        tab2_p_vla.add_tools(TapTool(renderers=[r_vla_quadx, r_vla_quady]))
        SRC_vla_quadx.on_change('selected', tab2_vla_region_select)
        SRC_vla_quady.on_change('selected', tab2_vla_region_select)

        tab2_r_vla_multi_line = tab2_p_vla.multi_line(xs='xs', ys='ys', line_color='line_color',
                                                      source=tab2_SRC_vlamap_contour, alpha=0.7, line_width=2)

        tab2_r_vla_ImgRgn_patch = tab2_p_vla.patch('xx', 'yy', source=tab2_SRC_ImgRgn_Patch,
                                                   fill_color=None, fill_alpha=0.5, line_color="white",
                                                   line_alpha=1, line_width=1)
        SRC_vla_lines = ColumnDataSource({'xs': [], 'ys': []})
        r_vla_line = tab2_p_vla.multi_line('xs', 'ys', source=SRC_vla_lines, line_color='cyan', line_width=1, alpha=0.8)

        tab2_LinkImg_HGHT = config_main['plot_config']['tab_FSview_base']['vla_hght']
        tab2_LinkImg_WDTH = config_main['plot_config']['tab_FSview_base']['vla_wdth']

        tab2_Div_LinkImg_plot = Div(text=""" """,
                                    width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_Slider_time_LinkImg = Slider(start=0, end=tab2_ntim - 1, value=0, step=1, title="time idx",
                                          width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'] * 2,
                                          callback_policy='mouseup')
        tab2_Slider_freq_LinkImg = Slider(start=0, end=tab2_nfreq - 1, value=0, step=1, title="freq idx",
                                          width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'] * 2,
                                          callback_policy='mouseup')

        pols = DButil.polsfromfitsheader(hdu.header)
        # pols = ['RR', 'LL', 'I', 'V']
        SRL = set(['RR', 'LL'])
        SXY = set(['XX', 'YY', 'XY', 'YX'])
        Spol = set(pols)
        if hdu.header['NAXIS4'] == 2 and len(SRL.intersection(Spol)) == 2:
            pols = pols + ['I', 'V']
        if hdu.header['NAXIS4'] == 4 and len(SXY.intersection(Spol)) == 4:
            pols = pols + ['I', 'V']

        tab2_Select_vla_pol = Select(title="Polarization:", value=pols[0], options=pols,
                                     width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_CTRLs_LinkImg = [tab2_Slider_time_LinkImg, tab2_Slider_freq_LinkImg, tab2_Select_vla_pol]
        for ctrl in tab2_CTRLs_LinkImg:
            ctrl.on_change('value', tab3_slider_LinkImg_update)

        tab2_BUT_XCorr = Button(label='XCorr',
                                width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                button_type='warning')
        tab2_BUT_XCorr.on_click(tab2_panel_XCorr_update)

        tab2_panel2_BUT_exit = Button(label='Exit FSview',
                                      width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                      button_type='danger')
        tab2_panel2_BUT_exit.on_click(tab2_panel_exit)

        # tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)
        # tab2_vla_square_selected = None
        # tab2_SRC_vla_square.on_change('selected', tab2_vla_region_select)
        tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                    options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                             "2048x2048"],
                                    width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_Select_MapRES.on_change('value', tab2_update_MapRES)
        rgnfitsfile = CleanID_dir + "CASA_imfit_region_fits.rgn"
        tab2_BUT_SavRgn = Button(label='Save Region',
                                 width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                 button_type='primary')
        tab2_BUT_SavRgn.on_click(tab2_save_region)

        tab2_input_tImfit = TextInput(value="Input the param here", title="pimfit task parameters:",
                                      width=config_main['plot_config']['tab_FSviewPrep']['input_tCLN_wdth'])
        tab2_Div_tImfit = Div(text='', width=config_main['plot_config']['tab_FSviewPrep']['tab2_Div_tImfit_wdth'])
        tab2_Div_tImfit2 = Div(text='', width=config_main['plot_config']['tab_FSviewPrep']['input_tCLN_wdth'])
        tab2_maxfit_checkbox = CheckboxGroup(labels=["parabolic fit"], active=[0])
        tab2_maxfit_checkbox.on_click(Domaxfit)

        tab2_BUT_tImfit_param_default()
        tab2_BUT_tImfit_param_ADD = Button(label='Add to Param',
                                           width=config_main['plot_config']['tab_FSviewPrep']['button_wdth'],
                                           button_type='primary')
        tab2_BUT_tImfit_param_ADD.on_click(tab2_BUT_tImfit_param_add)
        tab2_BUT_tImfit_param_DEL = Button(label='Delete Param',
                                           width=config_main['plot_config']['tab_FSviewPrep']['button_wdth'],
                                           button_type='warning')
        tab2_BUT_tImfit_param_DEL.on_click(tab2_BUT_tImfit_param_delete)
        tab2_BUT_tImfit_param_Default = Button(label='Default Param',
                                               width=config_main['plot_config']['tab_FSviewPrep']['button_wdth'])
        tab2_BUT_tImfit_param_Default.on_click(tab2_BUT_tImfit_param_default)
        tab2_SPCR_LFT_BUT_tImfit_param_DEL = Spacer(
            width=config_main['plot_config']['tab_FSviewPrep']['space_wdth10'])
        tab2_SPCR_LFT_BUT_tImfit_param_RELOAD = Spacer(
            width=config_main['plot_config']['tab_FSviewPrep']['space_wdth10'])
        tab2_SPCR_LFT_BUT_tImfit_param_DEFAULT = Spacer(
            width=config_main['plot_config']['tab_FSviewPrep']['space_wdth10'])

        tab2_BUT_tImfit_param_RELOAD = Button(label='reload Param',
                                              width=config_main['plot_config']['tab_FSviewPrep']['button_wdth'])
        tab2_BUT_tImfit_param_RELOAD.on_click(tab2_BUT_timfit_param_reload)
        tab2_SPCR_LFT_BUT_tImfit_param_reload = Spacer(
            width=config_main['plot_config']['tab_FSviewPrep']['space_wdth20'])

        tab2_BUT_tImfit = Button(label='imfit',
                                 width=config_main['plot_config']['tab_FSviewPrep']['button_wdth'],
                                 button_type='success')
        tab2_BUT_tImfit.on_click(tab2_BUT_tImfit_update)
        tab2_SPCR_ABV_BUT_timfit = Spacer(width=config_main['plot_config']['tab_FSviewPrep']['space_wdth10'])

        tab2_SPCR_LFT_Div_tImfit = Spacer(width=config_main['plot_config']['tab_FSviewPrep']['space_wdth10'])

        But_ClickMode = Button(label="ClickMode: Double",
                               width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                               button_type="primary")
        But_ClickMode.on_click(tab2_ClickMode_handler)

        lout2_1_1 = row(gridplot([[tab2_p_aia, tab2_p_hmi, tab2_p_vla]], toolbar_location='right'),
                        widgetbox(tab2_Select_MapRES, tab2_Select_vla_pol, tab2_Slider_time_LinkImg,
                                  tab2_Slider_freq_LinkImg, tab2_BUT_vdspec, tab2_BUT_SavRgn, tab2_Div_LinkImg_plot,
                                  width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'] * 2))
        lout2_1_2 = row(column(row(tab2_p_dspec, tab2_p_dspec_yPro),
                               tab2_p_dspec_xPro),
                        widgetbox(tab2_Select_pol, tab2_Select_bl,
                                  tab2_Select_colorspace, But_ClickMode, tab2_BUT_XCorr,
                                  tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                                  width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth']))
        lout2_1 = column(lout2_1_1, lout2_1_2)
        lout2_2 = tab2_SPCR_LFT_Div_tImfit
        lout2_3 = column(widgetbox(tab2_input_tImfit, tab2_maxfit_checkbox, tab2_Div_tImfit2,
                                   width=config_main['plot_config']['tab_FSview_base']['widgetbox_wdth'] * 2),
                         row(tab2_BUT_tImfit_param_ADD, tab2_SPCR_LFT_BUT_tImfit_param_DEL,
                             tab2_BUT_tImfit_param_DEL, tab2_SPCR_LFT_BUT_tImfit_param_DEFAULT,
                             tab2_BUT_tImfit_param_Default, tab2_SPCR_LFT_BUT_tImfit_param_RELOAD,
                             tab2_BUT_tImfit_param_RELOAD, tab2_SPCR_ABV_BUT_timfit, tab2_BUT_tImfit),
                         tab2_Div_tImfit)
        lout = row(lout2_1, lout2_2, lout2_3)

        curdoc().add_root(lout)
        curdoc().title = "FSview"
    else:
        tab_panel_Div_info = Div(
            text="""<p><b>Warning</b>: Synthesis images not found!!!</p>""",
            width=config_main['plot_config']['tab_FSview_base']['dspec_wdth'])
        lout = tab_panel_Div_info
        curdoc().add_root(lout)
        curdoc().title = "FSview"
else:
    raise SystemExit
