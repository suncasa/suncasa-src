import json
import os, sys
import pickle
import time
from collections import OrderedDict
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.misc import bytescale
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import scipy.ndimage as sn
import sunpy.map
from math import radians, cos, sin
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Slider, Button, TextInput, RadioButtonGroup, CheckboxGroup,
                          BoxSelectTool, LassoSelectTool, HoverTool, Spacer, LabelSet, Div)
from bokeh.models.mappers import LinearColorMapper
from bokeh.models.widgets import Select
from bokeh.palettes import Spectral11
from bokeh.plotting import figure, curdoc
import glob
from astropy.time import Time
from QLook_util import get_contour_data
from puffin import PuffinMap

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
'''load config file'''
with open(suncasa_dir + 'DataBrowser/config.json', 'r') as fp:
    config_plot = json.load(fp)
database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/'
spec_rs_tmax = config_plot['plot_config']['tab_FSview_base']['spec_rs_tmax']
spec_rs_fmax = config_plot['plot_config']['tab_FSview_base']['spec_rs_fmax']
spec_square_rs_ratio = config_plot['plot_config']['tab_FSview_base']['spec_square_rs_ratio']
dspec_fs_tmax = config_plot['plot_config']['tab_FSview_base']['dspec_fs_tmax']
dspec_fs_fmax = config_plot['plot_config']['tab_FSview_base']['dspec_fs_fmax']
with open('{}config_EvtID_curr.json'.format(database_dir), 'r') as fp:
    config_EvtID = json.load(fp)

do_spec_regrid = False

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
    return hdu


def goodchan(hdu):
    ndx = hdu.header["NAXIS1"]
    ndy = hdu.header["NAXIS2"]
    xc = ndx / 2
    yc = ndy / 2
    hdu_goodchan = \
        np.where(np.nanmean(hdu.data[0, :, yc - ndy / 16:yc + ndy / 16, xc - ndx / 16:xc + ndx / 16], axis=(-1, -2)))[0]
    return hdu_goodchan


def sdo_aia_scale(image=None, wavelength=None):
    if wavelength == '94':
        image[image > 30] = 30
        image[image < 0.5] = 0.5
        image = np.log10(image)
    elif wavelength == '131':
        image[image > 200] = 200
        image[image < 2] = 2
        image = np.log10(image)
    elif wavelength == '171':
        image[image > 2000] = 2000
        image[image < 50] = 50
        image = np.log10(image)
    return bytescale(image)


def sdomapfromlocalfile(wavelength=None, jdtime=None):
    aiafitspath = glob.glob(database_dir + event_id + '/AIA/aia_lev1_{}a*.fits'.format(wavelength))
    aiafits = [ll.split('/')[-1] for ll in aiafitspath]
    aiatimeline = [ll.replace('aia_lev1_{}a_'.format(wavelength), '').replace('.fits', '').replace('z_image_lev1', '')
                   for
                   ll in aiafits]
    aiatimeline = [ll.split('t')[0].replace('_', '-') + ' ' + ll.split('t')[1].replace('_', ':') for ll in
                   aiatimeline]
    aiatimeline = [ll[0:ll.rindex(':')] + '.' + ll[(ll.rindex(':') + 1):] for ll in aiatimeline]
    aiatimeline = Time(aiatimeline, format='iso', scale='utc')
    idxaia = np.argmin(np.abs(aiatimeline.jd - jdtime))
    filepath = aiafitspath[idxaia]
    aiamap = sunpy.map.Map(filepath)
    aiamap.data = sdo_aia_scale(image=aiamap.data / aiamap.exposure_time.value, wavelength=wavelength)
    return aiamap


def rebin_specdata(tab2_spec, bl_index, select_pol):
    global tab2_tim_square_rs, tab2_freq_square_rs, tab2_ntim_square_rs, tab2_nfreq_square_rs
    tab2_spec_sz = tab2_spec.shape
    spec_sz2, spec_sz1 = 10, 10
    if tab2_spec_sz[3] > spec_rs_tmax:
        spec_sz2 = next(i for i in xrange(1, 11) if i / 10. * tab2_spec_sz[3] > spec_rs_tmax)
    if tab2_spec_sz[2] > spec_rs_fmax:
        spec_sz1 = next(i for i in xrange(1, 11) if i / 10. * tab2_spec_sz[2] > spec_rs_fmax)
    spec_sz1_rs, spec_sz2_rs = spec_sz1 / 10.0, spec_sz2 / 10.0
    tab2_tim_square_rs = sn.interpolation.zoom(tab2_tim, spec_sz2_rs / spec_square_rs_ratio, order=1)
    tab2_freq_square_rs = sn.interpolation.zoom(tab2_freq, spec_sz1_rs / spec_square_rs_ratio, order=1)
    tab2_ntim_square_rs = len(tab2_tim_square_rs)
    tab2_nfreq_square_rs = len(tab2_freq_square_rs)


def tab2_vdspec_update():
    global tab2_Select_pol_opt, spec_pol_dict
    select_pol = tab2_Select_pol.value
    tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
    if tab2_BUT_vdspec.label == "VEC Dyn Spec":
        if tab2_vla_square_selected:
            tab2_Div_LinkImg_plot.text = '<p><b>Vector dynamic spectrum in calculating...</b></p>'
            tab2_Select_pol_opt = tab2_Select_pol.options
            tab2_Select_pol.options = pols
            tab2_BUT_vdspec.label = "Dyn Spec"
            idxmax = max(tab2_vla_square_selected)
            idxmin = min(tab2_vla_square_selected)
            x0pix, x1pix = idxmin % mapvlasize[0], idxmax % mapvlasize[0]
            y0pix, y1pix = idxmin / mapvlasize[0], idxmax / mapvlasize[0]
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
                        idxfreq = freq.index(freq_ref)
                        vla_l = hdu.data[0, :, y0pix:y1pix + 1, x0pix:x1pix + 1]
                        vla_r = hdu.data[1, :, y0pix:y1pix + 1, x0pix:x1pix + 1]
                        spec_plt_R[idxfreq:idxfreq + nfreq_hdu, ll] = \
                            np.nanmean(vla_l, axis=(-1, -2))[hdu_goodchan[0]:hdu_goodchan[-1] + 1]
                        spec_plt_L[idxfreq:idxfreq + nfreq_hdu, ll] = \
                            np.nanmean(vla_r, axis=(-1, -2))[hdu_goodchan[0]:hdu_goodchan[-1] + 1]
                spec_plt_R[spec_plt_R < 0] = 0
                spec_plt_L[spec_plt_L < 0] = 0
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
                spec_plt_R[spec_plt_R < 0] = 0
                spec_plt_L = spec_plt_R
            tab2_Div_LinkImg_plot.text = '<p><b>Vector dynamic spectrum calculated.</b></p>'
            spec_pol_dict = make_spec_plt(spec_plt_R, spec_plt_L)
            tab2_r_square_rs_selection_change(None, select_pol)
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
        spec_plt_pol = spec_pol_dict['spec']
        tab2_r_square_rs_selection_change(bl_index, select_pol)
        tab2_p_dspec.title.text = "Dynamic spectrum"
        tab2_Div_LinkImg_plot.text = ''


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


def tab2_update_dspec_image(attrname, old, new):
    global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
    select_pol = tab2_Select_pol.value
    select_bl = tab2_Select_bl.value
    bl_index = tab2_bl.index(select_bl)
    if tab2_BUT_vdspec.label == "VEC Dyn Spec":
        tab2_r_square_rs_selection_change(bl_index, select_pol)


def tab2_update_dspec_rs_image(attrname, old, new):
    select_pol = tab2_Select_pol.value
    select_bl = tab2_Select_bl.value
    bl_index = tab2_bl.index(select_bl)
    tab2_r_square_rs_selection_change(bl_index, select_pol)


def tab2_r_square_rs_selection_change(bl_index, select_pol):
    global dspecDF_frac, spec_pol_dict, dspecDF_rs0
    global tab2_dtim_fs, tab2_freq_fs, tab2_tim_ind0, tab2_tim_ind1, tab2_freq_ind0, tab2_freq_ind1
    global tab2_SRC_dspec_image, tab2_SRC_dspec_square
    tab2_SRC_dspec_image_rs.data = {'data': [spec_pol_dict['spec'][select_pol]], 'xx': [tab2_dtim],
                                    'yy': [tab2_freq]}
    tab2_r_square_rs_selected = tab2_SRC_dspec_square_rs.selected['1d']['indices']
    if tab2_r_square_rs_selected:
        dspecDF_rs = dspecDF_rs0.iloc[tab2_r_square_rs_selected, :]
        trs0, trs1 = dspecDF_rs['xx'].min(), dspecDF_rs['xx'].max()
        frs0, frs1 = dspecDF_rs['yy'].min(), dspecDF_rs['yy'].max()
        if trs1 > trs0 + dspec_fs_tmax * tab2_dt:
            trs1 = trs0 + dspec_fs_tmax * tab2_dt
        if frs1 > frs0 + dspec_fs_fmax * tab2_df:
            frs1 = frs0 + dspec_fs_fmax * tab2_df
        tab2_r_square_rs_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [trs0, trs1, trs1, trs0], 'yy': [frs0, frs0, frs1, frs1]})).data
        dspecDF_frac = dspecDF0[dspecDF0.time < trs1][dspecDF0.time >= trs0][dspecDF0.freq >= frs0][
            dspecDF0.freq < frs1]
        tab2_dtim_fs = pd.Series.unique(dspecDF_frac['time'])
        tab2_freq_fs = pd.Series.unique(dspecDF_frac['freq'])
        tab2_ntim_fs = len(tab2_dtim_fs)
        tab2_nfreq_fs = len(tab2_freq_fs)
        tab2_tim_ind0 = np.where(abs(tab2_dtim - tab2_dtim_fs[0]) < tab2_dt / 2.0)[0][0]
        tab2_tim_ind1 = np.where(abs(tab2_dtim - tab2_dtim_fs[-1]) < tab2_dt / 2.0)[0][0]
        tab2_freq_ind0 = np.where(abs(tab2_freq - tab2_freq_fs[0]) < tab2_df / 2.0)[0][0]
        tab2_freq_ind1 = np.where(abs(tab2_freq - tab2_freq_fs[-1]) < tab2_df / 2.0)[0][0]
        print tab2_tim_ind0, tab2_tim_ind1, tab2_freq_ind0, tab2_freq_ind1
        global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
        if tab2_BUT_vdspec.label == "VEC Dyn Spec":
            spec_plt_R_frac = tab2_spec[0, bl_index, tab2_freq_ind0:(tab2_freq_ind1 + 1),
                              tab2_tim_ind0:(tab2_tim_ind1 + 1)]
            spec_plt_L_frac = tab2_spec[1, bl_index, tab2_freq_ind0:(tab2_freq_ind1 + 1),
                              tab2_tim_ind0:(tab2_tim_ind1 + 1)]
        else:
            spec_plt_R_frac = spec_pol_dict['spec'][select_pol][tab2_freq_ind0:(tab2_freq_ind1 + 1),
                              tab2_tim_ind0:(tab2_tim_ind1 + 1)]
            spec_plt_L_frac = spec_pol_dict['spec'][select_pol][tab2_freq_ind0:(tab2_freq_ind1 + 1),
                              tab2_tim_ind0:(tab2_tim_ind1 + 1)]

        spec_pol_dict_frac = make_spec_plt(spec_plt_R_frac, spec_plt_L_frac)
        if select_pol == 'V':
            tab2_Select_colorspace.value = 'linear'
        if tab2_Select_colorspace.value == 'log' and select_pol != 'V':
            tab2_SRC_dspec_image.data = {'data': [np.log(spec_pol_dict_frac['spec'][select_pol])], 'xx': [tab2_dtim],
                                         'yy': [tab2_freq]}
        else:
            tab2_SRC_dspec_image.data = {'data': [spec_pol_dict_frac['spec'][select_pol]], 'xx': [tab2_dtim],
                                         'yy': [tab2_freq]}
        tab2_SRC_dspec_square.data['dspec'] = spec_pol_dict_frac['spec'][select_pol].flatten()
        tab2_p_dspec_xPro.y_range.start = spec_pol_dict_frac['min'][select_pol]
        tab2_p_dspec_xPro.y_range.end = spec_pol_dict_frac['max'][select_pol]
        tab2_p_dspec_yPro.x_range.start = spec_pol_dict_frac['min'][select_pol]
        tab2_p_dspec_yPro.x_range.end = spec_pol_dict_frac['max'][select_pol]
        tab2_SRC_dspec_square = ColumnDataSource(dspecDF_frac)
        tab2_p_dspec = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                              plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_wdth'],
                              plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_hght'],
                              x_range=(dspecDF_frac['time'].min(), dspecDF_frac['time'].max()),
                              y_range=(dspecDF_frac['freq'].min(), dspecDF_frac['freq'].max()),
                              toolbar_location="above")
        tab2_p_dspec.x_range.start =dspecDF_frac['time'].min()
        tab2_p_dspec.x_range.end = dspecDF_frac['time'].max()
        tab2_p_dspec.y_range.start =dspecDF_frac['freq'].min()
        tab2_p_dspec.y_range.end = dspecDF_frac['freq'].max()
        # tab2_SRC_dspec_image = ColumnDataSource(
        #     data={'data': [spec_pol_dict['spec'][select_pol][tab2_freq_ind0:(tab2_freq_ind1 + 1),
        #                    tab2_tim_ind0:(tab2_tim_ind1 + 1)]], 'xx': [tab2_dtim_fs], 'yy': [tab2_freq_fs]})
        tab2_p_dspec.image(image="data", x=tab2_dtim_fs[0], y=tab2_freq_fs[0],
                           dw=tab2_dtim_fs[-1] - tab2_dtim_fs[0],
                           dh=tab2_freq_fs[-1] - tab2_freq_fs[0],
                           source=tab2_SRC_dspec_image, palette=bokehpalette_jet)
        tab2_r_square = tab2_p_dspec.square('time', 'freq', source=tab2_SRC_dspec_square, fill_color=colors_dspec,
                                            fill_alpha=0.0,
                                            line_color=None, line_alpha=0.0, selection_fill_alpha=0.1,
                                            selection_fill_color='black',
                                            nonselection_fill_alpha=0.0,
                                            selection_line_alpha=0.2, selection_line_color='white',
                                            nonselection_line_alpha=0.0,
                                            size=max(
                                                config_plot['plot_config']['tab_FSview_base'][
                                                    'dspec_wdth'] / tab2_ntim_fs,
                                                config_plot['plot_config']['tab_FSview_base'][
                                                    'dspec_hght'] / tab2_nfreq_fs))
    else:
        tab2_r_square_rs_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [], 'yy': []})).data


# initial the source of maxfit centroid
def tab2_SRC_maxfit_centroid_init(dspecDFsel):
    start_timestamp = time.time()
    global SRC_maxfit_centroid
    SRC_maxfit_centroid = {}
    for ll in np.unique(dspecDFsel['time']):
        df_tmp = pd.DataFrame(
            {'freq': [], 'shape_longitude': [], 'shape_latitude': [], 'shape_majoraxis': [],
             'shape_minoraxis': [], 'peak': [],
             'shape_positionangle': []})
        SRC_maxfit_centroid[np.where(abs(tab2_dtim - ll) < 0.02)[0].tolist()[0]] = ColumnDataSource(df_tmp)
    print("---tab2_SRC_maxfit_centroid_init -- %s seconds ---" % (time.time() - start_timestamp))


def aia_submap_wavelength_selection(attrname, old, new):
    global tab3_r_aia_submap
    select_wave = tab2_Select_aia_wave.value
    print 'wavelength {} selected'.format(select_wave)
    aiamap = sdomapfromlocalfile(wavelength=select_wave, jdtime=xx[0] / 3600. / 24.)
    lengthx = vla_local_pfmap.dw[0] * u.arcsec
    lengthy = vla_local_pfmap.dh[0] * u.arcsec
    x0 = vla_local_pfmap.smap.center.x
    y0 = vla_local_pfmap.smap.center.y
    aiamap_submap = aiamap.submap(u.Quantity([x0 - lengthx / 2, x0 + lengthx / 2]),
                                  u.Quantity([y0 - lengthy / 2, y0 + lengthy / 2]))
    aia_submap_pfmap = PuffinMap(smap=aiamap_submap,
                                 plot_height=config_plot['plot_config']['tab_FSview_FitANLYS'][
                                     'aia_submap_hght'],
                                 plot_width=config_plot['plot_config']['tab_FSview_FitANLYS'][
                                     'aia_submap_wdth'],
                                 webgl=config_plot['plot_config']['WebGL'])
    tab3_r_aia_submap.data_source.data['data'] = aia_submap_pfmap.ImageSource().data['data']


def tab3_slider_LinkImg_update(attrname, old, new):
    global hdu
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
    hdufile = fits_LOCL_dir + dspecDF_frac.loc[tidx, :]['fits_local']
    if os.path.exists(hdufile):
        hdu = read_fits(hdufile)
        hdu_goodchan = goodchan(hdu)
        freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
        freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
        idxfreq = freq.index(freq_ref)
        fidx_hdu = fidx - idxfreq
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
                              plot_width=tab2_LinkImg_WDTH, webgl=config_plot['plot_config']['WebGL'])
            SRC_Img = pfmap.ImageSource()
            tab2_r_vla.data_source.data['data'] = SRC_Img.data['data']
            mapx, mapy = pfmap.meshgrid()
            mapx, mapy = mapx.value, mapy.value
            SRC_contour = get_contour_data(mapx, mapy, pfmap.smap.data)
            tab2_r_vla_multi_line.data_source.data = SRC_contour.data
            tab2_Div_LinkImg_plot.text = '<p><b>{}</b> loaded.</p>'.format(
                dspecDF_frac.loc[tidx, :]['fits_local'])
        else:
            tab2_Div_LinkImg_plot.text = '<p><b>freq idx</b> out of range.</p>'
    else:
        tab2_Div_LinkImg_plot.text = '<p><b>{}</b> not found.</p>'.format(
            dspecDF_frac.loc[tidx, :]['fits_local'])


def tab2_SRC_maxfit_centroid_update(dspecDFsel):
    start_timestamp = time.time()
    global SRC_maxfit_centroid, timebin
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        SRC_maxfit_centroid = {}
        for ll in np.unique(dspecDFsel['time']):
            dftmp = dspecDFsel[dspecDFsel.time == ll]
            dftmp = dftmp.dropna(how='any')
            df_tmp = pd.concat(
                [dftmp.loc[:, 'freq'], dftmp.loc[:, 'shape_longitude'], dftmp.loc[:, 'shape_latitude'],
                 dftmp.loc[:, 'shape_majoraxis'],
                 dftmp.loc[:, 'shape_minoraxis'], dftmp.loc[:, 'peak'],
                 dftmp.loc[:, 'shape_positionangle'] - np.pi / 2], axis=1)
            SRC_maxfit_centroid[np.where(abs(tab2_dtim - ll) < 0.02)[0].tolist()[0]] = ColumnDataSource(df_tmp)
    else:
        time_dspec = np.unique(dspecDFsel['time'])
        ntime_dspec = len(time_dspec)
        if timebin != 1:
            tidx = np.arange(0, ntime_dspec + 1, timebin)
            time_seq = time_dspec[0:0 + timebin]
            dftmp = dspecDFsel[dspecDFsel['time'].isin(time_seq)]
            dftmp = dftmp.dropna(how='any')
            dftmp_concat = pd.DataFrame(dict(dftmp.mean()), index=[0, ])
            for ll in tidx[1:]:
                time_seq = time_dspec[ll:ll + timebin]
                dftmp = dspecDFsel[dspecDFsel['time'].isin(time_seq)]
                dftmp = dftmp.dropna(how='any')
                dftmp_concat = dftmp_concat.append(pd.DataFrame(dict(dftmp.mean()), index=[0, ]),
                                                   ignore_index=True)
            SRC_maxfit_centroid = ColumnDataSource(
                dftmp_concat[
                    ['freq', 'shape_longitude', 'shape_latitude', 'shape_majoraxis', 'shape_minoraxis', 'peak',
                     'shape_positionangle']].dropna(
                    how='any'))
        else:
            dftmp = dspecDFsel.copy()
            dftmp = dftmp.dropna(how='any')
            df_tmp = pd.concat(
                [dftmp.loc[:, 'freq'], dftmp.loc[:, 'shape_longitude'], dftmp.loc[:, 'shape_latitude'],
                 dftmp.loc[:, 'shape_majoraxis'],
                 dftmp.loc[:, 'shape_minoraxis'], dftmp.loc[:, 'peak'],
                 dftmp.loc[:, 'shape_positionangle'] - np.pi / 2], axis=1)
            SRC_maxfit_centroid = ColumnDataSource(df_tmp)
    print("--- tab2_SRC_maxfit_centroid_update -- %s seconds ---" % (time.time() - start_timestamp))


def tab2_panel_exit():
    tab2_panel2_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    tab2_panel3_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


def tab2_aia_submap_square_selection_change(attrname, old, new):
    global tab3_SRC_dspec_vectorx, tab3_SRC_dspec_vectory
    global vmax_vx, vmax_vy, vmin_vx, vmin_vy, mean_vx, mean_vy
    global VdspecDF
    tab2_aia_submap_square_selected = tab2_SRC_aia_submap_square.selected['1d']['indices']
    if tab2_aia_submap_square_selected:
        ImgDF = ImgDF0.iloc[tab2_aia_submap_square_selected, :]
        xa0, xa1 = ImgDF['xx'].min(), ImgDF['xx'].max()
        ya0, ya1 = ImgDF['yy'].min(), ImgDF['yy'].max()
        mean_vx = (xa0 + xa1) / 2
        mean_vy = (ya0 + ya1) / 2
        tab3_r_aia_submap_rect.data_source.data['x'] = [mean_vx]
        tab3_r_aia_submap_rect.data_source.data['y'] = [mean_vy]
        tab3_r_aia_submap_rect.data_source.data['width'] = [(xa1 - xa0)]
        tab3_r_aia_submap_rect.data_source.data['height'] = [(ya1 - ya0)]
        vx = (VdspecDF['shape_longitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
        vmax_vx, vmin_vx = xa1, xa0
        vx[vx > vmax_vx] = vmax_vx
        vx[vx < vmin_vx] = vmin_vx
        tab3_SRC_dspec_vectorx.data['data'] = [vx]
        vy = (VdspecDF['shape_latitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
        vmax_vy, vmin_vy = ya1, ya0
        vy[vy > vmax_vy] = vmax_vy
        vy[vy < vmin_vy] = vmin_vy
        tab3_SRC_dspec_vectory.data['data'] = [vy]
        tab3_dspec_small_CTRLs_OPT['vmax_values_last'][1] = xa1
        tab3_dspec_small_CTRLs_OPT['vmax_values_last'][2] = ya1
        tab3_dspec_small_CTRLs_OPT['vmin_values_last'][1] = xa0
        tab3_dspec_small_CTRLs_OPT['vmin_values_last'][2] = ya0
    else:
        tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
        tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
        tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
        tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]


def VdspecDF_init():
    global VdspecDF, dspecDF_frac
    VdspecDF = pd.DataFrame()
    nrows_dspecDF = len(dspecDF_frac.index)
    VdspecDF['peak'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF_frac.index)
    VdspecDF['shape_longitude'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF_frac.index)
    VdspecDF['shape_latitude'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF_frac.index)


def VdspecDF_update(selected=None):
    global VdspecDF
    if selected:
        VdspecDF.loc[selected, 'shape_longitude'] = dspecDF_frac.loc[selected, 'shape_longitude']
        VdspecDF.loc[selected, 'shape_latitude'] = dspecDF_frac.loc[selected, 'shape_latitude']
        VdspecDF.loc[selected, 'peak'] = dspecDF_frac.loc[selected, 'peak']
    else:
        VdspecDF.loc[:, 'shape_longitude'] = dspecDF_frac.loc[:, 'shape_longitude']
        VdspecDF.loc[:, 'shape_latitude'] = dspecDF_frac.loc[:, 'shape_latitude']
        VdspecDF.loc[:, 'peak'] = dspecDF_frac.loc[:, 'peak']


def tab3_SRC_dspec_vector_init():
    global tab3_SRC_dspec_vector, tab3_SRC_dspec_vectorx, tab3_SRC_dspec_vectory
    global mean_amp_g, mean_vx, mean_vy, drange_amp_g, drange_vx, drange_vy
    global vmax_amp_g, vmax_vx, vmax_vy, vmin_amp_g, vmin_vx, vmin_vy
    start_timestamp = time.time()
    amp_g = (dspecDF_frac['peak'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_amp_g = np.nanmean(amp_g)
    drange_amp_g = 40.
    vmax_amp_g, vmin_amp_g = mean_amp_g + drange_amp_g * np.asarray([1., -1.])
    amp_g[amp_g > vmax_amp_g] = vmax_amp_g
    amp_g[amp_g < vmin_amp_g] = vmin_amp_g
    tab3_SRC_dspec_vector = ColumnDataSource(data={'data': [amp_g], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
    vx = (dspecDF_frac['shape_longitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_vx = np.nanmean(vx)
    drange_vx = 40.
    vmax_vx, vmin_vx = mean_vx + drange_vx * np.asarray([1., -1.])
    vx[vx > vmax_vx] = vmax_vx
    vx[vx < vmin_vx] = vmin_vx
    tab3_SRC_dspec_vectorx = ColumnDataSource(data={'data': [vx], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
    vy = (dspecDF_frac['shape_latitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_vy = np.nanmean(vy)
    drange_vy = 40.
    vmax_vy, vmin_vy = mean_vy + drange_vy * np.asarray([1., -1.])
    vy[vy > vmax_vy] = vmax_vy
    vy[vy < vmin_vy] = vmin_vy
    tab3_SRC_dspec_vectory = ColumnDataSource(data={'data': [vy], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
    tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
    tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
    tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
    tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]
    print("--- tab3_SRC_dspec_small_init -- %s seconds ---" % (time.time() - start_timestamp))


def tab3_SRC_dspec_vector_update():
    global tab3_SRC_dspec_vector, tab3_SRC_dspec_vectorx, tab3_SRC_dspec_vectory
    global mean_amp_g, mean_vx, mean_vy, drange_amp_g, drange_vx, drange_vy
    global vmax_amp_g, vmax_vx, vmax_vy, vmin_amp_g, vmin_vx, vmin_vy
    global VdspecDF
    start_timestamp = time.time()
    amp_g = (VdspecDF['peak'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_amp_g = np.nanmean(amp_g)
    drange_amp_g = 40.
    vmax_amp_g, vmin_amp_g = mean_amp_g + drange_amp_g * np.asarray([1., -1.])
    amp_g[amp_g > vmax_amp_g] = vmax_amp_g
    amp_g[amp_g < vmin_amp_g] = vmin_amp_g
    tab3_SRC_dspec_vector.data['data'] = [amp_g]
    # dspecDFselect = dspecDFselect.where(dspecDFselect['peak']>vmin_amp_g)
    # tab2_SRC_maxfit_centroid_update(dspecDFselect)
    # todo add threshold selection to the vector dynamic spectrum
    vx = (VdspecDF['shape_longitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_vx = np.nanmean(vx)
    drange_vx = 40.
    vmax_vx, vmin_vx = mean_vx + drange_vx * np.asarray([1., -1.])
    vx[vx > vmax_vx] = vmax_vx
    vx[vx < vmin_vx] = vmin_vx
    tab3_SRC_dspec_vectorx.data['data'] = [vx]
    vy = (VdspecDF['shape_latitude'].copy()).reshape(tab2_nfreq, tab2_ntim)
    mean_vy = np.nanmean(vy)
    drange_vy = 40.
    vmax_vy, vmin_vy = mean_vy + drange_vy * np.asarray([1., -1.])
    vy[vy > vmax_vy] = vmax_vy
    vy[vy < vmin_vy] = vmin_vy
    tab3_SRC_dspec_vectory.data['data'] = [vy]
    print("--- tab3_SRC_dspec_small_update -- %s seconds ---" % (time.time() - start_timestamp))


def tab2_dspec_selection_change(attrname, old, new):
    global tab2_dspec_selected
    tab2_dspec_selected = tab2_SRC_dspec_square.selected['1d']['indices']
    if tab2_dspec_selected:
        global dspecDF_select
        dspecDF_select = dspecDF_frac.iloc[tab2_dspec_selected, :]
        idx_selected = dspecDF_select.index[len(dspecDF_select) / 2]
        tidx = int(['{:.3f}'.format(ll) for ll in tab2_dtim].index(
            '{:.3f}'.format(dspecDF_frac.loc[idx_selected, :]['time'])))
        fidx = int(['{:.3f}'.format(ll) for ll in tab2_freq].index(
            '{:.3f}'.format(dspecDF_frac.loc[idx_selected, :]['freq'])))
        tab2_Slider_time_LinkImg.value = tidx
        tab2_Slider_freq_LinkImg.value = fidx


def tab2_vla_square_selection_change(attrname, old, new):
    global x0, x1, y0, y1
    tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
    if tab2_vla_square_selected:
        ImgDF = ImgDF0.iloc[tab2_vla_square_selected, :]
        x0, x1 = ImgDF['xx'].min(), ImgDF['xx'].max()
        y0, y1 = ImgDF['yy'].min(), ImgDF['yy'].max()
        tab2_r_vla_ImgRgn_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
    else:
        tab2_r_vla_ImgRgn_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [], 'yy': []})).data


def tab2_update_MapRES(attrname, old, new):
    start_timestamp = time.time()
    select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
    dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
    aia_resampled_map = aiamap.resample(dimensions)
    aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_AIA = aia_resampled_pfmap.ImageSource()
    tab2_r_aia.data_source.data['data'] = SRC_AIA.data['data']
    hmi_resampled_map = hmimap.resample(dimensions)
    hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_HMI = hmi_resampled_pfmap.ImageSource()
    tab2_r_hmi.data_source.data['data'] = SRC_HMI.data['data']
    print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


def tab2_save_region():
    tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
    if tab2_vla_square_selected:
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


def tab2_dspec_vector_selection_change(attrname, old, new):
    global tab2_dspec_vector_selected
    tab2_dspec_vector_selected = tab2_SRC_dspec_vector_square.selected['1d']['indices']
    if tab2_dspec_vector_selected:
        global dspecDF_select
        dspecDF_select = dspecDF_frac.iloc[tab2_dspec_vector_selected, :]
        VdspecDF_init()
        VdspecDF_update(selected=tab2_dspec_vector_selected)
        # tab3_SRC_dspec_vector_update(VdspecDF)
        tab2_SRC_maxfit_centroid_update(dspecDF_select)
        if tab3_BUT_animate_ONOFF.label == 'Animate OFF & Go':
            tab3_r_aia_submap_cross.visible = True
            tab3_r_dspec_vector_line.visible = False
            tab3_r_dspec_vectorx_line.visible = False
            tab3_r_dspec_vectory_line.visible = False
            tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid.data


def tab3_RBG_dspec_small_update(attrname, old, new):
    idx_p_dspec_small = tab3_RBG_dspec_small.active
    global tab3_dspec_small_CTRLs_OPT
    tab3_dspec_small_CTRLs_OPT['idx_p_dspec_small'] = idx_p_dspec_small
    tab3_dspec_small_CTRLs_OPT['radio_button_group_dspec_small_update_flag'] = True
    mean_values = tab3_dspec_small_CTRLs_OPT['mean_values']
    drange_values = tab3_dspec_small_CTRLs_OPT['drange_values']
    vmax_values_last = tab3_dspec_small_CTRLs_OPT['vmax_values_last']
    vmin_values_last = tab3_dspec_small_CTRLs_OPT['vmin_values_last']
    tab3_Slider_dspec_small_dmax.start = mean_values[idx_p_dspec_small] - drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmax.end = mean_values[idx_p_dspec_small] + 2 * drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmax.value = vmax_values_last[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.start = mean_values[idx_p_dspec_small] - 2 * drange_values[
        idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.end = mean_values[idx_p_dspec_small] + drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.value = vmin_values_last[idx_p_dspec_small]
    tab3_dspec_small_CTRLs_OPT['radio_button_group_dspec_small_update_flag'] = False


def tab3_BUT_dspec_small_reset_update():
    global VdspecDF, tab2_nfreq, tab2_ntim, tab3_dspec_small_CTRLs_OPT
    global tab3_SRC_dspec_vector, tab3_SRC_dspec_vectorx, tab3_SRC_dspec_vectory
    items_dspec_small = tab3_dspec_small_CTRLs_OPT['items_dspec_small']
    mean_values = tab3_dspec_small_CTRLs_OPT['mean_values']
    drange_values = tab3_dspec_small_CTRLs_OPT['drange_values']
    vmax_values = tab3_dspec_small_CTRLs_OPT['vmax_values']
    vmin_values = tab3_dspec_small_CTRLs_OPT['vmin_values']
    source_list = [tab3_SRC_dspec_vector, tab3_SRC_dspec_vectorx, tab3_SRC_dspec_vectory]
    for ll, item in enumerate(items_dspec_small):
        TmpData = (VdspecDF[item].copy()).reshape(tab2_nfreq, tab2_ntim)
        TmpData[TmpData > vmax_values[ll]] = vmax_values[ll]
        TmpData[TmpData < vmin_values[ll]] = vmin_values[ll]
        source_list[ll].data['data'] = [TmpData]
    idx_p_dspec_small = 0
    tab3_dspec_small_CTRLs_OPT['idx_p_dspec_small'] = idx_p_dspec_small
    tab3_RBG_dspec_small.active = idx_p_dspec_small
    tab3_Slider_dspec_small_dmax.start = mean_values[idx_p_dspec_small] - drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmax.end = mean_values[idx_p_dspec_small] + 2 * drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmax.value = vmax_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.start = mean_values[idx_p_dspec_small] - 2 * drange_values[
        idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.end = mean_values[idx_p_dspec_small] + drange_values[idx_p_dspec_small]
    tab3_Slider_dspec_small_dmin.value = vmin_values[idx_p_dspec_small]
    tab3_dspec_small_CTRLs_OPT['vmax_values_last'] = [ll for ll in vmax_values]
    tab3_dspec_small_CTRLs_OPT['vmin_values_last'] = [ll for ll in vmin_values]
    vmax_vx, vmax_vy = tab3_dspec_small_CTRLs_OPT['vmax_values_last'][1:]
    vmin_vx, vmin_vy = tab3_dspec_small_CTRLs_OPT['vmin_values_last'][1:]
    tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
    tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
    tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
    tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]


def tab3_BUT_dspec_small_resetall_update():
    VdspecDF_update()
    tab3_BUT_dspec_small_reset_update()
    print 'reset all'


def tab3_slider_dspec_small_update(attrname, old, new):
    global VdspecDF, tab2_nfreq, tab2_ntim, tab3_dspec_small_CTRLs_OPT
    items_dspec_small = tab3_dspec_small_CTRLs_OPT['items_dspec_small']
    idx_p_dspec_small = tab3_dspec_small_CTRLs_OPT['idx_p_dspec_small']
    dmax = tab3_Slider_dspec_small_dmax.value
    dmin = tab3_Slider_dspec_small_dmin.value
    if not tab3_dspec_small_CTRLs_OPT['radio_button_group_dspec_small_update_flag']:
        tab3_dspec_small_CTRLs_OPT['vmax_values_last'][idx_p_dspec_small] = dmax
        tab3_dspec_small_CTRLs_OPT['vmin_values_last'][idx_p_dspec_small] = dmin
    TmpData = (VdspecDF[items_dspec_small[idx_p_dspec_small]].copy()).reshape(tab2_nfreq, tab2_ntim)
    TmpData[TmpData > dmax] = dmax
    TmpData[TmpData < dmin] = dmin
    if idx_p_dspec_small == 0:
        global tab3_SRC_dspec_vector
        tab3_SRC_dspec_vector.data['data'] = [TmpData]
    elif idx_p_dspec_small == 1:
        global tab3_SRC_dspec_vectorx
        tab3_SRC_dspec_vectorx.data['data'] = [TmpData]
    elif idx_p_dspec_small == 2:
        global tab3_SRC_dspec_vectory
        tab3_SRC_dspec_vectory.data['data'] = [TmpData]
    vmax_vx, vmax_vy = tab3_dspec_small_CTRLs_OPT['vmax_values_last'][1:]
    vmin_vx, vmin_vy = tab3_dspec_small_CTRLs_OPT['vmin_values_last'][1:]
    tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
    tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
    tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
    tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]


def tab3_slider_ANLYS_idx_update(attrname, old, new):
    global tab2_dtim, tab2_freq, tab2_ntim, SRC_maxfit_centroid
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        tab3_Slider_ANLYS_idx.start = next(
            i for i in xrange(tab2_ntim) if tab2_dtim[i] >= tab3_p_dspec_vector.x_range.start)
        tab3_Slider_ANLYS_idx.end = next(
            i for i in xrange(tab2_ntim - 1, -1, -1) if tab2_dtim[i] <= tab3_p_dspec_vector.x_range.end) + 1
        indices_time = tab3_Slider_ANLYS_idx.value
        tab3_r_dspec_vector_line.visible = True
        tab3_r_dspec_vector_line.data_source.data = ColumnDataSource(
            pd.DataFrame({'time': [tab2_dtim[indices_time], tab2_dtim[indices_time]],
                          'freq': [tab2_freq[0], tab2_freq[-1]]})).data
        try:
            tab3_r_aia_submap_cross.visible = True
            tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid[indices_time].data
        except:
            tab3_r_aia_submap_cross.visible = False
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


def tab3_BUT_plot_xargs_default():
    global tab3_plot_xargs_dict
    tab3_plot_xargs_dict = OrderedDict()
    tab3_plot_xargs_dict['timebin'] = "1"
    tab3_plot_xargs_dict['timeline'] = "False"
    tab3_Div_plot_xargs_text = '<p>' + ';'.join(
        "<b>{}</b> = {}".format(key, val) for (key, val) in tab3_plot_xargs_dict.items()) + '</p>'
    tab3_Div_plot_xargs.text = tab3_Div_plot_xargs_text
    tab3_Div_Tb.text = '<p><b>Default xargs Restored.</b></p>'


def tab3_animate_update():
    global tab3_animate_step, tab2_dspec_vector_selected
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        if tab2_dspec_vector_selected:
            indices_time = tab3_Slider_ANLYS_idx.value + tab3_animate_step
            if (tab3_animate_step == timebin) and (indices_time > tab3_Slider_ANLYS_idx.end):
                indices_time = tab3_Slider_ANLYS_idx.start
            if (tab3_animate_step == -timebin) and (indices_time < tab3_Slider_ANLYS_idx.start):
                indices_time = tab3_Slider_ANLYS_idx.end
            tab3_Slider_ANLYS_idx.value = indices_time
            tab3_Div_Tb.text = """ """
        else:
            tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


def tab3_animate():
    global tab2_dspec_vector_selected
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        if tab3_BUT_PlayCTRL.label == 'Play':
            if tab2_dspec_vector_selected:
                tab3_BUT_PlayCTRL.label = 'Pause'
                tab3_BUT_PlayCTRL.button_type = 'danger'
                curdoc().add_periodic_callback(tab3_animate_update, 125)
                tab3_Div_Tb.text = """ """
            else:
                tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""
        else:
            tab3_BUT_PlayCTRL.label = 'Play'
            tab3_BUT_PlayCTRL.button_type = 'success'
            curdoc().remove_periodic_callback(tab3_animate_update)
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


def tab3_animate_step_CTRL():
    global tab3_animate_step, tab2_dspec_vector_selected
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        if tab2_dspec_vector_selected:
            if tab3_BUT_PlayCTRL.label == 'Pause':
                tab3_BUT_PlayCTRL.label = 'Play'
                tab3_BUT_PlayCTRL.button_type = 'success'
                curdoc().remove_periodic_callback(tab3_animate_update)
            idx = tab3_Slider_ANLYS_idx.value + tab3_animate_step
            if (tab3_animate_step == timebin) and (idx > tab3_Slider_ANLYS_idx.end):
                idx = tab3_Slider_ANLYS_idx.start
            elif (tab3_animate_step == -timebin) and (idx < tab3_Slider_ANLYS_idx.start):
                idx = tab3_Slider_ANLYS_idx.end
            tab3_Slider_ANLYS_idx.value = idx
            tab3_Div_Tb.text = """ """
        else:
            tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


def tab3_animate_FRWD_REVS():
    global tab3_animate_step
    if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
        if tab2_dspec_vector_selected:
            if tab3_animate_step == timebin:
                tab3_BUT_FRWD_REVS_CTRL.label = 'Reverse'
                tab3_animate_step = -timebin
            else:
                tab3_BUT_FRWD_REVS_CTRL.label = 'Forward'
                tab3_animate_step = timebin
            tab3_Div_Tb.text = """ """
        else:
            tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


def tab3_animate_onoff():
    if tab2_dspec_vector_selected:
        global tab3_plot_xargs_dict, timebin, timeline
        if not 'timebin' in tab3_plot_xargs_dict.keys():
            tab3_plot_xargs_dict['timebin'] = '1'
        if not 'timeline' in tab3_plot_xargs_dict.keys():
            tab3_plot_xargs_dict['timeline'] = 'False'
        txts = tab3_input_plot_xargs.value.strip()
        txts = txts.split(';')
        for txt in txts:
            txt = txt.strip()
            txt = txt.split('=')
            if len(txt) == 2:
                key, val = txt
                key, val = key.strip(), val.strip()
                if key == 'timebin':
                    if not (0 <= int(val) <= tab2_ntim - 1):
                        val = '1'
                    timebin = int(val)
                if key == 'timeline':
                    if val not in ['True', 'False']:
                        val = 'False'
                    timeline = json.loads(val.lower())
                tab3_plot_xargs_dict[key.strip()] = val.strip()
                if key not in ['timebin', 'timeline']:
                    tab3_plot_xargs_dict.pop(key, None)
            else:
                tab3_Div_plot_xargs.text = '<p>Input syntax: <b>timebin</b>=1; <b>linesytle</b>=False;' \
                                           'Any spaces will be ignored.</p>'

        tab3_Div_plot_xargs_text = '<p>' + ';'.join(
            "<b>{}</b> = {}".format(key, val) for (key, val) in tab3_plot_xargs_dict.items()) + '</p>'
        tab3_Div_plot_xargs.text = tab3_Div_plot_xargs_text
        tab3_animate_step = timebin
        tab3_Slider_ANLYS_idx.step = timebin
        if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
            tab3_BUT_animate_ONOFF.label = 'Animate OFF & Go'
            tab3_r_aia_submap_cross.visible = True
            tab3_r_aia_submap_line.visible = timeline
            tab3_r_dspec_vector_line.visible = False
            tab3_r_dspec_vectorx_line.visible = False
            tab3_r_dspec_vectory_line.visible = False
            tab2_SRC_maxfit_centroid_update(dspecDF_select)
            tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid.data
        else:
            tab3_BUT_animate_ONOFF.label = 'Animate ON & Go'
            tab3_r_aia_submap_cross.visible = True
            tab3_r_aia_submap_line.visible = False
            tab3_r_dspec_vector_line.visible = True
            tab3_r_dspec_vectorx_line.visible = True
            tab3_r_dspec_vectory_line.visible = True
            tab2_SRC_maxfit_centroid_update(dspecDF_select)
            indices_time = tab3_Slider_ANLYS_idx.value
            tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid[indices_time].data
            tab3_Div_Tb.text = """ """
    else:
        tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""


def tab2_update_MapRES(attrname, old, new):
    start_timestamp = time.time()
    select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
    dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
    aia_resampled_map = aiamap.resample(dimensions)
    aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_AIA = aia_resampled_pfmap.ImageSource()
    tab2_r_aia.data_source.data['data'] = SRC_AIA.data['data']
    hmi_resampled_map = hmimap.resample(dimensions)
    hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_HMI = hmi_resampled_pfmap.ImageSource()
    tab2_r_hmi.data_source.data['data'] = SRC_HMI.data['data']
    print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


def tab3_slider_LinkImg_update(attrname, old, new):
    global hdu
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
    hdufile = fits_LOCL_dir + dspecDF_frac.loc[tidx, :]['fits_local']
    if os.path.exists(hdufile):
        hdu = read_fits(hdufile)
        hdu_goodchan = goodchan(hdu)
        freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
        freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
        idxfreq = freq.index(freq_ref)
        fidx_hdu = fidx - idxfreq
        if hdu_goodchan[0] <= fidx_hdu <= hdu_goodchan[-1]:
            if select_vla_pol == 'RR':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :]
            elif select_vla_pol == 'LL':
                vladata = hdu.data[pols.index('LL'), fidx_hdu, :, :]
            elif select_vla_pol == 'I':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :] + hdu.data[pols.index('1'), fidx_hdu,
                                                                       :, :]
            elif select_vla_pol == 'V':
                vladata = hdu.data[pols.index('RR'), fidx_hdu, :, :] - hdu.data[pols.index('1'), fidx_hdu,
                                                                       :, :]
            pfmap = PuffinMap(vladata, hdu.header, plot_height=tab2_LinkImg_HGHT,
                              plot_width=tab2_LinkImg_WDTH, webgl=config_plot['plot_config']['WebGL'])
            SRC_Img = pfmap.ImageSource()
            tab2_r_vla.data_source.data['data'] = SRC_Img.data['data']
            mapx, mapy = pfmap.meshgrid()
            mapx, mapy = mapx.value, mapy.value
            SRC_contour = get_contour_data(mapx, mapy, pfmap.smap.data)
            tab2_r_vla_multi_line.data_source.data = SRC_contour.data
            tab2_Div_LinkImg_plot.text = '<p><b>{}</b> loaded.</p>'.format(
                dspecDF_frac.loc[tidx, :]['fits_local'])
        else:
            tab2_Div_LinkImg_plot.text = '<p><b>freq idx</b> out of range.</p>'
    else:
        tab2_Div_LinkImg_plot.text = '<p><b>{}</b> not found.</p>'.format(
            dspecDF_frac.loc[tidx, :]['fits_local'])


def tab2_panel_exit():
    tab2_panel2_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


def tab2_vla_square_selection_change(attrname, old, new):
    global x0, x1, y0, y1
    tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
    if tab2_vla_square_selected:
        ImgDF = ImgDF0.iloc[tab2_vla_square_selected, :]
        x0, x1 = ImgDF['xx'].min(), ImgDF['xx'].max()
        y0, y1 = ImgDF['yy'].min(), ImgDF['yy'].max()
        tab2_r_vla_ImgRgn_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
        idxmax = max(tab2_vla_square_selected)
        idxmin = min(tab2_vla_square_selected)
        x0pix, x1pix = idxmin % mapvlasize[0], idxmax % mapvlasize[0]
        y0pix, y1pix = idxmin / mapvlasize[0], idxmax / mapvlasize[0]
        tab2_tImfit_Param_dict['box'] = "'{},{},{},{}'".format(x0pix, y0pix, x1pix, y1pix)
        tab2_Div_tImfit_text = '<p><b>#  pimfit :: Fit one or more elliptical Gaussian components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
        tab2_Div_tImfit.text = tab2_Div_tImfit_text
    else:
        tab2_r_vla_ImgRgn_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [], 'yy': []})).data


def tab2_update_MapRES(attrname, old, new):
    start_timestamp = time.time()
    select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
    dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
    aia_resampled_map = aiamap.resample(dimensions)
    aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_AIA = aia_resampled_pfmap.ImageSource()
    tab2_r_aia.data_source.data['data'] = SRC_AIA.data['data']
    hmi_resampled_map = hmimap.resample(dimensions)
    hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
    SRC_HMI = hmi_resampled_pfmap.ImageSource()
    tab2_r_hmi.data_source.data['data'] = SRC_HMI.data['data']
    print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


def tab2_save_region():
    tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
    if tab2_vla_square_selected:
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
    tab2_Div_tImfit_text = '<p><b>#  pimfit :: Fit one or more elliptical Gaussian components \
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
            tab2_tImfit_Param_dict.pop(key)
        except:
            tab2_Div_tImfit2.text = '<p>Input syntax: <b>stokes</b>; <b>ncpu</b>; ' \
                                    '<b>region</b>. Any spaces will be ignored.</p>'
    tab2_Div_tImfit_text = '<p><b>#  pimfit :: Fit one or more elliptical Gaussian components \
    on an image region(s)</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    tab2_Div_tImfit.text = tab2_Div_tImfit_text


def tab2_BUT_tImfit_param_default():
    global tab2_tImfit_Param_dict
    tab2_tImfit_Param_dict = OrderedDict()
    vlafileliststr = "'" + "','".join(vlafile) + "'"
    tab2_tImfit_Param_dict['event_id'] = "'{}'".format(event_id.replace("/", ""))
    tab2_tImfit_Param_dict['struct_id'] = "'{}'".format(struct_id.replace("/", ""))
    tab2_tImfit_Param_dict['ncpu'] = "10"
    tab2_tImfit_Param_dict['box'] = "''"
    tab2_tImfit_Param_dict['stokes'] = "'{}'".format(tab2_Select_vla_pol.value)
    tab2_tImfit_Param_dict['mask'] = "''"
    tab2_tImfit_Param_dict['imagefiles'] = "[{}]".format(vlafileliststr)
    tab2_Div_tImfit_text = '<p><b>#  pimfit :: Fit one or more elliptical Gaussian components \
    on an image region(s)</b></p>' + ' '.join(
        "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tImfit_Param_dict.items())
    tab2_Div_tImfit.text = tab2_Div_tImfit_text
    tab2_Div_tImfit2.text = '<p><b>Default parameter Restored.</b></p>'


def tab2_BUT_timfit_param_reload():
    global tab2_tImfit_Param_dict
    try:
        with open(database_dir + event_id + struct_id + 'CASA_imfit_args.json', 'r') as fp:
            tab2_timfit_Param_dict = json.load(fp)
        tab2_Div_timfit_text = '<p><b>#  pimfit :: Fit one or more elliptical Gaussian components \
        on an image region(s)</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_timfit_Param_dict.items())
        tab2_Div_tImfit.text = tab2_Div_timfit_text
        tab2_Div_tImfit2.text = '<p>CASA pimfit arguments reload from config file in <b>{}</b>.</p>'.format(
            database_dir + event_id + struct_id)
    except:
        tab2_Div_tImfit2.text = '<p>' + database_dir + event_id + struct_id + \
                                '{}CASA_imfit_args.json not found!!!</p>'.format(
                                    database_dir + event_id + struct_id)


def tab2_BUT_tImfit_update():
    with open(database_dir + event_id + struct_id + 'CASA_imfit_args.json', 'w') as fp:
        json.dump(tab2_tImfit_Param_dict, fp)
    os.system('cp {}/DataBrowser/FSview/script_imfit.py {}'.format(suncasa_dir,
                                                                   database_dir + event_id + struct_id))
    tab2_Div_tImfit2.text = '<p>CASA pimfit script and arguments config\
     file saved to <b>{}</b>.</p>'.format(database_dir + event_id + struct_id)
    cwd = os.getcwd()
    try:
        tab2_Div_tImfit2.text = '<p>CASA imfit script and arguments config file saved to <b>{}.</b></p>\
        <p>CASA imfit is in processing.</p>'.format(
            database_dir + event_id + struct_id)
        os.chdir(database_dir + event_id + struct_id)
        suncasapy47 = config_plot['core']['casapy47']
        suncasapy47 = os.path.expandvars(suncasapy47)
        os.system('{} -c script_imfit.py'.format(suncasapy47))
        tab2_Div_tImfit2.text = '<p>imfit finished, go back to <b>QLook</b> \
        window, select StrID <b>{}</b> and click <b>FSview</b> button again.</p>'.format(
            struct_id[0:-1])
    except:
        tab2_Div_tImfit2.text = '<p>CASA imfit script and arguments config file \
        saved to <b>{}.</b></p><p>Do imfit with CASA manually.</p>'.format(
            database_dir + event_id + struct_id) + '<p>When finished, go back to <b>QLook</b> window, \
            select StrID <b>{}</b> and click <b>FSview</b> button again.</p>'.format(struct_id[0:-1])
    os.chdir(cwd)


event_id = config_EvtID['datadir']['event_id']
try:
    with open(database_dir + event_id + 'CurrFS.json', 'r') as fp:
        FS_config = json.load(fp)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']
FS_specfile = FS_config['datadir']['FS_specfile']
tab2_specdata = np.load(FS_specfile)
tab2_spec = tab2_specdata['spec']
tab2_npol = tab2_specdata['npol']
tab2_nbl = tab2_specdata['nbl']
tab2_ntim = tab2_specdata['ntim']
tab2_nfreq = tab2_specdata['nfreq']
tab2_tim = tab2_specdata['tim']
tab2_dt = np.median(np.diff(tab2_tim))
tab2_freq = tab2_specdata['freq'] / 1e9
tab2_freq = [float('{:.03f}'.format(ll)) for ll in tab2_freq]
tab2_df = np.median(np.diff(tab2_freq))

tab2_bl = tab2_specdata['bl'].item().split(';')
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
tim_map = ((np.tile(tab2_tim, tab2_nfreq).reshape(tab2_nfreq, tab2_ntim) / 3600. / 24. + 2400000.5)) * 86400.
freq_map = np.tile(tab2_freq, tab2_ntim).reshape(tab2_ntim, tab2_nfreq).swapaxes(0, 1)
xx = tim_map.flatten()
yy = freq_map.flatten()
fits_LOCL = config_EvtID['datadir']['fits_LOCL']
fits_GLOB = config_EvtID['datadir']['fits_GLOB']
fits_LOCL_dir = database_dir + event_id + struct_id + fits_LOCL
fits_GLOB_dir = database_dir + event_id + struct_id + fits_GLOB

FS_dspecDF = database_dir + event_id + struct_id + config_EvtID['datadir']['dspecDF']
if os.path.exists(FS_dspecDF):
    with open(FS_dspecDF, 'rb') as f:
        dspecDF0 = pickle.load(f)
    # dspecDF0 = pd.merge(dspecDF0.copy(), pd.DataFrame(
    #     {'time': xx - xx[0], 'freq': yy,
    #      'dspecR': tab2_spec_plt_R.ravel(),
    #      'dspecL': tab2_spec_plt_L.ravel(),
    #      'dspecV': tab2_spec_plt_V.ravel()}), how='outer', on=['time', 'freq'])
    dspecDF_frac = dspecDF0.copy()
    dspecDF_select = dspecDF0.copy()
    itemset1 = set(['shape_longitude', 'shape_latitude'])
    itemset2 = set(dspecDF0.columns.tolist())
    if len(itemset2.intersection(itemset1)) == 2:
        '''
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        #################################### FS_view base ######################################
        #######################┏━━━┓┏━━━┓╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋┏┓#####################################
        #######################┃┏━━┛┃┏━┓┃╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋┃┃#####################################
        #######################┃┗━━┓┃┗━━┓┏┓┏┓┏┓┏━━┓┏┓┏┓┏┓┃┗━┓┏━━┓┏━━┓┏━━┓#######################
        #######################┃┏━━┛┗━━┓┃┃┗┛┃┣┫┃┃━┫┃┗┛┗┛┃┃┏┓┃┃┏┓┃┃━━┫┃┃━┫#######################
        #######################┃┃╋╋╋┃┗━┛┃┗┓┏┛┃┃┃┃━┫┗┓┏┓┏┛┃┗┛┃┃┏┓┃┣━━┃┃┃━┫#######################
        #######################┗┛╋╋╋┗━━━┛╋┗┛╋┗┛┗━━┛╋┗┛┗┛╋┗━━┛┗┛┗┛┗━━┛┗━━┛#######################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        '''
        vlafile = glob.glob(fits_LOCL_dir + '*.fits')
        tab2_panel2_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>
                                first before closing the tab</p></b>""",
                                   width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_panel3_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>
                                first before closing the tab</p></b>""",
                                   width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()
        colors_dspec = [colors.rgb2hex(m) for m in colormap_jet((tab2_spec_plt.flatten() - rmin) / (rmax - rmin))]

        '''create the dynamic spectrum plot'''
        TOOLS = "crosshair,pan,wheel_zoom,tap,box_zoom,reset,save"
        tab2_SRC_dspec_square = ColumnDataSource(dspecDF_select)
        tab2_p_dspec = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                              plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_wdth'],
                              plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_hght'],
                              x_range=(tab2_dtim[0], tab2_dtim[-1]), y_range=(tab2_freq[0], tab2_freq[-1]),
                              toolbar_location="above")
        tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        tab2_p_dspec.axis.visible = True
        tab2_p_dspec.title.text = "Dynamic spectrum"
        tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
        tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
        tab2_SRC_dspec_image = ColumnDataSource(
            data={'data': [tab2_spec_plt], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
        tab2_p_dspec.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                           dh=tab2_freq[-1] - tab2_freq[0],
                           source=tab2_SRC_dspec_image, palette=bokehpalette_jet)

        # make the dspec data source selectable
        tab2_r_square = tab2_p_dspec.square('time', 'freq', source=tab2_SRC_dspec_square, fill_color=colors_dspec,
                                            fill_alpha=0.0,
                                            line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                            selection_fill_color='black',
                                            nonselection_fill_alpha=0.0,
                                            selection_line_alpha=0.0, selection_line_color='white',
                                            nonselection_line_alpha=0.0,
                                            size=max(
                                                config_plot['plot_config']['tab_FSview_base']['dspec_wdth'] / tab2_ntim,
                                                config_plot['plot_config']['tab_FSview_base'][
                                                    'dspec_hght'] / tab2_nfreq))

        # tab2_p_dspec.border_fill_color = "silver"
        tab2_p_dspec.border_fill_alpha = 0.4
        tab2_p_dspec.axis.major_tick_out = 0
        tab2_p_dspec.axis.major_tick_in = 5
        tab2_p_dspec.axis.minor_tick_out = 0
        tab2_p_dspec.axis.minor_tick_in = 3
        tab2_p_dspec.axis.major_tick_line_color = "white"
        tab2_p_dspec.axis.minor_tick_line_color = "white"

        tab2_Select_pol = Select(title="Polarization:", value='I', options=['RR', 'LL', 'I', 'V'],
                                 width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_Select_bl = Select(title="Baseline:", value=tab2_bl[0], options=tab2_bl,
                                width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                                        width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_p_dspec_xPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                                   plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_xPro_wdth'],
                                   plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_xPro_hght'],
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
        tab2_l_dspec_xPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                           source=tab2_SRC_dspec_xPro_hover,
                                           render_mode='canvas')
        tab2_l_dspec_xPro_hover.text_font_size = '5pt'
        tab2_p_dspec_xPro.add_layout(tab2_l_dspec_xPro_hover)
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

        tab2_p_dspec_yPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                                   plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_yPro_wdth'],
                                   plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_yPro_hght'],
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
        l_dspec_yPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                      source=tab2_SRC_dspec_yPro_hover,
                                      render_mode='canvas')
        l_dspec_yPro_hover.text_font_size = '5pt'
        tab2_p_dspec_yPro.add_layout(l_dspec_yPro_hover)
        tab2_p_dspec_yPro.title.text_font_size = '6pt'
        tab2_p_dspec_yPro.yaxis.visible = False
        tab2_p_dspec_yPro.background_fill_color = "beige"
        tab2_p_dspec_yPro.background_fill_alpha = 0.4
        tab2_p_dspec_yPro.xaxis.axis_label = 'Intensity [sfu]'
        tab2_p_dspec_yPro.yaxis.axis_label_text_font_size = '5px'
        # tab2_p_dspec_yPro.border_fill_color = "silver"
        tab2_p_dspec_yPro.border_fill_alpha = 0.4
        tab2_p_dspec_yPro.min_border_bottom = 0
        tab2_p_dspec_yPro.min_border_left = 0
        tab2_p_dspec_yPro.axis.major_tick_out = 0
        tab2_p_dspec_yPro.axis.major_tick_in = 5
        tab2_p_dspec_yPro.axis.minor_tick_out = 0
        tab2_p_dspec_yPro.axis.minor_tick_in = 3
        tab2_p_dspec_yPro.axis.major_tick_line_color = "black"
        tab2_p_dspec_yPro.axis.minor_tick_line_color = "black"

        tab2_BUT_vdspec = Button(label="VEC Dyn Spec",
                                 width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                 button_type="success")

        tab2_Select_pol_opt = ['RR', 'LL', 'I', 'V']

        tab2_BUT_vdspec.on_click(tab2_vdspec_update)

        tab2_ctrls = [tab2_Select_bl, tab2_Select_pol, tab2_Select_colorspace]
        for ctrl in tab2_ctrls:
            ctrl.on_change('value', tab2_update_dspec_image)

        # # Add a hover tool
        tooltips = None

        hover_JScode = """
            var nx = %d;
            var ny = %d;
            var data = {'x': [], 'y': []};
            var cdata = rs.get('data');
            var indices = cb_data.index['1d'].indices;
            var idx_offset = indices[0] - (indices[0] %% nx);
            for (i=0; i < nx; i++) {
                data['x'].push(cdata.time[i+idx_offset]);
                data['y'].push(cdata.dspec[i+idx_offset]);
            }
            rdx.set('data', data);
            idx_offset = indices[0] %% nx;
            data = {'x': [], 'y': []};
            for (i=0; i < ny; i++) {
                data['x'].push(cdata.dspec[i*nx+idx_offset]);
                data['y'].push(cdata.freq[i*nx+idx_offset]);
            }
            rdy.set('data', data);
            var time = cdata.timestr[indices[0]]+' '
            var freq = cdata.freq[indices[0]].toFixed(3)+'[GHz] '
            var dspec = cdata.dspec[indices[0]].toFixed(3)+ '[sfu]'
            var tooltips = freq + time + dspec
            data = {'x': [], 'y': [], 'tooltips': []};
            data['x'].push(cdata.time[indices[0]]);
            data['y'].push(cdata.dspec[indices[0]]);
            data['tooltips'].push(tooltips);
            rdx_hover.set('data', data);
            tooltips = time + freq + dspec
            data = {'x': [], 'y': [], 'tooltips': []};
            data['x'].push(cdata.dspec[indices[0]]);
            data['y'].push(cdata.freq[indices[0]]);
            data['tooltips'].push(tooltips);
            rdy_hover.set('data', data);
            """ % (tab2_ntim, tab2_nfreq)

        tab2_p_dspec_hover_callback = CustomJS(
            args={'rs': tab2_r_square.data_source, 'rdx': r_dspec_xPro.data_source, 'rdy': r_dspec_yPro.data_source,
                  'rdx_hover': r_dspec_xPro_hover.data_source,
                  'rdy_hover': r_dspec_yPro_hover.data_source}, code=hover_JScode)
        tab2_p_dspec_hover = HoverTool(tooltips=tooltips, callback=tab2_p_dspec_hover_callback,
                                       renderers=[tab2_r_square])
        tab2_p_dspec.add_tools(tab2_p_dspec_hover)

        # initial the VLA map contour source
        tab2_SRC_vlamap_contour = ColumnDataSource(
            data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})
        tab2_SRC_vlamap_peak = ColumnDataSource(
            data={'dspec': [], 'shape_longitude': [], 'shape_latitude': [], 'peak': []})
        tab2_SRC_maxfit_centroid_init(dspecDF_select)

        # initial the VLA map contour source
        tab2_SRC_vlamap_contour = ColumnDataSource(
            data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})
        tab2_SRC_vlamap_peak = ColumnDataSource(
            data={'dspec': [], 'shape_longitude': [], 'shape_latitude': [], 'peak': []})

        # import the vla image
        hdu = read_fits(vlafile[0])
        hdu_goodchan = goodchan(hdu)
        vla_local_pfmap = PuffinMap(hdu.data[0, hdu_goodchan[0], :, :], hdu.header,
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                    webgl=config_plot['plot_config']['WebGL'])
        # plot the contour of vla image
        mapx, mapy = vla_local_pfmap.meshgrid()
        mapx, mapy = mapx.value, mapy.value
        mapvlasize = mapy.shape
        ImgDF0 = pd.DataFrame({'xx': mapx.ravel(), 'yy': mapy.ravel()})
        tab2_SRC_vla_square = ColumnDataSource(ImgDF0)
        tab2_SRC_vlamap_contour = get_contour_data(mapx, mapy, vla_local_pfmap.smap.data)
        colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
        bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        tab2_SRC_ImgRgn_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))

        aiamap = sdomapfromlocalfile(wavelength='171', jdtime=xx[0] / 3600. / 24.)
        colormap = cm.get_cmap("sdoaia171")  # choose any matplotlib colormap here
        bokehpalette_sdoaia171 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        colormap = cm.get_cmap("sdoaia94")  # choose any matplotlib colormap here
        bokehpalette_sdoaia94 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        colormap = cm.get_cmap("sdoaia131")  # choose any matplotlib colormap here
        bokehpalette_sdoaia131 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        colormap = cm.get_cmap("gray")  # choose any matplotlib colormap here
        bokehpalette_gray = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        lengthx = vla_local_pfmap.dw[0] * u.arcsec
        lengthy = vla_local_pfmap.dh[0] * u.arcsec
        x0 = vla_local_pfmap.smap.center.x
        y0 = vla_local_pfmap.smap.center.y
        aiamap_submap = aiamap.submap(u.Quantity([x0 - lengthx / 2, x0 + lengthx / 2]),
                                      u.Quantity([y0 - lengthy / 2, y0 + lengthy / 2]))
        MapRES = 256
        dimensions = u.Quantity([MapRES, MapRES], u.pixel)
        aia_resampled_map = aiamap.resample(dimensions)

        # plot the global AIA image

        aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'],
                                        webgl=config_plot['plot_config']['WebGL'])

        tab2_p_aia, tab2_r_aia = aia_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                             palette=bokehpalette_sdoaia171)
        tab2_p_aia.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_aia.circle(x='shape_longitude', y='shape_latitude',  # size=10.*dspecDFselect.loc[76,:]['peak']/50.,
                          radius=3, radius_units='data', source=tab2_SRC_vlamap_peak, fill_alpha=0.8,
                          fill_color='#7c7e71',
                          line_color='#7c7e71')
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

        # plot the detail AIA image
        aia_submap_pfmap = PuffinMap(smap=aiamap_submap,
                                     plot_height=config_plot['plot_config']['tab_FSview_FitANLYS']['aia_submap_hght'],
                                     plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['aia_submap_wdth'],
                                     webgl=config_plot['plot_config']['WebGL'])

        tab2_SRC_aia_submap_square = ColumnDataSource(ImgDF0)
        tab3_p_aia_submap, tab3_r_aia_submap = aia_submap_pfmap.PlotMap(DrawLimb=True, DrawGrid=True,
                                                                        grid_spacing=20 * u.deg,
                                                                        title='EM sources centroid map',
                                                                        palette=bokehpalette_gray)
        tab2_r_aia_submap_square = tab3_p_aia_submap.square('xx', 'yy', source=tab2_SRC_aia_submap_square,
                                                            fill_alpha=0.0, fill_color=None,
                                                            line_color=None, line_alpha=0.0, selection_fill_alpha=0.5,
                                                            selection_fill_color=None,
                                                            nonselection_fill_alpha=0.0,
                                                            selection_line_alpha=0.0, selection_line_color=None,
                                                            nonselection_line_alpha=0.0,
                                                            size=4)
        tab3_p_aia_submap.add_tools(BoxSelectTool(renderers=[tab2_r_aia_submap_square]))

        # tab3_p_aia_submap.border_fill_color = "silver"
        tab3_p_aia_submap.border_fill_alpha = 0.4
        tab3_p_aia_submap.axis.major_tick_out = 0
        tab3_p_aia_submap.axis.major_tick_in = 5
        tab3_p_aia_submap.axis.minor_tick_out = 0
        tab3_p_aia_submap.axis.minor_tick_in = 3
        tab3_p_aia_submap.axis.major_tick_line_color = "white"
        tab3_p_aia_submap.axis.minor_tick_line_color = "white"
        color_mapper = LinearColorMapper(Spectral11)

        tab3_r_aia_submap_cross = tab3_p_aia_submap.cross(x='shape_longitude', y='shape_latitude', size=15,
                                                          color={'field': 'freq', 'transform': color_mapper},
                                                          line_width=3,
                                                          source=SRC_maxfit_centroid[tab2_dtim[0]], line_alpha=0.8)
        tab3_r_aia_submap_line = tab3_p_aia_submap.line(x='shape_longitude', y='shape_latitude', line_width=3,
                                                        line_color='black',
                                                        line_alpha=0.5,
                                                        source=SRC_maxfit_centroid[tab2_dtim[0]])
        tab3_r_aia_submap_line.visible = False
        tab3_SRC_aia_submap_rect = ColumnDataSource({'x': [], 'y': [], 'width': [], 'height': []})
        tab3_r_aia_submap_rect = tab3_p_aia_submap.rect(x='x', y='y', width='width', height='height', fill_alpha=0.1,
                                                        line_color='black', fill_color='black',
                                                        source=tab3_SRC_aia_submap_rect)

        tab2_Select_aia_wave = Select(title="Wavelenght:", value='171', options=['94', '131', '171'],
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_Select_aia_wave.on_change('value', aia_submap_wavelength_selection)
        colormap = cm.get_cmap("gray")  # choose any matplotlib colormap here
        bokehpalette_sdohmimag = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        hmimap = aiamap
        # todo fix this bug
        hmi_resampled_map = hmimap.resample(dimensions)
        hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                        webgl=config_plot['plot_config']['WebGL'])

        tab2_p_hmi, tab2_r_hmi = hmi_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                             x_range=tab2_p_aia.x_range,
                                                             y_range=tab2_p_aia.y_range,
                                                             palette=bokehpalette_sdohmimag)
        tab2_p_hmi.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_hmi.circle(x='shape_longitude', y='shape_latitude', radius=3, radius_units='data',
                          source=tab2_SRC_vlamap_peak,
                          fill_alpha=0.8,
                          fill_color='#7c7e71', line_color='#7c7e71')
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
        tab2_r_vla_square = tab2_p_vla.square('xx', 'yy', source=tab2_SRC_vla_square,
                                              fill_alpha=0.0, fill_color=None,
                                              line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                              selection_fill_color=None,
                                              nonselection_fill_alpha=0.0,
                                              selection_line_alpha=0.0, selection_line_color=None,
                                              nonselection_line_alpha=0.0,
                                              size=4)
        tab2_p_vla.add_tools(BoxSelectTool(renderers=[tab2_r_vla_square]))
        tab2_r_vla_multi_line = tab2_p_vla.multi_line(xs='xs', ys='ys', line_color='line_color',
                                                      source=tab2_SRC_vlamap_contour, alpha=0.7, line_width=2)

        tab2_r_vla_ImgRgn_patch = tab2_p_vla.patch('xx', 'yy', source=tab2_SRC_ImgRgn_Patch,
                                                   fill_color=None, fill_alpha=0.5, line_color="white",
                                                   line_alpha=1, line_width=1)

        tab2_LinkImg_HGHT = config_plot['plot_config']['tab_FSview_base']['vla_hght']
        tab2_LinkImg_WDTH = config_plot['plot_config']['tab_FSview_base']['vla_wdth']

        tab2_Div_LinkImg_plot = Div(text=""" """,
                                    width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_Slider_time_LinkImg = Slider(start=0, end=tab2_ntim - 1, value=0, step=1, title="time idx",
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                          callback_policy='mouseup')
        tab2_Slider_freq_LinkImg = Slider(start=0, end=tab2_nfreq - 1, value=0, step=1, title="freq idx",
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                          callback_policy='mouseup')

        stokeslist = ['{}'.format(int(ll)) for ll in
                      (hdu.header["CRVAL4"] + np.arange(hdu.header["NAXIS4"]) * hdu.header["CDELT4"])]
        stokesdict = {'1': 'I', '2': 'Q', '3': 'U', '4': 'V', '-1': 'RR', '-2': 'LL', '-3': 'RL', '-4': 'LR',
                      '-5': 'XX', '-6': 'YY', '-7': 'XY', '-8': 'YX'}
        pols = map(lambda x: stokesdict[x], stokeslist)
        # pols = ['RR', 'LL', 'I', 'V']
        SRL = set(['RR', 'LL'])
        SXY = set(['XX', 'YY', 'XY', 'YX'])
        Spol = set(pols)
        if hdu.header['NAXIS4'] == 2 and len(SRL.intersection(Spol)) == 2:
            pols = pols + ['I', 'V']
        if hdu.header['NAXIS4'] == 4 and len(SXY.intersection(Spol)) == 4:
            pols = pols + ['I', 'V']

        tab2_Select_vla_pol = Select(title="Polarization:", value=pols[0], options=pols,
                                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                line_color='white', source=tab2_source_idx_line_x)
        tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                line_color='white', source=tab2_source_idx_line_y)

        tab2_CTRLs_LinkImg = [tab2_Slider_time_LinkImg, tab2_Slider_freq_LinkImg, tab2_Select_vla_pol]
        for ctrl in tab2_CTRLs_LinkImg:
            ctrl.on_change('value', tab3_slider_LinkImg_update)

        tab2_LinkImg_HGHT = config_plot['plot_config']['tab_FSview_base']['vla_hght']
        tab2_LinkImg_WDTH = config_plot['plot_config']['tab_FSview_base']['vla_wdth']

        tab2_panel2_BUT_exit = Button(label='Exit FSview',
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                      button_type='danger')
        tab2_panel2_BUT_exit.on_click(tab2_panel_exit)

        tab2_panel3_BUT_exit = Button(label='Exit FSview',
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                      button_type='danger')
        tab2_panel3_BUT_exit.on_click(tab2_panel_exit)

        tab3_p_dspec_vector = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                     plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                     plot_height=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_hght'],
                                     x_range=(tab2_dtim[0], tab2_dtim[-1]),
                                     y_range=(tab2_freq[0], tab2_freq[-1]), toolbar_location='above')
        tab3_p_dspec_vectorx = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                      plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                      plot_height=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_hght'],
                                      x_range=tab3_p_dspec_vector.x_range,
                                      y_range=tab3_p_dspec_vector.y_range, toolbar_location='above')
        tab3_p_dspec_vectory = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                      plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                      plot_height=config_plot['plot_config']['tab_FSview_FitANLYS'][
                                                      'dspec_small_hght'] + 40,
                                      x_range=tab3_p_dspec_vector.x_range,
                                      y_range=tab3_p_dspec_vector.y_range, toolbar_location='above')
        tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        tab3_p_dspec_vector.xaxis.visible = False
        tab3_p_dspec_vectorx.xaxis.visible = False
        tab3_p_dspec_vector.title.text = "Vector Dynamic spectrum (Intensity)"
        tab3_p_dspec_vectorx.title.text = "Vector Dynamic spectrum (Vx)"
        tab3_p_dspec_vectory.title.text = "Vector Dynamic spectrum (Vy)"
        tab3_p_dspec_vectory.xaxis.axis_label = 'Seconds since ' + tim0_char
        tab3_p_dspec_vector.yaxis.axis_label = 'Frequency [GHz]'
        tab3_p_dspec_vectorx.yaxis.axis_label = 'Frequency [GHz]'
        tab3_p_dspec_vectory.yaxis.axis_label = 'Frequency [GHz]'
        # tab3_p_dspec_vector.border_fill_color = "silver"
        tab3_p_dspec_vector.border_fill_alpha = 0.4
        tab3_p_dspec_vector.axis.major_tick_out = 0
        tab3_p_dspec_vector.axis.major_tick_in = 5
        tab3_p_dspec_vector.axis.minor_tick_out = 0
        tab3_p_dspec_vector.axis.minor_tick_in = 3
        tab3_p_dspec_vector.axis.major_tick_line_color = "white"
        tab3_p_dspec_vector.axis.minor_tick_line_color = "white"
        # tab3_p_dspec_vectorx.border_fill_color = "silver"
        tab3_p_dspec_vectorx.border_fill_alpha = 0.4
        tab3_p_dspec_vectorx.axis.major_tick_out = 0
        tab3_p_dspec_vectorx.axis.major_tick_in = 5
        tab3_p_dspec_vectorx.axis.minor_tick_out = 0
        tab3_p_dspec_vectorx.axis.minor_tick_in = 3
        tab3_p_dspec_vectorx.axis.major_tick_line_color = "white"
        tab3_p_dspec_vectorx.axis.minor_tick_line_color = "white"
        # tab3_p_dspec_vectory.border_fill_color = "silver"
        tab3_p_dspec_vectory.border_fill_alpha = 0.4
        tab3_p_dspec_vectory.axis.major_tick_out = 0
        tab3_p_dspec_vectory.axis.major_tick_in = 5
        tab3_p_dspec_vectory.axis.minor_tick_out = 0
        tab3_p_dspec_vectory.axis.minor_tick_in = 3
        tab3_p_dspec_vectory.axis.major_tick_line_color = "white"
        tab3_p_dspec_vectory.axis.minor_tick_line_color = "white"
        tab3_p_dspec_vector.add_tools(BoxSelectTool())
        tab3_p_dspec_vector.add_tools(LassoSelectTool())
        tab3_p_dspec_vector.select(BoxSelectTool).select_every_mousemove = False
        tab3_p_dspec_vector.select(LassoSelectTool).select_every_mousemove = False
        tab2_SRC_aia_submap_square.on_change('selected', tab2_aia_submap_square_selection_change)

        VdspecDF_init()
        VdspecDF_update()
        tab3_SRC_dspec_vector_init()
        tab3_p_dspec_vector.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                  dh=tab2_freq[-1] - tab2_freq[0],
                                  palette=bokehpalette_jet, source=tab3_SRC_dspec_vector)
        tab3_p_dspec_vectorx.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                   dh=tab2_freq[-1] - tab2_freq[0],
                                   palette=bokehpalette_jet, source=tab3_SRC_dspec_vectorx)
        tab3_p_dspec_vectory.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                   dh=tab2_freq[-1] - tab2_freq[0],
                                   palette=bokehpalette_jet, source=tab3_SRC_dspec_vectory)
        tab3_source_idx_line = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab3_r_dspec_vector_line = tab3_p_dspec_vector.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                            line_color='white', source=tab3_source_idx_line)
        tab3_r_dspec_vectorx_line = tab3_p_dspec_vectorx.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                              line_color='white',
                                                              source=tab3_source_idx_line)
        tab3_r_dspec_vectory_line = tab3_p_dspec_vectory.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                              line_color='white',
                                                              source=tab3_source_idx_line)
        tab2_SRC_dspec_vector_square = ColumnDataSource(dspecDF_frac)
        tab2_r_dspec_vector_square = tab3_p_dspec_vector.square('time', 'freq', source=tab2_SRC_dspec_vector_square,
                                                                fill_color=colors_dspec,
                                                                fill_alpha=0.0,
                                                                line_color=None, line_alpha=0.0,
                                                                selection_fill_alpha=0.2,
                                                                selection_fill_color='black',
                                                                nonselection_fill_alpha=0.0,
                                                                selection_line_alpha=0.0, selection_line_color='white',
                                                                nonselection_line_alpha=0.0,
                                                                size=max(
                                                                    config_plot['plot_config']['tab_FSview_FitANLYS'][
                                                                        'dspec_small_wdth'] / tab2_ntim,
                                                                    config_plot['plot_config']['tab_FSview_FitANLYS'][
                                                                        'dspec_small_hght'] / tab2_nfreq))

        tab2_dspec_selected = None
        tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)
        tab2_vla_square_selected = None
        tab2_SRC_vla_square.on_change('selected', tab2_vla_square_selection_change)
        tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                    options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                             "2048x2048"],
                                    width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_Select_MapRES.on_change('value', tab2_update_MapRES)
        rgnfitsfile = database_dir + event_id + struct_id + "CASA_imfit_region_fits.rgn"

        tab2_BUT_SavRgn = Button(label='Save Region',
                                 width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                 button_type='primary')
        tab2_BUT_SavRgn.on_click(tab2_save_region)
        tab2_dspec_vector_selected = None
        tab2_SRC_dspec_vector_square.on_change('selected', tab2_dspec_vector_selection_change)
        tab3_dspec_small_CTRLs_OPT = dict(mean_values=[mean_amp_g, mean_vx, mean_vy],
                                          drange_values=[drange_amp_g, drange_vx, drange_vy],
                                          vmax_values=[vmax_amp_g, vmax_vx, vmax_vy],
                                          vmin_values=[vmin_amp_g, vmin_vx, vmin_vy],
                                          vmax_values_last=[vmax_amp_g, vmax_vx, vmax_vy],
                                          vmin_values_last=[vmin_amp_g, vmin_vx, vmin_vy],
                                          items_dspec_small=['peak', 'shape_longitude', 'shape_latitude'],
                                          labels_dspec_small=["Flux", "X-pos", "Y-pos"], idx_p_dspec_small=0,
                                          radio_button_group_dspec_small_update_flag=False)

        tab3_RBG_dspec_small = RadioButtonGroup(labels=tab3_dspec_small_CTRLs_OPT['labels_dspec_small'], active=0)
        tab3_BUT_dspec_small_reset = Button(label='Reset DRange',
                                            width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab3_Slider_dspec_small_dmax = Slider(start=mean_amp_g, end=mean_amp_g + 2 * drange_amp_g, value=vmax_amp_g,
                                              step=1,
                                              title='dmax', callback_throttle=250)
        tab3_Slider_dspec_small_dmin = Slider(start=mean_amp_g - 2 * drange_amp_g, end=mean_amp_g, value=vmin_amp_g,
                                              step=1,
                                              title='dmin', callback_throttle=250)
        tab3_RBG_dspec_small.on_change('active', tab3_RBG_dspec_small_update)
        tab3_BUT_dspec_small_reset.on_click(tab3_BUT_dspec_small_reset_update)
        tab3_BUT_dspec_small_resetall = Button(label='Reset All',
                                               width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab3_BUT_dspec_small_resetall.on_click(tab3_BUT_dspec_small_resetall_update)
        tab3_CTRLs_dspec_small = [tab3_Slider_dspec_small_dmax, tab3_Slider_dspec_small_dmin]
        for ctrl in tab3_CTRLs_dspec_small:
            ctrl.on_change('value', tab3_slider_dspec_small_update)

        tab3_RBG_TimeFreq = RadioButtonGroup(labels=["time", "freq"], active=0)
        tab3_Slider_ANLYS_idx = Slider(start=0, end=tab2_ntim - 1, value=0, step=1, title="time idx", width=450)
        tab3_Slider_ANLYS_idx.on_change('value', tab3_slider_ANLYS_idx_update)
        tab3_Div_Tb = Div(text=""" """, width=400)
        timebin = 1
        timeline = False
        tab3_animate_step = timebin
        tab3_BUT_PlayCTRL = Button(label='Play', width=60, button_type='success')
        tab3_BUT_PlayCTRL.on_click(tab3_animate)
        tab3_BUT_StepCTRL = Button(label='Step', width=60, button_type='primary')
        tab3_BUT_StepCTRL.on_click(tab3_animate_step_CTRL)
        tab3_BUT_FRWD_REVS_CTRL = Button(label='Forward', width=60, button_type='warning')
        tab3_BUT_FRWD_REVS_CTRL.on_click(tab3_animate_FRWD_REVS)
        tab3_BUT_animate_ONOFF = Button(label='Animate ON & Go', width=80)
        tab3_BUT_animate_ONOFF.on_click(tab3_animate_onoff)
        tab3_Div_plot_xargs = Div(text='', width=300)
        tab3_BUT_plot_xargs_default()
        tab3_SPCR_LFT_BUT_Step = Spacer(width=10, height=10)
        tab3_SPCR_LFT_BUT_REVS_CTRL = Spacer(width=10, height=10)
        tab3_SPCR_LFT_BUT_animate_ONOFF = Spacer(width=20, height=10)
        tab3_input_plot_xargs = TextInput(value='Input the param here', title="Plot parameters:", width=300)
        # todo add RCP LCP check box
        tab3_CheckboxGroup_pol = CheckboxGroup(labels=["RCP", "LCP"], active=[0, 1])

        tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                    options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                             "2048x2048"],
                                    width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_Select_MapRES.on_change('value', tab2_update_MapRES)

        lout2_1 = row(gridplot([[tab2_p_aia, tab2_p_hmi, tab2_p_vla]], toolbar_location='right'),
                      widgetbox(tab2_Select_MapRES, tab2_Select_vla_pol, tab2_Slider_time_LinkImg,
                                tab2_Slider_freq_LinkImg, tab2_BUT_vdspec, tab2_BUT_SavRgn, tab2_Div_LinkImg_plot,
                                width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth']))
        if do_spec_regrid:
            # lout2_2_1 = column(tab2_p_dspec_rs, row(tab2_p_dspec, tab2_p_dspec_yPro), tab2_p_dspec_xPro)
            pass
        else:
            lout2_2_1 = column(row(tab2_p_dspec, tab2_p_dspec_yPro), tab2_p_dspec_xPro)
        lout2_2_2 = widgetbox(tab2_Select_pol, tab2_Select_bl,
                              tab2_Select_colorspace,
                              tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                              width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        lout2_2 = row(lout2_2_1, lout2_2_2)
        panel2 = column(lout2_1, lout2_2)

        lout3_1 = column(tab3_p_aia_submap, tab3_Slider_ANLYS_idx,
                         row(tab3_BUT_PlayCTRL, tab3_SPCR_LFT_BUT_Step, tab3_BUT_StepCTRL,
                             tab3_SPCR_LFT_BUT_REVS_CTRL,
                             tab3_BUT_FRWD_REVS_CTRL, tab3_SPCR_LFT_BUT_animate_ONOFF,
                             tab3_BUT_animate_ONOFF), tab3_input_plot_xargs, tab3_Div_plot_xargs)
        lout3_2 = column(gridplot([tab3_p_dspec_vector], [tab3_p_dspec_vectorx], [tab3_p_dspec_vectory],
                                  toolbar_location='right'), tab3_Div_Tb)
        lout3_3 = widgetbox(tab3_RBG_dspec_small, tab3_Slider_dspec_small_dmax, tab3_Slider_dspec_small_dmin,
                            tab3_BUT_dspec_small_reset, tab3_BUT_dspec_small_resetall, tab2_Select_aia_wave,
                            tab2_panel3_BUT_exit,
                            tab2_panel3_Div_exit,
                            width=200)
        panel3 = row(lout3_1, lout3_2, lout3_3)
        # tab2 = Panel(child=panel2, title="FS View")
        # tab3 = Panel(child=panel3, title="FitANLYS")
        #
        # tabs_top = Tabs(tabs=[tab2, tab3])

        lout = column(panel2, panel3)

        # def timeout_callback():
        #     print 'timeout'
        #     raise SystemExit


        curdoc().add_root(lout)
        curdoc().title = "FSview"
    else:
        '''
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        #################################### FS_view_prep ######################################
        #######################┏━━━┓┏━━━┓#######################################################
        #######################┃┏━━┛┃┏━┓┃#######################################################
        #######################┃┗━━┓┃┗━━┓┏┓┏┓┏┓┏━━┓┏┓┏┓┏┓┏━━┓┏━┓┏━━┓┏━━┓########################
        #######################┃┏━━┛┗━━┓┃┃┗┛┃┣┫┃┃━┫┃┗┛┗┛┃┃┏┓┃┃┏┛┃┃━┫┃┏┓┃########################
        #######################┃┃╋╋╋┃┗━┛┃┗┓┏┛┃┃┃┃━┫┗┓┏┓┏┛┃┗┛┃┃┃╋┃┃━┫┃┗┛┃########################
        #######################┗┛╋╋╋┗━━━┛╋┗┛╋┗┛┗━━┛╋┗┛┗┛╋┃┏━┛┗┛╋┗━━┛┃┏━┛########################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        ########################################################################################
        '''
        vlafile = glob.glob(fits_LOCL_dir + '*.fits')
        if len(vlafile) > 0:
            with open(database_dir + event_id + struct_id + 'CASA_CLN_args.json', 'rb') as fp:
                tab2_tCLN_Param_dict = json.load(fp)
            tab2_panel2_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>\
                                    first before closing the tab</p></b>""",
                                       width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
            rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()
            colors_dspec = [colors.rgb2hex(m) for m in colormap_jet((tab2_spec_plt.flatten() - rmin) / (rmax - rmin))]

            '''create the regridded dynamic spectrum plot'''
            TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
            if tab2_ntim > spec_rs_tmax or tab2_nfreq > spec_rs_fmax:
                do_spec_regrid = True
                tab2_dspec_fs_ntim = spec_rs_tmax if spec_rs_tmax < tab2_ntim else tab2_ntim
                tab2_dspec_fs_nfreq = spec_rs_fmax if spec_rs_fmax < tab2_nfreq else tab2_nfreq
                dspecDF_frac = \
                    dspecDF0[dspecDF0.time < tab2_dtim[0] + tab2_dspec_fs_ntim * tab2_dt][
                        dspecDF0.time >= tab2_dtim[0]][
                        dspecDF0.freq >= tab2_freq[0]][dspecDF0.freq < tab2_freq[0] + tab2_dspec_fs_nfreq * tab2_df]
                rebin_specdata(tab2_spec, bl_index, tab2_pol)
                tab2_p_dspec_rs = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                                         plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_rs_wdth'],
                                         plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_rs_hght'],
                                         x_range=(tab2_dtim[0], tab2_dtim[-1]), y_range=(tab2_freq[0], tab2_freq[-1]),
                                         toolbar_location="above")
                tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
                tab2_p_dspec_rs.axis.visible = True
                tab2_p_dspec_rs.title.text = "Dynamic spectrum"
                tab2_p_dspec_rs.xaxis.axis_label = 'Seconds since ' + tim0_char
                tab2_p_dspec_rs.yaxis.axis_label = 'Frequency [GHz]'
                tab2_SRC_dspec_image_rs = ColumnDataSource(
                    data={'data': [tab2_spec_plt], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
                tab2_p_dspec_rs.image(image="data", x=tab2_dtim[0], y=tab2_freq[0],
                                      dw=tab2_dtim[-1] - tab2_dtim[0],
                                      dh=tab2_freq[-1] - tab2_freq[0],
                                      source=tab2_SRC_dspec_image_rs, palette=bokehpalette_jet)
                tim_map_square_rs = (np.tile(tab2_tim_square_rs, tab2_nfreq_square_rs).reshape(tab2_nfreq_square_rs, \
                                                                                               tab2_ntim_square_rs) / 3600. / 24. + 2400000.5) * 86400.
                freq_map_square_rs = np.tile(tab2_freq_square_rs, tab2_ntim_square_rs).reshape(tab2_ntim_square_rs, \
                                                                                               tab2_nfreq_square_rs).swapaxes(
                    0, 1)
                dspecDF_rs0 = pd.DataFrame({'xx': tim_map_square_rs.ravel() - xx[0], 'yy': freq_map_square_rs.ravel()})
                tab2_SRC_dspec_square_rs = ColumnDataSource(dspecDF_rs0)
                tab2_r_square_rs = tab2_p_dspec_rs.square('xx', 'yy', source=tab2_SRC_dspec_square_rs, fill_color=None,
                                                          fill_alpha=0.0,
                                                          line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                                          selection_fill_color='black',
                                                          nonselection_fill_alpha=0.0,
                                                          selection_line_alpha=0.0, selection_line_color='white',
                                                          nonselection_line_alpha=0.0,
                                                          size=min(
                                                              config_plot['plot_config']['tab_FSview_base'][
                                                                  'dspec_wdth'] / tab2_ntim_square_rs,
                                                              config_plot['plot_config']['tab_FSview_base'][
                                                                  'dspec_hght'] / tab2_nfreq_square_rs))
                tab2_SRC_r_square_rs_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
                tab2_r_square_rs_patch = tab2_p_dspec_rs.patch('xx', 'yy', source=tab2_SRC_r_square_rs_Patch,
                                                               fill_color=None, fill_alpha=0.5, line_color="white",
                                                               line_alpha=1, line_width=1)
                tab2_p_dspec_rs.add_tools(BoxSelectTool(renderers=[tab2_r_square_rs]))
                tab2_SRC_dspec_square_rs.on_change('selected', tab2_update_dspec_rs_image)

            tab2_dtim_fs = pd.Series.unique(dspecDF_frac['time'])
            tab2_freq_fs = pd.Series.unique(dspecDF_frac['freq'])
            tab2_ntim_fs = len(tab2_dtim_fs)
            tab2_nfreq_fs = len(tab2_freq_fs)
            tab2_tim_ind0 = np.where(abs(tab2_dtim - tab2_dtim_fs[0]) < tab2_dt / 2.0)[0][0]
            tab2_tim_ind1 = np.where(abs(tab2_dtim - tab2_dtim_fs[-1]) < tab2_dt / 2.0)[0][0]
            tab2_freq_ind0 = np.where(abs(tab2_freq - tab2_freq_fs[0]) < tab2_df / 2.0)[0][0]
            tab2_freq_ind1 = np.where(abs(tab2_freq - tab2_freq_fs[-1]) < tab2_df / 2.0)[0][0]
            '''create the dynamic spectrum plot'''
            TOOLS = "crosshair,pan,wheel_zoom,tap,box_zoom,reset,save"
            tab2_SRC_dspec_square = ColumnDataSource(dspecDF_frac)
            tab2_p_dspec = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                                  plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_wdth'],
                                  plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_hght'],
                                  x_range=(dspecDF_frac['time'].min(), dspecDF_frac['time'].max()),
                                  y_range=(dspecDF_frac['freq'].min(), dspecDF_frac['freq'].max()),
                                  toolbar_location="above")
            tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
            tab2_p_dspec.axis.visible = True
            tab2_p_dspec.title.text = "Dynamic spectrum"
            tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
            tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
            tab2_SRC_dspec_image = ColumnDataSource(
                data={'data': [tab2_spec_plt[tab2_freq_ind0:tab2_freq_ind1 + 1, tab2_tim_ind0:tab2_tim_ind1 + 1]],
                      'xx': [tab2_dtim_fs], 'yy': [tab2_freq_fs]})
            tab2_p_dspec.image(image="data", x=tab2_dtim_fs[0], y=tab2_freq_fs[0],
                               dw=tab2_dtim_fs[-1] - tab2_dtim_fs[0],
                               dh=tab2_freq_fs[-1] - tab2_freq_fs[0],
                               source=tab2_SRC_dspec_image, palette=bokehpalette_jet)

            # make the dspec data source selectable
            tab2_r_square = tab2_p_dspec.square('time', 'freq', source=tab2_SRC_dspec_square, fill_color=colors_dspec,
                                                fill_alpha=0.0,
                                                line_color=None, line_alpha=0.0, selection_fill_alpha=0.1,
                                                selection_fill_color='black',
                                                nonselection_fill_alpha=0.0,
                                                selection_line_alpha=0.2, selection_line_color='white',
                                                nonselection_line_alpha=0.0,
                                                size=max(
                                                    config_plot['plot_config']['tab_FSview_base'][
                                                        'dspec_wdth'] / tab2_ntim_fs,
                                                    config_plot['plot_config']['tab_FSview_base'][
                                                        'dspec_hght'] / tab2_nfreq_fs))

            # tab2_p_dspec.border_fill_color = "silver"
            tab2_p_dspec.border_fill_alpha = 0.4
            tab2_p_dspec.axis.major_tick_out = 0
            tab2_p_dspec.axis.major_tick_in = 5
            tab2_p_dspec.axis.minor_tick_out = 0
            tab2_p_dspec.axis.minor_tick_in = 3
            tab2_p_dspec.axis.major_tick_line_color = "white"
            tab2_p_dspec.axis.minor_tick_line_color = "white"

            tab2_Select_pol = Select(title="Polarization:", value='I', options=['RR', 'LL', 'I', 'V'],
                                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
            tab2_Select_bl = Select(title="Baseline:", value=tab2_bl[0], options=tab2_bl,
                                    width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
            tab2_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                                            width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_p_dspec_xPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                                       plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_xPro_wdth'],
                                       plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_xPro_hght'],
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
            tab2_l_dspec_xPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                               source=tab2_SRC_dspec_xPro_hover,
                                               render_mode='canvas')
            tab2_l_dspec_xPro_hover.text_font_size = '5pt'
            tab2_p_dspec_xPro.add_layout(tab2_l_dspec_xPro_hover)
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

            tab2_p_dspec_yPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                                       plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_yPro_wdth'],
                                       plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_yPro_hght'],
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
            l_dspec_yPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                          source=tab2_SRC_dspec_yPro_hover,
                                          render_mode='canvas')
            l_dspec_yPro_hover.text_font_size = '5pt'
            tab2_p_dspec_yPro.add_layout(l_dspec_yPro_hover)
            tab2_p_dspec_yPro.title.text_font_size = '6pt'
            tab2_p_dspec_yPro.yaxis.visible = False
            tab2_p_dspec_yPro.background_fill_color = "beige"
            tab2_p_dspec_yPro.background_fill_alpha = 0.4
            tab2_p_dspec_yPro.xaxis.axis_label = 'Intensity [sfu]'
            tab2_p_dspec_yPro.yaxis.axis_label_text_font_size = '5px'
            # tab2_p_dspec_yPro.border_fill_color = "silver"
            tab2_p_dspec_yPro.border_fill_alpha = 0.4
            tab2_p_dspec_yPro.min_border_bottom = 0
            tab2_p_dspec_yPro.min_border_left = 0
            tab2_p_dspec_yPro.axis.major_tick_out = 0
            tab2_p_dspec_yPro.axis.major_tick_in = 5
            tab2_p_dspec_yPro.axis.minor_tick_out = 0
            tab2_p_dspec_yPro.axis.minor_tick_in = 3
            tab2_p_dspec_yPro.axis.major_tick_line_color = "black"
            tab2_p_dspec_yPro.axis.minor_tick_line_color = "black"

            tab2_BUT_vdspec = Button(label="VEC Dyn Spec",
                                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                     button_type="success")
            tab2_Select_pol_opt = ['RR', 'LL', 'I', 'V']
            tab2_BUT_vdspec.on_click(tab2_vdspec_update)
            tab2_ctrls = [tab2_Select_bl, tab2_Select_pol, tab2_Select_colorspace]
            for ctrl in tab2_ctrls:
                ctrl.on_change('value', tab2_update_dspec_image)

            # # Add a hover tool
            tooltips = None

            hover_JScode = """
                var nx = %d;
                var ny = %d;
                var data = {'x': [], 'y': []};
                var cdata = rs.get('data');
                var indices = cb_data.index['1d'].indices;
                var idx_offset = indices[0] - (indices[0] %% nx);
                for (i=0; i < nx; i++) {
                    data['x'].push(cdata.time[i+idx_offset]);
                    data['y'].push(cdata.dspec[i+idx_offset]);
                }
                rdx.set('data', data);
                idx_offset = indices[0] %% nx;
                data = {'x': [], 'y': []};
                for (i=0; i < ny; i++) {
                    data['x'].push(cdata.dspec[i*nx+idx_offset]);
                    data['y'].push(cdata.freq[i*nx+idx_offset]);
                }
                rdy.set('data', data);
                var time = cdata.timestr[indices[0]]+' '
                var freq = cdata.freq[indices[0]].toFixed(3)+'[GHz] '
                var dspec = cdata.dspec[indices[0]].toFixed(3)+ '[sfu]'
                var tooltips = freq + time + dspec
                data = {'x': [], 'y': [], 'tooltips': []};
                data['x'].push(cdata.time[indices[0]]);
                data['y'].push(cdata.dspec[indices[0]]);
                data['tooltips'].push(tooltips);
                rdx_hover.set('data', data);
                tooltips = time + freq + dspec
                data = {'x': [], 'y': [], 'tooltips': []};
                data['x'].push(cdata.dspec[indices[0]]);
                data['y'].push(cdata.freq[indices[0]]);
                data['tooltips'].push(tooltips);
                rdy_hover.set('data', data);
                """ % (tab2_ntim, tab2_nfreq)

            tab2_p_dspec_hover_callback = CustomJS(
                args={'rs': tab2_r_square.data_source, 'rdx': r_dspec_xPro.data_source, 'rdy': r_dspec_yPro.data_source,
                      'rdx_hover': r_dspec_xPro_hover.data_source,
                      'rdy_hover': r_dspec_yPro_hover.data_source}, code=hover_JScode)
            tab2_p_dspec_hover = HoverTool(tooltips=tooltips, callback=tab2_p_dspec_hover_callback,
                                           renderers=[tab2_r_square])
            tab2_p_dspec.add_tools(tab2_p_dspec_hover)

            # initial the VLA map contour source
            tab2_SRC_vlamap_contour = ColumnDataSource(
                data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})
            tab2_SRC_vlamap_peak = ColumnDataSource(
                data={'dspec': [], 'shape_longitude': [], 'shape_latitude': [], 'peak': []})

            # import the vla image
            hdu = read_fits(vlafile[0])
            hdu_goodchan = goodchan(hdu)
            vla_local_pfmap = PuffinMap(hdu.data[0, hdu_goodchan[0], :, :], hdu.header,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                        webgl=config_plot['plot_config']['WebGL'])
            # plot the contour of vla image
            mapx, mapy = vla_local_pfmap.meshgrid()
            mapx, mapy = mapx.value, mapy.value
            mapvlasize = mapy.shape
            ImgDF0 = pd.DataFrame({'xx': mapx.ravel(), 'yy': mapy.ravel()})
            tab2_SRC_vla_square = ColumnDataSource(ImgDF0)
            tab2_SRC_vlamap_contour = get_contour_data(mapx, mapy, vla_local_pfmap.smap.data)
            colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
            bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            tab2_SRC_ImgRgn_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))

            colormap = cm.get_cmap("sdoaia171")  # choose any matplotlib colormap here
            bokehpalette_sdoaia171 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            # aiamap = sunpy.map.Map(filepath)
            aiamap = sdomapfromlocalfile(wavelength='171', jdtime=xx[0] / 3600. / 24.)
            MapRES = 256
            dimensions = u.Quantity([MapRES, MapRES], u.pixel)
            aia_resampled_map = aiamap.resample(dimensions)

            # plot the global AIA image
            aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'],
                                            webgl=config_plot['plot_config']['WebGL'])

            tab2_p_aia, tab2_r_aia = aia_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                                 palette=bokehpalette_sdoaia171)
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
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                            webgl=config_plot['plot_config']['WebGL'])

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
            tab2_r_vla_square = tab2_p_vla.square('xx', 'yy', source=tab2_SRC_vla_square,
                                                  fill_alpha=0.0, fill_color=None,
                                                  line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                                  selection_fill_color=None,
                                                  nonselection_fill_alpha=0.0,
                                                  selection_line_alpha=0.0, selection_line_color=None,
                                                  nonselection_line_alpha=0.0,
                                                  size=4)
            tab2_p_vla.add_tools(BoxSelectTool(renderers=[tab2_r_vla_square]))
            tab2_r_vla_multi_line = tab2_p_vla.multi_line(xs='xs', ys='ys', line_color='line_color',
                                                          source=tab2_SRC_vlamap_contour, alpha=0.7, line_width=2)

            tab2_r_vla_ImgRgn_patch = tab2_p_vla.patch('xx', 'yy', source=tab2_SRC_ImgRgn_Patch,
                                                       fill_color=None, fill_alpha=0.5, line_color="white",
                                                       line_alpha=1, line_width=1)

            tab2_LinkImg_HGHT = config_plot['plot_config']['tab_FSview_base']['vla_hght']
            tab2_LinkImg_WDTH = config_plot['plot_config']['tab_FSview_base']['vla_wdth']

            tab2_Div_LinkImg_plot = Div(text=""" """,
                                        width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_Slider_time_LinkImg = Slider(start=0, end=tab2_ntim - 1, value=0, step=1, title="time idx",
                                              width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                              callback_policy='mouseup')
            tab2_Slider_freq_LinkImg = Slider(start=0, end=tab2_nfreq - 1, value=0, step=1, title="freq idx",
                                              width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                              callback_policy='mouseup')

            stokeslist = ['{}'.format(int(ll)) for ll in
                          (hdu.header["CRVAL4"] + np.arange(hdu.header["NAXIS4"]) * hdu.header["CDELT4"])]
            stokesdict = {'1': 'I', '2': 'Q', '3': 'U', '4': 'V', '-1': 'RR', '-2': 'LL', '-3': 'RL', '-4': 'LR',
                          '-5': 'XX', '-6': 'YY', '-7': 'XY', '-8': 'YX'}
            pols = map(lambda x: stokesdict[x], stokeslist)
            # pols = ['RR', 'LL', 'I', 'V']
            SRL = set(['RR', 'LL'])
            SXY = set(['XX', 'YY', 'XY', 'YX'])
            Spol = set(pols)
            if hdu.header['NAXIS4'] == 2 and len(SRL.intersection(Spol)) == 2:
                pols = pols + ['I', 'V']
            if hdu.header['NAXIS4'] == 4 and len(SXY.intersection(Spol)) == 4:
                pols = pols + ['I', 'V']

            tab2_Select_vla_pol = Select(title="Polarization:", value=pols[0], options=pols,
                                         width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
            tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                    line_color='white', source=tab2_source_idx_line_x)
            tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
            tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                    line_color='white', source=tab2_source_idx_line_y)

            tab2_CTRLs_LinkImg = [tab2_Slider_time_LinkImg, tab2_Slider_freq_LinkImg, tab2_Select_vla_pol]
            for ctrl in tab2_CTRLs_LinkImg:
                ctrl.on_change('value', tab3_slider_LinkImg_update)

            tab2_panel2_BUT_exit = Button(label='Exit FSview',
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                          button_type='danger')
            tab2_panel2_BUT_exit.on_click(tab2_panel_exit)
            tab2_dspec_selected = None

            tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)
            tab2_vla_square_selected = None
            tab2_SRC_vla_square.on_change('selected', tab2_vla_square_selection_change)
            tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                        options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                                 "2048x2048"],
                                        width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_Select_MapRES.on_change('value', tab2_update_MapRES)
            rgnfitsfile = database_dir + event_id + struct_id + "CASA_imfit_region_fits.rgn"
            tab2_BUT_SavRgn = Button(label='Save Region',
                                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                     button_type='primary')
            tab2_BUT_SavRgn.on_click(tab2_save_region)
            tab2_input_tImfit = TextInput(value="Input the param here", title="pimfit task parameters:",
                                          width=config_plot['plot_config']['tab_FSviewPrep']['input_tCLN_wdth'])
            tab2_Div_tImfit = Div(text='', width=config_plot['plot_config']['tab_FSviewPrep']['tab2_Div_tImfit_wdth'])
            tab2_Div_tImfit2 = Div(text='', width=config_plot['plot_config']['tab_FSviewPrep']['input_tCLN_wdth'])

            tab2_BUT_tImfit_param_default()
            tab2_BUT_tImfit_param_ADD = Button(label='Add to Param',
                                               width=config_plot['plot_config']['tab_FSviewPrep']['button_wdth'],
                                               button_type='primary')
            tab2_BUT_tImfit_param_ADD.on_click(tab2_BUT_tImfit_param_add)
            tab2_BUT_tImfit_param_DEL = Button(label='Delete Param',
                                               width=config_plot['plot_config']['tab_FSviewPrep']['button_wdth'],
                                               button_type='warning')
            tab2_BUT_tImfit_param_DEL.on_click(tab2_BUT_tImfit_param_delete)
            tab2_BUT_tImfit_param_Default = Button(label='Default Param',
                                                   width=config_plot['plot_config']['tab_FSviewPrep']['button_wdth'])
            tab2_BUT_tImfit_param_Default.on_click(tab2_BUT_tImfit_param_default)
            tab2_SPCR_LFT_BUT_tImfit_param_DEL = Spacer(
                width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth10'])
            tab2_SPCR_LFT_BUT_tImfit_param_RELOAD = Spacer(
                width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth10'])
            tab2_SPCR_LFT_BUT_tImfit_param_DEFAULT = Spacer(
                width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth10'])

            tab2_BUT_tImfit_param_RELOAD = Button(label='reload Param',
                                                  width=config_plot['plot_config']['tab_FSviewPrep']['button_wdth'])
            tab2_BUT_tImfit_param_RELOAD.on_click(tab2_BUT_timfit_param_reload)
            tab2_SPCR_LFT_BUT_tImfit_param_reload = Spacer(
                width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth20'])

            tab2_BUT_tImfit = Button(label='imfit',
                                     width=config_plot['plot_config']['tab_FSviewPrep']['button_wdth'],
                                     button_type='success')
            tab2_BUT_tImfit.on_click(tab2_BUT_tImfit_update)
            tab2_SPCR_ABV_BUT_timfit = Spacer(width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth10'])

            tab2_SPCR_LFT_Div_tImfit = Spacer(width=config_plot['plot_config']['tab_FSviewPrep']['space_wdth10'])

            lout2_1_1 = row(gridplot([[tab2_p_aia, tab2_p_hmi, tab2_p_vla]], toolbar_location='right'),
                            widgetbox(tab2_Select_MapRES, tab2_Select_vla_pol, tab2_Slider_time_LinkImg,
                                      tab2_Slider_freq_LinkImg, tab2_BUT_vdspec, tab2_BUT_SavRgn, tab2_Div_LinkImg_plot,
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth']))
            if do_spec_regrid:
                lout2_1_2 = row(column(tab2_p_dspec_rs, row(tab2_p_dspec, tab2_p_dspec_yPro),
                                       tab2_p_dspec_xPro),
                                widgetbox(tab2_Select_pol, tab2_Select_bl,
                                          tab2_Select_colorspace,
                                          tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth']))
            else:
                lout2_1_2 = row(column(row(tab2_p_dspec, tab2_p_dspec_yPro),
                                       tab2_p_dspec_xPro),
                                widgetbox(tab2_Select_pol, tab2_Select_bl,
                                          tab2_Select_colorspace,
                                          tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth']))
            lout2_1 = column(lout2_1_1, lout2_1_2)
            lout2_2 = tab2_SPCR_LFT_Div_tImfit
            lout2_3 = column(tab2_input_tImfit, tab2_Div_tImfit2,
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
                width=config_plot['plot_config']['tab_FSview_base']['dspec_wdth'])
            lout = tab_panel_Div_info
            curdoc().add_root(lout)
            curdoc().title = "FSview"
else:
    '''
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    #################################### FSview2CASA #######################################
    ##################┏━━━┓┏━━━┓╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋┏━━━┓┏━━━┓┏━━━┓┏━━━┓┏━━━┓###################
    ##################┃┏━━┛┃┏━┓┃╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋╋┃┏━┓┃┃┏━┓┃┃┏━┓┃┃┏━┓┃┃┏━┓┃###################
    ##################┃┗━━┓┃┗━━┓┏┓┏┓┏┓┏━━┓┏┓┏┓┏┓┗┛┏┛┃┃┃╋┗┛┃┃╋┃┃┃┗━━┓┃┃╋┃┃###################
    ##################┃┏━━┛┗━━┓┃┃┗┛┃┣┫┃┃━┫┃┗┛┗┛┃┏━┛┏┛┃┃╋┏┓┃┗━┛┃┗━━┓┃┃┗━┛┃###################
    ##################┃┃╋╋╋┃┗━┛┃┗┓┏┛┃┃┃┃━┫┗┓┏┓┏┛┃┃┗━┓┃┗━┛┃┃┏━┓┃┃┗━┛┃┃┏━┓┃###################
    ##################┗┛╋╋╋┗━━━┛╋┗┛╋┗┛┗━━┛╋┗┛┗┛╋┗━━━┛┗━━━┛┗┛╋┗┛┗━━━┛┗┛╋┗┛###################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    ########################################################################################
    '''
    tab2_panel_Div_exit = Div(
        text="""<p><b>Warning</b>: Click the <b>Exit FSview2CASA</b>
                first before closing the tab</p></b>""",
        width=config_plot['plot_config']['tab_FSview2CASA']['widgetbox_wdth1'])
    timestrs = []
    for ii in range(len(xx)):
        t0 = xx[ii]  # -0.5*t_int
        timestr0 = Time(t0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        timestrs.append(timestr0.split(' ')[1])
    dspecDF0 = pd.DataFrame({'time': xx - xx[0],
                             'freq': yy,
                             'dspec': tab2_spec_plt.flatten(),
                             'timestr': timestrs})

    rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()
    colors_dspec = [colors.rgb2hex(m) for m in colormap_jet((tab2_spec_plt.flatten() - rmin) / (rmax - rmin))]

    TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"

    tab2_SRC_dspec_square = ColumnDataSource(dspecDF0)

    '''create the dynamic spectrum plot'''
    tab2_p_dspec = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                          plot_width=config_plot['plot_config']['tab_FSview2CASA']['dspec_wdth'],
                          plot_height=config_plot['plot_config']['tab_FSview2CASA']['dspec_hght'],
                          x_range=(tab2_dtim[0], tab2_dtim[-1]), y_range=(tab2_freq[0], tab2_freq[-1]),
                          toolbar_location="above")
    tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
    tab2_p_dspec.axis.visible = True
    tab2_p_dspec.title.text = "Dynamic spectrum"
    tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
    tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
    tab2_SRC_dspec_image = ColumnDataSource(data={'data': [tab2_spec_plt], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
    tab2_p_dspec.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                       dh=tab2_freq[-1] - tab2_freq[0],
                       source=tab2_SRC_dspec_image, palette=bokehpalette_jet)

    # make the dspec data source selectable
    tab2_r_square = tab2_p_dspec.square('time', 'freq', source=tab2_SRC_dspec_square, fill_color=None, fill_alpha=0.0,
                                        line_color=None, line_alpha=0.0, selection_fill_color=None,
                                        selection_fill_alpha=0.0, nonselection_fill_alpha=0.0,
                                        selection_line_alpha=0.0, nonselection_line_alpha=0.0,
                                        size=max(
                                            config_plot['plot_config']['tab_FSview2CASA']['dspec_wdth'] / float(
                                                tab2_ntim),
                                            config_plot['plot_config']['tab_FSview2CASA']['dspec_hght'] / float(
                                                tab2_nfreq)))
    tab2_SRC_dspec_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
    tab2_r_dspec_patch = tab2_p_dspec.patch('xx', 'yy', source=tab2_SRC_dspec_Patch,
                                            fill_color=None, fill_alpha=0.5, line_color="Magenta",
                                            line_alpha=1.0, line_width=1)
    tab2_p_dspec.add_tools(BoxSelectTool())
    tab2_p_dspec.select(BoxSelectTool).select_every_mousemove = False
    # tab2_p_dspec.border_fill_color = "silver"
    tab2_p_dspec.border_fill_alpha = 0.4
    tab2_p_dspec.axis.major_tick_out = 0
    tab2_p_dspec.axis.major_tick_in = 5
    tab2_p_dspec.axis.minor_tick_out = 0
    tab2_p_dspec.axis.minor_tick_in = 3
    tab2_p_dspec.axis.major_tick_line_color = "white"
    tab2_p_dspec.axis.minor_tick_line_color = "white"

    tab2_dspec_selected = None

    timestart = xx[0]


    def tab2_dspec_selection_change(attrname, old, new):
        global tab2_dspec_selected
        tab2_dspec_selected = tab2_SRC_dspec_square.selected['1d']['indices']

        if tab2_dspec_selected:
            global dspecDF_select, tab2_tCLN_Param_dict
            dspecDF_select = dspecDF0.copy()
            dspecDF_select = dspecDF_select.iloc[tab2_dspec_selected, :]
            x0, x1 = dspecDF_select['time'].min(), dspecDF_select['time'].max()
            y0, y1 = dspecDF_select['freq'].min(), dspecDF_select['freq'].max()
            tab2_r_dspec_patch.data_source.data = ColumnDataSource(
                pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
            time0, time1 = dspecDF_select['time'].min() + timestart, dspecDF_select['time'].max() + timestart
            t0_char = Time(time0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
            date0_char = t0_char.split(' ')[0].replace('-', '/')
            time0_char = t0_char.split(' ')[1]
            t1_char = Time(time1 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
            date1_char = t1_char.split(' ')[0].replace('-', '/')
            time1_char = t1_char.split(' ')[1]
            tab2_tCLN_Param_dict['timerange'] = "'{}/{}~{}/{}'".format(date0_char, time0_char, date1_char, time1_char)
            freq0, freq1 = dspecDF_select['freq'].min(), dspecDF_select['freq'].max()
            freqrange = "'{:.3f}~{:.3f} GHz'".format(freq0, freq1)
            tab2_tCLN_Param_dict['freqrange'] = freqrange
            tab2_Div_tCLN_text = ' '.join(
                "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
            tab2_Div_tCLN.text = tab2_Div_tCLN_text
        else:
            tab2_r_dspec_patch.data_source.data = ColumnDataSource(
                pd.DataFrame({'xx': [], 'yy': []})).data


    tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)

    tab2_Select_pol = Select(title="Polarization:", value='I', options=['RR', 'LL', 'I', 'V'],
                             width=config_plot['plot_config']['tab_FSview2CASA']['widgetbox_wdth1'])
    tab2_Select_bl = Select(title="Baseline:", value=tab2_bl[0], options=tab2_bl,
                            width=config_plot['plot_config']['tab_FSview2CASA']['widgetbox_wdth1'])
    tab2_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                                    width=config_plot['plot_config']['tab_FSview2CASA']['widgetbox_wdth1'])

    tab2_p_dspec_xPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                               plot_width=config_plot['plot_config']['tab_FSview2CASA']['dspec_xPro_wdth'],
                               plot_height=config_plot['plot_config']['tab_FSview2CASA']['dspec_xPro_hght'],
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
    tab2_l_dspec_xPro_hover = LabelSet(x='x', y='y', text='tooltips', level='glyph',
                                       source=tab2_SRC_dspec_xPro_hover,
                                       render_mode='canvas')
    tab2_l_dspec_xPro_hover.text_font_size = '5pt'
    tab2_p_dspec_xPro.add_layout(tab2_l_dspec_xPro_hover)
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

    tab2_p_dspec_yPro = figure(tools='', webgl=config_plot['plot_config']['WebGL'],
                               plot_width=config_plot['plot_config']['tab_FSview2CASA']['dspec_yPro_wdth'],
                               plot_height=config_plot['plot_config']['tab_FSview2CASA']['dspec_yPro_hght'],
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
    l_dspec_yPro_hover.text_font_size = '5pt'
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
            tab2_SRC_dspec_image.data = {'data': [np.log(spec_plt)], 'xx': [tab2_dtim], 'yy': [tab2_freq]}
        else:
            tab2_SRC_dspec_image.data = {'data': [spec_plt], 'xx': [tab2_dtim], 'yy': [tab2_freq]}
        tab2_SRC_dspec_square.data['dspec'] = spec_plt.flatten()
        tab2_p_dspec_xPro.y_range.start = spec_plt_min
        tab2_p_dspec_xPro.y_range.end = spec_plt_max
        tab2_p_dspec_yPro.x_range.start = spec_plt_min
        tab2_p_dspec_yPro.x_range.end = spec_plt_max


    tab2_ctrls = [tab2_Select_bl, tab2_Select_pol, tab2_Select_colorspace]
    for ctrl in tab2_ctrls:
        ctrl.on_change('value', tab2_update_dspec_image)

    # # Add a hover tool
    tooltips = None

    hover_JScode = """
        var nx = %d;
        var ny = %d;
        var data = {'x': [], 'y': []};
        var cdata = rs.get('data');
        var indices = cb_data.index['1d'].indices;
        var idx_offset = indices[0] - (indices[0] %% nx);
        for (i=0; i < nx; i++) {
            data['x'].push(cdata.time[i+idx_offset]);
            data['y'].push(cdata.dspec[i+idx_offset]);
        }
        rdx.set('data', data);
        idx_offset = indices[0] %% nx;
        data = {'x': [], 'y': []};
        for (i=0; i < ny; i++) {
            data['x'].push(cdata.dspec[i*nx+idx_offset]);
            data['y'].push(cdata.freq[i*nx+idx_offset]);
        }
        rdy.set('data', data);
        var time = cdata.timestr[indices[0]]+' '
        var freq = cdata.freq[indices[0]].toFixed(3)+'[GHz] '
        var dspec = cdata.dspec[indices[0]].toFixed(3)+ '[sfu]'
        var tooltips = freq + time + dspec
        data = {'x': [], 'y': [], 'tooltips': []};
        data['x'].push(cdata.time[indices[0]]);
        data['y'].push(cdata.dspec[indices[0]]);
        data['tooltips'].push(tooltips);
        rdx_hover.set('data', data);
        tooltips = time + freq + dspec
        data = {'x': [], 'y': [], 'tooltips': []};
        data['x'].push(cdata.dspec[indices[0]]);
        data['y'].push(cdata.freq[indices[0]]);
        data['tooltips'].push(tooltips);
        rdy_hover.set('data', data);
        """ % (tab2_ntim, tab2_nfreq)

    tab2_p_dspec_hover_callback = CustomJS(
        args={'rs': tab2_r_square.data_source, 'rdx': r_dspec_xPro.data_source, 'rdy': r_dspec_yPro.data_source,
              'rdx_hover': r_dspec_xPro_hover.data_source,
              'rdy_hover': r_dspec_yPro_hover.data_source}, code=hover_JScode)
    tab2_p_dspec_hover = HoverTool(tooltips=tooltips, callback=tab2_p_dspec_hover_callback,
                                   renderers=[tab2_r_square])
    tab2_p_dspec.add_tools(tab2_p_dspec_hover)

    tab2_input_tCLN = TextInput(value="Input the param here", title="Clean task parameters:",
                                width=config_plot['plot_config']['tab_FSview2CASA']['input_tCLN_wdth'])
    tab2_Div_tCLN = Div(text='', width=config_plot['plot_config']['tab_FSview2CASA']['tab2_Div_tCLN_wdth'])
    tab2_Div_tCLN2 = Div(text='', width=config_plot['plot_config']['tab_FSview2CASA']['input_tCLN_wdth'])


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
        global tab2_tCLN_Param_dict
        tab2_tCLN_Param_dict = OrderedDict()
        tab2_tCLN_Param_dict['workdir'] = "'./'"
        tab2_tCLN_Param_dict['vis'] = "''"
        tab2_tCLN_Param_dict['imageprefix'] = "'slfcal/{}'".format(struct_id)
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
        tab2_tCLN_Param_dict['mode'] = "'channel'"
        tab2_tCLN_Param_dict['imagermode'] = "'csclean'"
        tab2_tCLN_Param_dict['weighting'] = "'natural'"
        tab2_tCLN_Param_dict['gain'] = '0.1'
        tab2_tCLN_Param_dict['psfmode'] = "'clark'"
        tab2_tCLN_Param_dict['imsize'] = ""'[128, 128]'""
        tab2_tCLN_Param_dict['cell'] = "['5.0arcsec', '5.0arcsec']"
        tab2_tCLN_Param_dict['phasecenter'] = "'J2000 14h26m22.7351 -14d29m29.801'"
        tab2_tCLN_Param_dict['mask'] = "''"
        tab2_tCLN_Param_dict['stokes'] = "'RRLL'"
        tab2_tCLN_Param_dict['uvtaper'] = 'True'
        tab2_tCLN_Param_dict['outertaper'] = "'15.0arcsec'"
        tab2_tCLN_Param_dict['uvrange'] = "''"
        tab2_tCLN_Param_dict['niter'] = "200"
        tab2_tCLN_Param_dict['usescratch'] = "False"
        tab2_tCLN_Param_dict['interactive'] = "False"
        tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
        tab2_Div_tCLN.text = tab2_Div_tCLN_text
        tab2_Div_tCLN2.text = '<p><b>Default parameter Restored.</b></p>'


    tab2_BUT_tCLN_param_default()

    tab2_BUT_tCLN_param_ADD = Button(label='Add to Param',
                                     width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                                     button_type='primary')
    tab2_BUT_tCLN_param_ADD.on_click(tab2_BUT_tCLN_param_add)
    tab2_BUT_tCLN_param_DEL = Button(label='Delete Param',
                                     width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                                     button_type='warning')
    tab2_BUT_tCLN_param_DEL.on_click(tab2_BUT_tCLN_param_delete)
    tab2_BUT_tCLN_param_Default = Button(label='Default Param',
                                         width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'])
    tab2_BUT_tCLN_param_Default.on_click(tab2_BUT_tCLN_param_default)
    tab2_SPCR_LFT_BUT_tCLN_param_DEL = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_BUT_tCLN_param_RELOAD = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_BUT_tCLN_param_SAVE = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_BUT_tCLN_param_DEFAULT = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_BUT_tCLN_param_FSVIEW = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_Div_tCLN2 = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth50'],
                                     height=config_plot['plot_config']['tab_FSview2CASA']['space_hght10'])


    def tab2_BUT_tCLN_param_reload():
        global tab2_tCLN_Param_dict
        with open(database_dir + event_id + struct_id + 'CASA_CLN_args.json', 'r') as fp:
            tab2_tCLN_Param_dict = json.load(fp)
        tab2_Div_tCLN_text = '<p><b>#  ptclean :: Parallelized clean in consecutive time steps</b></p>' + ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
        tab2_Div_tCLN.text = tab2_Div_tCLN_text
        tab2_Div_tCLN2.text = '<p>CASA arguments reload from config file in <b>{}</b>.</p>'.format(
            database_dir + event_id + struct_id)


    def tab2_BUT_tCLN_param_save():
        with open(database_dir + event_id + struct_id + 'CASA_CLN_args.json', 'w') as fp:
            json.dump(tab2_tCLN_Param_dict, fp)
        os.system(
            'cp {}/DataBrowser/FSview/script_clean.py {}'.format(suncasa_dir, database_dir + event_id + struct_id))
        tab2_Div_tCLN2.text = '<p>CASA script and arguments config file saved to <b>{}</b>.</p>'.format(
            database_dir + event_id + struct_id)
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
        os.system('cp {} {}'.format(database_dir + event_id + struct_id + 'CASA_CLN_args.json', workdir))
        os.system('cp {}/DataBrowser/FSview/script_clean.py {}'.format(suncasa_dir, workdir))

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
        with open(database_dir + event_id + struct_id + 'dspecDF-save', 'wb') as fp:
            pickle.dump(dspecDFout, fp)
        tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}</b>. '.format(
            database_dir + event_id + struct_id) + 'Click the <b>clean</b> button to clean. When finished, \
            go back to <b>QLook</b> window, select StrID <b>{}</b> and \
            click <b>FSview</b> button again.</p>'.format(
            database_dir + event_id + struct_id, struct_id[0:-1])


    def tab2_BUT_tCLN_clean():
        cwd = os.getcwd()
        try:
            tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}.</b></p>\
            <p>CASA clean is in processing.</p>'.format(database_dir + event_id + struct_id)
            os.chdir(database_dir + event_id + struct_id)
            suncasapy46 = config_plot['core']['casapy46']
            suncasapy46 = os.path.expandvars(suncasapy46)
            os.system('{} -c script_clean.py'.format(suncasapy46))
            tab2_Div_tCLN2.text = '<p>Clean finished, go back to <b>QLook</b> window, select StrID <b>{}</b> and \
                click <b>FSview</b> button again.</p>'.format(
                struct_id[0:-1])
        except:
            tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}.</b></p>\
            <p>Do image clean with CASA manually.</p>'.format(database_dir + event_id + struct_id) + '<p>When finished,\
             go back to <b>QLook</b> window, select StrID <b>{}</b>\
              and click <b>FSview</b> button again.</p>'.format(struct_id[0:-1])
        os.chdir(cwd)


    tab2_BUT_tCLN_param_RELOAD = Button(label='reload Param',
                                        width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'])
    tab2_BUT_tCLN_param_RELOAD.on_click(tab2_BUT_tCLN_param_reload)
    tab2_BUT_tCLN_param_SAVE = Button(label='save Param',
                                      width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                                      button_type='success')
    tab2_BUT_tCLN_param_SAVE.on_click(tab2_BUT_tCLN_param_save)
    tab2_BUT_tCLN_CLEAN = Button(label='clean',
                                 width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                                 button_type='success')
    tab2_SPCR_ABV_BUT_tCLN = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'],
                                    height=config_plot['plot_config']['tab_FSview2CASA']['space_hght18'])
    tab2_BUT_tCLN_CLEAN.on_click(tab2_BUT_tCLN_clean)
    tab2_SPCR_LFT_BUT_CLEAN = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_BUT_FS_view = Button(label='FS view', width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                              button_type='primary')


    def tab2_panel2_exit():
        tab2_panel_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
        raise SystemExit


    tab2_panel2_BUT_exit = Button(label='Exit FSview2CASA',
                                  width=config_plot['plot_config']['tab_FSview2CASA']['widgetbox_wdth1'],
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
                                                 tab2_panel2_BUT_exit,
                                                 tab2_panel_Div_exit,
                                                 width=config_plot['plot_config']['tab_FSview2CASA'][
                                                     'widgetbox_wdth2'])), tab2_Div_tCLN2))), tab2_Div_tCLN)

    curdoc().add_root(panel2)
    curdoc().title = "FSview2CASA"
