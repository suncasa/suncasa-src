# The plot server must be running
# Go to http://localhost:5006/bokeh to view this plot
import json
import os,sys
import pickle
import time
from collections import OrderedDict
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import sunpy.map
from math import radians, cos, sin
from astropy.io import fits
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Slider, Button, TextInput, RadioButtonGroup, CheckboxGroup,
                          BoxSelectTool, LassoSelectTool, HoverTool, Spacer, LabelSet, Div)
from bokeh.models.mappers import LinearColorMapper
from bokeh.models.widgets import Select
from bokeh.palettes import Spectral11
from bokeh.plotting import figure, curdoc
import glob
from astropy.time import Time
from QLook_util import get_contour_data, twoD_Gaussian
from puffin import PuffinMap

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

'''load config file'''
with open('../config.json', 'r') as fp:
    config_plot = json.load(fp)
database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir)+'/'
with open('{}config_EvtID_curr.json'.format(database_dir), 'r') as fp:
    config_EvtID = json.load(fp)

'''define the colormaps'''
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]
colormap = cm.get_cmap("cubehelix")  # choose any matplotlib colormap here
bokehpalette_SynthesisImg = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]

'''
-------------------------- panel 2,3   --------------------------
'''

event_id = config_EvtID['datadir']['event_id']
try:
    with open(database_dir + event_id + 'CurrFS.json', 'r') as fp:
        FS_config = json.load(fp)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']
FS_specfile = FS_config['datadir']['FS_specfile']
# try:
tab2_specdata = np.load(FS_specfile)
tab2_spec = tab2_specdata['spec']
tab2_npol = tab2_specdata['npol']
tab2_nbl = tab2_specdata['nbl']
tab2_ntim = tab2_specdata['ntim']
tab2_nfreq = tab2_specdata['nfreq']
tab2_tim = tab2_specdata['tim']
tab2_freq = tab2_specdata['freq'] / 1e9
tab2_bl = tab2_specdata['bl'].item().split(';')
bl_index = 0
tab2_pol = 'I'
sz_spec = tab2_spec.shape
tab2_spec_plt_R = tab2_spec[0, bl_index, :, :]
tab2_spec_plt_L = tab2_spec[1, bl_index, :, :]
tab2_spec_plt_I = (tab2_spec[0, bl_index, :, :] + tab2_spec[1, bl_index, :, :]) / 2.
tab2_spec_plt_V = (tab2_spec[0, bl_index, :, :] - tab2_spec[1, bl_index, :, :]) / 2.
spec_plt_max_IRL = (int(max(tab2_spec_plt_R.max(), tab2_spec_plt_L.max(), tab2_spec_plt_I.max())) / 100 + 1) * 100
spec_plt_min_IRL = (int(min(tab2_spec_plt_R.min(), tab2_spec_plt_L.min(), tab2_spec_plt_I.min())) / 100) * 100
spec_plt_max_V = (max(abs(int(tab2_spec_plt_V.max())), abs(int(tab2_spec_plt_V.min()))) / 100 + 1) * 100
spec_plt_min_V = -spec_plt_max_V
if tab2_pol == 'RR':
    tab2_spec_plt = tab2_spec[0, bl_index, :, :]
    spec_plt_max = spec_plt_max_IRL
    spec_plt_min = spec_plt_min_IRL
elif tab2_pol == 'LL':
    tab2_spec_plt = tab2_spec[1, bl_index, :, :]
    spec_plt_max = spec_plt_max_IRL
    spec_plt_min = spec_plt_min_IRL
elif tab2_pol == 'I':
    tab2_spec_plt = (tab2_spec[0, bl_index, :, :] + tab2_spec[1, bl_index, :, :]) / 2.
    spec_plt_max = spec_plt_max_IRL
    spec_plt_min = spec_plt_min_IRL
elif tab2_pol == 'V':
    tab2_spec_plt = (tab2_spec[0, bl_index, :, :] - tab2_spec[1, bl_index, :, :]) / 2.
    spec_plt_max = spec_plt_max_V
    spec_plt_min = spec_plt_min_V

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
        dspecDF = dspecDF0.copy()
    itemset1 = set(['x_pos', 'y_pos'])
    itemset2 = set(dspecDF.columns.tolist())
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
        tab2_panel2_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>
                                first before closing the tab</p></b>""",
                                   width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_panel3_Div_exit = Div(text="""<p><b>Warning</b>: Click the <b>Exit FSview</b>
                                first before closing the tab</p></b>""",
                                   width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        rmax, rmin = tab2_spec_plt.max(), tab2_spec_plt.min()
        colors_dspec = [colors.rgb2hex(m) for m in colormap_jet((tab2_spec_plt.flatten() - rmin) / (rmax - rmin))]

        TOOLS = "crosshair,pan,wheel_zoom,tap,box_zoom,reset,save"

        tab2_SRC_dspec_square = ColumnDataSource(dspecDF)

        '''create the dynamic spectrum plot'''
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
                                            line_color=None, line_alpha=0.0, selection_fill_alpha=0.1,
                                            selection_fill_color='black',
                                            nonselection_fill_alpha=0.0,
                                            selection_line_alpha=0.2, selection_line_color='white',
                                            nonselection_line_alpha=0.0,
                                            size=max(
                                                config_plot['plot_config']['tab_FSview_base']['dspec_wdth'] / tab2_ntim,
                                                config_plot['plot_config']['tab_FSview_base'][
                                                    'dspec_hght'] / tab2_nfreq))

        tab2_p_dspec.add_tools(BoxSelectTool())
        tab2_p_dspec.add_tools(LassoSelectTool())
        tab2_p_dspec.select(BoxSelectTool).select_every_mousemove = False
        tab2_p_dspec.select(LassoSelectTool).select_every_mousemove = False
        tab2_p_dspec.border_fill_color = "silver"
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

        map = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

        tab2_p_dspec_xPro = figure(tools='',
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
        tab2_p_dspec_xPro.border_fill_color = "silver"
        tab2_p_dspec_xPro.border_fill_alpha = 0.4
        tab2_p_dspec_xPro.axis.major_tick_out = 0
        tab2_p_dspec_xPro.axis.major_tick_in = 5
        tab2_p_dspec_xPro.axis.minor_tick_out = 0
        tab2_p_dspec_xPro.axis.minor_tick_in = 3
        tab2_p_dspec_xPro.axis.major_tick_line_color = "black"
        tab2_p_dspec_xPro.axis.minor_tick_line_color = "black"

        tab2_p_dspec_yPro = figure(tools='',
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
        tab2_p_dspec_yPro.border_fill_color = "silver"
        tab2_p_dspec_yPro.border_fill_alpha = 0.4
        tab2_p_dspec_yPro.min_border_bottom = 0
        tab2_p_dspec_yPro.min_border_left = 0
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
            if select_pol == 'RR':
                spec_plt = spec_plt_R
                spec_plt_max = spec_plt_max_IRL
                spec_plt_min = spec_plt_min_IRL
            elif select_pol == 'LL':
                spec_plt = spec_plt_L
                spec_plt_max = spec_plt_max_IRL
                spec_plt_min = spec_plt_min_IRL
            elif select_pol == 'I':
                spec_plt = spec_plt_I
                spec_plt_max = spec_plt_max_IRL
                spec_plt_min = spec_plt_min_IRL
            elif select_pol == 'V':
                spec_plt = spec_plt_V
                spec_plt_max = spec_plt_max_V
                spec_plt_min = spec_plt_min_V
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

        url = ""
        tab2_SRC_p_dspec_thumb = ColumnDataSource(
            dict(url=[url], x1=[0], y1=[0], w1=[config_plot['plot_config']['tab_FSview_base']['dspec_thumb_wdth']],
                 h1=[config_plot['plot_config']['tab_FSview_base']['dspec_thumb_hght']], ))
        tab2_p_dspec_thumb = figure(tools='',
                                    plot_width=config_plot['plot_config']['tab_FSview_base']['dspec_thumb_wdth'],
                                    plot_height=config_plot['plot_config']['tab_FSview_base']['dspec_thumb_hght'],
                                    x_range=(0, config_plot['plot_config']['tab_FSview_base']['dspec_thumb_wdth']),
                                    y_range=(0, config_plot['plot_config']['tab_FSview_base']['dspec_thumb_hght']),
                                    title="EVLA thumbnail",
                                    toolbar_location=None)
        r_dspec_thumb = tab2_p_dspec_thumb.image_url(url="url", x="x1", y="y1", w="w1", h="h1",
                                                     source=tab2_SRC_p_dspec_thumb, anchor='bottom_left')
        tab2_p_dspec_thumb.xaxis.visible = False
        tab2_p_dspec_thumb.yaxis.visible = False
        tab2_p_dspec_thumb.title.text_font_size = '6pt'
        tab2_p_dspec_thumb.border_fill_color = "silver"
        tab2_p_dspec_thumb.border_fill_alpha = 0.4

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
            rdt.data['url'] = []
            rdt.data['url'].push(cdata.thumbnail[indices[0]])
            rdt.trigger('change');
            """ % (tab2_ntim, tab2_nfreq)

        tab2_p_dspec_hover_callback = CustomJS(
            args={'rs': tab2_r_square.data_source, 'rdx': r_dspec_xPro.data_source, 'rdy': r_dspec_yPro.data_source,
                  'rdt': r_dspec_thumb.data_source, 'rdx_hover': r_dspec_xPro_hover.data_source,
                  'rdy_hover': r_dspec_yPro_hover.data_source}, code=hover_JScode)
        tab2_p_dspec_hover = HoverTool(tooltips=tooltips, callback=tab2_p_dspec_hover_callback,
                                       renderers=[tab2_r_square])
        tab2_p_dspec.add_tools(tab2_p_dspec_hover)

        # initial the VLA map contour source
        tab2_SRC_vlamap_contour = ColumnDataSource(
            data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})
        tab2_SRC_vlamap_peak = ColumnDataSource(data={'dspec': [], 'x_pos': [], 'y_pos': [], 'amp_gaus': []})


        # initial the source of maxfit centroid
        def tab2_SRC_maxfit_centroid_init(dspecDF):
            start_timestamp = time.time()
            global SRC_maxfit_centroid
            SRC_maxfit_centroid = {}
            for ll in np.unique(dspecDF['time']):
                df_tmp = pd.DataFrame(
                    {'freq': [], 'x_pos': [], 'y_pos': [], 'x_width': [], 'y_width': [], 'amp_gaus': [],
                     'theta': [],
                     'amp_offset': []})
                SRC_maxfit_centroid[np.where(abs(tab2_dtim - ll) < 0.02)[0].tolist()[0]] = ColumnDataSource(df_tmp)
            print("---tab2_SRC_maxfit_centroid_init -- %s seconds ---" % (time.time() - start_timestamp))


        tab2_SRC_maxfit_centroid_init(dspecDF)

        # import the vla image
        if dspecDF.loc[76, :]['fits_exist']:
            hdulist = fits.open(fits_GLOB_dir + dspecDF.loc[76, :]['fits_global'])
            hdu = hdulist[0]
            vla_global_pfmap = PuffinMap(hdu.data[0, 0, :, :], hdu.header,
                                         plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                         plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'],
                                         palette=bokehpalette_SynthesisImg)
            hdulist = fits.open(fits_LOCL_dir + dspecDF.loc[76, :]['fits_local'])
            hdu = hdulist[0]
            vla_local_pfmap = PuffinMap(hdu.data[0, 0, :, :], hdu.header)
            # plot the contour of vla image
            popt = [dspecDF.loc[76, :]['amp_gaus'], dspecDF.loc[76, :]['x_pos'], dspecDF.loc[76, :]['y_pos'],
                    dspecDF.loc[76, :]['x_width'], dspecDF.loc[76, :]['y_width'], dspecDF.loc[76, :]['theta'],
                    dspecDF.loc[76, :]['amp_offset']]
            mapx, mapy = vla_local_pfmap.meshgrid()
            mapx, mapy = mapx.value, mapy.value
            vlamap_fitted = twoD_Gaussian((mapx, mapy), *popt).reshape(vla_local_pfmap.smap.data.shape)
            tab2_SRC_vlamap_contour = get_contour_data(mapx, mapy, vlamap_fitted)
            tab2_SRC_vlamap_peak = ColumnDataSource(
                data={'dspec': [dspecDF.loc[76, :]['dspec']], 'x_pos': [dspecDF.loc[76, :]['x_pos']],
                      'y_pos': [dspecDF.loc[76, :]['y_pos']], 'amp_gaus': [dspecDF.loc[76, :]['amp_gaus']]})

        # import the aia image
        # from sunpy.net.helioviewer import HelioviewerClient
        #
        # hv = HelioviewerClient()
        # filepath = hv.download_jp2(jdutil.jd_to_datetime(xx[0] / 3600. / 24.), observatory='SDO', instrument='AIA',
        #                            detector='AIA', measurement='171',
        #                            directory=database_dir + event_id + struct_id + config_EvtID['datadir']['J2000'],
        #                            overwrite=True)
        filepath = database_dir + event_id + struct_id + config_EvtID['datadir'][
            'J2000'] + '2014_11_01__16_45_59_34__SDO_AIA_AIA_171.jp2'
        colormap = cm.get_cmap("sdoaia171")  # choose any matplotlib colormap here
        bokehpalette_sdoaia171 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        aiamap = sunpy.map.Map(filepath)
        lengthx = vla_local_pfmap.dw[0] * u.arcsec
        lengthy = vla_local_pfmap.dh[0] * u.arcsec
        x0 = vla_local_pfmap.smap.center.x
        y0 = vla_local_pfmap.smap.center.y
        aiamap_submap = aiamap.submap(u.Quantity([x0 - lengthx, x0 + lengthx]),
                                      u.Quantity([y0 - lengthy, y0 + lengthy]))
        MapRES = 256
        dimensions = u.Quantity([MapRES, MapRES], u.pixel)
        aia_resampled_map = aiamap.resample(dimensions)

        # plot the global AIA image

        aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'])

        tab2_p_aia, tab2_r_aia = aia_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                             palette=bokehpalette_sdoaia171)
        tab2_p_aia.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_aia.circle(x='x_pos', y='y_pos',  # size=10.*dspecDF.loc[76,:]['amp_gaus']/50.,
                          radius=3, radius_units='data', source=tab2_SRC_vlamap_peak, fill_alpha=0.8,
                          fill_color='#7c7e71',
                          line_color='#7c7e71')
        tab2_p_aia.title.text_font_size = '6pt'
        tab2_p_aia.border_fill_color = "silver"
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
                                     plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['aia_submap_wdth'])

        tab3_p_aia_submap, tab3_r_aia_submap = aia_submap_pfmap.PlotMap(DrawLimb=True, DrawGrid=True,
                                                                        grid_spacing=20 * u.deg,
                                                                        title='EM sources centroid map',
                                                                        palette=bokehpalette_sdoaia171)

        # tab3_p_aia_submap.add_tools(ResizeTool())
        tab3_p_aia_submap.border_fill_color = "silver"
        tab3_p_aia_submap.border_fill_alpha = 0.4
        tab3_p_aia_submap.axis.major_tick_out = 0
        tab3_p_aia_submap.axis.major_tick_in = 5
        tab3_p_aia_submap.axis.minor_tick_out = 0
        tab3_p_aia_submap.axis.minor_tick_in = 3
        tab3_p_aia_submap.axis.major_tick_line_color = "white"
        tab3_p_aia_submap.axis.minor_tick_line_color = "white"
        color_mapper = LinearColorMapper(Spectral11)

        tab3_r_aia_submap_cross = tab3_p_aia_submap.cross(x='x_pos', y='y_pos', size=15,
                                                          color={'field': 'freq', 'transform': color_mapper},
                                                          line_width=3,
                                                          source=SRC_maxfit_centroid[tab2_dtim[0]], line_alpha=0.8)
        tab3_r_aia_submap_line = tab3_p_aia_submap.line(x='x_pos', y='y_pos', line_width=3, line_color='black',
                                                        line_alpha=0.5,
                                                        source=SRC_maxfit_centroid[tab2_dtim[0]])
        tab3_r_aia_submap_line.visible = False
        tab3_SRC_aia_submap_rect = ColumnDataSource({'x': [], 'y': [], 'width': [], 'height': []})
        tab3_r_aia_submap_rect = tab3_p_aia_submap.rect(x='x', y='y', width='width', height='height', fill_alpha=0.1,
                                                        line_color='black', fill_color='black',
                                                        source=tab3_SRC_aia_submap_rect)

        # plot the global HMI image
        # filepath = hv.download_jp2(jdutil.jd_to_datetime(xx[0] / 3600. / 24.), observatory='SDO', instrument='HMI',
        #                            detector='HMI', measurement='magnetogram',
        #                            directory=database_dir + event_id + struct_id + config_EvtID['datadir']['J2000'],
        #                            overwrite=True)
        filepath = database_dir + event_id + struct_id + config_EvtID['datadir'][
            'J2000'] + '2014_11_01__16_46_58_605__SDO_HMI_HMI_magnetogram.jp2'
        colormap = cm.get_cmap("gray")  # choose any matplotlib colormap here
        bokehpalette_sdohmimag = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        hmimap = sunpy.map.Map(filepath)
        # plot the global HMI image
        hmi_resampled_map = hmimap.resample(dimensions)
        hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'])

        tab2_p_hmi, tab2_r_hmi = hmi_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                             x_range=tab2_p_aia.x_range,
                                                             y_range=tab2_p_aia.y_range,
                                                             palette=bokehpalette_sdohmimag)
        tab2_p_hmi.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                              line_width=2)
        tab2_p_hmi.circle(x='x_pos', y='y_pos', radius=3, radius_units='data', source=tab2_SRC_vlamap_peak,
                          fill_alpha=0.8,
                          fill_color='#7c7e71', line_color='#7c7e71')
        tab2_p_hmi.yaxis.visible = False
        tab2_p_hmi.border_fill_color = "silver"
        tab2_p_hmi.border_fill_alpha = 0.4
        tab2_p_hmi.axis.major_tick_out = 0
        tab2_p_hmi.axis.major_tick_in = 5
        tab2_p_hmi.axis.minor_tick_out = 0
        tab2_p_hmi.axis.minor_tick_in = 3
        tab2_p_hmi.axis.major_tick_line_color = "white"
        tab2_p_hmi.axis.minor_tick_line_color = "white"
        tab2_p_hmi.grid.grid_line_color = None
        tab2_p_hmi.background_fill_color = "black"

        # plot the global vla image
        if dspecDF.loc[76, :]['fits_exist']:
            tab2_p_vla, tab2_r_vla = vla_global_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                              palette=bokehpalette_SynthesisImg,
                                                              x_range=tab2_p_aia.x_range,
                                                              y_range=tab2_p_aia.y_range)
            tab2_p_vla.title.text_font_size = '6pt'
            tab2_p_vla.yaxis.visible = False
            tab2_p_vla.border_fill_color = "silver"
            tab2_p_vla.border_fill_alpha = 0.4
            tab2_p_vla.axis.major_tick_out = 0
            tab2_p_vla.axis.major_tick_in = 5
            tab2_p_vla.axis.minor_tick_out = 0
            tab2_p_vla.axis.minor_tick_in = 3
            tab2_p_vla.axis.major_tick_line_color = "white"
            tab2_p_vla.axis.minor_tick_line_color = "white"
            tab2_p_vla.grid.grid_line_color = None
            tab2_p_vla.background_fill_color = "black"

            tab2_r_vla_multi_line = tab2_p_vla.multi_line(xs='xs', ys='ys', line_color='line_color',
                                                          source=tab2_SRC_vlamap_contour, alpha=0.7, line_width=2)
            tab2_r_vla_circle = tab2_p_vla.circle(x='x_pos', y='y_pos',
                                                  # size=10.*dspecDF.loc[76,:]['amp_gaus']/50.,
                                                  radius=3, radius_units='data', source=tab2_SRC_vlamap_peak,
                                                  fill_alpha=0.8, fill_color='#7c7e71',
                                                  line_color='#7c7e71')


        def tab2_SRC_maxfit_centroid_update(dspecDF):
            start_timestamp = time.time()
            global SRC_maxfit_centroid, timebin
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                SRC_maxfit_centroid = {}
                for ll in np.unique(dspecDF['time']):
                    dftmp = dspecDF[dspecDF.time == ll]
                    dftmp = dftmp.dropna(how='any')
                    df_tmp = pd.concat(
                        [dftmp.loc[:, 'freq'], dftmp.loc[:, 'x_pos'], dftmp.loc[:, 'y_pos'], dftmp.loc[:, 'x_width'],
                         dftmp.loc[:, 'y_width'], dftmp.loc[:, 'amp_gaus'], dftmp.loc[:, 'theta'] - np.pi / 2,
                         dftmp.loc[:, 'amp_offset']], axis=1)
                    SRC_maxfit_centroid[np.where(abs(tab2_dtim - ll) < 0.02)[0].tolist()[0]] = ColumnDataSource(df_tmp)
            else:
                time_dspec = np.unique(dspecDF['time'])
                ntime_dspec = len(time_dspec)
                if timebin != 1:
                    tidx = np.arange(0, ntime_dspec + 1, timebin)
                    time_seq = time_dspec[0:0 + timebin]
                    dftmp = dspecDF[dspecDF['time'].isin(time_seq)]
                    dftmp = dftmp.dropna(how='any')
                    dftmp_concat = pd.DataFrame(dict(dftmp.mean()), index=[0, ])
                    for ll in tidx[1:]:
                        time_seq = time_dspec[ll:ll + timebin]
                        dftmp = dspecDF[dspecDF['time'].isin(time_seq)]
                        dftmp = dftmp.dropna(how='any')
                        dftmp_concat = dftmp_concat.append(pd.DataFrame(dict(dftmp.mean()), index=[0, ]),
                                                           ignore_index=True)
                    SRC_maxfit_centroid = ColumnDataSource(
                        dftmp_concat[
                            ['freq', 'x_pos', 'y_pos', 'x_width', 'y_width', 'amp_gaus', 'theta', 'amp_offset']].dropna(
                            how='any'))
                else:
                    dftmp = dspecDF.copy()
                    dftmp = dftmp.dropna(how='any')
                    df_tmp = pd.concat(
                        [dftmp.loc[:, 'freq'], dftmp.loc[:, 'x_pos'], dftmp.loc[:, 'y_pos'], dftmp.loc[:, 'x_width'],
                         dftmp.loc[:, 'y_width'], dftmp.loc[:, 'amp_gaus'], dftmp.loc[:, 'theta'] - np.pi / 2,
                         dftmp.loc[:, 'amp_offset']], axis=1)
                    SRC_maxfit_centroid = ColumnDataSource(df_tmp)
            print("--- tab2_SRC_maxfit_centroid_update -- %s seconds ---" % (time.time() - start_timestamp))


        tab2_LinkImg_HGHT = config_plot['plot_config']['tab_FSview_base']['vla_hght']
        tab2_LinkImg_WDTH = config_plot['plot_config']['tab_FSview_base']['vla_wdth']


        def tab2_LinkImg_replot_update():
            idx_selected = dspecDF.index[len(dspecDF) / 2]
            if dspecDF.loc[idx_selected, :]['fits_exist']:
                hdulist = fits.open(fits_GLOB_dir + dspecDF.loc[idx_selected, :]['fits_global'])
                hdu = hdulist[0]
                pfmap = PuffinMap(hdu.data[0, 0, :, :], hdu.header, plot_height=tab2_LinkImg_HGHT,
                                  plot_width=tab2_LinkImg_WDTH)
                SRC_Img = pfmap.ImageSource()
                tab2_r_vla.data_source.data['data'] = SRC_Img.data['data']
                popt = [dspecDF.loc[idx_selected, :]['amp_gaus'], dspecDF.loc[idx_selected, :]['x_pos'],
                        dspecDF.loc[idx_selected, :]['y_pos'], dspecDF.loc[idx_selected, :]['x_width'],
                        dspecDF.loc[idx_selected, :]['y_width'], dspecDF.loc[idx_selected, :]['theta'],
                        dspecDF.loc[idx_selected, :]['amp_offset']]
                hdulist = fits.open(fits_LOCL_dir + dspecDF.loc[idx_selected, :]['fits_local'])
                hdu = hdulist[0]
                pfmap_local = PuffinMap(hdu.data[0, 0, :, :], hdu.header)
                mapx, mapy = pfmap_local.meshgrid()
                mapx, mapy = mapx.value, mapy.value
                vlamap_fitted = twoD_Gaussian((mapx, mapy), *popt).reshape(pfmap_local.smap.data.shape)
                SRC_contour = get_contour_data(mapx, mapy, vlamap_fitted)
                tab2_r_vla_multi_line.data_source.data = SRC_contour.data
                SRC_peak = ColumnDataSource(data={'dspec': [dspecDF.loc[idx_selected, :]['dspec']],
                                                  'x_pos': [dspecDF.loc[idx_selected, :]['x_pos']],
                                                  'y_pos': [dspecDF.loc[idx_selected, :]['y_pos']],
                                                  'amp_gaus': [dspecDF.loc[idx_selected, :]['amp_gaus']]})
                tab2_r_vla_circle.data_source.data = SRC_peak.data


        tab2_BUT_LinkImg_replot = Button(label='replot',
                                         width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])
        tab2_BUT_LinkImg_replot.on_click(tab2_LinkImg_replot_update)


        def tab2_panel_exit():
            tab2_panel2_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
            tab2_panel3_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
            raise SystemExit


        tab2_panel2_BUT_exit = Button(label='Exit FSview',
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                      button_type='danger')
        tab2_panel2_BUT_exit.on_click(tab2_panel_exit)

        tab2_panel3_BUT_exit = Button(label='Exit FSview',
                                      width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                      button_type='danger')
        tab2_panel3_BUT_exit.on_click(tab2_panel_exit)

        tab3_p_dspec_small = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                    plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                    plot_height=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_hght'],
                                    x_range=(tab2_dtim[0], tab2_dtim[-1]),
                                    y_range=(tab2_freq[0], tab2_freq[-1]), toolbar_location='above')
        tab3_p_dspecvx_small = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                      plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                      plot_height=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_hght'],
                                      x_range=tab3_p_dspec_small.x_range,
                                      y_range=tab3_p_dspec_small.y_range, toolbar_location='above')
        tab3_p_dspecvy_small = figure(tools='pan,wheel_zoom,box_zoom,save,reset',
                                      plot_width=config_plot['plot_config']['tab_FSview_FitANLYS']['dspec_small_wdth'],
                                      plot_height=config_plot['plot_config']['tab_FSview_FitANLYS'][
                                                      'dspec_small_hght'] + 40,
                                      x_range=tab3_p_dspec_small.x_range,
                                      y_range=tab3_p_dspec_small.y_range, toolbar_location='above')
        tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        tab3_p_dspec_small.xaxis.visible = False
        tab3_p_dspecvx_small.xaxis.visible = False
        tab3_p_dspec_small.title.text = "Vector Dynamic spectrum (Intensity)"
        tab3_p_dspecvx_small.title.text = "Vector Dynamic spectrum (Vx)"
        tab3_p_dspecvy_small.title.text = "Vector Dynamic spectrum (Vy)"
        tab3_p_dspecvy_small.xaxis.axis_label = 'Seconds since ' + tim0_char
        tab3_p_dspec_small.yaxis.axis_label = 'Frequency [GHz]'
        tab3_p_dspecvx_small.yaxis.axis_label = 'Frequency [GHz]'
        tab3_p_dspecvy_small.yaxis.axis_label = 'Frequency [GHz]'
        tab3_p_dspec_small.border_fill_color = "silver"
        tab3_p_dspec_small.border_fill_alpha = 0.4
        tab3_p_dspec_small.axis.major_tick_out = 0
        tab3_p_dspec_small.axis.major_tick_in = 5
        tab3_p_dspec_small.axis.minor_tick_out = 0
        tab3_p_dspec_small.axis.minor_tick_in = 3
        tab3_p_dspec_small.axis.major_tick_line_color = "white"
        tab3_p_dspec_small.axis.minor_tick_line_color = "white"
        tab3_p_dspecvx_small.border_fill_color = "silver"
        tab3_p_dspecvx_small.border_fill_alpha = 0.4
        tab3_p_dspecvx_small.axis.major_tick_out = 0
        tab3_p_dspecvx_small.axis.major_tick_in = 5
        tab3_p_dspecvx_small.axis.minor_tick_out = 0
        tab3_p_dspecvx_small.axis.minor_tick_in = 3
        tab3_p_dspecvx_small.axis.major_tick_line_color = "white"
        tab3_p_dspecvx_small.axis.minor_tick_line_color = "white"
        tab3_p_dspecvy_small.border_fill_color = "silver"
        tab3_p_dspecvy_small.border_fill_alpha = 0.4
        tab3_p_dspecvy_small.axis.major_tick_out = 0
        tab3_p_dspecvy_small.axis.major_tick_in = 5
        tab3_p_dspecvy_small.axis.minor_tick_out = 0
        tab3_p_dspecvy_small.axis.minor_tick_in = 3
        tab3_p_dspecvy_small.axis.major_tick_line_color = "white"
        tab3_p_dspecvy_small.axis.minor_tick_line_color = "white"

        dspecDFtmp = pd.DataFrame()
        nrows_dspecDF = len(dspecDF0.index)
        dspecDFtmp['x_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
        dspecDFtmp['y_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
        dspecDFtmp['amp_gaus'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
        dspecDFtmp.loc[:, 'x_pos'] = dspecDF0.loc[:, 'x_pos']
        dspecDFtmp.loc[:, 'y_pos'] = dspecDF0.loc[:, 'y_pos']
        dspecDFtmp.loc[:, 'amp_gaus'] = dspecDF0.loc[:, 'amp_gaus']


        def tab3_SRC_dspec_small_init():
            global tab3_SRC_dspec_small, tab3_SRC_dspecvx_small, tab3_SRC_dspecvy_small
            global mean_amp_g, mean_vx, mean_vy, drange_amp_g, drange_vx, drange_vy
            global vmax_amp_g, vmax_vx, vmax_vy, vmin_amp_g, vmin_vx, vmin_vy
            start_timestamp = time.time()
            amp_g = (dspecDF0['amp_gaus'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_amp_g = np.nanmean(amp_g)
            drange_amp_g = 40.
            vmax_amp_g, vmin_amp_g = mean_amp_g + drange_amp_g * np.asarray([1., -1.])
            amp_g[amp_g > vmax_amp_g] = vmax_amp_g
            amp_g[amp_g < vmin_amp_g] = vmin_amp_g
            tab3_SRC_dspec_small = ColumnDataSource(data={'data': [amp_g], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
            vx = (dspecDF0['x_pos'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_vx = np.nanmean(vx)
            drange_vx = 40.
            vmax_vx, vmin_vx = mean_vx + drange_vx * np.asarray([1., -1.])
            vx[vx > vmax_vx] = vmax_vx
            vx[vx < vmin_vx] = vmin_vx
            tab3_SRC_dspecvx_small = ColumnDataSource(data={'data': [vx], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
            vy = (dspecDF0['y_pos'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_vy = np.nanmean(vy)
            drange_vy = 40.
            vmax_vy, vmin_vy = mean_vy + drange_vy * np.asarray([1., -1.])
            vy[vy > vmax_vy] = vmax_vy
            vy[vy < vmin_vy] = vmin_vy
            tab3_SRC_dspecvy_small = ColumnDataSource(data={'data': [vy], 'xx': [tab2_dtim], 'yy': [tab2_freq]})
            tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
            tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
            tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
            tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]
            print("--- tab3_SRC_dspec_small_init -- %s seconds ---" % (time.time() - start_timestamp))


        def tab3_SRC_dspec_small_update(dspecDFtmp):
            global tab3_SRC_dspec_small, tab3_SRC_dspecvx_small, tab3_SRC_dspecvy_small
            global mean_amp_g, mean_vx, mean_vy, drange_amp_g, drange_vx, drange_vy
            global vmax_amp_g, vmax_vx, vmax_vy, vmin_amp_g, vmin_vx, vmin_vy
            start_timestamp = time.time()
            amp_g = (dspecDFtmp['amp_gaus'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_amp_g = np.nanmean(amp_g)
            drange_amp_g = 40.
            vmax_amp_g, vmin_amp_g = mean_amp_g + drange_amp_g * np.asarray([1., -1.])
            amp_g[amp_g > vmax_amp_g] = vmax_amp_g
            amp_g[amp_g < vmin_amp_g] = vmin_amp_g
            tab3_SRC_dspec_small.data['data'] = [amp_g]
            vx = (dspecDFtmp['x_pos'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_vx = np.nanmean(vx)
            drange_vx = 40.
            vmax_vx, vmin_vx = mean_vx + drange_vx * np.asarray([1., -1.])
            vx[vx > vmax_vx] = vmax_vx
            vx[vx < vmin_vx] = vmin_vx
            tab3_SRC_dspecvx_small.data['data'] = [vx]
            vy = (dspecDFtmp['y_pos'].copy()).reshape(tab2_nfreq, tab2_ntim)
            mean_vy = np.nanmean(vy)
            drange_vy = 40.
            vmax_vy, vmin_vy = mean_vy + drange_vy * np.asarray([1., -1.])
            vy[vy > vmax_vy] = vmax_vy
            vy[vy < vmin_vy] = vmin_vy
            tab3_SRC_dspecvy_small.data['data'] = [vy]
            print("--- tab3_SRC_dspec_small_update -- %s seconds ---" % (time.time() - start_timestamp))


        tab3_SRC_dspec_small_init()

        tab3_p_dspec_small.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                 dh=tab2_freq[-1] - tab2_freq[0],
                                 palette=bokehpalette_jet, source=tab3_SRC_dspec_small)
        tab3_p_dspecvx_small.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                   dh=tab2_freq[-1] - tab2_freq[0],
                                   palette=bokehpalette_jet, source=tab3_SRC_dspecvx_small)
        tab3_p_dspecvy_small.image(image="data", x=tab2_dtim[0], y=tab2_freq[0], dw=tab2_dtim[-1] - tab2_dtim[0],
                                   dh=tab2_freq[-1] - tab2_freq[0],
                                   palette=bokehpalette_jet, source=tab3_SRC_dspecvy_small)
        tab3_source_idx_line = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
        tab3_r_dspec_small_line = tab3_p_dspec_small.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                          line_color='white', source=tab3_source_idx_line)
        tab3_r_dspecvx_small_line = tab3_p_dspecvx_small.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                              line_color='white',
                                                              source=tab3_source_idx_line)
        tab3_r_dspecvy_small_line = tab3_p_dspecvy_small.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                              line_color='white',
                                                              source=tab3_source_idx_line)
        tab2_dspec_selected = None


        def dspecDFtmp_init():
            global dspecDFtmp
            dspecDFtmp = pd.DataFrame()
            nrows_dspecDF = len(dspecDF0.index)
            dspecDFtmp['amp_gaus'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
            dspecDFtmp['x_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
            dspecDFtmp['y_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)


        def tab2_dspec_selection_change(attrname, old, new):
            global tab2_dspec_selected
            tab2_dspec_selected = tab2_SRC_dspec_square.selected['1d']['indices']
            if tab2_dspec_selected:
                global dspecDF
                dspecDF = dspecDF0.iloc[tab2_dspec_selected, :]
                dspecDFtmp_init()
                dspecDFtmp.loc[tab2_dspec_selected, 'x_pos'] = dspecDF0.loc[tab2_dspec_selected, 'x_pos']
                dspecDFtmp.loc[tab2_dspec_selected, 'y_pos'] = dspecDF0.loc[tab2_dspec_selected, 'y_pos']
                dspecDFtmp.loc[tab2_dspec_selected, 'amp_gaus'] = dspecDF0.loc[tab2_dspec_selected, 'amp_gaus']
                tab3_SRC_dspec_small_update(dspecDFtmp)
                tab2_SRC_maxfit_centroid_update(dspecDF)
                if tab3_BUT_animate_ONOFF.label == 'Animate OFF & Go':
                    tab3_r_aia_submap_cross.visible = True
                    tab3_r_dspec_small_line.visible = False
                    tab3_r_dspecvx_small_line.visible = False
                    tab3_r_dspecvy_small_line.visible = False
                    tab2_SRC_maxfit_centroid_update(dspecDF)
                    tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid.data


        tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)

        tab3_dspec_small_CTRLs_OPT = dict(mean_values=[mean_amp_g, mean_vx, mean_vy],
                                          drange_values=[drange_amp_g, drange_vx, drange_vy],
                                          vmax_values=[vmax_amp_g, vmax_vx, vmax_vy],
                                          vmin_values=[vmin_amp_g, vmin_vx, vmin_vy],
                                          vmax_values_last=[vmax_amp_g, vmax_vx, vmax_vy],
                                          vmin_values_last=[vmin_amp_g, vmin_vx, vmin_vy],
                                          items_dspec_small=['amp_gaus', 'x_pos', 'y_pos'],
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


        tab3_RBG_dspec_small.on_change('active', tab3_RBG_dspec_small_update)


        def tab3_BUT_dspec_small_reset_update():
            global dspecDFtmp, tab2_nfreq, tab2_ntim, tab3_dspec_small_CTRLs_OPT
            global tab3_SRC_dspec_small, tab3_SRC_dspecvx_small, tab3_SRC_dspecvy_small
            items_dspec_small = tab3_dspec_small_CTRLs_OPT['items_dspec_small']
            mean_values = tab3_dspec_small_CTRLs_OPT['mean_values']
            drange_values = tab3_dspec_small_CTRLs_OPT['drange_values']
            vmax_values = tab3_dspec_small_CTRLs_OPT['vmax_values']
            vmin_values = tab3_dspec_small_CTRLs_OPT['vmin_values']
            source_list = [tab3_SRC_dspec_small, tab3_SRC_dspecvx_small, tab3_SRC_dspecvy_small]
            for ll, item in enumerate(items_dspec_small):
                TmpData = (dspecDFtmp[item].copy()).reshape(tab2_nfreq, tab2_ntim)
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


        tab3_BUT_dspec_small_reset.on_click(tab3_BUT_dspec_small_reset_update)


        def tab3_slider_dspec_small_update(attrname, old, new):
            global dspecDFtmp, tab2_nfreq, tab2_ntim, tab3_dspec_small_CTRLs_OPT
            items_dspec_small = tab3_dspec_small_CTRLs_OPT['items_dspec_small']
            idx_p_dspec_small = tab3_dspec_small_CTRLs_OPT['idx_p_dspec_small']
            dmax = tab3_Slider_dspec_small_dmax.value
            dmin = tab3_Slider_dspec_small_dmin.value
            if not tab3_dspec_small_CTRLs_OPT['radio_button_group_dspec_small_update_flag']:
                tab3_dspec_small_CTRLs_OPT['vmax_values_last'][idx_p_dspec_small] = dmax
                tab3_dspec_small_CTRLs_OPT['vmin_values_last'][idx_p_dspec_small] = dmin
            TmpData = (dspecDFtmp[items_dspec_small[idx_p_dspec_small]].copy()).reshape(tab2_nfreq, tab2_ntim)
            TmpData[TmpData > dmax] = dmax
            TmpData[TmpData < dmin] = dmin
            if idx_p_dspec_small == 0:
                global tab3_SRC_dspec_small
                tab3_SRC_dspec_small.data['data'] = [TmpData]
            elif idx_p_dspec_small == 1:
                global tab3_SRC_dspecvx_small
                tab3_SRC_dspecvx_small.data['data'] = [TmpData]
            elif idx_p_dspec_small == 2:
                global tab3_SRC_dspecvy_small
                tab3_SRC_dspecvy_small.data['data'] = [TmpData]
            vmax_vx, vmax_vy = tab3_dspec_small_CTRLs_OPT['vmax_values_last'][1:]
            vmin_vx, vmin_vy = tab3_dspec_small_CTRLs_OPT['vmin_values_last'][1:]
            tab3_r_aia_submap_rect.data_source.data['x'] = [(vmax_vx + vmin_vx) / 2]
            tab3_r_aia_submap_rect.data_source.data['y'] = [(vmax_vy + vmin_vy) / 2]
            tab3_r_aia_submap_rect.data_source.data['width'] = [(vmax_vx - vmin_vx)]
            tab3_r_aia_submap_rect.data_source.data['height'] = [(vmax_vy - vmin_vy)]


        tab3_CTRLs_dspec_small = [tab3_Slider_dspec_small_dmax, tab3_Slider_dspec_small_dmin]
        for ctrl in tab3_CTRLs_dspec_small:
            ctrl.on_change('value', tab3_slider_dspec_small_update)

        # todo add the time/freq selection
        tab3_RBG_TimeFreq = RadioButtonGroup(labels=["time", "freq"], active=0)
        tab3_Slider_ANLYS_idx = Slider(start=0, end=tab2_ntim - 1, value=0, step=1, title="time idx", width=450)


        def tab3_slider_ANLYS_idx_update(attrname, old, new):
            global tab2_dtim, tab2_freq, tab2_ntim, SRC_maxfit_centroid
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                tab3_Slider_ANLYS_idx.start = next(
                    i for i in xrange(tab2_ntim) if tab2_dtim[i] >= tab3_p_dspec_small.x_range.start)
                tab3_Slider_ANLYS_idx.end = next(
                    i for i in xrange(tab2_ntim - 1, -1, -1) if tab2_dtim[i] <= tab3_p_dspec_small.x_range.end) + 1
                indices_time = tab3_Slider_ANLYS_idx.value
                tab3_r_dspec_small_line.visible = True
                tab3_r_dspec_small_line.data_source.data = ColumnDataSource(
                    pd.DataFrame({'time': [tab2_dtim[indices_time], tab2_dtim[indices_time]],
                                  'freq': [tab2_freq[0], tab2_freq[-1]]})).data
                try:
                    tab3_r_aia_submap_cross.visible = True
                    tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid[indices_time].data
                except:
                    tab3_r_aia_submap_cross.visible = False
            else:
                tab3_Div_Tb.text = """<p><b>Warning: Animate is OFF!!!</b></p>"""


        tab3_Slider_ANLYS_idx.on_change('value', tab3_slider_ANLYS_idx_update)

        tab3_Div_Tb = Div(text=""" """, width=400)


        def tab3_BUT_plot_xargs_default():
            global tab3_plot_xargs_dict
            tab3_plot_xargs_dict = OrderedDict()
            tab3_plot_xargs_dict['timebin'] = "1"
            tab3_plot_xargs_dict['timeline'] = "False"
            tab3_Div_plot_xargs_text = '<p>' + ';'.join(
                "<b>{}</b> = {}".format(key, val) for (key, val) in tab3_plot_xargs_dict.items()) + '</p>'
            tab3_Div_plot_xargs.text = tab3_Div_plot_xargs_text
            tab3_Div_Tb.text = '<p><b>Default xargs Restored.</b></p>'


        timebin = 1
        timeline = False
        tab3_animate_step = timebin


        def tab3_animate_update():
            global tab3_animate_step, tab2_dspec_selected
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                if tab2_dspec_selected:
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
            global tab2_dspec_selected
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                if tab3_BUT_PlayCTRL.label == 'Play':
                    if tab2_dspec_selected:
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


        tab3_BUT_PlayCTRL = Button(label='Play', width=60, button_type='success')
        tab3_BUT_PlayCTRL.on_click(tab3_animate)


        def tab3_animate_step_CTRL():
            global tab3_animate_step, tab2_dspec_selected
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                if tab2_dspec_selected:
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


        tab3_BUT_StepCTRL = Button(label='Step', width=60, button_type='primary')
        tab3_BUT_StepCTRL.on_click(tab3_animate_step_CTRL)


        def tab3_animate_FRWD_REVS():
            global tab3_animate_step
            if tab3_BUT_animate_ONOFF.label == 'Animate ON & Go':
                if tab2_dspec_selected:
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


        tab3_BUT_FRWD_REVS_CTRL = Button(label='Forward', width=60, button_type='warning')
        tab3_BUT_FRWD_REVS_CTRL.on_click(tab3_animate_FRWD_REVS)


        def tab3_animate_onoff():
            if tab2_dspec_selected:
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
                    tab3_r_dspec_small_line.visible = False
                    tab3_r_dspecvx_small_line.visible = False
                    tab3_r_dspecvy_small_line.visible = False
                    tab2_SRC_maxfit_centroid_update(dspecDF)
                    tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid.data
                else:
                    tab3_BUT_animate_ONOFF.label = 'Animate ON & Go'
                    tab3_r_aia_submap_cross.visible = True
                    tab3_r_aia_submap_line.visible = False
                    tab3_r_dspec_small_line.visible = True
                    tab3_r_dspecvx_small_line.visible = True
                    tab3_r_dspecvy_small_line.visible = True
                    tab2_SRC_maxfit_centroid_update(dspecDF)
                    indices_time = tab3_Slider_ANLYS_idx.value
                    tab3_r_aia_submap_cross.data_source.data = SRC_maxfit_centroid[indices_time].data
                    tab3_Div_Tb.text = """ """
            else:
                tab3_Div_Tb.text = """<p><b>Warning: Select time and frequency from the Dynamic Spectrum first!!!</b></p>"""


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

        # todo add AIA & HMI resolution selection
        tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                    options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                             "2048x2048"],
                                    width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])


        def tab2_update_MapRES(attrname, old, new):
            start_timestamp = time.time()
            select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
            dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
            aia_resampled_map = aiamap.resample(dimensions)
            aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'])
            SRC_AIA = aia_resampled_pfmap.ImageSource()
            tab2_r_aia.data_source.data['data'] = SRC_AIA.data['data']
            hmi_resampled_map = hmimap.resample(dimensions)
            hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'])
            SRC_HMI = hmi_resampled_pfmap.ImageSource()
            tab2_r_hmi.data_source.data['data'] = SRC_HMI.data['data']
            print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


        tab2_Select_MapRES.on_change('value', tab2_update_MapRES)
        # todo add the threshold selection (overplot another gylph)
        # todo add the dmax dmin and reset

        panel2 = column(
            row(gridplot([[tab2_p_aia, tab2_p_hmi, tab2_p_vla]], toolbar_location='right'),
                widgetbox(tab2_Select_MapRES, tab2_BUT_LinkImg_replot,
                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])),
            row(column(row(tab2_p_dspec, tab2_p_dspec_yPro),
                       row(tab2_p_dspec_xPro, tab2_p_dspec_thumb)),
                widgetbox(tab2_Select_pol, tab2_Select_bl,
                          tab2_Select_colorspace,
                          tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])))

        panel3 = row(column(tab3_p_aia_submap, tab3_Slider_ANLYS_idx,
                            row(tab3_BUT_PlayCTRL, tab3_SPCR_LFT_BUT_Step, tab3_BUT_StepCTRL,
                                tab3_SPCR_LFT_BUT_REVS_CTRL,
                                tab3_BUT_FRWD_REVS_CTRL, tab3_SPCR_LFT_BUT_animate_ONOFF,
                                tab3_BUT_animate_ONOFF), tab3_input_plot_xargs, tab3_Div_plot_xargs),
                     column(gridplot([tab3_p_dspec_small], [tab3_p_dspecvx_small], [tab3_p_dspecvy_small],
                                     toolbar_location='right'), tab3_Div_Tb),
                     widgetbox(tab3_RBG_dspec_small, tab3_Slider_dspec_small_dmax, tab3_Slider_dspec_small_dmin,
                               tab3_BUT_dspec_small_reset, tab2_panel3_BUT_exit, tab2_panel3_Div_exit,
                               width=200))
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

            TOOLS = "crosshair,pan,wheel_zoom,tap,box_zoom,reset,save"

            tab2_SRC_dspec_square = ColumnDataSource(dspecDF)

            '''create the dynamic spectrum plot'''
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
                                                line_color=None, line_alpha=0.0, selection_fill_alpha=0.1,
                                                selection_fill_color='black',
                                                nonselection_fill_alpha=0.0,
                                                selection_line_alpha=0.2, selection_line_color='white',
                                                nonselection_line_alpha=0.0,
                                                size=max(
                                                    config_plot['plot_config']['tab_FSview_base'][
                                                        'dspec_wdth'] / tab2_ntim,
                                                    config_plot['plot_config']['tab_FSview_base'][
                                                        'dspec_hght'] / tab2_nfreq))

            tab2_p_dspec.border_fill_color = "silver"
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

            map = Select(title="ColorSpace:", value="linear", options=["linear", "log"],
                         width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_p_dspec_xPro = figure(tools='',
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
            tab2_p_dspec_xPro.border_fill_color = "silver"
            tab2_p_dspec_xPro.border_fill_alpha = 0.4
            tab2_p_dspec_xPro.axis.major_tick_out = 0
            tab2_p_dspec_xPro.axis.major_tick_in = 5
            tab2_p_dspec_xPro.axis.minor_tick_out = 0
            tab2_p_dspec_xPro.axis.minor_tick_in = 3
            tab2_p_dspec_xPro.axis.major_tick_line_color = "black"
            tab2_p_dspec_xPro.axis.minor_tick_line_color = "black"

            tab2_p_dspec_yPro = figure(tools='',
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
            tab2_p_dspec_yPro.border_fill_color = "silver"
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


            def tab2_dspec_image_plt(select_pol):
                global spec_plt_R, spec_plt_L, spec_plt_I, spec_plt_V
                spec_plt_max_IRL = int(
                    max(spec_plt_R.max(), spec_plt_L.max(), spec_plt_I.max())) * 1.2
                spec_plt_min_IRL = (int(min(spec_plt_R.min(), spec_plt_L.min(), spec_plt_I.min())) / 10) * 10
                spec_plt_max_V = max(abs(int(spec_plt_V.max())), abs(int(spec_plt_V.min()))) * 1.2
                spec_plt_min_V = -spec_plt_max_V
                if select_pol == 'RR':
                    spec_plt = spec_plt_R
                    spec_plt_max = spec_plt_max_IRL
                    spec_plt_min = spec_plt_min_IRL
                elif select_pol == 'LL':
                    spec_plt = spec_plt_L
                    spec_plt_max = spec_plt_max_IRL
                    spec_plt_min = spec_plt_min_IRL
                elif select_pol == 'I':
                    spec_plt = spec_plt_I
                    spec_plt_max = spec_plt_max_IRL
                    spec_plt_min = spec_plt_min_IRL
                elif select_pol == 'V':
                    spec_plt = spec_plt_V
                    spec_plt_max = spec_plt_max_V
                    spec_plt_min = spec_plt_min_V
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


            def tab2_vdspec_update():
                global spec_plt_R, spec_plt_L, spec_plt_I, spec_plt_V
                select_pol = tab2_Select_pol.value
                tab2_vla_square_selected = tab2_SRC_vla_square.selected['1d']['indices']
                if tab2_BUT_vdspec.label == "VEC Dyn Spec":
                    if tab2_vla_square_selected:
                        tab2_BUT_vdspec.label = "Dyn Spec"
                        idxmax = max(tab2_vla_square_selected)
                        idxmin = min(tab2_vla_square_selected)
                        x0, x1 = idxmin % mapvlasize[0], idxmax % mapvlasize[0]
                        y0, y1 = idxmin / mapvlasize[0], idxmax / mapvlasize[0]
                        print x0, x1, y0, y1
                        spec_plt_R = np.zeros((tab2_nfreq, tab2_ntim))
                        spec_plt_L = np.zeros((tab2_nfreq, tab2_ntim))
                        spec_plt_I = np.zeros((tab2_nfreq, tab2_ntim))
                        spec_plt_V = np.zeros((tab2_nfreq, tab2_ntim))
                        for ll in xrange(tab2_ntim):
                            hdufile = fits_LOCL_dir + dspecDF0.loc[ll, :]['fits_local']
                            if os.path.exists(hdufile):
                                hdulist = fits.open(hdufile)
                                hdu = hdulist[0]
                                nfreq_hdu = hdu.header['NAXIS3']
                                freq_ref = '{:.3f}'.format(hdu.header['CRVAL3'] / 1e9)
                                freq = ['{:.3f}'.format(fq) for fq in tab2_freq]
                                idxfreq_ref = freq.index(freq_ref)
                                # select_pol == 'RR':
                                vladata = hdu.data[0, :, y0:y1 + 1, x0:x1 + 1]
                                spec_plt_R[idxfreq_ref:idxfreq_ref + nfreq_hdu, ll] = np.nanmean(vladata, axis=(-1, -2))
                                spec_plt_R[spec_plt_R < 0] = 0
                                # select_pol == 'LL':
                                vladata = hdu.data[1, :, y0:y1 + 1, x0:x1 + 1]
                                spec_plt_L[idxfreq_ref:idxfreq_ref + nfreq_hdu, ll] = np.nanmean(vladata, axis=(-1, -2))
                                spec_plt_L[spec_plt_L < 0] = 0
                                # select_pol == 'I':
                                vladata = hdu.data[0, :, y0:y1 + 1, x0:x1 + 1] + hdu.data[1, :, y0:y1 + 1, x0:x1 + 1]
                                spec_plt_I[idxfreq_ref:idxfreq_ref + nfreq_hdu, ll] = np.nanmean(vladata, axis=(-1, -2))
                                spec_plt_I[spec_plt_I < 0] = 0
                                # select_pol == 'V':
                                vladata = hdu.data[0, :, y0:y1 + 1, x0:x1 + 1] - hdu.data[1, :, y0:y1 + 1, x0:x1 + 1]
                                spec_plt_V[idxfreq_ref:idxfreq_ref + nfreq_hdu, ll] = np.nanmean(vladata, axis=(-1, -2))
                                spec_plt_V[spec_plt_V < 0] = 0
                        tab2_dspec_image_plt(select_pol)
                        tab2_p_dspec.title.text = "Vector Dynamic spectrum"
                    else:
                        tab2_Div_LinkImg_plot.text = '<p><b>Warning:</b> select a region first.</p>'
                else:
                    select_bl = tab2_Select_bl.value
                    tab2_BUT_vdspec.label = "VEC Dyn Spec"
                    bl_index = tab2_bl.index(select_bl)
                    spec_plt_R = tab2_spec[0, bl_index, :, :]
                    spec_plt_L = tab2_spec[1, bl_index, :, :]
                    spec_plt_I = (tab2_spec[0, bl_index, :, :] + tab2_spec[1, bl_index, :, :]) / 2.
                    spec_plt_V = (tab2_spec[0, bl_index, :, :] - tab2_spec[1, bl_index, :, :]) / 2.
                    tab2_dspec_image_plt(select_pol)
                    tab2_p_dspec.title.text = "Dynamic spectrum"


            tab2_BUT_vdspec.on_click(tab2_vdspec_update)


            def tab2_update_dspec_image(attrname, old, new):
                global spec_plt_R, spec_plt_L, spec_plt_I, spec_plt_V
                global tab2_spec, tab2_dtim, tab2_freq, tab2_bl
                select_pol = tab2_Select_pol.value
                select_bl = tab2_Select_bl.value
                bl_index = tab2_bl.index(select_bl)
                if tab2_BUT_vdspec.label == "VEC Dyn Spec":
                    spec_plt_R = tab2_spec[0, bl_index, :, :]
                    spec_plt_L = tab2_spec[1, bl_index, :, :]
                    spec_plt_I = (tab2_spec[0, bl_index, :, :] + tab2_spec[1, bl_index, :, :]) / 2.
                    spec_plt_V = (tab2_spec[0, bl_index, :, :] - tab2_spec[1, bl_index, :, :]) / 2.
                tab2_dspec_image_plt(select_pol)


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
            tab2_SRC_vlamap_peak = ColumnDataSource(data={'dspec': [], 'x_pos': [], 'y_pos': [], 'amp_gaus': []})

            # import the vla image
            hdulist = fits.open(vlafile[0])
            hdu = hdulist[0]
            vla_local_pfmap = PuffinMap(hdu.data[0, 0, :, :], hdu.header,
                                        plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                        plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'])
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
            # import the aia image
            # from sunpy.net.helioviewer import HelioviewerClient
            #
            # hv = HelioviewerClient()
            # filepath = hv.download_jp2(jdutil.jd_to_datetime(xx[0] / 3600. / 24.), observatory='SDO', instrument='AIA',
            #                            detector='AIA', measurement='171',
            #                            directory=database_dir + event_id + struct_id + config_EvtID['datadir']['J2000'],
            #                            overwrite=True)
            filepath = database_dir + event_id + struct_id + config_EvtID['datadir'][
                'J2000'] + '2014_11_01__16_45_59_34__SDO_AIA_AIA_171.jp2'
            colormap = cm.get_cmap("sdoaia171")  # choose any matplotlib colormap here
            bokehpalette_sdoaia171 = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            aiamap = sunpy.map.Map(filepath)
            MapRES = 256
            dimensions = u.Quantity([MapRES, MapRES], u.pixel)
            aia_resampled_map = aiamap.resample(dimensions)

            # plot the global AIA image

            aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'])

            tab2_p_aia, tab2_r_aia = aia_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                                 palette=bokehpalette_sdoaia171)
            tab2_p_aia.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                                  line_width=2)
            tab2_p_aia.title.text_font_size = '6pt'
            tab2_p_aia.border_fill_color = "silver"
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

            # plot the global HMI image
            # filepath = hv.download_jp2(jdutil.jd_to_datetime(xx[0] / 3600. / 24.), observatory='SDO', instrument='HMI',
            #                            detector='HMI', measurement='magnetogram',
            #                            directory=database_dir + event_id + struct_id + config_EvtID['datadir']['J2000'],
            #                            overwrite=True)
            filepath = database_dir + event_id + struct_id + config_EvtID['datadir'][
                'J2000'] + '2014_11_01__16_46_58_605__SDO_HMI_HMI_magnetogram.jp2'
            colormap = cm.get_cmap("gray")  # choose any matplotlib colormap here
            bokehpalette_sdohmimag = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            hmimap = sunpy.map.Map(filepath)
            # plot the global HMI image
            hmi_resampled_map = hmimap.resample(dimensions)
            hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                            plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                            plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'])

            tab2_p_hmi, tab2_r_hmi = hmi_resampled_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                                                 x_range=tab2_p_aia.x_range,
                                                                 y_range=tab2_p_aia.y_range,
                                                                 palette=bokehpalette_sdohmimag)
            tab2_p_hmi.multi_line(xs='xs', ys='ys', line_color='line_color', source=tab2_SRC_vlamap_contour, alpha=0.7,
                                  line_width=2)
            tab2_p_hmi.yaxis.visible = False
            tab2_p_hmi.border_fill_color = "silver"
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
            tab2_p_vla.border_fill_color = "silver"
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

            tab2_Select_vla_pol = Select(title="Polarization:", value='RR', options=['RR', 'LL', 'I', 'V'],
                                         width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])

            tab2_source_idx_line_x = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
            tab2_r_dspec_line_x = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                    line_color='white', source=tab2_source_idx_line_x)
            tab2_source_idx_line_y = ColumnDataSource(pd.DataFrame({'time': [], 'freq': []}))
            tab2_r_dspec_line_y = tab2_p_dspec.line(x='time', y='freq', line_width=1.5, line_alpha=0.8,
                                                    line_color='white', source=tab2_source_idx_line_y)


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
                hdufile = fits_LOCL_dir + dspecDF0.loc[tidx, :]['fits_local']
                if os.path.exists(hdufile):
                    hdulist = fits.open(hdufile)
                    hdu = hdulist[0]
                    if select_vla_pol == 'RR':
                        vladata = hdu.data[0, fidx, :, :]
                    elif select_vla_pol == 'LL':
                        vladata = hdu.data[1, fidx, :, :]
                    elif select_vla_pol == 'I':
                        vladata = hdu.data[0, fidx, :, :] + hdu.data[1, fidx, :, :]
                    elif select_vla_pol == 'V':
                        vladata = hdu.data[0, fidx, :, :] - hdu.data[1, fidx, :, :]
                    pfmap = PuffinMap(vladata, hdu.header, plot_height=tab2_LinkImg_HGHT,
                                      plot_width=tab2_LinkImg_WDTH)
                    SRC_Img = pfmap.ImageSource()
                    tab2_r_vla.data_source.data['data'] = SRC_Img.data['data']
                    mapx, mapy = pfmap.meshgrid()
                    mapx, mapy = mapx.value, mapy.value
                    SRC_contour = get_contour_data(mapx, mapy, pfmap.smap.data)
                    tab2_r_vla_multi_line.data_source.data = SRC_contour.data
                    tab2_Div_LinkImg_plot.text = '<p><b>{}</b> loaded.</p>'.format(
                        dspecDF0.loc[tidx, :]['fits_local'])
                else:
                    tab2_Div_LinkImg_plot.text = '<p><b>{}</b> not found.</p>'.format(
                        dspecDF0.loc[tidx, :]['fits_local'])


            tab2_CTRLs_LinkImg = [tab2_Slider_time_LinkImg, tab2_Slider_freq_LinkImg, tab2_Select_vla_pol]
            for ctrl in tab2_CTRLs_LinkImg:
                ctrl.on_change('value', tab3_slider_LinkImg_update)


            def tab2_panel_exit():
                tab2_panel2_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
                raise SystemExit


            tab2_panel2_BUT_exit = Button(label='Exit FSview',
                                          width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                          button_type='danger')
            tab2_panel2_BUT_exit.on_click(tab2_panel_exit)

            tab2_dspec_selected = None


            def dspecDFtmp_init():
                global dspecDFtmp
                dspecDFtmp = pd.DataFrame()
                nrows_dspecDF = len(dspecDF0.index)
                dspecDFtmp['amp_gaus'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
                dspecDFtmp['x_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)
                dspecDFtmp['y_pos'] = pd.Series([np.nan] * nrows_dspecDF, index=dspecDF0.index)


            def tab2_dspec_selection_change(attrname, old, new):
                global tab2_dspec_selected
                tab2_dspec_selected = tab2_SRC_dspec_square.selected['1d']['indices']
                if tab2_dspec_selected:
                    global dspecDF
                    dspecDF = dspecDF0.iloc[tab2_dspec_selected, :]
                    idx_selected = dspecDF.index[len(dspecDF) / 2]
                    # print ['{:.3f}'.format(ll) for ll in tab2_dtim], '{:.3f}'.format(
                    #     dspecDF0.loc[idx_selected, :]['time'])
                    # print ['{:.3f}'.format(ll) for ll in tab2_freq], '{:.3f}'.format(
                    #     dspecDF0.loc[idx_selected, :]['freq'])
                    tidx = int(['{:.3f}'.format(ll) for ll in tab2_dtim].index(
                        '{:.3f}'.format(dspecDF0.loc[idx_selected, :]['time'])))
                    fidx = int(['{:.3f}'.format(ll) for ll in tab2_freq].index(
                        '{:.3f}'.format(dspecDF0.loc[idx_selected, :]['freq'])))
                    tab2_Slider_time_LinkImg.value = tidx
                    tab2_Slider_freq_LinkImg.value = fidx


            tab2_SRC_dspec_square.on_change('selected', tab2_dspec_selection_change)
            tab2_vla_square_selected = None


            def tab2_vla_square_selection_change(attrname, old, new):
                global x0, x1, y0, y1, x0pix, x1pix, y0pix, y1pix
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


            tab2_SRC_vla_square.on_change('selected', tab2_vla_square_selection_change)

            # todo add AIA & HMI resolution selection
            tab2_Select_MapRES = Select(title="Img resolution:", value='{}x{}'.format(MapRES, MapRES),
                                        options=["32x32", "64x64", "128x128", "256x256", "512x512", "1024x1024",
                                                 "2048x2048"],
                                        width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])


            def tab2_update_MapRES(attrname, old, new):
                start_timestamp = time.time()
                select_MapRES = int(tab2_Select_MapRES.value.split('x')[0])
                dimensions = u.Quantity([select_MapRES, select_MapRES], u.pixel)
                aia_resampled_map = aiamap.resample(dimensions)
                aia_resampled_pfmap = PuffinMap(smap=aia_resampled_map,
                                                plot_height=config_plot['plot_config']['tab_FSview_base']['aia_hght'],
                                                plot_width=config_plot['plot_config']['tab_FSview_base']['aia_wdth'])
                SRC_AIA = aia_resampled_pfmap.ImageSource()
                tab2_r_aia.data_source.data['data'] = SRC_AIA.data['data']
                hmi_resampled_map = hmimap.resample(dimensions)
                hmi_resampled_pfmap = PuffinMap(smap=hmi_resampled_map,
                                                plot_height=config_plot['plot_config']['tab_FSview_base']['vla_hght'],
                                                plot_width=config_plot['plot_config']['tab_FSview_base']['vla_wdth'])
                SRC_HMI = hmi_resampled_pfmap.ImageSource()
                tab2_r_hmi.data_source.data['data'] = SRC_HMI.data['data']
                print("---tab2_update_MapRES -- %s seconds ---" % (time.time() - start_timestamp))


            tab2_Select_MapRES.on_change('value', tab2_update_MapRES)


            # todo add the threshold selection (overplot another gylph)


            def tab2_save_region():
                pangle = hdu.header['p_angle']
                x0Deg, x1Deg, y0Deg, y1Deg = (x0 - hdu.header['CRVAL1']) / 3600., (x1 - hdu.header['CRVAL1']) / 3600., (
                    y0 - hdu.header['CRVAL2']) / 3600., (y1 - hdu.header['CRVAL2']) / 3600.
                # todo find p0 in fits header
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
                rgnfits = '#CRTFv0 CASA Region Text Format version 0\nbox [[{}], [{}]] coord=J2000, linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font="DejaVu Sans", fontsize=11, fontstyle=normal, usetex=false'.format(
                    ', '.join(c0fits.to_string('hmsdms').split(' ')), ', '.join(c1fits.to_string('hmsdms').split(' ')))
                print rgnfits
                rgnfitsfile = database_dir + event_id + struct_id + "CASA_imfit_region_fits.rgn"
                with open(rgnfitsfile, "w") as fp:
                    fp.write(rgnfits)
                tab2_Div_LinkImg_plot.text = '<p>region saved to <b>{}</b>.</p>'.format(rgnfitsfile)


            tab2_BUT_SavRgn = Button(label='Save Region',
                                     width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                                     button_type='primary')
            tab2_BUT_SavRgn.on_click(tab2_save_region)

            panel2 = column(
                row(gridplot([[tab2_p_aia, tab2_p_hmi, tab2_p_vla]], toolbar_location='right'),
                    widgetbox(tab2_Select_MapRES, tab2_Select_vla_pol, tab2_Slider_time_LinkImg,
                              tab2_Slider_freq_LinkImg, tab2_BUT_vdspec, tab2_BUT_SavRgn, tab2_Div_LinkImg_plot,
                              width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])),
                row(column(row(tab2_p_dspec, tab2_p_dspec_yPro),
                           tab2_p_dspec_xPro),
                    widgetbox(tab2_Select_pol, tab2_Select_bl,
                              tab2_Select_colorspace,
                              tab2_panel2_BUT_exit, tab2_panel2_Div_exit,
                              width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'])))

            lout = panel2

            # def timeout_callback():
            #     print 'timeout'
            #     raise SystemExit


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
    tab2_p_dspec.border_fill_color = "silver"
    tab2_p_dspec.border_fill_alpha = 0.4
    tab2_p_dspec.axis.major_tick_out = 0
    tab2_p_dspec.axis.major_tick_in = 5
    tab2_p_dspec.axis.minor_tick_out = 0
    tab2_p_dspec.axis.minor_tick_in = 3
    tab2_p_dspec.axis.major_tick_line_color = "white"
    tab2_p_dspec.axis.minor_tick_line_color = "white"

    tab2_dspec_selected = None


    def tab2_dspec_selection_change(attrname, old, new):
        global tab2_dspec_selected
        tab2_dspec_selected = tab2_SRC_dspec_square.selected['1d']['indices']
        if tab2_dspec_selected:
            global dspecDF
            dspecDF = dspecDF0.copy()
            dspecDF = dspecDF.iloc[tab2_dspec_selected, :]
            x0, x1 = dspecDF['time'].min(), dspecDF['time'].max()
            y0, y1 = dspecDF['freq'].min(), dspecDF['freq'].max()
            tab2_r_dspec_patch.data_source.data = ColumnDataSource(
                pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
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

    tab2_p_dspec_xPro = figure(tools='', plot_width=config_plot['plot_config']['tab_FSview2CASA']['dspec_xPro_wdth'],
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
    tab2_p_dspec_xPro.border_fill_color = "silver"
    tab2_p_dspec_xPro.border_fill_alpha = 0.4
    tab2_p_dspec_xPro.axis.major_tick_out = 0
    tab2_p_dspec_xPro.axis.major_tick_in = 5
    tab2_p_dspec_xPro.axis.minor_tick_out = 0
    tab2_p_dspec_xPro.axis.minor_tick_in = 3
    tab2_p_dspec_xPro.axis.major_tick_line_color = "black"
    tab2_p_dspec_xPro.axis.minor_tick_line_color = "black"

    tab2_p_dspec_yPro = figure(tools='', plot_width=config_plot['plot_config']['tab_FSview2CASA']['dspec_yPro_wdth'],
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
    tab2_p_dspec_yPro.border_fill_color = "silver"
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
        if select_pol == 'RR':
            spec_plt = spec_plt_R
            spec_plt_max = spec_plt_max_IRL
            spec_plt_min = spec_plt_min_IRL
        elif select_pol == 'LL':
            spec_plt = spec_plt_L
            spec_plt_max = spec_plt_max_IRL
            spec_plt_min = spec_plt_min_IRL
        elif select_pol == 'I':
            spec_plt = spec_plt_I
            spec_plt_max = spec_plt_max_IRL
            spec_plt_min = spec_plt_min_IRL
        elif select_pol == 'V':
            spec_plt = spec_plt_V
            spec_plt_max = spec_plt_max_V
            spec_plt_min = spec_plt_min_V
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

    timestart = xx[0]


    def tab2_BUT_tCLN_param_add():
        tab2_Div_tCLN2.text = ' '
        txts = tab2_input_tCLN.value.strip()
        txts = txts.split(';')
        for txt in txts:
            txt = txt.strip()
            # todo box selection to timerange and freqrange automaticlly for clean in FSview2CASA
            if txt == 'timerange':
                time0, time1 = dspecDF['time'].min() + timestart, dspecDF['time'].max() + timestart
                date_char = Time(timestart / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date').iso
                t0_char = Time(time0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
                t0_char = t0_char.split(' ')[1]
                t1_char = Time(time1 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
                t1_char = t1_char.split(' ')[1]
                tab2_tCLN_Param_dict['timerange'] = "'{}T{}~{}T{}'".format(date_char, t0_char, date_char, t1_char)
            elif txt == 'freqrange':
                freq0, freq1 = dspecDF['freq'].min(), dspecDF['freq'].max()
                freqrange = "'{:.3f}~{:.3f} GHz'".format(freq0, freq1)
                tab2_tCLN_Param_dict['freqrange'] = freqrange
            else:
                txt = txt.split('=')
                if len(txt) == 2:
                    key, val = txt
                    tab2_tCLN_Param_dict[key.strip()] = val.strip()
                else:
                    tab2_Div_tCLN2.text = '<p>Input syntax: <b>uvtaper</b>=True; <b>niter</b>=200; ' \
                                          '<b>cell</b>=["5.0arcsec", "5.0arcsec"]. Any spaces will be ignored.</p>'
        tab2_Div_tCLN_text = ' '.join(
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
        tab2_Div_tCLN_text = ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
        tab2_Div_tCLN.text = tab2_Div_tCLN_text


    def tab2_BUT_tCLN_param_default():
        global tab2_tCLN_Param_dict
        tab2_tCLN_Param_dict = OrderedDict()
        tab2_tCLN_Param_dict['vis'] = "''"
        tab2_tCLN_Param_dict['imageprefix'] = "'slfcal/{}'".format(struct_id)
        tab2_tCLN_Param_dict['ncpu'] = "10"
        tab2_tCLN_Param_dict['twidth'] = "1"
        tab2_tCLN_Param_dict['doreg'] = "True"
        tab2_tCLN_Param_dict['timerange'] = "''"
        tab2_tCLN_Param_dict['uvrange'] = "''"
        tab2_tCLN_Param_dict['antenna'] = "''"
        tab2_tCLN_Param_dict['ephemfile'] = "'horizons_sun_20141101.radecp'"
        tab2_tCLN_Param_dict['msinfofile'] = "'SUN01_20141101.T163940-164700.50ms.cal.msinfo.npz'"
        tab2_tCLN_Param_dict['struct_id'] = "'{}'".format(struct_id.replace("/", ""))
        tab2_tCLN_Param_dict['event_id'] = "'{}'".format(event_id.replace("/", ""))
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
        tab2_tCLN_Param_dict['outertaper'] = "[]"
        tab2_tCLN_Param_dict['uvrange'] = "''"
        tab2_tCLN_Param_dict['niter'] = "200"
        tab2_tCLN_Param_dict['usescratch'] = "False"
        tab2_tCLN_Param_dict['interactive'] = "False"
        tab2_Div_tCLN_text = ' '.join(
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
    tab2_SPCR_LFT_BUT_tCLN_param_FSVIEW = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_BUT_tCLN_param_DEFAULT = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
    tab2_SPCR_LFT_Div_tCLN2 = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth50'],
                                     height=config_plot['plot_config']['tab_FSview2CASA']['space_hght10'])


    def tab2_BUT_tCLN_param_reload():
        with open(database_dir + event_id + struct_id + 'CASA_CLN_args.json', 'r') as fp:
            tab2_tCLN_Param_dict = json.load(fp)
        tab2_Div_tCLN_text = ' '.join(
            "<p><b>{}</b> = {}</p>".format(key, val) for (key, val) in tab2_tCLN_Param_dict.items())
        tab2_Div_tCLN.text = tab2_Div_tCLN_text
        tab2_Div_tCLN2.text = '<p>CASA arguments reload from config file in <b>{}</b>.</p>'.format(
            database_dir + event_id + struct_id)


    def tab2_BUT_tCLN_clean():
        with open(database_dir + event_id + struct_id + 'CASA_CLN_args.json', 'w') as fp:
            json.dump(tab2_tCLN_Param_dict, fp)
        os.system('cp FSview/script_process.py {}'.format(database_dir + event_id + struct_id))
        # os.system('cp FSview/script_preprocess.py {}'.format(database_dir + event_id + struct_id))
        tab2_Div_tCLN2.text = '<p>CASA script and arguments config file saved to <b>{}</b>.</p>'.format(
            database_dir + event_id + struct_id)
        timestrs = []
        fits_local = []
        fits_global = []
        for ii in range(len(xx)):
            t0 = xx[ii]
            datestr = Time(t0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date').iso
            timestr0 = Time(t0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
            timestr0 = timestr0.split(' ')[1]
            timestr = datestr.replace("-", "") + 'T' + timestr0.replace(":", "")
            timestrs.append(timestr0)
            fits_local.append(timestr + '.fits')
            fits_global.append(timestr + '.fits')
        dspecDF_tmp = pd.DataFrame({'time': xx - xx[0],
                                    'freq': yy,
                                    'timestr': timestrs,
                                    'dspec': tab2_spec_plt.flatten(),
                                    'fits_local': fits_local,
                                    'fits_global': fits_global})
        with open(database_dir + event_id + struct_id + 'dspecDF-save', 'wb') as fp:
            pickle.dump(dspecDF_tmp, fp)
        tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}</b>. '.format(
            database_dir + event_id + struct_id) + 'Go back to <b>QLook</b> window, select StrID <b>{}</b> and \
            click <b>FSview</b> button again.</p>'.format(
            database_dir + event_id + struct_id, struct_id[0:-1])
        cwd = os.getcwd()
        try:
            os.chdir(database_dir + event_id + struct_id)
            suncasapy = config_plot['core']['casapy']
            suncasapy = os.path.expandvars(suncasapy)
            os.system('{} -c script_process.py'.format(suncasapy))
            tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}</b></p>. <p>CASA clean is in processing.</p>'.format(
                database_dir + event_id + struct_id) + '<p>When finished, go back to <b>QLook</b> window, select StrID <b>{}</b> and \
                click <b>FSview</b> button again.</p>'.format(
                database_dir + event_id + struct_id, struct_id[0:-1])
        except:
            tab2_Div_tCLN2.text = '<p>CASA script, arguments config file and dspecDF-save saved to <b>{}</b></p>. <p>Do image clean with CASA manually.</p>'.format(
                database_dir + event_id + struct_id) + '<p>When finished, go back to <b>QLook</b> window, select StrID <b>{}</b> and \
                click <b>FSview</b> button again.</p>'.format(
                database_dir + event_id + struct_id, struct_id[0:-1])
        os.chdir(cwd)

    tab2_BUT_tCLN_param_RELOAD = Button(label='reload Param',
                                        width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'])
    tab2_BUT_tCLN_param_RELOAD.on_click(tab2_BUT_tCLN_param_reload)
    tab2_SPCR_LFT_BUT_tCLN_param_reload = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth20'],
                                                 height=config_plot['plot_config']['tab_FSview2CASA']['space_hght10'])
    tab2_BUT_tCLN_CLEAN = Button(label='ToClean',
                                 width=config_plot['plot_config']['tab_FSview2CASA']['button_wdth'],
                                 button_type='success')
    tab2_SPCR_ABV_BUT_tCLN = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'],
                                    height=config_plot['plot_config']['tab_FSview2CASA']['space_hght18'])
    tab2_BUT_tCLN_CLEAN.on_click(tab2_BUT_tCLN_clean)
    tab2_SPCR_LFT_BUT_FS_view = Spacer(width=config_plot['plot_config']['tab_FSview2CASA']['space_wdth10'])
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
                                                  tab2_BUT_tCLN_param_RELOAD, tab2_SPCR_LFT_BUT_FS_view,
                                                  tab2_BUT_tCLN_CLEAN)),
                                       widgetbox(tab2_Select_pol, tab2_Select_bl, tab2_Select_colorspace,
                                                 tab2_panel2_BUT_exit,
                                                 tab2_panel_Div_exit,
                                                 width=config_plot['plot_config']['tab_FSview2CASA'][
                                                     'widgetbox_wdth2'])), tab2_Div_tCLN2))), tab2_Div_tCLN)

    curdoc().add_root(panel2)
    curdoc().title = "FSview2CASA"
