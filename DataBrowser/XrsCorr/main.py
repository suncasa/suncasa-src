import json
import os, sys
import pickle
import time
# import bokeh.palettes as bp
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Slider, Button, TextInput, RadioButtonGroup, CheckboxGroup,
                          BoxSelectTool, LassoSelectTool, HoverTool, Spacer, LabelSet, Div)
from bokeh.models import LinearColorMapper, ColorBar, LogColorMapper
from bokeh.plotting import figure, curdoc
from astropy.time import Time
from suncasa.utils import DButil

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
'''load config file'''
with open(suncasa_dir + 'DataBrowser/config.json', 'r') as fp:
    config_plot = json.load(fp)
database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/'
with open('{}config_EvtID_curr.json'.format(database_dir), 'r') as fp:
    config_EvtID = json.load(fp)

# do_spec_regrid = False

'''define the colormaps'''
colormap = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
colormap = cm.get_cmap("RdBu")  # choose any matplotlib colormap here
bokehpalette_RdBu = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
colormap = cm.get_cmap("Blues")  # choose any matplotlib colormap here
bokehpalette_Blues = [colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]


def exit():
    raise SystemExit


event_id = config_EvtID['datadir']['event_id']
try:
    with open(database_dir + event_id + 'CurrFS.json', 'r') as fp:
        FS_config = json.load(fp)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']

event_id = config_EvtID['datadir']['event_id']
try:
    with open(database_dir + event_id + 'CurrFS.json', 'r') as fp:
        FS_config = json.load(fp)
except:
    print 'Error: No CurrFS.json found!!!'
    raise SystemExit
struct_id = FS_config['datadir']['struct_id']
CC_save = database_dir + event_id + struct_id + 'CC_save.npz'
# try:
CC_savedata = np.load(CC_save)
spec = CC_savedata['spec']
ccmax = CC_savedata['ccmax']
ccpeak = CC_savedata['ccpeak']
tim = CC_savedata['tim']
dtim = tim - tim[0]
dt = np.median(np.diff(tim))
freq = CC_savedata['freq']
df = np.median(np.diff(freq))
ntim = len(tim)
nfreq = len(freq)
freqa = CC_savedata['freqa']
freqv = CC_savedata['freqv']
fidxa = CC_savedata['fidxa']
fidxv = CC_savedata['fidxv']
timfit = CC_savedata['timfit']
dtimfit = timfit - timfit[0]
ntimfit = CC_savedata['ntimfit']
specfit = CC_savedata['specfit']
dtfit = np.median(np.diff(timfit))

CCmaxDF = pd.DataFrame(
    {'ccmax': ccmax.ravel(), 'ccpeak': ccpeak.ravel() * dtfit * 1000.0, 'freqa': freqa.ravel(), 'freqv': freqv.ravel(),
     'fidxa': fidxa.ravel(), 'fidxv': fidxv.ravel()})
TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
p_dspec = figure(tools=TOOLS,
                 plot_width=config_plot['plot_config']['tab_XrsCorr']['dspec_wdth'] +
                            config_plot['plot_config']['tab_XrsCorr']['dspec_wdth_offset'],
                 plot_height=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                 x_range=(dtim[0], dtim[-1]), y_range=(freq[0], freq[-1]),
                 toolbar_location="above")
tim0_char = Time(tim[0] / 3600. / 24., format='mjd', scale='utc', precision=3, out_subfmt='date_hms').iso
p_dspec.axis.visible = True
p_dspec.title.text = "Dynamic spectrum"
p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
p_dspec.yaxis.axis_label = 'Frequency [GHz]'
p_dspec.border_fill_alpha = 0.4
p_dspec.axis.major_tick_out = 0
p_dspec.axis.major_tick_in = 5
p_dspec.axis.minor_tick_out = 0
p_dspec.axis.minor_tick_in = 3
p_dspec.axis.major_tick_line_color = "white"
p_dspec.axis.minor_tick_line_color = "white"
r_dspec = p_dspec.image(image=[spec], x=dtim[0], y=freq[0],
                        dw=dtim[-1] - dtim[0], dh=freq[-1] - freq[0],
                        color_mapper=LinearColorMapper(palette=bokehpalette_jet, low=np.amin(spec),
                                                       high=np.amax(spec)))

TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
p_dspecfit = figure(tools=TOOLS,
                    plot_width=config_plot['plot_config']['tab_XrsCorr']['dspec_wdth'] -
                               config_plot['plot_config']['tab_XrsCorr']['dspec_wdth_offset'],
                    plot_height=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                    x_range=p_dspec.x_range, y_range=p_dspec.y_range,
                    toolbar_location="above")
tim0_char = Time(tim[0] / 3600. / 24., format='mjd', scale='utc', precision=3, out_subfmt='date_hms').iso
p_dspecfit.axis.visible = True
# p_dspecfit.yaxis.visible = False
p_dspecfit.title.text = "Interpolated Dynamic spectrum"
p_dspecfit.xaxis.axis_label = 'Seconds since ' + tim0_char
p_dspecfit.yaxis.axis_label = 'Frequency [GHz]'
p_dspecfit.border_fill_alpha = 0.4
p_dspecfit.axis.major_tick_out = 0
p_dspecfit.axis.major_tick_in = 5
p_dspecfit.axis.minor_tick_out = 0
p_dspecfit.axis.minor_tick_in = 3
p_dspecfit.axis.major_tick_line_color = "white"
p_dspecfit.axis.minor_tick_line_color = "white"
r_dspecfit = p_dspecfit.image(image=[specfit], x=dtimfit[0], y=freq[0],
                              dw=dtimfit[-1] - dtimfit[0], dh=freq[-1] - freq[0],
                              color_mapper=LinearColorMapper(palette=bokehpalette_jet, low=np.amin(spec),
                                                             high=np.amax(spec)))

TOOLS = "save"
SRC_CCmax_square = ColumnDataSource(CCmaxDF)
p_CCmax = figure(tools=TOOLS,
                 plot_width=config_plot['plot_config']['tab_XrsCorr']['CCmap_wdth'],
                 plot_height=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                 x_range=(freq[0], freq[-1]), y_range=(freq[0], freq[-1]))
p_CCmax.axis.visible = True
p_CCmax.title.text = "Xross Correlation maximum"
p_CCmax.xaxis.axis_label = 'Frequency [GHz]'
p_CCmax.yaxis.axis_label = 'Frequency [GHz]'
ccmaxplt = ccmax
r_CCmax = p_CCmax.image(image=[ccmaxplt], x=freq[0], y=freq[0],
                        dw=freq[-1] - freq[0],
                        dh=freq[-1] - freq[0],
                        color_mapper=LinearColorMapper(palette=bokehpalette_Blues, low=0.5, high=1))
r_CCmax_square = p_CCmax.square('freqa', 'freqv', source=SRC_CCmax_square, fill_color=None,
                                fill_alpha=0.0,
                                line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                selection_fill_color=None,
                                nonselection_fill_alpha=0.0,
                                selection_line_alpha=0.0, selection_line_color=None,
                                nonselection_line_alpha=0.0,
                                size=min(config_plot['plot_config']['tab_XrsCorr']['dspec_wdth'] / (nfreq - 1),
                                         config_plot['plot_config']['tab_XrsCorr']['dspec_hght'] / (nfreq - 1)))
SRC_CCmax_line1 = ColumnDataSource({'x': [], 'y': []})
r_CCmax_line1 = p_CCmax.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange', source=SRC_CCmax_line1)
SRC_CCmax_line2 = ColumnDataSource({'x': [], 'y': []})
r_CCmax_line2 = p_CCmax.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red', source=SRC_CCmax_line2)

cm_CCmax = LinearColorMapper(palette=bokehpalette_Blues, low=0.5, high=1.0)
cb_CCmax = ColorBar(color_mapper=cm_CCmax, label_standoff=5, width=5, border_line_color=None, location=(0, 0))
p_CCmax.add_layout(cb_CCmax, 'right')

p_CCpeak = figure(tools=TOOLS,
                  plot_width=config_plot['plot_config']['tab_XrsCorr']['CCmap_wdth'],
                  plot_height=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                  x_range=p_CCmax.x_range, y_range=p_CCmax.y_range)
p_CCpeak.axis.visible = True
p_CCpeak.title.text = "Xcorr lag [ms]"
p_CCpeak.xaxis.axis_label = 'Frequency [GHz]'
p_CCpeak.yaxis.axis_label = 'Frequency [GHz]'

lagmaxd = (ntimfit - 1) / 2
lagmax = lagmaxd * dtfit * 1000.0
ccpeakplt = ccpeak * dtfit * 1000.0

r_CCpeak = p_CCpeak.image(image=[ccpeakplt], x=freq[0], y=freq[0],
                          dw=freq[-1] - freq[0],
                          dh=freq[-1] - freq[0],
                          color_mapper=LinearColorMapper(palette=bokehpalette_RdBu, low=-lagmax / 3, high=lagmax / 3))

r_CCpeak_square = p_CCpeak.square('freqa', 'freqv', source=SRC_CCmax_square, fill_color=None,
                                  fill_alpha=0.0,
                                  line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                  selection_fill_color=None,
                                  nonselection_fill_alpha=0.0,
                                  selection_line_alpha=0.0, selection_line_color=None,
                                  nonselection_line_alpha=0.0,
                                  size=min(config_plot['plot_config']['tab_XrsCorr']['dspec_wdth'] / (nfreq - 1),
                                           config_plot['plot_config']['tab_XrsCorr']['dspec_hght'] / (nfreq - 1)))
SRC_CCpeak_line1 = ColumnDataSource({'x': [], 'y': []})
r_CCpeak_line1 = p_CCpeak.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange', source=SRC_CCmax_line1)
SRC_CCpeak_line2 = ColumnDataSource({'x': [], 'y': []})
r_CCpeak_line2 = p_CCpeak.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red', source=SRC_CCmax_line2)

cm_CCpeak = LinearColorMapper(palette=bokehpalette_RdBu, low=-lagmax / 3, high=lagmax / 3)
cb_CCpeak = ColorBar(color_mapper=cm_CCpeak, label_standoff=5, width=5, border_line_color=None, location=(0, 0))
p_CCpeak.add_layout(cb_CCpeak, 'right')

SRC_freq_line1 = ColumnDataSource({'x': [], 'y': []})
r_freq_line1 = p_dspec.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange', source=SRC_freq_line1)
SRC_freq_line2 = ColumnDataSource({'x': [], 'y': []})
r_freq_line2 = p_dspec.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red', source=SRC_freq_line2)
r_fit_freq_line1 = p_dspecfit.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange',
                                   source=r_freq_line1.data_source)
r_fit_freq_line2 = p_dspecfit.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red',
                                   source=r_freq_line2.data_source)

TOOLS = "crosshair,save"
p_dspec_lines = figure(tools=TOOLS,
                       plot_width=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                       plot_height=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                       x_range=(dtim[0], dtim[-1]), y_range=(np.amin(spec), np.amax(spec)),
                       toolbar_location="above")
p_dspec_lines.axis.visible = True
p_dspec_lines.title.text = "Light curves"
p_dspec_lines.xaxis.axis_label = 'Seconds since ' + tim0_char
p_dspec_lines.yaxis.axis_label = 'Flux [sfu]'
p_dspec_lines.border_fill_alpha = 0.4
p_dspec_lines.axis.major_tick_out = 0
p_dspec_lines.axis.major_tick_in = 5
p_dspec_lines.axis.minor_tick_out = 0
p_dspec_lines.axis.minor_tick_in = 3
p_dspec_lines.axis.major_tick_line_color = "white"
p_dspec_lines.axis.minor_tick_line_color = "white"

SRC_dspec_prof1 = ColumnDataSource({'x': [], 'y': []})
r_dspec_prof1 = p_dspec_lines.line(x='x', y='y', alpha=0.3, line_width=2, line_color='orange', source=SRC_dspec_prof1)
SRC_dspec_prof2 = ColumnDataSource({'x': [], 'y': []})
r_dspec_prof2 = p_dspec_lines.line(x='x', y='y', alpha=0.3, line_width=2, line_color='red', source=SRC_dspec_prof2)
SRC_dspecfit_prof1 = ColumnDataSource({'x': [], 'y': []})
r_dspecfit_prof1 = p_dspec_lines.line(x='x', y='y', alpha=0.7, line_width=1, line_color='orange',
                                      source=SRC_dspecfit_prof1)
SRC_dspecfit_prof2 = ColumnDataSource({'x': [], 'y': []})
r_dspecfit_prof2 = p_dspec_lines.line(x='x', y='y', alpha=0.7, line_width=1, line_color='red',
                                      source=SRC_dspecfit_prof2)
tooltips = [("freqa, freqv", "@freqa, @freqv"), ("max, lag [ms]", "@ccmax, @ccpeak"), ]

hover_JScode = """
    var CCfdata = CC_SQR.data;
    var spec = specplt.get('data').data[0];
    var specfit = specplt.get('data').datafit[0];
    var dtim = time.get('data').data[0];
    var dtimfit = timefit.get('data').data[0];
    var t0 = %s;
    var t1 = %s;
    var nx = %d;
    var f0 = %s;
    var f1 = %s;
    var nxfit = %d;
    var indices = cb_data.index['1d'].indices;

    var data = {'x': [], 'y': []};
    data['x'].push(t0);
    data['x'].push(t1);
    for (i=0; i < 2; i++) {
        data['y'].push(CCfdata.freqa[indices[0]]);
    }
    r_freq1.set('data',data);

    data = {'x': [], 'y': []};
    for (i=0; i < 2; i++) {
        data['y'].push(CCfdata.freqv[indices[0]]);
    }
    data['x'].push(t0);
    data['x'].push(t1);
    r_freq2.set('data',data);

    var data = {'x': [], 'y': []};
    data['y'].push(f0);
    data['y'].push(f1);
    for (i=0; i < 2; i++) {
        data['x'].push(CCfdata.freqa[indices[0]]);
    }
    r_CCmaxl1.set('data',data);

    var data = {'x': [], 'y': []};
    data['x'].push(f0);
    data['x'].push(f1);
    for (i=0; i < 2; i++) {
        data['y'].push(CCfdata.freqv[indices[0]]);
    }
    r_CCmaxl2.set('data',data);

    data = {'x': [], 'y': []};
    for (i=0; i < nx; i++) {
        data['x'].push(dtim[i]);
        data['y'].push(spec[ CCfdata.fidxa[indices[0]]*nx+i]);
    }
    r_prof1.set('data',data);

    data = {'x': [], 'y': []};
    for (i=0; i < nx; i++) {
        data['x'].push(dtim[i]);
        data['y'].push(spec[CCfdata.fidxv[indices[0]]*nx+i]);
    }
    r_prof2.set('data',data);

    data = {'x': [], 'y': []};
    for (i=0; i < nxfit; i++) {
        data['x'].push(dtimfit[i]);
        data['y'].push(specfit[CCfdata.fidxa[indices[0]]*nxfit+i]);
    }
    r_fit_prof1.set('data',data);

    data = {'x': [], 'y': []};
    for (i=0; i < nxfit; i++) {
        data['x'].push(dtimfit[i]);
        data['y'].push(specfit[CCfdata.fidxv[indices[0]]*nxfit+i]);
    }
    r_fit_prof2.set('data',data);
    """ % (dtim[0], dtim[-1], ntim, freq[0], freq[-1], ntimfit)

CJSargs = {'CC_SQR': SRC_CCmax_square, 'r_freq1': r_freq_line1.data_source, 'r_freq2': r_freq_line2.data_source,
           'specplt': ColumnDataSource({'data': [spec.ravel()], 'datafit': [specfit.ravel()]}),
           'time': ColumnDataSource({'data': [dtim]}),
           'timefit': ColumnDataSource({'data': [dtimfit]}), 'r_prof1': r_dspec_prof1.data_source,
           'r_prof2': r_dspec_prof2.data_source, 'r_fit_prof1': r_dspecfit_prof1.data_source,
           'r_fit_prof2': r_dspecfit_prof2.data_source, 'r_CCmaxl1': r_CCmax_line1.data_source,
           'r_CCmaxl2': r_CCmax_line2.data_source}

p_CCmax_hover_callback = CustomJS(args=CJSargs, code=hover_JScode)
p_CCmax_hover = HoverTool(tooltips=tooltips, callback=p_CCmax_hover_callback,
                          renderers=[r_CCmax_square])
p_CCmax.add_tools(p_CCmax_hover)

p_CCpeak_hover_callback = CustomJS(args=CJSargs, code=hover_JScode)
p_CCpeak_hover = HoverTool(tooltips=tooltips, callback=p_CCpeak_hover_callback,
                           renderers=[r_CCpeak_square])
p_CCpeak.add_tools(p_CCpeak_hover)


def Slider_watershed_update(attrname, old, new):
    wtshdpercent = Slider_watershed.value / 100.0
    specfit = CC_savedata['specfit'].copy()
    spec = CC_savedata['spec'].copy()
    for fidx in xrange(nfreq):
        specfitslice = specfit[fidx, :]
        slicemax, slicemin = specfitslice.max(), specfitslice.min()
        thrshd = slicemin + (slicemax - slicemin) * wtshdpercent
        specfitslice[specfitslice < thrshd] = np.nan
        specfit[fidx, :] = specfitslice
    r_dspecfit.data_source.data['image'] = [specfit]
    # CC_dict = DButil.XrsCorrMap(specfit, timfit, freq, doxscale=False)
    # r_CCmax.data_source.data['image'] = [CC_dict['ccmax']]
    # r_CCpeak.data_source.data['image'] = [CC_dict['ccpeak'] * dtfit * 1000.0]



Slider_watershed = Slider(start=0, end=100, value=0, step=5, title="watershed",
                          width=config_plot['plot_config']['tab_XrsCorr']['widgetbox_wdth'],
                          callback_policy='mouseup')

# This data source is just used to communicate / trigger the real callback
Slider_watershed.on_change('value', Slider_watershed_update)


def Slider_threshold_update(attrname, old, new):
    thrshdpercent = Slider_threshold.value / 100.0
    specfit = CC_savedata['specfit'].copy()
    fidxflag = np.log(specfit.max(axis=1)) < thrshdpercent * np.log(specfit.max())
    ccmaxplt = ccmax
    ccpeakplt = ccpeak * dtfit * 1000.0
    for fflag in freq[fidxflag]:
        ccmaxplt[freqa == fflag] = np.nan
        ccmaxplt[freqv == fflag] = np.nan
        ccpeakplt[freqa == fflag] = np.nan
        ccpeakplt[freqv == fflag] = np.nan
    maspecfit = np.tile(fidxflag, ntimfit).reshape(ntimfit, nfreq).swapaxes(0, 1)
    specfit[maspecfit] = np.nan
    r_dspecfit.data_source.data['image'] = [specfit]
    r_CCmax.data_source.data['image'] = [ccmaxplt]
    r_CCpeak.data_source.data['image'] = [ccpeakplt]


Slider_threshold = Slider(start=0, end=100, value=0, step=1, title="threshold",
                          width=config_plot['plot_config']['tab_XrsCorr']['widgetbox_wdth'],
                          callback_policy='mouseup')

Slider_threshold.on_change('value', Slider_threshold_update)

BUT_exit = Button(label='Exit FSview',
                  width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                  button_type='danger')
BUT_exit.on_click(exit)

lout = column(row(gridplot([[p_dspec, p_dspecfit]], toolbar_location='right'), p_dspec_lines),
              row(gridplot([[p_CCmax, p_CCpeak]], toolbar_location='right'),
                  column(Slider_watershed, Slider_threshold, BUT_exit)))

curdoc().add_root(lout)
curdoc().title = "XrsCorr"

# except:
#     print 'Error: No CC_save found!!!'
#     raise SystemExit
