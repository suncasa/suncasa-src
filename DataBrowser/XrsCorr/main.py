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
from bokeh.models import LinearColorMapper, ColorBar
from bokeh.plotting import figure, curdoc
from astropy.time import Time

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

CCmaxDF = pd.DataFrame(
    {'ccmax': ccmax.ravel(), 'ccpeak': ccpeak.ravel(), 'freqa': freqa.ravel(), 'freqv': freqv.ravel(),
     'fidxa': fidxa.ravel(), 'fidxv': fidxv.ravel()})
TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
p_dspec = figure(tools=TOOLS,
                 plot_width=config_plot['plot_config']['tab_XrsCorr']['dspec_wdth'],
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
                        dw=dtim[-1] - dtim[0],
                        dh=freq[-1] - freq[0], palette=bokehpalette_jet)

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
r_CCmax = p_CCmax.image(image=[ccmax], x=freq[0], y=freq[0],
                        dw=freq[-1] - freq[0],
                        dh=freq[-1] - freq[0], palette=bokehpalette_Blues)
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


cm_CCmax = LinearColorMapper(palette=bokehpalette_Blues, low=0, high=1.0)
cb_CCmax = ColorBar(color_mapper=cm_CCmax, label_standoff=5, width=5, border_line_color=None, location=(0, 0))
p_CCmax.add_layout(cb_CCmax, 'right')

p_CCpeak = figure(tools=TOOLS,
                  plot_width=config_plot['plot_config']['tab_XrsCorr']['CCmap_wdth'],
                  plot_height=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                  x_range=p_CCmax.x_range, y_range=p_CCmax.y_range)
p_CCpeak.axis.visible = True
p_CCpeak.title.text = "Xross Correlation lag"
p_CCpeak.xaxis.axis_label = 'Frequency [GHz]'
p_CCpeak.yaxis.axis_label = 'Frequency [GHz]'
r_CCpeak = p_CCpeak.image(image=[ccpeak], x=freq[0], y=freq[0],
                          dw=freq[-1] - freq[0],
                          dh=freq[-1] - freq[0], palette=bokehpalette_RdBu)

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

lagmax = ntim / 2
cm_CCpeak = LinearColorMapper(palette=bokehpalette_RdBu, low=-lagmax, high=lagmax)
cb_CCpeak = ColorBar(color_mapper=cm_CCpeak, label_standoff=5, width=5, border_line_color=None, location=(0, 0))
p_CCpeak.add_layout(cb_CCpeak, 'right')

SRC_freq_line1 = ColumnDataSource({'x': [], 'y': []})
r_freq_line1 = p_dspec.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange', source=SRC_freq_line1)
SRC_freq_line2 = ColumnDataSource({'x': [], 'y': []})
r_freq_line2 = p_dspec.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red', source=SRC_freq_line2)

TOOLS = "crosshair,save"
p_dspec_lines = figure(tools=TOOLS,
                       plot_width=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                       plot_height=config_plot['plot_config']['tab_XrsCorr']['CCmap_hght'],
                       x_range=(dtim[0], dtim[-1]), y_range=(np.amin(spec), np.amax(spec)),
                       toolbar_location="above")
p_dspec_lines.axis.visible = True
p_dspec_lines.title.text = "Lgiht curves"
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
r_dspec_prof1 = p_dspec_lines.line(x='x', y='y', alpha=0.6, line_width=2, line_color='orange', source=SRC_dspec_prof1)
SRC_dspec_prof2 = ColumnDataSource({'x': [], 'y': []})
r_dspec_prof2 = p_dspec_lines.line(x='x', y='y', alpha=0.6, line_width=2, line_color='red',
                                   source=SRC_dspec_prof2)
tooltips = [("(freqa,freqv)", "(@freqa, @freqv)"), ("max, lag", "(@ccmax,@ccpeak)"), ]
# tooltips = """
# <div>
#     <div>
#         <img
#             src="@imgs" height="42" alt="@imgs" width="42"
#             style="float: left; margin: 0px 15px 15px 0px;"
#             border="2"
#         ></img>
#     </div>
#     <div>
#         <span style="font-size: 17px; font-weight: bold;">@desc</span>
#         <span style="font-size: 15px; color: #966;">[$index]</span>
#     </div>
#     <div>
#         <span>@fonts{safe}</span>
#     </div>
#     <div>
#         <span style="font-size: 15px;">Location</span>
#         <span style="font-size: 10px; color: #696;">($x, $y)</span>
#     </div>
# </div>
# """

hover_JScode = """
    var CCfdata = CC_SQR.data;
    var spec = specplt.get('data').data[0];
    var dtim = time.get('data').data[0];
    var t0 = %s;
    var t1 = %s;
    var nx = %d;
    var f0 = %s;
    var f1 = %s;
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
    """ % (dtim[0], dtim[-1], ntim, freq[0],freq[-1])

CJSargs = {'CC_SQR': SRC_CCmax_square, 'r_freq1': r_freq_line1.data_source, 'r_freq2': r_freq_line2.data_source,
           'specplt': ColumnDataSource({'data': [spec.ravel()]}), 'time': ColumnDataSource({'data': [dtim]}),
           'r_prof1': r_dspec_prof1.data_source,
           'r_prof2': r_dspec_prof2.data_source,
           'r_CCmaxl1':r_CCmax_line1.data_source,
           'r_CCmaxl2':r_CCmax_line2.data_source}

p_CCmax_hover_callback = CustomJS(args=CJSargs, code=hover_JScode)
p_CCmax_hover = HoverTool(tooltips=tooltips, callback=p_CCmax_hover_callback,
                          renderers=[r_CCmax_square])
p_CCmax.add_tools(p_CCmax_hover)

p_CCpeak_hover_callback = CustomJS(args=CJSargs, code=hover_JScode)
p_CCpeak_hover = HoverTool(tooltips=tooltips, callback=p_CCpeak_hover_callback,
                           renderers=[r_CCpeak_square])
p_CCpeak.add_tools(p_CCpeak_hover)

BUT_exit = Button(label='Exit FSview',
                  width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                  button_type='danger')
BUT_exit.on_click(exit)

lout = column(row(p_dspec, p_dspec_lines), row(gridplot([[p_CCmax, p_CCpeak]], toolbar_location='right'), BUT_exit))

curdoc().add_root(lout)
curdoc().title = "XrsCorr"

# except:
#     print 'Error: No CC_save found!!!'
#     raise SystemExit
