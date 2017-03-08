import astropy.units as u
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import glob
import json
import os
import bokeh.palettes as bp
from astropy.time import Time
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (ColumnDataSource, Slider, Button, TextInput, CheckboxGroup, RadioGroup,
                          BoxSelectTool, TapTool, Div, Spacer, Range1d)
from bokeh.models.mappers import LogColorMapper, LinearColorMapper
from bokeh.plotting import figure, curdoc
from suncasa.utils import DButil
from suncasa.utils.puffin import getAIApalette

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"


def exit_update():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
'''load config file'''
with open(suncasa_dir + 'aiaBrowser/config.json', 'r') as fp:
    config_plot = json.load(fp)

database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/aiaBrowserData/'
if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))
SDO_dir = database_dir + 'Download/'
if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))
infile = database_dir + 'MkPlot_args.json'
try:
    with open(infile, 'rb') as fp:
        MkPlot_args_dict = json.load(fp)
except:
    print 'Error: No MkPlot_args.json found!!!'
    raise SystemExit

PlotID = MkPlot_args_dict['PlotID']
try:
    stackplt_save = database_dir + 'stackplt-' + PlotID + '.npy'
    stackpltdict = np.load(stackplt_save).item()
except:
    print 'Error: No stackplt-xxxxxxxxxx.npz found!!!'
    raise SystemExit

intens = stackpltdict['intens']
tim = stackpltdict['time']
dist = stackpltdict['dist']
wavelngth = stackpltdict['wavelength']
dtim = (tim - tim[0]) * 24. * 3600
dt = np.median(np.diff(dtim))
ntim = len(tim)
ndist = len(dist)

TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save"
p_stackplt = figure(tools=TOOLS,
                    plot_width=config_plot['plot_config']['tab_StackPlt']['stackplt_wdth'],
                    plot_height=config_plot['plot_config']['tab_StackPlt']['stackplt_hght'],
                    x_range=(dtim[0], dtim[-1]), y_range=(dist[0], dist[-1]),
                    toolbar_location="above")
tim0_char = Time(tim[0], format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
p_stackplt.axis.visible = True
p_stackplt.title.text = "Stack plot"
p_stackplt.xaxis.axis_label = 'Seconds since ' + tim0_char
p_stackplt.yaxis.axis_label = 'distance [arcsec]'
p_stackplt.border_fill_alpha = 0.4
p_stackplt.axis.major_tick_out = 0
p_stackplt.axis.major_tick_in = 5
p_stackplt.axis.minor_tick_out = 0
p_stackplt.axis.minor_tick_in = 3
p_stackplt.axis.major_tick_line_color = "white"
p_stackplt.axis.minor_tick_line_color = "white"
if stackpltdict['observatory'] == 'SDO' and stackpltdict['instrument'] == 'AIA':
    palette = getAIApalette(wavelngth)
    clrange = DButil.sdo_aia_scale_dict(wavelength=wavelngth)
    colormapper = LogColorMapper(palette=palette, low=clrange['low'], high=clrange['high'])
else:
    palette = bp.viridis(256)
    colormapper = LinearColorMapper(palette=palette)

r_stackplt = p_stackplt.image(image=[intens], x=dtim[0], y=dist[0], dw=dtim[-1] - dtim[0], dh=dist[-1] - dist[0],
                              color_mapper=colormapper)
Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b>
            first before closing the tab</p></b>""", width=config_plot['plot_config']['tab_MkPlot']['button_wdth'])
BUT_exit = Button(label='Exit', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit_update)

lout = row(p_stackplt, column(BUT_exit, Div_info))

# def timeout_callback():
#     print 'timeout'
#     raise SystemExit


curdoc().add_root(lout)
curdoc().title = "StackPlt"
