import json
import os, sys
import pickle
import time
from collections import OrderedDict
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.misc import bytescale
import bokeh.palettes as bp
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
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]


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
fidx1 = CC_savedata['fidx1']
fidx2 = CC_savedata['fidx2']

CCmaxDF = pd.DataFrame({'fidx1': fidx1, 'fidx2': fidx2})
TOOLS = "crosshair,pan,wheel_zoom,tap,box_zoom,reset,save"
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
# SRC_dspec_image = ColumnDataSource(
#     data={'data': [spec_plt], 'xx': [dtim], 'yy': [freq]})
r_dspec = p_dspec.image(image=[spec], x=dtim[0], y=freq[0],
                        dw=dtim[-1] - dtim[0],
                        dh=freq[-1] - freq[0], palette=bokehpalette_jet)

p_CCmax = figure(tools=TOOLS,
                 plot_width=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                 plot_height=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                 x_range=(freq[0], freq[-1]), y_range=(freq[0], freq[-1]))
p_CCmax.axis.visible = True
p_CCmax.title.text = "Xross Correlation maximum"
p_CCmax.xaxis.axis_label = 'Frequency [GHz]'
p_CCmax.yaxis.axis_label = 'Frequency [GHz]'
r_CCmax = p_CCmax.image(image=[ccmax], x=freq[0], y=freq[0],
                        dw=freq[-1] - freq[0],
                        dh=freq[-1] - freq[0], palette=bp.RdBu[11])

p_CCpeak = figure(tools=TOOLS,
                  plot_width=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                  plot_height=config_plot['plot_config']['tab_XrsCorr']['dspec_hght'],
                  x_range=p_CCmax.x_range, y_range=p_CCmax.y_range)
p_CCpeak.axis.visible = True
p_CCpeak.title.text = "Xross Correlation lag"
p_CCpeak.xaxis.axis_label = 'Frequency [GHz]'
p_CCpeak.yaxis.axis_label = 'Frequency [GHz]'
r_CCpeak = p_CCpeak.image(image=[ccpeak], x=freq[0], y=freq[0],
                          dw=freq[-1] - freq[0],
                          dh=freq[-1] - freq[0], palette=bp.RdBu[11])

BUT_exit = Button(label='Exit FSview',
                  width=config_plot['plot_config']['tab_FSview_base']['widgetbox_wdth'],
                  button_type='danger')
BUT_exit.on_click(exit)

lout = column(row(p_dspec, BUT_exit), gridplot([[p_CCmax, p_CCpeak]], toolbar_location='right'))

curdoc().add_root(lout)
curdoc().title = "XrsCorr"
# except:
#     print 'Error: No CC_save found!!!'
#     raise SystemExit
