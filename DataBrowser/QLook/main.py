# The plot server must be running
# Go to http://localhost:5006/bokeh to view this plot
import json
import os
import os.path
import time
from sys import platform
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import scipy.ndimage as sn
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (ColumnDataSource, Button, TextInput, DataTable, TableColumn, BoxSelectTool, TapTool,
                          HoverTool, Spacer, Div)
from bokeh.models.widgets import Select
from bokeh.plotting import figure, curdoc
from astropy.time import Time


__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

'''Prepare the ports for FSview'''

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
'''load config file'''
with open(suncasa_dir + 'DataBrowser/config.json', 'r') as fp:
    config_plot = json.load(fp)
database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir)+'/'
with open('{}config_EvtID_curr.json'.format(database_dir), 'r') as fp:
    config_EvtID = json.load(fp)


'''define the colormaps'''
colormap_jet = cm.get_cmap("jet")  # choose any matplotlib colormap here
bokehpalette_jet = [colors.rgb2hex(m) for m in colormap_jet(np.arange(colormap_jet.N))]

'''
-------------------------- panel 1 --------------------------
'''

start_timestamp = time.time()
database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir)+'/'
event_id = config_EvtID['datadir']['event_id']
specfile = database_dir + event_id + config_EvtID['datadir']['event_specfile']


def load_specdata(specfile=None):
    global tab1_specdata, tab1_spec, tab1_tim, tab1_freq, tab1_dtim, tab1_spec_plt, tab1_bl
    tab1_specdata = np.load(specfile)
    tab1_bl = tab1_specdata['bl'].item().split(';')
    tab1_pol = 'I'
    bl_index = 0
    tab1_spec = tab1_specdata['spec'][:, :, :, :]
    tab1_tim = tab1_specdata['tim'][:]
    tab1_freq = tab1_specdata['freq'] / 1e9
    tab1_spec_sz = tab1_spec.shape
    spec_sz2, spec_sz1 = 10, 10
    if tab1_spec_sz[3] > 1750:
        spec_sz2 = next(i for i in xrange(1, 10) if i / 10. * tab1_spec_sz[3] > 1750)
    if tab1_spec_sz[2] > 250:
        spec_sz1 = next(i for i in xrange(1, 10) if i / 10. * tab1_spec_sz[2] > 250)
    tab1_spec = sn.interpolation.zoom(tab1_spec, [1, 1, spec_sz1 / 10.0, spec_sz2 / 10.0], order=1)
    tab1_tim = sn.interpolation.zoom(tab1_tim, spec_sz2 / 10.0, order=1)
    tab1_freq = sn.interpolation.zoom(tab1_freq, spec_sz1 / 10.0, order=1)

    if tab1_pol == 'RR':
        tab1_spec_plt = tab1_spec[0, bl_index, :, :]
    elif tab1_pol == 'LL':
        tab1_spec_plt = tab1_spec[1, bl_index, :, :]
    elif tab1_pol == 'I':
        tab1_spec_plt = (tab1_spec[0, bl_index, :, :] + tab1_spec[1, bl_index, :, :]) / 2.
    elif tab1_pol == 'V':
        tab1_spec_plt = (tab1_spec[0, bl_index, :, :] - tab1_spec[1, bl_index, :, :]) / 2.
    tab1_dtim = tab1_tim - tab1_tim[0]


load_specdata(specfile)
TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
'''create the dynamic spectrum plot'''
tab1_p_dspec = figure(tools=TOOLS, webgl=config_plot['plot_config']['WebGL'],
                      plot_width=config_plot['plot_config']['tab_QLook']['dspec_wdth'],
                      plot_height=config_plot['plot_config']['tab_QLook']['dspec_hght'],
                      x_range=(tab1_dtim[0], tab1_dtim[-1]),
                      y_range=(tab1_freq[0], tab1_freq[-1]), toolbar_location="above")
tim0_char = Time((tab1_tim[0] / 3600. / 24. + 2400000.5) * 86400. / 3600. / 24., format='jd', scale='utc',
                 precision=3).iso
tab1_p_dspec.axis.visible = True
tab1_p_dspec.title.text = "Dynamic spectrum"
tab1_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
tab1_p_dspec.yaxis.axis_label = 'Frequency [GHz]'
tab1_p_dspec.border_fill_color = "silver"
tab1_p_dspec.border_fill_alpha = 0.4
tab1_p_dspec.axis.major_tick_out = 0
tab1_p_dspec.axis.major_tick_in = 5
tab1_p_dspec.axis.minor_tick_out = 0
tab1_p_dspec.axis.minor_tick_in = 3
tab1_p_dspec.axis.major_tick_line_color = "white"
tab1_p_dspec.axis.minor_tick_line_color = "white"

tab1_SRC_dspec = ColumnDataSource(data={'data': [tab1_spec_plt], 'xx': [tab1_dtim], 'yy': [tab1_freq]})
tab1_r_dspec = tab1_p_dspec.image(image="data", x=tab1_dtim[0], y=tab1_freq[0], dw=tab1_dtim[-1] - tab1_dtim[0],
                                  dh=tab1_freq[-1] - tab1_freq[0], source=tab1_SRC_dspec, palette=bokehpalette_jet)

tab1_spec_sz = tab1_spec.shape
ratio_spec_sz2, ratio_spec_sz1 = 10, 10
if tab1_spec_sz[3] > 520:
    ratio_spec_sz2 = next(i for i in xrange(1, 10) if i / 10. * tab1_spec_sz[3] > 520)
if tab1_spec_sz[2] > 75:
    ratio_spec_sz1 = next(i for i in xrange(1, 10) if i / 10. * tab1_spec_sz[2] > 75)
tab1_tim_square = sn.interpolation.zoom(tab1_tim, ratio_spec_sz2 / 10., order=1)
tab1_freq_square = sn.interpolation.zoom(tab1_freq, ratio_spec_sz1 / 10., order=1)
tab1_nfreq_square, = tab1_freq_square.shape
tab1_ntim_square, = tab1_tim_square.shape
tim_map = ((np.tile(tab1_tim_square, tab1_nfreq_square).reshape(tab1_nfreq_square,
                                                                tab1_ntim_square) / 3600. / 24. + 2400000.5)) * 86400.
freq_map = np.tile(tab1_freq_square, tab1_ntim_square).reshape(tab1_ntim_square, tab1_nfreq_square).swapaxes(0, 1)
xx = tim_map.flatten()
yy = freq_map.flatten()
tab1_dspecDF_square = pd.DataFrame({'time': xx - xx[0], 'freq': yy})
tab1_SRC_dspec_square = ColumnDataSource(tab1_dspecDF_square)
tab1_r_square = tab1_p_dspec.square('time', 'freq', source=tab1_SRC_dspec_square, fill_color=None, fill_alpha=0.0,
                                    line_color=None, line_alpha=0.0, selection_fill_color=None,
                                    selection_fill_alpha=0.0, nonselection_fill_alpha=0.0,
                                    selection_line_alpha=0.0, nonselection_line_alpha=0.0)
tab2_SRC_dspec_Patch = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
tab2_r_dspec_patch = tab1_p_dspec.patch('xx', 'yy', source=tab2_SRC_dspec_Patch,
                                        fill_color=None, fill_alpha=0.5, line_color="Magenta",
                                        line_alpha=1.0, line_width=1)

## ----------------baseline & polarization selection------------------------
tab1_Select_bl = Select(title="Baseline:", value=tab1_bl[0], options=tab1_bl, width=150)
tab1_Select_pol = Select(title="Polarization:", value="I", options=["RR", "LL", "I", "V"], width=150)
tab1_Select_colorspace = Select(title="ColorSpace:", value="linear", options=["linear", "log"], width=150)


def tab1_update_dspec(attrname, old, new):
    global tab1_spec, tab1_dtim, tab1_freq, tab1_bl
    select_pol = tab1_Select_pol.value
    select_bl = tab1_Select_bl.value
    bl_index = tab1_bl.index(select_bl)
    if select_pol == 'RR':
        spec_plt = tab1_spec[0, bl_index, :, :]
    elif select_pol == 'LL':
        spec_plt = tab1_spec[1, bl_index, :, :]
    elif select_pol == 'I':
        spec_plt = (tab1_spec[0, bl_index, :, :] + tab1_spec[1, bl_index, :, :]) / 2.
    elif select_pol == 'V':
        spec_plt = (tab1_spec[0, bl_index, :, :] - tab1_spec[1, bl_index, :, :]) / 2.
        tab1_Select_colorspace.value = 'linear'
    if tab1_Select_colorspace.value == 'log' and select_pol != 'V':
        spec_plt = np.log(spec_plt)
    tab1_SRC_dspec.data = {'data': [spec_plt], 'xx': [tab1_dtim], 'yy': [tab1_freq]}


tab1_ctrls = [tab1_Select_bl, tab1_Select_pol, tab1_Select_colorspace]
for ctrl in tab1_ctrls:
    ctrl.on_change('value', tab1_update_dspec)
try:
    os.system('cp {}StrID_list.json {}StrID_list_tmp.json'.format(database_dir + event_id, database_dir + event_id))
    StrIDList = pd.read_json(database_dir + event_id + 'StrID_list_tmp.json')
    StrIDList = StrIDList.sort_values(by='timeran', ascending=1)
    StrIDList['time'] = [ll - tab1_tim[0] for ll in StrIDList['time']]
except:
    StrIDList = pd.DataFrame({'date': [], 'freq': [], 'freqran': [], 'str_id': [], 'time': [], 'timeran': []})
    StrIDList.to_json(database_dir + event_id + 'StrID_list.json')
    os.system('cp {}StrID_list.json {}StrID_list_tmp.json'.format(database_dir + event_id, database_dir + event_id))

tab1_SRC_StrIDPatch = ColumnDataSource(StrIDList)
tab1_render_patch = tab1_p_dspec.patches('time', 'freq', source=tab1_SRC_StrIDPatch, hover_fill_color="OrangeRed",
                                         hover_line_color='OrangeRed', hover_fill_alpha=0.2, hover_line_alpha=1.0,
                                         selection_fill_color=None,
                                         selection_fill_alpha=0.0, selection_line_color="OrangeRed",
                                         selection_line_alpha=1, nonselection_fill_color=None,
                                         nonselection_fill_alpha=0.0,
                                         nonselection_line_color="white", nonselection_line_alpha=0.6,
                                         line_width=1, fill_color=None, fill_alpha=0.0,
                                         line_color="white")

tooltips = [("StrID", "@str_id"), ("Date", "@date"), ("TimeRange", "@timeran"), ("FreqRange", "@freqran")]
tab1_HoverTool_dspec = HoverTool(tooltips=tooltips, anchor='top_left', renderers=[tab1_render_patch])
tab1_p_dspec.add_tools(tab1_HoverTool_dspec)
tab1_p_dspec.add_tools(BoxSelectTool())
tab1_p_dspec.add_tools(TapTool(renderers=[tab1_render_patch]))
tab1_p_dspec.select(BoxSelectTool).select_every_mousemove = False
tab1_p_dspec.axis.major_tick_out = 0
tab1_p_dspec.axis.major_tick_in = 5
tab1_p_dspec.axis.minor_tick_out = 0
tab1_p_dspec.axis.minor_tick_in = 3
tab1_p_dspec.axis.major_tick_line_color = "white"
tab1_p_dspec.axis.minor_tick_line_color = "white"

tab1_TbCols = [TableColumn(field="str_id", title="StrID"), TableColumn(field="timeran", title="Time Range"),
               TableColumn(field="freqran", title="Freq Range"), ]
tab1_DataTb_dspec = DataTable(source=tab1_render_patch.data_source, columns=tab1_TbCols,
                              width=config_plot['plot_config']['tab_QLook']['StrID_DataTb_wdth'],
                              height=config_plot['plot_config']['tab_QLook']['StrID_DataTb_hght'])  # , editable=True)

tab1_Div_Tb = Div(text=""" """, width=config_plot['plot_config']['tab_QLook']['StrID_DataTb_wdth'])
tab1_Div_exit = Div(text="""
<p><b>Warning</b>: 1. Click the <b>Exit QLook</b> first before closing the tab</p>
<p><b>Warning</b>: 2. <b>FSview</b> or <b>FSview2CASA</b> tabs will disconnect if <b>Exit QLook is clicked</b></p>""",
                    width=150)

tab1_selected_dspec_square = None


def tab1_SRC_dspec_square_select(attrname, old, new):
    global tab1_selected_dspec_square
    tab1_selected_dspec_square = tab1_SRC_dspec_square.selected['1d']['indices']
    if tab1_selected_dspec_square:
        global dftmp
        dftmp = tab1_dspecDF_square.iloc[tab1_selected_dspec_square, :]
        x0, x1 = dftmp['time'].min(), dftmp['time'].max()
        y0, y1 = dftmp['freq'].min(), dftmp['freq'].max()
        tab2_r_dspec_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
    else:
        tab2_r_dspec_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [], 'yy': []})).data


tab1_SRC_dspec_square.on_change('selected', tab1_SRC_dspec_square_select)

tab1_input_StrID = TextInput(value="Type in here", title="New StrID:",
                             width=config_plot['plot_config']['tab_QLook']['StrID_DataTb_BUT_wdth'])
timestart = xx[0]


def tab1_update_addStrID():
    global dftmp
    if tab1_selected_dspec_square:
        time0, time1 = dftmp['time'].min() + timestart, dftmp['time'].max() + timestart
        freq0, freq1 = dftmp['freq'].min(), dftmp['freq'].max()
        date_char = Time(timestart / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date').iso
        t0_char = Time(time0 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        t0_char = t0_char.split(' ')[1]
        t1_char = Time(time1 / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
        t1_char = t1_char.split(' ')[1]
        time0, time1 = (time0 / 86400. - 2400000.5) * 24. * 3600., (time1 / 86400. - 2400000.5) * 24. * 3600.
        StrID = pd.DataFrame(dict(time=[[time0, time1, time1, time0]], freq=[[freq0, freq0, freq1, freq1]],
                                  str_id=[[tab1_input_StrID.value]], date=[[date_char]],
                                  timeran=[[t0_char + '~' + t1_char]],
                                  freqran=[["{:.3f}~{:.3f} GHz".format(freq0, freq1)]]))
        StrIDList = pd.read_json(database_dir + event_id + 'StrID_list_tmp.json')
        StrIDList = pd.concat([StrIDList, StrID])
        StrIDList = StrIDList.sort_values(by='timeran', ascending=1)
        StrIDList.index = range(StrIDList.index.size)
        StrIDList.to_json(database_dir + event_id + 'StrID_list_tmp.json')
        StrIDList['time'] = [ll - tab1_tim_square[0] for ll in StrIDList['time']]
        tab1_render_patch.data_source.data = ColumnDataSource(StrIDList).data
        tab1_Div_Tb.text = """<p>added <b>""" + tab1_input_StrID.value + """</b>  to the list</p>"""
    else:
        tab1_Div_Tb.text = """<p><b>Warning: No time & freq range selected.
                                Use the box select tool to select time & freq
                                range in the dynamic spectrum first!!!</b></p>"""


tab1_selected_StrID_entry = None


def tab1_selection_StrID_entry(attrname, old, new):
    global tab1_selected_StrID_entry
    tab1_selected_StrID_entry = tab1_SRC_StrIDPatch.selected['1d']['indices']


tab1_SRC_StrIDPatch.on_change('selected', tab1_selection_StrID_entry)


def tab1_update_deleteStrID():
    global tab1_tim_square, tab1_selected_StrID_entry, database_dir, event_id
    if tab1_selected_StrID_entry:
        StrIDList = pd.read_json(database_dir + event_id + 'StrID_list_tmp.json')
        StrIDList = StrIDList.sort_values(by='timeran', ascending=1)
        StrIDList = StrIDList.drop(StrIDList.index[tab1_selected_StrID_entry])
        StrIDList.index = range(StrIDList.index.size)
        StrIDList.to_json(database_dir + event_id + 'StrID_list_tmp.json')
        StrIDList['time'] = [ll - tab1_tim_square[0] for ll in StrIDList['time']]
        tab1_render_patch.data_source.data = ColumnDataSource(StrIDList).data
        tab1_Div_Tb.text = """<p>removed <b>""" + StrIDList.iloc[tab1_selected_StrID_entry[0]]['str_id'][
            0] + """</b> from the list</p>"""
        tab1_selected_StrID_entry = None
    else:
        tab1_Div_Tb.text = """<p><b>Warning: No StrID selected. Select one StrID first!!!</b></p>"""


port = 5100
ports = []
ports.append(port)


def tab1_update_FSviewStrID():
    global dftmp
    global port, ports
    if tab1_selected_StrID_entry:
        StrIDList = pd.read_json(database_dir + event_id + 'StrID_list_tmp.json')
        StrIDList = StrIDList.sort_values(by='timeran', ascending=1)
        StrID = StrIDList.iloc[tab1_selected_StrID_entry[0]]
        out_json = database_dir + event_id + StrID['str_id'][0] + '.json'
        StrID.to_json(out_json)
        out_json = database_dir + event_id + 'CurrFS.json'
        struct_id = StrID['str_id'][0] + '/'
        FS_config = {'datadir': {'event_id': event_id, 'struct_id': struct_id,
                                 'FS_specfile': database_dir + event_id + struct_id + StrID['str_id'][0] + '_' +
                                                StrID['date'][0] + 'T' + str(StrID['timeran'][0]).translate(None,
                                                                                                            ':') + '.spec.npz'}}
        with open(out_json, 'w') as fp:
            json.dump(FS_config, fp)
        in_json = database_dir + event_id + 'CurrFS.json'
        with open(in_json, 'r') as fp:
            FS_config = json.load(fp)
        FS_specfile = FS_config['datadir']['FS_specfile']
        FS_dspecDF = database_dir + event_id + struct_id + 'dspecDF-save'
        time0, time1 = StrID['time'][0], StrID['time'][1]
        freq0, freq1 = StrID['freq'][0], StrID['freq'][-1]
        bl = tab1_specdata['bl']
        spec = tab1_specdata['spec']
        npol = tab1_specdata['npol']
        nbl = tab1_specdata['nbl']
        ntim = tab1_specdata['ntim']
        nfreq = tab1_specdata['nfreq']
        tim = tab1_specdata['tim'][:]
        freq = tab1_specdata['freq'] / 1e9
        timeidx0 = next(i for i in xrange(ntim) if tim[i] >= time0)
        timeidx1 = next(i for i in xrange(ntim - 1, -1, -1) if tim[i] <= time1) + 1
        freqidx0 = next(i for i in xrange(nfreq) if freq[i] >= freq0)
        freqidx1 = next(i for i in xrange(nfreq - 1, -1, -1) if freq[i] <= freq1) + 1
        spec = spec[:, :, freqidx0:(freqidx1 + 1), timeidx0:(timeidx1 + 1)]
        tim = tim[timeidx0:(timeidx1 + 1)]
        freq = freq[freqidx0:(freqidx1 + 1)] * 1.0e9
        ntim = len(tim)
        nfreq = len(freq)
        if not os.path.exists(FS_specfile):
            struct_id_dir = database_dir + event_id + struct_id
            if not os.path.exists(struct_id_dir):
                os.system('mkdir {}'.format(struct_id_dir))
            np.savez(FS_specfile, spec=spec, tim=tim, freq=freq, bl=bl, npol=npol, nbl=nbl, nfreq=nfreq, ntim=ntim)
            # tab1_Div_Tb.text = """<p><b>""" + FS_specfile + """</b> saved >>>>>> Click the <b>FS veiw button</b> again to make aperture synthesis images</p>"""

        print 'bokeh serve FSview --show --port {} &'.format(port)
        os.system('cd .. & bokeh serve {}DataBrowser/FSview --show --port {} &'.format(suncasa_dir,port))
        port += 1
        ports.append(port)
        if os.path.exists(FS_dspecDF):
            tab1_Div_Tb.text = """<p>sent StrID to <b>""" + database_dir + StrID['str_id'][0] + """.json</b></p>
                <p>sent FS_config to <b>""" + database_dir + event_id + """CurrFS.json</b></p>
                <p>Check the <b>FS_view</b> in the <b>new tab</b></p>"""
        else:
            tab1_Div_Tb.text = """<p>Check the <b>FS_clean </b> in the <b>new tab</b></p>"""
    else:
        tab1_Div_Tb.text = """<p><b>Warning: No StrID selected. Select one StrID first!!!</b></p>"""


def tab1_update_saveStrID():
    os.system('cp {}StrID_list_tmp.json {}StrID_list.json'.format(database_dir + event_id, database_dir + event_id))
    tab1_Div_Tb.text = """<p>StrID data saved to <b>""" + '{}StrID_list.json</b></p>'.format(database_dir + event_id)


def tab1_update_reloadStrID():
    os.system('cp {}StrID_list.json {}StrID_list_tmp.json'.format(database_dir + event_id, database_dir + event_id))
    StrIDList = pd.read_json(database_dir + event_id + 'StrID_list_tmp.json')
    StrIDList = StrIDList.sort_values(by='timeran', ascending=1)
    StrIDList['time'] = [ll - tab1_tim_square[0] for ll in StrIDList['time']]
    tab1_render_patch.data_source.data = ColumnDataSource(StrIDList).data
    tab1_Div_Tb.text = """<p>StrID data reloaded from <b>""" + '{}StrID_list.json</b></p>'.format(
        database_dir + event_id)


tab1_BUT_OPT = dict(width=config_plot['plot_config']['tab_QLook']['StrID_DataTb_BUT_wdth'])
tab1_BUT_addStrID = Button(label='Add to StrID', button_type='success', **tab1_BUT_OPT)
tab1_BUT_addStrID.on_click(tab1_update_addStrID)
tab1_BUT_deleteStrID = Button(label='Delete StrID', button_type='danger', **tab1_BUT_OPT)
tab1_BUT_deleteStrID.on_click(tab1_update_deleteStrID)
tab1_BUT_saveStrID = Button(label='Save StrID', button_type='primary', **tab1_BUT_OPT)
tab1_BUT_saveStrID.on_click(tab1_update_saveStrID)
tab1_BUT_reloadStrID = Button(label='Reload StrID', button_type='warning', **tab1_BUT_OPT)
tab1_BUT_reloadStrID.on_click(tab1_update_reloadStrID)
tab1_BUT_FSviewStrID = Button(label='FS View', button_type='primary', **tab1_BUT_OPT)
tab1_BUT_FSviewStrID.on_click(tab1_update_FSviewStrID)

tab1_SPCR_LFT_DataTb_dspec = Spacer(width=10, height=100)
tab1_SPCR_RGT_DataTb_dspec = Spacer(width=20, height=100)
tab1_SPCR_ABV_DataTb_dspec = Spacer(width=100, height=18)
tab1_SPCR_LFT_But = Spacer(width=10, height=25)
tab1_SPCR_LFT_DataTb_evt = Spacer(width=20, height=15)
tab1_SPCR_ABV_DataTb_evt = Spacer(width=100, height=18)


def tab1_exit():
    tab1_Div_exit.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    print 'You may close the tab anytime you like.'
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


tab1_BUT_exit = Button(label='Exit QLook', width=150, button_type='danger')
tab1_BUT_exit.on_click(tab1_exit)

panel1 = column(tab1_p_dspec,
                row(widgetbox(tab1_Select_bl, tab1_Select_pol, tab1_Select_colorspace, tab1_BUT_exit, tab1_Div_exit,
                              width=150),
                    tab1_SPCR_LFT_DataTb_evt, tab1_SPCR_LFT_DataTb_dspec, column(tab1_DataTb_dspec, tab1_Div_Tb),
                    tab1_SPCR_LFT_But,
                    widgetbox(tab1_BUT_FSviewStrID, tab1_input_StrID, tab1_BUT_addStrID, tab1_BUT_deleteStrID,
                              tab1_BUT_saveStrID,
                              tab1_BUT_reloadStrID)))

print("--- %s seconds ---" % (time.time() - start_timestamp))

lout = panel1

curdoc().add_root(lout)
curdoc().title = "QLook tool"
