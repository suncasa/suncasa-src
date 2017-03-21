import json
import os
import json
import pickle
import numpy as np
import pandas as pd
import sunpy.map
from sys import platform
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Slider, Button, RadioButtonGroup, CheckboxGroup,
                          BoxSelectTool, TapTool, LassoSelectTool, Tool, HoverTool, Spacer, LabelSet, Div)
from bokeh.models.widgets import Select, TextInput
from bokeh.plotting import figure, curdoc
from bokeh.core.properties import Instance
import glob
from astropy.time import Time
import calendar
import drms
from sunpy.lightcurve import GOESLightCurve
from sunpy.time import TimeRange
from suncasa.utils import DButil

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

if platform == "linux" or platform == "linux2":
    print 'Runing EvtBrowser in Linux platform'
    # for ll in xrange(5010, 5010 + 10):
    #     os.system('fuser -n tcp -k {}'.format(ll))
elif platform == "darwin":
    print 'Runing EvtBrowser in OS X platform'
    # for ll in xrange(5010, 5010 + 10):
    #     os.system(
    #         'port=($(lsof -i tcp:{}|grep python2.7 |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
    #     os.system(
    #         'port=($(lsof -i tcp:{}|grep Google |cut -f2 -d" ")); [[ -n "$port" ]] && kill -9 $port'.format(ll))
elif platform == "win32":
    print 'Runing EvtBrowser in Windows platform'

ports = []

'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
'''load config file'''
config_main = DButil.loadjsonfile(suncasa_dir + 'DataBrowser/config.json')

database_dir = os.path.expandvars(config_main['datadir']['database']) + '/aiaBrowserData/'
if not os.path.exists(database_dir):
    os.makedirs(database_dir)

SDOdir = DButil.getSDOdir(config_main, database_dir, suncasa_dir)
print SDOdir
if not os.path.exists(SDOdir):
    os.makedirs(SDOdir)

YY_select_st = Select(title="Start Time: Year", value="2014", options=['{:04d}'.format(ll) for ll in range(2010, 2025)],
                      width=config_main['plot_config']['tab_aiaBrowser']['YY_select_wdth'])
MM_select_st = Select(title="Month", value="01", options=['{:02d}'.format(ll) for ll in range(1, 13)],
                      width=config_main['plot_config']['tab_aiaBrowser']['MM_select_wdth'])
DD_select_st = Select(title="Day", value="01", options=['{:02d}'.format(ll) for ll in range(1, 32)],
                      width=config_main['plot_config']['tab_aiaBrowser']['DD_select_wdth'])
hh_select_st = Select(title="hh", value="00", options=['{:02d}'.format(ll) for ll in range(0, 25)],
                      width=config_main['plot_config']['tab_aiaBrowser']['hh_select_wdth'])
mm_select_st = Select(title="mm", value="00", options=['{:02d}'.format(ll) for ll in range(0, 60)],
                      width=config_main['plot_config']['tab_aiaBrowser']['mm_select_wdth'])
ss_select_st = Select(title="ss", value="00", options=['{:02d}'.format(ll) for ll in range(0, 60)],
                      width=config_main['plot_config']['tab_aiaBrowser']['ss_select_wdth'])

YY_select_ed = Select(title="End Time: Year", value="2014", options=['{:04d}'.format(ll) for ll in range(2010, 2025)],
                      width=config_main['plot_config']['tab_aiaBrowser']['YY_select_wdth'])
MM_select_ed = Select(title="Month", value="01", options=['{:02d}'.format(ll) for ll in range(1, 13)],
                      width=config_main['plot_config']['tab_aiaBrowser']['MM_select_wdth'])
DD_select_ed = Select(title="Day", value="01", options=['{:02d}'.format(ll) for ll in range(1, 32)],
                      width=config_main['plot_config']['tab_aiaBrowser']['DD_select_wdth'])
hh_select_ed = Select(title="hh", value="23", options=['{:02d}'.format(ll) for ll in range(0, 24)],
                      width=config_main['plot_config']['tab_aiaBrowser']['hh_select_wdth'])
mm_select_ed = Select(title="mm", value="59", options=['{:02d}'.format(ll) for ll in range(0, 60)],
                      width=config_main['plot_config']['tab_aiaBrowser']['mm_select_wdth'])
ss_select_ed = Select(title="ss", value="59", options=['{:02d}'.format(ll) for ll in range(0, 60)],
                      width=config_main['plot_config']['tab_aiaBrowser']['ss_select_wdth'])
Text_PlotID = TextInput(value='', title="Plot ID:", width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])

print database_dir + 'Trange_args.json'
try:
    with open(database_dir + 'Trange_args.json', 'r') as fp:
        trcfg = json.load(fp)
    print trcfg['aiaBrowser']['tst']
    YY_select_st.value, MM_select_st.value = trcfg['aiaBrowser']['tst'][0:4], trcfg['aiaBrowser']['tst'][5:7]
    DD_select_st.value, hh_select_st.value = trcfg['aiaBrowser']['tst'][8:10], trcfg['aiaBrowser']['tst'][11:13]
    mm_select_st.value, ss_select_st.value = trcfg['aiaBrowser']['tst'][14:16], trcfg['aiaBrowser']['tst'][17:19]
    YY_select_ed.value, MM_select_ed.value = trcfg['aiaBrowser']['ted'][0:4], trcfg['aiaBrowser']['ted'][5:7]
    DD_select_ed.value, hh_select_ed.value = trcfg['aiaBrowser']['ted'][8:10], trcfg['aiaBrowser']['ted'][11:13]
    mm_select_ed.value, ss_select_ed.value = trcfg['aiaBrowser']['ted'][14:16], trcfg['aiaBrowser']['ted'][17:19]
except:
    trcfg = {}


def YY_select_st_select_update(attrname, old, new):
    YY_select_ed.value = YY_select_st.value
    DD_select_st.options = ['{:02d}'.format(ll) for ll in
                            range(1, calendar.monthrange(int(YY_select_st.value), int(MM_select_st.value))[1] + 1)]
    DD_select_ed.options = DD_select_st.options
    Text_PlotID.value = getPlotID()
    selected_time['YY_select_st'] = YY_select_st.value


YY_select_st.on_change('value', YY_select_st_select_update)


def MM_select_st_select_update(attrname, old, new):
    MM_select_ed.value = MM_select_st.value
    DD_select_st.options = ['{:02d}'.format(ll) for ll in
                            range(1, calendar.monthrange(int(YY_select_st.value), int(MM_select_st.value))[1] + 1)]
    DD_select_ed.options = DD_select_st.options
    Text_PlotID.value = getPlotID()
    selected_time['MM_select_st'] = MM_select_st.value


MM_select_st.on_change('value', MM_select_st_select_update)


def DD_select_st_select_update(attrname, old, new):
    DD_select_ed.value = DD_select_st.value
    Text_PlotID.value = getPlotID()
    selected_time['DD_select_st'] = DD_select_st.value


DD_select_st.on_change('value', DD_select_st_select_update)


def hh_select_st_select_update(attrname, old, new):
    hh_select_ed.value = hh_select_st.value
    Text_PlotID.value = getPlotID()
    selected_time['hh_select_st'] = hh_select_st.value


hh_select_st.on_change('value', hh_select_st_select_update)


def mm_select_st_select_update(attrname, old, new):
    mm_select_ed.value = mm_select_st.value
    Text_PlotID.value = getPlotID()
    selected_time['mm_select_st'] = mm_select_st.value


mm_select_st.on_change('value', mm_select_st_select_update)


def ss_select_st_select_update(attrname, old, new):
    ss_select_ed.value = ss_select_st.value
    Text_PlotID.value = getPlotID()
    selected_time['ss_select_st'] = ss_select_st.value


ss_select_st.on_change('value', ss_select_st_select_update)


def YY_select_ed_select_update(attrname, old, new):
    if checktime():
        DD_select_ed.options = ['{:02d}'.format(ll) for ll in
                                range(1, calendar.monthrange(int(YY_select_ed.value), int(MM_select_ed.value))[1] + 1)]
        Text_PlotID.value = getPlotID()
        selected_time['YY_select_ed'] = YY_select_ed.value
    else:
        YY_select_ed.value = selected_time['YY_select_ed']


YY_select_ed.on_change('value', YY_select_ed_select_update)


def MM_select_ed_select_update(attrname, old, new):
    if checktime():
        DD_select_ed.options = ['{:02d}'.format(ll) for ll in
                                range(1, calendar.monthrange(int(YY_select_ed.value), int(MM_select_ed.value))[1] + 1)]
        Text_PlotID.value = getPlotID()
        selected_time['MM_select_ed'] = MM_select_ed.value
    else:
        MM_select_ed.value = selected_time['MM_select_ed']


MM_select_ed.on_change('value', MM_select_ed_select_update)


def DD_select_ed_select_update(attrname, old, new):
    if checktime():
        Text_PlotID.value = getPlotID()
        selected_time['DD_select_ed'] = DD_select_ed.value
    else:
        DD_select_ed.value = selected_time['DD_select_ed']


DD_select_ed.on_change('value', DD_select_ed_select_update)


def hh_select_ed_select_update(attrname, old, new):
    if checktime():
        Text_PlotID.value = getPlotID()
        selected_time['hh_select_ed'] = hh_select_ed.value
    else:
        hh_select_ed.value = selected_time['hh_select_ed']


hh_select_ed.on_change('value', hh_select_ed_select_update)


def mm_select_ed_select_update(attrname, old, new):
    if checktime():
        Text_PlotID.value = getPlotID()
        selected_time['mm_select_ed'] = mm_select_ed.value
    else:
        mm_select_ed.value = selected_time['mm_select_ed']


mm_select_ed.on_change('value', mm_select_ed_select_update)


def ss_select_ed_select_update(attrname, old, new):
    if checktime():
        Text_PlotID.value = getPlotID()
        selected_time['ss_select_ed'] = ss_select_ed.value
    else:
        ss_select_ed.value = selected_time['ss_select_ed']


ss_select_ed.on_change('value', ss_select_ed_select_update)

Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b> first before closing the tab</p></b>""",
               width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])
Div_JSOC_info = Div(text="""""",
                    width=config_main['plot_config']['tab_aiaBrowser']['divJSOCinfo_wdth'])

Text_Cadence = TextInput(value='12s', title="Cadence:",
                         width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])
Text_email = TextInput(value='', title="JSOC registered email:",
                       width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])
try:
    email = config_main['core']['JSOC_reg_email']
    Text_email.value = email
except:
    pass

global selected_time
selected_time = {'YY_select_st': YY_select_st.value, 'MM_select_st': MM_select_st.value,
                 'DD_select_st': DD_select_st.value, 'hh_select_st': hh_select_st.value,
                 'mm_select_st': mm_select_st.value, 'ss_select_st': ss_select_st.value,
                 'YY_select_ed': YY_select_ed.value, 'MM_select_ed': MM_select_ed.value,
                 'DD_select_ed': DD_select_ed.value, 'hh_select_ed': hh_select_ed.value,
                 'mm_select_ed': mm_select_ed.value, 'ss_select_ed': ss_select_ed.value}


def gettimestr(YY, MM, DD, hh, mm, ss):
    return '{}-{}-{} {}:{}:{}'.format(YY, MM, DD, hh, mm, ss)


def gettime():
    tststr = gettimestr(YY_select_st.value, MM_select_st.value, DD_select_st.value, hh_select_st.value,
                        mm_select_st.value, ss_select_st.value)
    tedstr = gettimestr(YY_select_ed.value, MM_select_ed.value, DD_select_ed.value, hh_select_ed.value,
                        mm_select_ed.value, ss_select_ed.value)
    return Time([tststr, tedstr], format='iso', scale='utc')


def checktime():
    [tst, ted] = gettime()
    if ted.mjd <= tst.mjd:
        return False
    else:
        return True


def getPlotID():
    [tst, ted] = gettime()
    tr = TimeRange(tst.iso, ted.iso)
    if tr.minutes.value < 60:
        dur = '{:.0f}'.format(tr.seconds).replace(' ', '')
    elif tr.hours.value < 24:
        dur = '{:.0f}'.format(tr.minutes).replace(' ', '')
    elif tr.days.value < 10:
        dur = '{:.0f}'.format(tr.hours).replace(' ', '')
    else:
        dur = '{:.0f}'.format(tr.days).replace(' ', '')
    return tst.iso.replace(' ', 'T').replace('-', '').replace(':', '').split('.')[0] + '_' + dur


Text_PlotID.value = getPlotID()

wavelngth_list = ["goes", "HMI_Magnetogram", "1700", "1600", "304", "171", "193", "211", "335", "94", "131", "All"]
Wavelngth_checkbox = CheckboxGroup(labels=wavelngth_list, active=[0],
                                   width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])
serieslist = {}
for ll in wavelngth_list:
    if ll == "HMI_Magnetogram":
        serieslist[ll] = 'hmi.M_45s'
    elif ll in ["1700", "1600"]:
        serieslist[ll] = 'aia.lev1_uv_24s'
    elif ll in ['goes', 'All', 'None']:
        serieslist[ll] = ''
    else:
        serieslist[ll] = 'aia.lev1_euv_12s'


def Wavelngth_checkbox_change(attrname, old, new):
    labelsactive = [Wavelngth_checkbox.labels[ll] for ll in Wavelngth_checkbox.active]
    if 'All' in labelsactive:
        Wavelngth_checkbox.active = range(len(Wavelngth_checkbox.labels) - 1)
        Wavelngth_checkbox.labels = wavelngth_list[0:-1] + ['None']
    elif 'None' in labelsactive:
        Wavelngth_checkbox.active = [0]
        Wavelngth_checkbox.labels = wavelngth_list


Wavelngth_checkbox.on_change('active', Wavelngth_checkbox_change)


def get_num(x):
    return int(''.join(ele for ele in x if ele.isdigit()))


global MkPlot_args_dict
MkPlot_args_dict = {}


def DownloadData():
    [tst, ted] = gettime()
    if ted.mjd <= tst.mjd:
        Div_JSOC_info.text = '''Error: start time must occur earlier than end time. please re-enter start time and end time!!!'''
    elif len(Wavelngth_checkbox.active) == 0:
        Div_JSOC_info.text = '''Error: at least choose one wavelength!!!'''
    else:
        Div_JSOC_info.text = ''''''
        c = drms.Client(verbose=True)
        export_protocol = 'fits'
        if not Text_email.value:
            Div_JSOC_info.text = '''Error: provide your JSOC registered email address!!!'''
        else:
            if not c.check_email(Text_email.value):
                raise RuntimeError('Email address is not valid or not registered.')

        labelsactive = [Wavelngth_checkbox.labels[ll] for ll in Wavelngth_checkbox.active]
        if 'goes' in labelsactive:
            global goes
            tr = TimeRange(tst.iso, ted.iso)
            goes = GOESLightCurve.create(tr)
            fout = database_dir + 'goes-' + Text_PlotID.value
            with open(fout, 'wb') as fp:
                pickle.dump(goes, fp)
            Div_JSOC_info.text = """<p>{} saved.</p>""".format(fout)
            MkPlot_args_dict['goesfile'] = os.path.basename(fout)
        for series in ['hmi.M_45s', 'aia.lev1_uv_24s', 'aia.lev1_euv_12s']:
            waves = []
            for ll in labelsactive:
                if serieslist[ll] == series:
                    waves.append(ll)
            if len(waves) > 0:
                try:
                    tsel = tst.iso.replace(' ', 'T') + '_TAI-' + ted.iso.replace(' ', 'T') + '_TAI'
                    wave = ','.join(waves)
                    cadence = Text_Cadence.value
                    if cadence[-1] == 's' and get_num(cadence) < 12:
                        cadence = '12s'
                        Text_Cadence.value = cadence
                    if series == 'aia.lev1_uv_24s':
                        if cadence[-1] == 's' and get_num(cadence) < 24:
                            cadence = '24s'
                    if series == 'hmi.M_45s':
                        wave = ''
                        if cadence[-1] == 's' and get_num(cadence) < 45:
                            cadence = '45s'
                    segments = 'image'
                    qstr = '%s[%s@%s][%s]{%s}' % (series, tsel, cadence, wave, segments)
                    print qstr
                    r = c.export(qstr, method='url', protocol=export_protocol, email=Text_email.value)
                    Div_JSOC_info.text = Div_JSOC_info.text + """<p>Submitting export request {}...</p>""".format(qstr)
                    Div_JSOC_info.text = Div_JSOC_info.text + """<p>Request URL: {}</p>""".format(r.request_url)
                    Div_JSOC_info.text = Div_JSOC_info.text + """<p>{:d} file(s) available for download.</p>""".format(
                        len(r.urls))
                    idx2download = DButil.FileNotInList(r.data['filename'],
                                                        DButil.readsdofile(datadir=SDOdir, wavelength=wave,
                                                                           jdtime=[tst.jd, ted.jd],
                                                                           isexists=True))
                    if len(idx2download) > 0:
                        r.download(SDOdir, index=idx2download)
                        # Div_JSOC_info.text = Div_JSOC_info.text + """<p>Target file(s) existed.</p>"""

                    filename = glob.glob(SDOdir + '*.fits')
                    dirs = DButil.getsdodir(filename)
                    for ll, dd in enumerate(dirs['dir']):
                        if not os.path.exists(SDOdir + dd):
                            os.makedirs(SDOdir + dd)
                        os.system('mv {}/*{}*.fits {}{}'.format(SDOdir, dirs['timstr'][ll], SDOdir, dd))

                    Div_JSOC_info.text = Div_JSOC_info.text + """<p>Download finished.</p>"""
                    Div_JSOC_info.text = Div_JSOC_info.text + """<p>Download directory: {}</p>""".format(
                        os.path.abspath(SDOdir))
                except:
                    print qstr + ' fail to export'


BUT_DownloadData = Button(label='Download Data', width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'],
                          button_type='primary')
BUT_DownloadData.on_click(DownloadData)


def MkPlot():
    [tst, ted] = gettime()
    if ted.mjd <= tst.mjd:
        Div_info.text = '''Error: start time must occur earlier than end time. please re-enter start time and end time!!!'''
    else:
        labelsactive = [Wavelngth_checkbox.labels[ll] for ll in Wavelngth_checkbox.active]
        Slabelsactive = set(labelsactive)
        Swavelength = set(["1700", "1600", "304", "171", "193", "211", "335", "94", "131"])
        if len(Swavelength.intersection(Slabelsactive)) == 1:
            MkPlot_args_dict['wavelength'] = list(Swavelength.intersection(Slabelsactive))[0]
            MkPlot_args_dict['tst'] = tst.iso
            MkPlot_args_dict['ted'] = ted.iso
            MkPlot_args_dict['PlotID'] = Text_PlotID.value
            fout = database_dir + 'MkPlot_args.json'
            with open(fout, 'w') as fp:
                json.dump(MkPlot_args_dict, fp)
            Div_info.text = '''<p><b>{} saved.</b></p>'''.format(fout)
            port = DButil.getfreeport()
            print 'bokeh serve {}aiaBrowser/MkPlot --show --port {} &'.format(suncasa_dir, port)
            os.system('bokeh serve {}aiaBrowser/MkPlot --show --port {} &'.format(suncasa_dir, port))
            ports.append(port)
            Div_info.text = Div_info.text + """<p>Check the <b>MkPlot</b> in the <b>new tab</b></p>"""
        else:
            Div_info.text = Div_info.text + """<p>Choose one <b>AIA</b> wavelength to <b>MkPlot!!</b></p>"""


BUT_MkPlot = Button(label='MkPlot', width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'],
                    button_type='success')
BUT_MkPlot.on_click(MkPlot)

Text_sdodir = TextInput(value=SDOdir, title="SDO Directory:",
                        width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])


def Buttonaskdir_handler():
    import Tkinter
    import tkFileDialog
    global SDOdir
    tkRoot = Tkinter.Tk()
    tkRoot.withdraw()  # Close the root window
    in_path = tkFileDialog.askdirectory(initialdir=SDOdir, parent=tkRoot) + '/'
    tkRoot.destroy()
    if in_path:
        Text_sdodir.value = in_path
        SDOdir = in_path
        config_main['datadir']['SDOdir'] = SDOdir
        fout = suncasa_dir + 'DataBrowser/config.json'
        DButil.updatejsonfile(fout, config_main)
        print in_path


But_sdodir = Button(label='Directory', width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'])
But_sdodir.on_click(Buttonaskdir_handler)

DiffImglabelsdict = config_main['plot_config']['tab_MkPlot']['imagetype']
DiffImglabels = []
DiffImglabelsdict_r = {}
for k, v in DiffImglabelsdict.items():
    DiffImglabels.append(v)
    DiffImglabelsdict_r[v] = k


def Select_DiffImg_update(attrname, old, new):
    global MkPlot_args_dict
    MkPlot_args_dict['imagetype'] = DiffImglabelsdict_r[Select_DiffImg.value]
    # outfile = database_dir + 'MkPlot_args.json'
    # DButil.updatejsonfile(outfile, MkPlot_args_dict)
    # Div_info.text = """<p><b>Refresh</b> the <b>web page</b> to load diff images</p>"""


Select_DiffImg = Select(title="Difference images:", value=DiffImglabelsdict['image'], options=DiffImglabels,
                        width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
Select_DiffImg.on_change('value', Select_DiffImg_update)


def exit_update():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


BUT_exit = Button(label='Exit', width=config_main['plot_config']['tab_aiaBrowser']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit_update)

SPCR_LFT_widgetbox = Spacer(width=50, height=10)
SPCR_RGT_widgetbox = Spacer(width=50, height=10)

lout = row(column(row(YY_select_st, MM_select_st, DD_select_st, hh_select_st, mm_select_st, ss_select_st),
                  row(YY_select_ed, MM_select_ed, DD_select_ed, hh_select_ed, mm_select_ed, ss_select_ed)),
           SPCR_LFT_widgetbox, Wavelngth_checkbox, SPCR_RGT_widgetbox,
           widgetbox(Text_Cadence, Text_email, Text_sdodir, But_sdodir, BUT_DownloadData,
                     Text_PlotID, Select_DiffImg,
                     BUT_MkPlot, BUT_exit, Div_info,
                     width=config_main['plot_config']['tab_aiaBrowser']['button_wdth']))
# def timeout_callback():
#     print 'timeout'
#     raise SystemExit


curdoc().add_root(lout)
curdoc().title = "AIAbrowser"
