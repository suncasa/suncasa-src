import glob
import json
import os
import astropy.units as u
import bokeh.palettes as bp
import numpy as np
import pandas as pd
import sunpy.map
from astropy.time import Time
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (ColumnDataSource, Slider, Button, TextInput, CheckboxGroup, RadioGroup,
                          BoxSelectTool, TapTool, Div)
from bokeh.plotting import figure, curdoc
from scipy.interpolate import splev, splrep, splprep
from suncasa.utils.puffin import PuffinMap
from suncasa.utils import DButil

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"

'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
'''load config file'''
with open(suncasa_dir + 'aiaBrowser/config.json', 'r') as fp:
    config_plot = json.load(fp)

database_dir = config_plot['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/aiaBrowserData/'
if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))
download_dir = database_dir + 'Download/'
if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))
fin = database_dir + 'MkPlot_args.json'
with open(fin, 'rb') as fp:
    MkPlot_args_dict = json.load(fp)

aiafile = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/database/sun_20141101_t191020-191040.50ms/AIA/aia_lev1_171a_2014_11_01t19_10_23_34z_image_lev1.fits.fits'
aiamap = sunpy.map.Map(aiafile)
MapRES = config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_RES']
dimensions = u.Quantity([MapRES, MapRES], u.pixel)
aia_RSPmap = aiamap.resample(dimensions)
aia_RSPmap.data = DButil.sdo_aia_scale(image=aia_RSPmap.data / aia_RSPmap.exposure_time.value,
                                       wavelength='171')
aia_RSP_pfmap = PuffinMap(smap=aia_RSPmap,
                          plot_height=config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_hght'],
                          plot_width=config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
p_aiamap, r_aiamap = aia_RSP_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                           palette=bp.viridis(256))
mapx_aia_RSPmap, mapy_aia_RSPmap = aia_RSP_pfmap.meshgrid(rescale=0.2)
mapx_aia_RSPmap, mapy_aia_RSPmap = mapx_aia_RSPmap.value, mapy_aia_RSPmap.value
ImgDF0_aia_RSPmap = pd.DataFrame({'xx': mapx_aia_RSPmap.ravel(), 'yy': mapy_aia_RSPmap.ravel()})
SRC_aia_RSPmap_square = ColumnDataSource(ImgDF0_aia_RSPmap)
r_aia_RSPmap_square = p_aiamap.square('xx', 'yy', source=SRC_aia_RSPmap_square,
                                      fill_alpha=0.0, fill_color=None,
                                      line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                      selection_fill_color=None,
                                      nonselection_fill_alpha=0.0,
                                      selection_line_alpha=0.0, selection_line_color=None,
                                      nonselection_line_alpha=0.0,
                                      size=max(config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_hght'] /
                                               config_plot['plot_config']['tab_MkPlot']['aia_RSP_squareny'],
                                               config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'] /
                                               config_plot['plot_config']['tab_MkPlot']['aia_RSP_squarenx']))

x0 = aia_RSP_pfmap.smap.center.x
y0 = aia_RSP_pfmap.smap.center.y
lengthx = config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'] * u.arcsec
lengthy = config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'] * u.arcsec
x0, x1 = x0 - lengthx / 2, x0 + lengthx / 2
y0, y1 = y0 - lengthy / 2, y0 + lengthy / 2
aia_submap = aiamap.submap(u.Quantity([x0, x1]), u.Quantity([y0, y1]))
aia_submap.data = DButil.sdo_aia_scale(image=aia_submap.data / aia_submap.exposure_time.value,
                                       wavelength='171')
xaxis_aia_submap, yaxis_aia_submap = aia_submap.data.shape
xaxis_aia_submap_canvas, yaxis_aia_submap_canvas = aia_submap.data.shape
Canvas_scaleX, Canvas_scaleY = float(xaxis_aia_submap) / xaxis_aia_submap_canvas, float(
    yaxis_aia_submap) / yaxis_aia_submap_canvas
Canvas_scale = np.sqrt(Canvas_scaleX ** 2 + Canvas_scaleY ** 2)
SRC_aia_RSPmap_patch = ColumnDataSource(
    pd.DataFrame({'xx': [x0.value, x1.value, x1.value, x0.value], 'yy': [y0.value, y0.value, y1.value, y1.value]}))
r_aia_RSPmap_square_patch = p_aiamap.patch('xx', 'yy', source=SRC_aia_RSPmap_patch,
                                           fill_color=None, fill_alpha=0.5, line_color="white",
                                           line_alpha=1, line_width=1)
p_aiamap.add_tools(BoxSelectTool(renderers=[r_aia_RSPmap_square]))


def update_aia_submap_image(attrname, old, new):
    global ImgDF0_aia_RSPmap, aia_submap
    global xaxis_aia_submap, yaxis_aia_submap
    global Canvas_scale
    r_aia_RSPmap_square_selected = SRC_aia_RSPmap_square.selected['1d']['indices']
    if r_aia_RSPmap_square_selected:
        ImgDF_aia_RSPmap = ImgDF0_aia_RSPmap.iloc[r_aia_RSPmap_square_selected, :]
        x0, x1 = ImgDF_aia_RSPmap['xx'].min(), ImgDF_aia_RSPmap['xx'].max()
        y0, y1 = ImgDF_aia_RSPmap['yy'].min(), ImgDF_aia_RSPmap['yy'].max()
        patch_size = min(x1 - x0, y1 - y0, config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'])
        x1, y0 = x0 + patch_size, y1 - patch_size
        r_aia_RSPmap_square_patch.data_source.data = ColumnDataSource(
            pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
        aia_submap = aiamap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                   u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        aia_submap.data = DButil.sdo_aia_scale(image=aia_submap.data / aia_submap.exposure_time.value,
                                               wavelength='171')
        xaxis_aia_submap, yaxis_aia_submap = aia_submap.data.shape
        Canvas_scaleX, Canvas_scaleY = float(xaxis_aia_submap) / xaxis_aia_submap_canvas, float(
            yaxis_aia_submap) / yaxis_aia_submap_canvas
        Canvas_scale = np.sqrt(Canvas_scaleX ** 2 + Canvas_scaleY ** 2)
        r_aia_submap.data_source.data['image'] = [aia_submap.data]
        Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_aia_submap ** 2 + yaxis_aia_submap ** 2) / 4)
        ClearDraw()


SRC_aia_RSPmap_square.on_change('selected', update_aia_submap_image)

aia_submap_pfmap = PuffinMap(smap=aia_submap,
                             plot_height=config_plot['plot_config']['tab_MkPlot']['aia_submap_hght'],
                             plot_width=config_plot['plot_config']['tab_MkPlot']['aia_submap_wdth'])
p_aia_submap, r_aia_submap = aia_submap_pfmap.PlotMap(DrawLimb=False, DrawGrid=False,
                                                      palette=bp.viridis(256), ignore_coord=True)
mapx_aia_submap, mapy_aia_submap = aia_submap_pfmap.meshgridpix(rescale=0.125)
mapx_aia_submap, mapy_aia_submap = mapx_aia_submap.value, mapy_aia_submap.value
ImgDF0_aia_submap = pd.DataFrame({'xx': mapx_aia_submap.ravel(), 'yy': mapy_aia_submap.ravel()})
pointDF_aia_submap = pd.DataFrame({'xx': [], 'yy': []})  ## the selected point to fit in aia submap
SRC_aia_submap_square = ColumnDataSource(ImgDF0_aia_submap)
r_aia_submap_square = p_aia_submap.square('xx', 'yy', source=SRC_aia_submap_square,
                                          fill_alpha=0.0, fill_color=None,
                                          line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                          selection_fill_color=None,
                                          nonselection_fill_alpha=0.0,
                                          selection_line_alpha=0.0, selection_line_color=None,
                                          nonselection_line_alpha=0.0,
                                          size=max(config_plot['plot_config']['tab_MkPlot']['aia_submap_hght'] /
                                                   config_plot['plot_config']['tab_MkPlot']['aia_RSP_squareny'],
                                                   config_plot['plot_config']['tab_MkPlot']['aia_submap_wdth'] /
                                                   config_plot['plot_config']['tab_MkPlot']['aia_RSP_squarenx']))
p_aia_submap.add_tools(TapTool(renderers=[r_aia_submap_square]))
SRC_aia_submap_line = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
SRC_aia_submap_line0 = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
SRC_aia_submap_line1 = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
r_aia_submap_line = p_aia_submap.line('xx', 'yy', source=SRC_aia_submap_line, line_color='white', line_width=2,
                                      alpha=0.7)
r_aia_submap_line0 = p_aia_submap.line('xx', 'yy', source=SRC_aia_submap_line0, line_color='white', line_width=1.5,
                                       alpha=0.7, line_dash='dashed')
r_aia_submap_line1 = p_aia_submap.line('xx', 'yy', source=SRC_aia_submap_line1, line_color='white', line_width=1.5,
                                       alpha=0.7, line_dash='dashed')
r_aia_submap_cross = p_aia_submap.cross(x='xx', y='yy', size=15, line_width=2, alpha=0.7, line_color='firebrick',
                                        source=ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})))

TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
tab2_p_dspec = figure(tools=TOOLS,
                      plot_width=config_plot['plot_config']['tab_MkPlot']['time_plot_wdth'],
                      plot_height=config_plot['plot_config']['tab_MkPlot']['time_plot_hght'],
                      # x_range=(tab2_dtim[0], tab2_dtim[-1]), y_range=(tab2_freq[0], tab2_freq[-1]),
                      toolbar_location="above")
# tim0_char = Time(xx[0] / 3600. / 24., format='jd', scale='utc', precision=3, out_subfmt='date_hms').iso
tab2_p_dspec.axis.visible = True
tab2_p_dspec.title.text = "light Curve"
# tab2_p_dspec.xaxis.axis_label = 'Seconds since ' + tim0_char
# tab2_p_dspec.yaxis.axis_label = 'Frequency [GHz]'




BUT_ClearDraw = Button(label='ClearDraw', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                       button_type='primary')
BUT_UndoDraw = Button(label='Undo', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                      button_type='warning')
BUT_MkStackplt = Button(label='MkStackPlt', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                        button_type='success')
BUT_next = Button(label='→', width=config_plot['plot_config']['tab_MkPlot']['button_wdth_short'],
                  button_type='success')
BUT_prev = Button(label='←', width=config_plot['plot_config']['tab_MkPlot']['button_wdth_short'],
                  button_type='warning')
BUT_last = Button(label='⇥', width=config_plot['plot_config']['tab_MkPlot']['button_wdth_short'],
                  button_type='warning')
BUT_first = Button(label='⇤', width=config_plot['plot_config']['tab_MkPlot']['button_wdth_short'],
                   button_type='warning')

Text_Cutwdth = TextInput(value='5.0', title="Cut Width (pix):")
Text_CutAng = TextInput(value='10.0', title="Cut Angle (deg):")
Text_SlitLgth = TextInput(value='{:.0f}'.format(np.sqrt(xaxis_aia_submap ** 2 + yaxis_aia_submap ** 2) / 4),
                          title="Slit Length (pix):")

Hideitemlist = ["Hide Slit", "Hide Points"]
Hideitem_checkbox = CheckboxGroup(labels=Hideitemlist, active=[])
fitmethoddict = {'0': "Polyfit", '1': "Spline", '2': "Param_Spline"}
fitMeth_radiogroup = RadioGroup(labels=["Polyfit", "Spline", "Parametric Spline"], active=0)
fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
ascending = True
clearpoint = False


def aia_submap_square_selection_change(attrname, old, new):
    global pointDF_aia_submap, ascending, cutslit
    aia_submap_square_selected = SRC_aia_submap_square.selected['1d']['indices']
    if aia_submap_square_selected:
        ImgDF_aia_submap = ImgDF0_aia_submap.iloc[aia_submap_square_selected, :]
        DFidx_selected = ImgDF_aia_submap.index[0]
        pointDF_aia_submap = pointDF_aia_submap.append(ImgDF0_aia_submap.loc[DFidx_selected, :], ignore_index=True)
        if len(pointDF_aia_submap.index) == 2:
            if pointDF_aia_submap.loc[pointDF_aia_submap.index[0], 'xx'] > pointDF_aia_submap.loc[
                pointDF_aia_submap.index[1], 'xx']:
                ascending = False
        cutslit = MakeSlit(pointDF_aia_submap)
        pushslit2plt(clearpoint=clearpoint)


SRC_aia_submap_square.on_change('selected', aia_submap_square_selection_change)
Slider_smoothing_factor = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title='smoothing factor')


def FitSlit(xx, yy, cutwidth, cutang, cutlength, s=None, method='Polyfit'):
    if len(xx) <= 3 or method == 'Polyfit':
        '''polynomial fit'''
        out = DButil.polyfit(xx, yy, cutlength, len(xx) - 1 if len(xx) <= 3 else 2)
        xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    else:
        if method == 'Param_Spline':
            '''parametic spline fit'''
            out = DButil.paramspline(xx, yy, cutlength, s=s)
            xs, ys, posangs = out['xs'], out['ys'], out['posangs']
        else:
            '''spline fit'''
            out = DButil.spline(xx, yy, cutlength, s=s)
            xs, ys, posangs = out['xs'], out['ys'], out['posangs']
    if not ascending and fitmethod != 'Param_Spline':
        xs, ys = xs[::-1], ys[::-1]
        posangs = posangs[::-1]
    dists = DButil.findDist(xs, ys)
    posangs2 = posangs + np.pi / 2
    cutwidths = np.cumsum(dists) * np.tan(cutang) + cutwidth
    xs0 = xs - cutwidths / 2. * np.cos(posangs2)
    ys0 = ys - cutwidths / 2. * np.sin(posangs2)
    xs1 = xs + cutwidths / 2. * np.cos(posangs2)
    ys1 = ys + cutwidths / 2. * np.sin(posangs2)
    return {'xcen': xs, 'ycen': ys, 'xs0': xs0, 'ys0': ys0, 'xs1': xs1, 'ys1': ys1, 'cutwidth': cutwidths,
            'posang': posangs, 'posangs': posangs2, 'dist': dists}


def MakeSlit(pointDF):
    global smoth_factor0, smoth_factor1, fitmethod
    if fitmethod != 'Param_Spline':
        pointDFtmp = pointDF.sort_values(['xx'], ascending=[True])
    else:
        pointDFtmp = pointDF
    xx = pointDFtmp.loc[:, 'xx']
    yy = pointDFtmp.loc[:, 'yy']
    if len(pointDFtmp.index) == 1:
        r_aia_submap_cross.data_source.data = ColumnDataSource(pd.DataFrame(pointDFtmp)).data
        cutslit = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                   'posang': [], 'posangs': [], 'dist': []}
    else:
        if len(pointDFtmp.index) <= 3:
            cutslit = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                              float(Text_SlitLgth.value), method=fitmethod)
        else:
            smoth_factor0 = len(yy) - np.sqrt(2 * len(yy))
            smoth_factor1 = np.var(yy) * (len(yy)) * 3
            cutslit = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                              float(Text_SlitLgth.value), s=np.exp(
                    Slider_smoothing_factor.value * (np.log(smoth_factor1) - np.log(smoth_factor0))) * smoth_factor0,
                              method=fitmethod)
    return cutslit


def pushslit2plt(clearpoint=False):
    r_aia_submap_line.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslit['xcen'], 'yy': cutslit['ycen']})).data
    r_aia_submap_line0.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslit['xs0'], 'yy': cutslit['ys0']})).data
    r_aia_submap_line1.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslit['xs1'], 'yy': cutslit['ys1']})).data
    if clearpoint:
        r_aia_submap_cross.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    else:
        r_aia_submap_cross.data_source.data = ColumnDataSource(pd.DataFrame(pointDF_aia_submap)).data


def clearslitplt():
    r_aia_submap_line.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_aia_submap_line0.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_aia_submap_line1.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_aia_submap_cross.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data


def slider_smoothing_factor_update(attrname, old, new):
    global pointDF_aia_submap, smoth_factor0, smoth_factor1, cutslit
    # if len(pointDF_aia_submap.index) > 3 and fitmethod == 'Polyfit':
    #     fitMeth_radiogroup.active = 1
    if len(pointDF_aia_submap.index) > 3:
        if fitmethod == 'Polyfit':
            fitMeth_radiogroup.active = 1
        if fitmethod == 'Spline':
            pointDFtmp = pointDF_aia_submap.sort_values(['xx'], ascending=[True])
        else:
            pointDFtmp = pointDF_aia_submap
        xx = pointDFtmp.loc[:, 'xx']
        yy = pointDFtmp.loc[:, 'yy']
        cutslit = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                          float(Text_SlitLgth.value), s=np.exp(
                Slider_smoothing_factor.value * (np.log(smoth_factor1) - np.log(smoth_factor0))) * smoth_factor0,
                          method=fitmethod)
        pushslit2plt(clearpoint=clearpoint)


Slider_smoothing_factor.on_change('value', slider_smoothing_factor_update)


def UndoDraw():
    global pointDF_aia_submap, smoth_factor0, smoth_factor1, cutslit
    if len(pointDF_aia_submap.index) > 1:
        pointDF_aia_submap.drop(len(pointDF_aia_submap.index) - 1, inplace=True)
        cutslit = MakeSlit(pointDF_aia_submap)
        pushslit2plt(clearpoint=clearpoint)
    else:
        ClearDraw()


BUT_UndoDraw.on_click(UndoDraw)


def ClearDraw():
    global pointDF_aia_submap, ascending
    Slider_smoothing_factor.value = 1.0
    ascending = True
    pointDF_aia_submap = pd.DataFrame({'xx': [], 'yy': []})
    clearslitplt()


BUT_ClearDraw.on_click(ClearDraw)


def MkStackplt():
    pass


BUT_MkStackplt.on_click(MkStackplt)


def HideItemUpdate(new):
    global clearpoint
    clearpoint = len(Hideitem_checkbox.active) == 1 and Hideitemlist[Hideitem_checkbox.active[0]] == "Hide Points"
    if len(Hideitem_checkbox.active) <= 1:
        pushslit2plt(clearpoint=clearpoint)
    else:
        clearslitplt()


Hideitem_checkbox.on_click(HideItemUpdate)


def fitMethUpdate(new):
    global fitmethod, cutslit
    fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
    cutslit = MakeSlit(pointDF_aia_submap)
    pushslit2plt(clearpoint=clearpoint)


fitMeth_radiogroup.on_click(fitMethUpdate)

Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b>
            first before closing the tab</p></b>""", width=config_plot['plot_config']['tab_MkPlot']['button_wdth'])


def default_fitparam():
    Text_Cutwdth.value = '5.0'
    Text_CutAng.value = '10.0'
    Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_aia_submap ** 2 + yaxis_aia_submap ** 2) / 4)


BUT_default_fitparam = Button(label='Default', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                              button_type='success')
BUT_default_fitparam.on_click(default_fitparam)


def exit():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


BUT_exit = Button(label='Exit', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit)

lbutt_play = row(BUT_first, BUT_prev, BUT_next, BUT_last)
lout = row(column(p_aiamap, lbutt_play), p_aia_submap,
           widgetbox(Text_SlitLgth, Text_Cutwdth, Text_CutAng, Slider_smoothing_factor,
                     Hideitem_checkbox, fitMeth_radiogroup,
                     BUT_MkStackplt, BUT_UndoDraw,
                     BUT_ClearDraw, BUT_default_fitparam,
                     BUT_exit, Div_info,
                     width=config_plot['plot_config']['tab_MkPlot']['button_wdth']))

# def timeout_callback():
#     print 'timeout'
#     raise SystemExit


curdoc().add_root(lout)
curdoc().title = "MkPlot"
