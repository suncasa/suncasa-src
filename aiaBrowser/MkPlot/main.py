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
                          BoxSelectTool, TapTool, Div, Spacer, Range1d)
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

sdofile = '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/database/aiaBrowserData/Download/aia.lev1_euv_12s.2014-01-01T190936Z.171.image_lev1.fits'
sdomap = sunpy.map.Map(sdofile)
MapRES = config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_RES']
dimensions = u.Quantity([MapRES, MapRES], u.pixel)
sdo_RSPmap = sdomap.resample(dimensions)
# sdo_RSPmap.data = DButil.sdo_aia_scale(image=sdo_RSPmap.data / sdo_RSPmap.exposure_time.value,
#                                        wavelength='171')
sdo_RSP_pfmap = PuffinMap(smap=sdo_RSPmap,
                          plot_height=config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_hght'],
                          plot_width=config_plot['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
p_sdomap, r_sdomap = sdo_RSP_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg)
mapx_sdo_RSPmap, mapy_sdo_RSPmap = sdo_RSP_pfmap.meshgrid(rescale=0.2)
mapx_sdo_RSPmap, mapy_sdo_RSPmap = mapx_sdo_RSPmap.value, mapy_sdo_RSPmap.value
ImgDF0_sdo_RSPmap = pd.DataFrame({'xx': mapx_sdo_RSPmap.ravel(), 'yy': mapy_sdo_RSPmap.ravel()})
SRC_sdo_RSPmap_square = ColumnDataSource(ImgDF0_sdo_RSPmap)
r_sdo_RSPmap_square = p_sdomap.square('xx', 'yy', source=SRC_sdo_RSPmap_square,
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

x0 = sdo_RSP_pfmap.smap.center.x
y0 = sdo_RSP_pfmap.smap.center.y
lengthx = config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'] * u.arcsec
lengthy = config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'] * u.arcsec
x0, x1 = x0 - lengthx / 2, x0 + lengthx / 2
y0, y1 = y0 - lengthy / 2, y0 + lengthy / 2
sdosubmp = sdomap.submap(u.Quantity([x0, x1]), u.Quantity([y0, y1]))
# sdosubmp.data = DButil.sdo_aia_scale(image=sdosubmp.data / sdosubmp.exposure_time.value,
#                                      wavelength='171')
xaxis_sdosubmp, yaxis_sdosubmp = sdosubmp.data.shape
xaxis_sdosubmp_canvas, yaxis_sdosubmp_canvas = sdosubmp.data.shape
Canvas_scaleX, Canvas_scaleY = float(xaxis_sdosubmp) / xaxis_sdosubmp_canvas, float(
    yaxis_sdosubmp) / yaxis_sdosubmp_canvas
# Canvas_scale = np.sqrt(Canvas_scaleX ** 2 + Canvas_scaleY ** 2)
Canvas_scale = (Canvas_scaleX + Canvas_scaleY) / 2.0
SRC_sdo_RSPmap_patch = ColumnDataSource(
    pd.DataFrame({'xx': [x0.value, x1.value, x1.value, x0.value], 'yy': [y0.value, y0.value, y1.value, y1.value]}))
r_sdo_RSPmap_square_patch = p_sdomap.patch('xx', 'yy', source=SRC_sdo_RSPmap_patch,
                                           fill_color=None, fill_alpha=0.5, line_color="white",
                                           line_alpha=1, line_width=1)
p_sdomap.add_tools(BoxSelectTool(renderers=[r_sdo_RSPmap_square]))


def update_sdosubmp_image(attrname, old, new):
    global ImgDF0_sdo_RSPmap, sdosubmp
    global xaxis_sdosubmp, yaxis_sdosubmp
    global Canvas_scale
    r_sdo_RSPmap_square_selected = SRC_sdo_RSPmap_square.selected['1d']['indices']
    if r_sdo_RSPmap_square_selected:
        ImgDF_sdo_RSPmap = ImgDF0_sdo_RSPmap.iloc[r_sdo_RSPmap_square_selected, :]
        x0, x1 = ImgDF_sdo_RSPmap['xx'].min(), ImgDF_sdo_RSPmap['xx'].max()
        y0, y1 = ImgDF_sdo_RSPmap['yy'].min(), ImgDF_sdo_RSPmap['yy'].max()
        print x1 > x0 + mapx_sdo_RSPmap[0, 1] - mapx_sdo_RSPmap[0, 0] and y1 > y0 + mapy_sdo_RSPmap[0, 1] - \
                                                                               mapy_sdo_RSPmap[0, 0]
        if x1 > x0 + mapx_sdo_RSPmap[0, 1] - mapx_sdo_RSPmap[0, 0] and y1 > y0 + mapy_sdo_RSPmap[0, 1] - \
                mapy_sdo_RSPmap[0, 0]:
            ## select at least 4 points
            patch_size = min(x1 - x0, y1 - y0, config_plot['plot_config']['tab_MkPlot']['aia_submap_sz_max'])
            x1, y0 = x0 + patch_size, y1 - patch_size
            r_sdo_RSPmap_square_patch.data_source.data = ColumnDataSource(
                pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
            sdosubmp = sdomap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                     u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
            xaxis_sdosubmp, yaxis_sdosubmp = sdosubmp.data.shape
            Canvas_scaleX, Canvas_scaleY = float(xaxis_sdosubmp) / xaxis_sdosubmp_canvas, float(
                yaxis_sdosubmp) / yaxis_sdosubmp_canvas
            # Canvas_scale = np.sqrt(Canvas_scaleX ** 2 + Canvas_scaleY ** 2)
            Canvas_scale = (Canvas_scaleX + Canvas_scaleY) / 2.0
            sdosubmp.data[~np.isnan(sdosubmp.data)] = sdosubmp.data[~np.isnan(
                sdosubmp.data)] / sdosubmp.exposure_time.value * 0.2
            r_sdosubmp.data_source.data['image'] = [sdosubmp.data]
            Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)
            ClearDraw()


SRC_sdo_RSPmap_square.on_change('selected', update_sdosubmp_image)

sdosubmp_pfmap = PuffinMap(smap=sdosubmp,
                           plot_height=config_plot['plot_config']['tab_MkPlot']['aia_submap_hght'],
                           plot_width=config_plot['plot_config']['tab_MkPlot']['aia_submap_wdth'])
p_sdosubmp, r_sdosubmp = sdosubmp_pfmap.PlotMap(DrawLimb=False, DrawGrid=False, ignore_coord=True)
mapx_sdosubmp, mapy_sdosubmp = sdosubmp_pfmap.meshgridpix(rescale=0.125)
mapx_sdosubmp, mapy_sdosubmp = mapx_sdosubmp.value, mapy_sdosubmp.value
ImgDF0_sdosubmp_canvas = pd.DataFrame({'xx': mapx_sdosubmp.ravel(), 'yy': mapy_sdosubmp.ravel()})
tapPointDF_sdosubmp_canvas = pd.DataFrame({'xx': [], 'yy': []})  ## the selected point to fit in aia submap
SRC_sdosubmp_canvas_square = ColumnDataSource(ImgDF0_sdosubmp_canvas)
r_sdosubmp_square = p_sdosubmp.square('xx', 'yy', source=SRC_sdosubmp_canvas_square,
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
p_sdosubmp.add_tools(TapTool(renderers=[r_sdosubmp_square]))
SRC_sdosubmp_line = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
SRC_sdosubmp_line0 = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
SRC_sdosubmp_line1 = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []}))
r_sdosubmp_line = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line, line_color='white', line_width=2,
                                  alpha=0.7)
r_sdosubmp_line0 = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line0, line_color='white', line_width=1.5,
                                   alpha=0.7, line_dash='dashed')
r_sdosubmp_line1 = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line1, line_color='white', line_width=1.5,
                                   alpha=0.7, line_dash='dashed')
r_sdosubmp_cross = p_sdosubmp.cross(x='xx', y='yy', size=15, line_width=2, alpha=0.7, line_color='firebrick',
                                    source=ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})))

TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
p_lighcurve = figure(tools=TOOLS,
                     plot_width=config_plot['plot_config']['tab_MkPlot']['time_plot_wdth'],
                     plot_height=config_plot['plot_config']['tab_MkPlot']['time_plot_hght'],
                     toolbar_location="above")
p_lighcurve.axis.visible = True
p_lighcurve.title.text = "light Curve"
p_lighcurve.xaxis.axis_label = 'pixel'
p_lighcurve.yaxis.axis_label = 'Flux'
p_lighcurve.border_fill_alpha = 0.4
p_lighcurve.axis.major_tick_out = 0
p_lighcurve.axis.major_tick_in = 5
p_lighcurve.axis.minor_tick_out = 0
p_lighcurve.axis.minor_tick_in = 3
SRC_lightcurve = ColumnDataSource({'x': [], 'y': []})
r_lighcurve = p_lighcurve.line(x='x', y='y', alpha=0.8, line_width=1, line_color='black', source=SRC_lightcurve)

BUT_ClearDraw = Button(label='ClearDraw', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                       button_type='primary')
BUT_UndoDraw = Button(label='Undo', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                      button_type='warning')
BUT_MkStackplt = Button(label='MkStackPlt', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                        button_type='success')

buttons_play = DButil.buttons_play(plot_width=config_plot['plot_config']['tab_MkPlot']['button_wdth_short'])
for ctrl in buttons_play.buttons:
    ctrl.on_click(buttons_play.play)

Text_Cutwdth = TextInput(value='5.0', title="Cut Width (pix):")
Text_CutAng = TextInput(value='10.0', title="Cut Angle (deg):")
Text_SlitLgth = TextInput(value='{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4),
                          title="Slit Length (pix):")

Hideitemlist = ["Hide Slit&Points", "Hide Points only"]
Hideitem_checkbox = CheckboxGroup(labels=Hideitemlist, active=[])
fitmethoddict = {'0': "Polyfit", '1': "Spline", '2': "Param_Spline"}
fitMeth_radiogroup = RadioGroup(labels=["Polyfit", "Spline", "Parametric Spline"], active=0)
fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
ascending = True
clearpoint = False


def sdosubmp_square_selection_change(attrname, old, new):
    global tapPointDF_sdosubmp_canvas, ascending, cutslitplt
    sdosubmp_square_selected = SRC_sdosubmp_canvas_square.selected['1d']['indices']
    if sdosubmp_square_selected:
        ImgDF_sdosubmp = ImgDF0_sdosubmp_canvas.iloc[sdosubmp_square_selected, :]
        DFidx_selected = ImgDF_sdosubmp.index[0]
        tapPointDF_sdosubmp_canvas = tapPointDF_sdosubmp_canvas.append(ImgDF0_sdosubmp_canvas.loc[DFidx_selected, :],
                                                                       ignore_index=True)
        if len(tapPointDF_sdosubmp_canvas.index) == 2:
            if tapPointDF_sdosubmp_canvas.loc[tapPointDF_sdosubmp_canvas.index[0], 'xx'] > \
                    tapPointDF_sdosubmp_canvas.loc[
                        tapPointDF_sdosubmp_canvas.index[1], 'xx']:
                ascending = False
        cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
        pushslit2plt(clearpoint=clearpoint)


SRC_sdosubmp_canvas_square.on_change('selected', sdosubmp_square_selection_change)
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
    if not ascending and (fitmethod != 'Param_Spline' or len(xx) <= 3):
        xs, ys = xs[::-1], ys[::-1]
        posangs = posangs[::-1]
    dist = DButil.findDist(xs, ys)
    dists = np.cumsum(dist)
    posangs2 = posangs + np.pi / 2
    cutwidths = dists * np.tan(cutang) + cutwidth
    xs0 = xs - cutwidths / 2. * np.cos(posangs2)
    ys0 = ys - cutwidths / 2. * np.sin(posangs2)
    xs1 = xs + cutwidths / 2. * np.cos(posangs2)
    ys1 = ys + cutwidths / 2. * np.sin(posangs2)
    return {'xcen': xs, 'ycen': ys, 'xs0': xs0, 'ys0': ys0, 'xs1': xs1, 'ys1': ys1, 'cutwidth': cutwidths,
            'posang': posangs, 'posangs': posangs2, 'dist': dists}


def MakeSlit(pointDF):
    global smoth_factor0, smoth_factor1, fitmethod
    if fitmethod != 'Param_Spline' or len(pointDF.index) <= 3:
        pointDFtmp = pointDF.sort_values(['xx'], ascending=[True])
    else:
        pointDFtmp = pointDF
    xx = pointDFtmp.loc[:, 'xx'].values
    yy = pointDFtmp.loc[:, 'yy'].values
    if len(pointDFtmp.index) == 1:
        r_sdosubmp_cross.data_source.data = ColumnDataSource(pd.DataFrame(pointDFtmp)).data
        cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                      'posang': [], 'posangs': [], 'dist': []}
    else:
        if len(pointDFtmp.index) <= 3:
            cutslitplt = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                                 float(Text_SlitLgth.value), method=fitmethod)
        else:
            smoth_factor0 = len(yy) - np.sqrt(2 * len(yy))
            smoth_factor1 = np.var(yy) * (len(yy)) * 3
            cutslitplt = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                                 float(Text_SlitLgth.value), s=np.exp(
                    Slider_smoothing_factor.value * (np.log(smoth_factor1) - np.log(smoth_factor0))) * smoth_factor0,
                                 method=fitmethod)
    return cutslitplt


def pushslit2plt(clearpoint=False):
    r_sdosubmp_line.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslitplt['xcen'], 'yy': cutslitplt['ycen']})).data
    r_sdosubmp_line0.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslitplt['xs0'], 'yy': cutslitplt['ys0']})).data
    r_sdosubmp_line1.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': cutslitplt['xs1'], 'yy': cutslitplt['ys1']})).data
    if clearpoint:
        r_sdosubmp_cross.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    else:
        r_sdosubmp_cross.data_source.data = ColumnDataSource(pd.DataFrame(tapPointDF_sdosubmp_canvas)).data


def clearslitplt():
    r_sdosubmp_line.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_sdosubmp_line0.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_sdosubmp_line1.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data
    r_sdosubmp_cross.data_source.data = ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})).data


def slider_smoothing_factor_update(attrname, old, new):
    global tapPointDF_sdosubmp_canvas, smoth_factor0, smoth_factor1, cutslitplt
    # if len(pointDF_sdosubmp.index) > 3 and fitmethod == 'Polyfit':
    #     fitMeth_radiogroup.active = 1
    if len(tapPointDF_sdosubmp_canvas.index) > 3:
        if fitmethod == 'Polyfit':
            fitMeth_radiogroup.active = 1
        if fitmethod == 'Spline':
            pointDFtmp = tapPointDF_sdosubmp_canvas.sort_values(['xx'], ascending=[True])
        else:
            pointDFtmp = tapPointDF_sdosubmp_canvas
        xx = pointDFtmp.loc[:, 'xx']
        yy = pointDFtmp.loc[:, 'yy']
        cutslitplt = FitSlit(xx, yy, float(Text_Cutwdth.value) / Canvas_scale, np.radians(float(Text_CutAng.value)),
                             float(Text_SlitLgth.value), s=np.exp(
                Slider_smoothing_factor.value * (np.log(smoth_factor1) - np.log(smoth_factor0))) * smoth_factor0,
                             method=fitmethod)
        pushslit2plt(clearpoint=clearpoint)


Slider_smoothing_factor.on_change('value', slider_smoothing_factor_update)


def UndoDraw():
    global tapPointDF_sdosubmp_canvas, smoth_factor0, smoth_factor1, cutslitplt
    if len(tapPointDF_sdosubmp_canvas.index) > 1:
        tapPointDF_sdosubmp_canvas.drop(len(tapPointDF_sdosubmp_canvas.index) - 1, inplace=True)
        cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
        pushslit2plt(clearpoint=clearpoint)
    else:
        ClearDraw()


BUT_UndoDraw.on_click(UndoDraw)


def ClearDraw():
    global tapPointDF_sdosubmp_canvas, ascending
    Slider_smoothing_factor.value = 1.0
    ascending = True
    tapPointDF_sdosubmp_canvas = pd.DataFrame({'xx': [], 'yy': []})
    clearslitplt()


BUT_ClearDraw.on_click(ClearDraw)


def MkStackplt():
    num = len(cutslitplt['xcen'])
    intens = np.zeros(num)
    for ll in xrange(num):
        xs0, xs1 = cutslitplt['xs0'][ll], cutslitplt['xs1'][ll]
        ys0, ys1 = cutslitplt['ys0'][ll], cutslitplt['ys1'][ll]
        x, y = DButil.canvaspix_to_data(sdosubmp, np.array([xs0, xs1]) * Canvas_scale,
                                        np.array([ys0, ys1]) * Canvas_scale)
        xpix, ypix = DButil.data_to_mappixel(sdosubmp, x, y)
        inten = DButil.improfile(sdosubmp.data, xpix, ypix, interp='nearest')
        intens[ll] = np.mean(inten)
    r_lighcurve.data_source.data = {'x': cutslitplt['dist'], 'y': intens}
    p_lighcurve.x_range = Range1d(0, cutslitplt['dist'][-1])
    p_lighcurve.y_range = Range1d(0, np.max(intens))



BUT_MkStackplt.on_click(MkStackplt)


def HideItemUpdate(new):
    global clearpoint
    clearpoint = len(Hideitem_checkbox.active) == 1 and Hideitemlist[Hideitem_checkbox.active[0]] == "Hide Points only"
    if len(Hideitem_checkbox.active) <= 1:
        if len(Hideitem_checkbox.active) == 1 and Hideitemlist[Hideitem_checkbox.active[0]] == "Hide Slit&Points":
            clearslitplt()
        else:
            pushslit2plt(clearpoint=clearpoint)
    else:
        clearslitplt()


Hideitem_checkbox.on_click(HideItemUpdate)


def fitMethUpdate(new):
    global fitmethod, cutslitplt
    fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
    cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
    pushslit2plt(clearpoint=clearpoint)


fitMeth_radiogroup.on_click(fitMethUpdate)

Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b>
            first before closing the tab</p></b>""", width=config_plot['plot_config']['tab_MkPlot']['button_wdth'])


def default_fitparam():
    Text_Cutwdth.value = '5.0'
    Text_CutAng.value = '10.0'
    Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)


BUT_default_fitparam = Button(label='Default', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'],
                              button_type='success')
BUT_default_fitparam.on_click(default_fitparam)


def exit():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


BUT_exit = Button(label='Exit', width=config_plot['plot_config']['tab_MkPlot']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit)

# SPCR_LFT_lbutt_play = Spacer(width=100)
# lbutt_play = row(BUT_first, BUT_prev, BUT_play, BUT_next, BUT_last)
lbutt_play = row(buttons_play.buttons[0], buttons_play.buttons[1], buttons_play.buttons[2], buttons_play.buttons[3],
                 buttons_play.buttons[4])
lout = row(column(p_sdomap, lbutt_play, p_lighcurve), p_sdosubmp,
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
