import os
import astropy.units as u
import bokeh.palettes as bp
import pickle
import numpy as np
import pandas as pd
import sunpy.map
from astropy.time import Time
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (ColumnDataSource, Slider, Button, TextInput, CheckboxGroup, CheckboxButtonGroup, RadioGroup,
                          BoxSelectTool, TapTool, Div, Spacer, Range1d)
from bokeh.models.widgets import Dropdown, RangeSlider, Select, Panel, Tabs
from bokeh.plotting import figure, curdoc
from sunpy.time import TimeRange
from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
from suncasa.utils.puffin import PuffinMap
from suncasa.utils import DButil
import Tkinter
import tkFileDialog

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"


def update_sdosubmp_image(fileidx):
    global sdomap, sdosubmp
    if not sdosubmpdict:
        sdomap = sunpy.map.Map(sdofile[fileidx])
        sdomap = DButil.normalize_aiamap(sdomap)
        sdosubmp = sdomap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                 u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
    else:
        sdosubmp = sdosubmpdict['submc'].maps[fileidx]
    r_sdosubmp.data_source.data['image'] = [sdosubmp.data]
    p_sdosubmp.title.text = sdosubmp.name
    try:
        getimprofile(sdosubmp, cutslitplt)
    except:
        ValueError('make a cut slit first!')


def aiamap_wavelength_selection(attrname, old, new):
    global select_wave, sdofile
    select_wave = Select_aia_wave.value
    MkPlot_args_dict['wavelength'] = select_wave
    outfile = database_dir + 'MkPlot_args.json'
    DButil.updatejsonfile(outfile, MkPlot_args_dict)
    try:
        SaveSlit(slitfile)
    except:
        pass
    Div_info.text = """<p><b>Refresh</b> the <b>web page</b> to load new wavelength of AIA</p>"""


def Running_diff_img(sdompdict, ratio=False):
    global datastep
    datastep = Slider_datadt.value / AIAcadence
    maps = sdosubmpdict['submc'].maps
    for sidx, smap in enumerate(maps):
        if sidx < datastep:
            sdompdict['mask'][sidx] = False
        else:
            sdompdict['mask'][sidx] = True
            if ratio:
                sdompdict['submc'].maps[sidx].data = smap.data / maps[sidx - datastep].data
            else:
                sdompdict['submc'].maps[sidx].data = smap.data - maps[sidx - datastep].data
        Div_info.text = """<p>{}</p>""".format(
            DButil.ProgressBar(sidx + 1, nsdofile, suffix='Update', decimals=0, length=14, empfill='=', fill='#'))
    return sdompdict


def Base_diff_img(sdompdict, ratio=False):
    maps = sdosubmpdict['submc'].maps
    for sidx, smap in enumerate(maps):
        if sidx == 0:
            sdompdict['mask'][sidx] = False
        else:
            sdompdict['mask'][sidx] = True
            if ratio:
                sdompdict['submc'].maps[sidx].data = smap.data / maps[0].data
            else:
                sdompdict['submc'].maps[sidx].data = smap.data - maps[0].data
        Div_info.text = """<p>{}</p>""".format(
            DButil.ProgressBar(sidx + 1, nsdofile, suffix='Update', decimals=0, length=14, empfill='=', fill='#'))
    return sdompdict


def DiffImg_update():
    global sdosubmpdict
    if sdosubmpdict:
        if Select_DiffImg.value == 'No diff images':
            LoadChunk_handler()
        else:
            UpdateSubChunk()
            if Select_DiffImg.value == 'Running diff':
                sdosubmpdict = Running_diff_img(sdosubmpdict, ratio=False)
            elif Select_DiffImg.value == 'Base diff':
                sdosubmpdict = Base_diff_img(sdosubmpdict, ratio=False)
            elif Select_DiffImg.value == 'Running diff/Ratio':
                sdosubmpdict = Running_diff_img(sdosubmpdict, ratio=True)
            elif Select_DiffImg.value == 'Base diff/Ratio':
                sdosubmpdict = Base_diff_img(sdosubmpdict, ratio=True)
            trange_update(updatemask=False)
    else:
        Div_info.text = """<p>load the <b>chunk</b> first!!!</p>"""


def Select_DiffImg_update(attrname, old, new):
    global imagetype
    imagetype = DiffImglabelsdict_r[Select_DiffImg.value]
    if imagetype == 'image':
        LoadChunk_handler()
    else:
        DiffImg_update()
    update_sdosubmp_image(Slider_sdoidx.value - 1)
    MkPlot_args_dict['imagetype'] = DiffImglabelsdict_r[Select_DiffImg.value]
    outfile = database_dir + 'MkPlot_args.json'
    DButil.updatejsonfile(outfile, MkPlot_args_dict)
    Div_info.text = """<p><b>Refresh</b> the <b>web page</b> to load diff images</p>"""


def Slider_datadt_update(attrname, old, new):
    if Select_DiffImg.value != 'No diff images':
        if sdosubmpdict:
            BUT_loadchunk.label = 'UpdateChunk'
    else:
        Slider_datadt.value = Slider_datadt.start
        # DiffImg_update()
        # update_sdosubmp_image(Slider_sdoidx.value - 1)


def update_sdosubmp_region(x0, x1, y0, y1):
    global sdosubmp
    global xaxis_sdosubmp, yaxis_sdosubmp
    global Canvas_scale, Canvas_scaleX, Canvas_scaleY
    r_sdo_RSPmap_square_patch.data_source.data = ColumnDataSource(
        pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})).data
    sdosubmp = sdomap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                             u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
    xaxis_sdosubmp, yaxis_sdosubmp = sdosubmp.data.shape
    Canvas_scaleX, Canvas_scaleY = float(xaxis_sdosubmp) / xaxis_sdosubmp_canvas, float(
        yaxis_sdosubmp) / yaxis_sdosubmp_canvas
    Canvas_scale = (Canvas_scaleX + Canvas_scaleY) / 2.0
    r_sdosubmp.data_source.data['image'] = [sdosubmp.data]
    try:
        if x0 >= sdosubmpdict['FOV'][0] and x1 <= sdosubmpdict['FOV'][2] and y0 >= sdosubmpdict['FOV'][
            1] and y1 <= \
                sdosubmpdict['FOV'][3]:
            UpdateSubChunk()
            trange_update()
            BUT_loadchunk.label = 'UpdateChunk'
        else:
            BUT_loadchunk.label = 'LoadChunk'
    except:
        pass


def sdosubmp_region_select(attrname, old, new):
    global tapPointDF_sdo_RSPmap, sdo_RSPmap_quadselround
    global sdosubmp, x0, x1, y0, y1
    global xaxis_sdosubmp, yaxis_sdosubmp
    global Canvas_scale, Canvas_scaleX, Canvas_scaleY
    sdo_RSPmap_quadselround += 1
    if sdo_RSPmap_quadselround == 2:
        sdo_RSPmap_quadx_selected = SRC_sdo_RSPmap_quadx.selected['1d']['indices']
        sdo_RSPmap_quady_selected = SRC_sdo_RSPmap_quady.selected['1d']['indices']
        SRC_sdo_RSPmap_quadx.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        SRC_sdo_RSPmap_quady.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        if len(sdo_RSPmap_quadx_selected) > 0 and len(sdo_RSPmap_quady_selected) > 0:
            tapPointDF_sdo_RSPmap = tapPointDF_sdo_RSPmap.append(
                pd.Series({'xx': mapx_sdo_RSPmap[sdo_RSPmap_quady_selected[0], sdo_RSPmap_quadx_selected[0]],
                           'yy': mapy_sdo_RSPmap[sdo_RSPmap_quady_selected[0], sdo_RSPmap_quadx_selected[0]]}),
                ignore_index=True)
            if len(tapPointDF_sdo_RSPmap.index) == 1:
                r_sdo_RSPmap_line.data_source.data = {
                    'xs': [[mapx_sdo_RSPmap[sdo_RSPmap_quady_selected[0], sdo_RSPmap_quadx_selected[0]]] * 2,
                           [mapx_sdo_RSPmap[sdo_RSPmap_quady_selected[0], 0],
                            mapx_sdo_RSPmap[sdo_RSPmap_quady_selected[0], -1]]],
                    'ys': [[mapy_sdo_RSPmap[0, sdo_RSPmap_quadx_selected[0]],
                            mapy_sdo_RSPmap[-1, sdo_RSPmap_quadx_selected[0]]],
                           [mapy_sdo_RSPmap[sdo_RSPmap_quady_selected[0], sdo_RSPmap_quadx_selected[0]]] * 2]}
            elif len(tapPointDF_sdo_RSPmap.index) == 2:
                x0, x1 = tapPointDF_sdo_RSPmap['xx'].min(), tapPointDF_sdo_RSPmap['xx'].max()
                y0, y1 = tapPointDF_sdo_RSPmap['yy'].min(), tapPointDF_sdo_RSPmap['yy'].max()
                if x1 > x0 + mapx_sdo_RSPmap[0, 1] - mapx_sdo_RSPmap[0, 0] and y1 > y0 + mapy_sdo_RSPmap[0, 1] - \
                        mapy_sdo_RSPmap[0, 0]:
                    ## select at least 4 points
                    patch_size = min(max(x1 - x0, y1 - y0),
                                     config_main['plot_config']['tab_MkPlot']['aia_submap_sz_max'])
                    x1, y0 = x0 + patch_size, y1 - patch_size
                    update_sdosubmp_region(x0, x1, y0, y1)
                    Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)
                    ClearDraw()
                tapPointDF_sdo_RSPmap = pd.DataFrame({'xx': [], 'yy': []})
                r_sdo_RSPmap_line.data_source.data = {'xs': [], 'ys': []}
        sdo_RSPmap_quadselround = 0


def ImgOption_checkbox_handler(new):
    LoadChunk_handler()


def LoadChunk():
    sdosubmplist = []
    timestamps = []
    for sidx, sfile in enumerate(sdofile):
        sdomaptmp = DButil.normalize_aiamap(sunpy.map.Map(sfile))
        sdosubmptmp = sdomaptmp.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                       u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        sdosubmplist.append(sdosubmptmp)
        timestamps.append(Time(sdosubmptmp.meta['date-obs'].replace('T', ' '), format='iso', scale='utc').jd)
        Div_info.text = """<p>{}</p>""".format(
            DButil.ProgressBar(sidx + 1, nsdofile + 1, suffix='Load', decimals=0, length=16, empfill='=', fill='#'))
    mc = sunpy.map.Map(sdosubmplist, cube=True)
    sdompdict = {'FOV': [x0, y0, x1, y1], 'subFOV': [x0, y0, x1, y1], 'mc': mc, 'time': np.array(timestamps)}
    sdompdict['submc'] = LoadSubChunk(sdompdict)
    mask = np.logical_and(sdompdict['time'] >= trange.jd[0] + Slider_trange.range[0] / 24. / 3600,
                          sdompdict['time'] <= trange.jd[0] + Slider_trange.range[1] / 24. / 3600)
    sdompdict['mask'] = mask
    Div_info.text = """<p>{}</p>""".format(
        DButil.ProgressBar(nsdofile + 1, nsdofile + 1, suffix='Load', decimals=0, length=16, empfill='=', fill='#'))
    return sdompdict


def LoadSubChunk(sdompdict):
    sdosubmplist = []
    for sidx, smap in enumerate(sdompdict['mc'].maps):
        sdosubmptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                  u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        if ImgOption_checkbox.active == [1]:
            if sidx == 0:
                grid = DButil.smapmeshgrid2(sdosubmptmp)
            sdosubmptmp = DButil.smapradialfilter(sdosubmptmp, grid=grid)
        sdosubmplist.append(sdosubmptmp)
        Div_info.text = """<p>{}</p>""".format(
            DButil.ProgressBar(sidx + 1, nsdofile + 1, suffix='Update', decimals=0, length=14, empfill='=', fill='#'))
    if 0 in ImgOption_checkbox.active:
        submc = mapcube_solar_derotate(sunpy.map.Map(sdosubmplist, cube=True))
    else:
        submc = sunpy.map.Map(sdosubmplist, cube=True)
    # if imagetype != 'image':
    #     DiffImg_update()
    Div_info.text = """<p>{}</p>""".format(
        DButil.ProgressBar(nsdofile + 1, nsdofile + 1, suffix='Update', decimals=0, length=14, empfill='=', fill='#'))
    return submc


def UpdateSubChunk():
    global sdosubmpdict
    sdosubmpdict['submc'] = LoadSubChunk(sdosubmpdict)
    sdosubmpdict['subFOV'] = [x0, y0, x1, y1]
    sdosubmpdict['mask'] = np.logical_and(
        sdosubmpdict['time'] >= trange.jd[0] + Slider_trange.range[0] / 24. / 3600,
        sdosubmpdict['time'] <= trange.jd[0] + Slider_trange.range[1] / 24. / 3600)


def LoadChunk_handler():
    global sdosubmpdict
    if sdosubmpdict:
        if x0 >= sdosubmpdict['FOV'][0] and x1 <= sdosubmpdict['FOV'][2] and y0 >= sdosubmpdict['FOV'][
            1] and y1 <= \
                sdosubmpdict['FOV'][3]:
            UpdateSubChunk()
        else:
            sdosubmpdict = LoadChunk()
    else:
        sdosubmpdict = LoadChunk()
    print imagetype
    trange_reset()
    if imagetype != 'image':
        DiffImg_update()
        trange_update(updatemask=False)
    else:
        trange_update()
    update_sdosubmp_image(Slider_sdoidx.value - 1)


def Slider_sdoidxstep_update(attrname, old, new):
    global sdofileidxstep
    sdofileidxstep = Slider_sdoidxstep.value
    Slider_sdoidx.step = sdofileidxstep


def Slider_sdoidx_update(attrname, old, new):
    global sdofileidx
    update_sdosubmp_image(Slider_sdoidx.value - 1)
    sdofileidx = Slider_sdoidx.value - 1


def ButtonNext_handler():
    global sdofileidx
    sdofileidx += sdofileidxstep
    if sdofileidx > sdofileinbound[1] or sdofileidx < sdofileinbound[0]:
        sdofileidx = sdofileinbound[0]
    sdofileidx = sdofileidx % nsdofile
    Slider_sdoidx.value = sdofileidx + 1


def ButtonPrev_handler():
    global sdofileidx
    sdofileidx -= sdofileidxstep
    if sdofileidx < sdofileinbound[0] or sdofileidx > sdofileinbound[1]:
        sdofileidx = sdofileinbound[1]
    sdofileidx = sdofileidx % nsdofile
    Slider_sdoidx.value = sdofileidx + 1


def ButtonLast_handler():
    global sdofileidx
    sdofileidx = sdofileinbound[1]
    sdofileidx = sdofileidx % nsdofile
    Slider_sdoidx.value = sdofileidx + 1


def ButtonFirst_handler():
    global sdofileidx
    sdofileidx = sdofileinbound[0]
    sdofileidx = sdofileidx % nsdofile
    Slider_sdoidx.value = sdofileidx + 1


def ButtonPlay_handler():
    global sdofileidx
    if ButtonsPlayCTRL.buttons[2].label == '>':
        if sdosubmpdict:
            ButtonsPlayCTRL.buttons[2].label = '||'
            curdoc().add_periodic_callback(sdosubmp_animate_update, 150)
        else:
            Div_info.text = """<p>load the <b>chunk</b> first!!!</p>"""
    else:
        ButtonsPlayCTRL.buttons[2].label = '>'
        curdoc().remove_periodic_callback(sdosubmp_animate_update)


def sdosubmp_animate_update():
    if sdosubmpdict['submc']:
        ButtonNext_handler()


def sdosubmp_quad_selection_change(attrname, old, new):
    global tapPointDF_sdosubmp_canvas, ascending, cutslitplt, sdosubmp_quadselround
    sdosubmp_quadselround += 1
    if sdosubmp_quadselround == 2:
        sdosubmp_quadx_selected = SRC_sdosubmp_canvas_quadx.selected['1d']['indices']
        sdosubmp_quady_selected = SRC_sdosubmp_canvas_quady.selected['1d']['indices']
        SRC_sdosubmp_canvas_quadx.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        SRC_sdosubmp_canvas_quady.selected = {'0d': {'glyph': None, 'indices': []}, '1d': {'indices': []}, '2d': {}}
        if len(sdosubmp_quadx_selected) > 0 and len(sdosubmp_quady_selected) > 0:
            tapPointDF_sdosubmp_canvas = tapPointDF_sdosubmp_canvas.append(
                pd.Series({'xx': mapx_sdosubmp[sdosubmp_quady_selected[0], sdosubmp_quadx_selected[0]],
                           'yy': mapy_sdosubmp[sdosubmp_quady_selected[0], sdosubmp_quadx_selected[0]]}),
                ignore_index=True)
            if len(tapPointDF_sdosubmp_canvas.index) == 2:
                if tapPointDF_sdosubmp_canvas.loc[tapPointDF_sdosubmp_canvas.index[0], 'xx'] > \
                        tapPointDF_sdosubmp_canvas.loc[
                            tapPointDF_sdosubmp_canvas.index[1], 'xx']:
                    ascending = False
            cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
            getimprofile(sdosubmp, cutslitplt)
            pushslit2plt(cutslitplt, clearpoint=clearpoint)
        sdosubmp_quadselround = 0


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
            'posangs': posangs, 'posangs2': posangs2, 'dist': dists}


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
                      'posangs': [], 'posangs2': [], 'dist': []}
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


def pushslit2plt(cutslit, clearpoint=False):
    r_sdosubmp_line.data_source.data = {'xx': cutslit['xcen'], 'yy': cutslit['ycen']}
    r_sdosubmp_line0.data_source.data = {'xx': cutslit['xs0'], 'yy': cutslit['ys0']}
    r_sdosubmp_line1.data_source.data = {'xx': cutslit['xs1'], 'yy': cutslit['ys1']}
    if clearpoint:
        r_sdosubmp_cross.data_source.data = {'xx': [], 'yy': []}
    else:
        r_sdosubmp_cross.data_source.data = ColumnDataSource(tapPointDF_sdosubmp_canvas).data


def clearslitplt():
    r_sdosubmp_line.data_source.data = {'xx': [], 'yy': []}
    r_sdosubmp_line0.data_source.data = {'xx': [], 'yy': []}
    r_sdosubmp_line1.data_source.data = {'xx': [], 'yy': []}
    r_sdosubmp_cross.data_source.data = {'xx': [], 'yy': []}


def slider_smoothing_factor_update(attrname, old, new):
    global tapPointDF_sdosubmp_canvas, smoth_factor0, smoth_factor1, cutslitplt
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
        pushslit2plt(cutslitplt, clearpoint=clearpoint)


def UndoDraw():
    global tapPointDF_sdosubmp_canvas, smoth_factor0, smoth_factor1, cutslitplt
    if len(tapPointDF_sdosubmp_canvas.index) > 1:
        tapPointDF_sdosubmp_canvas.drop(len(tapPointDF_sdosubmp_canvas.index) - 1, inplace=True)
        cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
        getimprofile(sdosubmp, cutslitplt)
        pushslit2plt(cutslitplt, clearpoint=clearpoint)
        if len(tapPointDF_sdosubmp_canvas.index) == 1:
            r_lighcurve.data_source.data = {'x': [], 'y': []}
    else:
        ClearDraw()


def ClearDraw():
    global tapPointDF_sdosubmp_canvas, ascending, cutslitplt
    Slider_smoothing_factor.value = 1.0
    ascending = True
    tapPointDF_sdosubmp_canvas = pd.DataFrame({'xx': [], 'yy': []})
    clearslitplt()
    cutslitplt = {'xcen': [], 'ycen': [], 'xs0': [], 'ys0': [], 'xs1': [], 'ys1': [], 'cutwidth': [],
                  'posangs': [], 'posangs2': [], 'dist': []}
    r_lighcurve.data_source.data = {'x': [], 'y': []}


def getimprofile(smap, cutslit, plot=True):
    num = len(cutslit['xcen'])
    if num > 1:
        ctslit = cutslit.copy()
        for key in ctslit.keys():
            if not key in ['posangs', 'posangs2']:
                if 'x' in key:
                    ctslit[key] = ctslit[key] * Canvas_scaleX
                elif 'y' in key:
                    ctslit[key] = ctslit[key] * Canvas_scaleY
                else:
                    ctslit[key] = ctslit[key] * Canvas_scale
        intens = np.zeros(num)
        for ll in xrange(num):
            # x, y = DButil.canvaspix_to_data(sdosubmp, np.array([xs0, xs1]),
            #                                 np.array([ys0, ys1]))
            # xpix, ypix = DButil.data_to_mappixel(sdosubmp, x, y)
            inten = DButil.improfile(smap.data, [ctslit['xs0'][ll], ctslit['xs1'][ll]],
                                     [ctslit['ys0'][ll], ctslit['ys1'][ll]], interp='nearest')
            intens[ll] = np.mean(inten)
        intens = intens / smap.exposure_time.value
        intensdist = {'x': ctslit['dist'] * np.mean([smap.scale[0].value, smap.scale[1].value]), 'y': intens}
        if plot:
            r_lighcurve.data_source.data = intensdist
        return intensdist
        # p_lighcurve.x_range = Range1d(0, cutslit['dist'][-1])
        # p_lighcurve.y_range = Range1d(0, np.max(intens))


def HideItemUpdate(new):
    global clearpoint
    clearpoint = len(Hideitem_checkbox.active) == 1 and Hideitemlist[Hideitem_checkbox.active[0]] == "Hide Points only"
    if len(Hideitem_checkbox.active) <= 1:
        if len(Hideitem_checkbox.active) == 1 and Hideitemlist[Hideitem_checkbox.active[0]] == "Hide Slit&Points":
            clearslitplt()
        else:
            pushslit2plt(cutslitplt, clearpoint=clearpoint)
    else:
        clearslitplt()


def fitMethUpdate(new):
    global fitmethod, cutslitplt
    fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
    cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
    getimprofile(sdosubmp, cutslitplt)
    pushslit2plt(cutslitplt, clearpoint=clearpoint)


def default_fitparam():
    Text_Cutwdth.value = '5.0'
    Text_CutAng.value = '0.0'
    Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)


def LoadSlit(fin):
    global cutslitplt, tapPointDF_sdosubmp_canvas, ascending
    global fitmethod, smoth_factor1, smoth_factor2
    global x0, x1, y0, y1
    if os.path.exists(fin):
        cutslitdict = pickle.load(open(fin, 'rb'))
        tapPointDF_sdosubmp_canvas = cutslitdict['tapPointDF_sdosubmp_canvas']
        ascending = cutslitdict['ascending']
        fitmethod = cutslitdict['fitmethod']
        smoth_factor1 = cutslitdict['smoth_factor1']
        smoth_factor1 = cutslitdict['smoth_factor2']
        smoth_factor = cutslitdict['smoth_factor']
        Cutwdth = cutslitdict['Cutwdth']
        CutAng = cutslitdict['CutAng']
        SlitLgth = cutslitdict['SlitLgth']
        x0 = cutslitdict['x0']
        x1 = cutslitdict['x1']
        y0 = cutslitdict['y0']
        y1 = cutslitdict['y1']
        # for key, val in cutslitdict.items():
        #     exec (key + '=val')
        Text_Cutwdth.value = Cutwdth
        Text_CutAng.value = CutAng
        Text_SlitLgth.value = SlitLgth
        fitMeth_radiogroup.active = int(fitmethoddict.keys()[fitmethoddict.values().index(fitmethod)])
        Slider_smoothing_factor.value = smoth_factor
        update_sdosubmp_region(x0, x1, y0, y1)
        cutslitplt = MakeSlit(tapPointDF_sdosubmp_canvas)
        getimprofile(sdosubmp, cutslitplt)
        pushslit2plt(cutslitplt, clearpoint=clearpoint)
        Div_info.text = """<p><b>cutslit settings</b> in {} <b>loaded</b>.</p>""".format(fin)
    else:
        Div_info.text = """<p>No cutslit file found!!</p>"""


def SaveSlit(fout):
    if cutslitplt:
        cutslitdict = {'tapPointDF_sdosubmp_canvas': tapPointDF_sdosubmp_canvas,
                       'ascending': ascending, 'fitmethod': fitmethod, 'smoth_factor1': smoth_factor1,
                       'smoth_factor2': smoth_factor1, 'smoth_factor': Slider_smoothing_factor.value,
                       'Cutwdth': Text_Cutwdth.value, 'CutAng': Text_CutAng.value,
                       'SlitLgth': Text_SlitLgth.value, 'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1}
        with open(fout, 'wb') as fp:
            pickle.dump(cutslitdict, fp)
        Div_info.text = """<p><b>cutslit settings saved</b> to {}.</p>""".format(fout)
    else:
        Div_info.text = """<p>make a slit first!!!</p>"""


def MkStackplt():
    global stackpltdict
    if sdosubmpdict:
        if len(cutslitplt['xcen']) > 0:
            stackplt = []
            for sidx, smap in enumerate(sdosubmpdict['submc'].maps):
                if sdosubmpdict['mask'][sidx]:
                    intens = getimprofile(smap, cutslitplt, plot=False)
                    stackplt.append(intens['y'])
            if len(stackplt) > 1:
                stackplt = np.vstack(stackplt)
                stackplt = stackplt.transpose()
                stackpltdict = {'zz': stackplt, 'x': sdosubmpdict['time'][sdosubmpdict['mask']], 'y': intens['x'],
                                'wavelength': '{:.0f}'.format(smap.wavelength.value), 'observatory': smap.observatory,
                                'instrument': smap.instrument[0:3]}
                Div_info.text = """<p>click <b>Stackplt</b> to view the stack plot</p>"""
                return True
            else:
                Div_info.text = """<p>less than two frame in selected time range!!!</p>"""
                return False
        else:
            Div_info.text = """<p>make a <b>cut slit</b> first!!!</p>"""
            return False
    else:
        Div_info.text = """<p>load the <b>chunk</b> first!!!</p>"""
        return False


def ViewStackplt():
    if MkStackplt():
        stackplt_save = database_dir + 'stackplt-' + PlotID + '.npy'
        np.save(stackplt_save, stackpltdict)
        port = DButil.getfreeport()
        print 'bokeh serve {}aiaBrowser/StackPlt --show --port {} &'.format(suncasa_dir, port)
        os.system('bokeh serve {}aiaBrowser/StackPlt --show --port {} &'.format(suncasa_dir, port))
        ports.append(port)
        Div_info.text = """<p><p>Check the <b>StackPlt</b> in the <b>new tab</b>.</p>"""


def trange_reset():
    global trangesec
    dur = (sdosubmpdict['time'][-1] - sdosubmpdict['time'][0]) * 24 * 3600
    deltat = float('{:.1f}'.format((sdosubmpdict['time'][0] - trange.jd[0]) * 24. * 3600.))
    trangesec = (deltat, float('{:.1f}'.format(dur)) + deltat)
    Slider_trange.start = trangesec[0]
    Slider_trange.end = trangesec[1]
    Slider_trange.range = trangesec
    # Slider_trange.title = 'Time range [seconds since {}]'.format(
    #     Time(sdosubmpdict['time'][0], format='jd', scale='utc').iso)
    # update of the title doesn't work


def trange_update(updatemask=True):
    global sdofileinbound, nsdofileinbound
    if sdosubmpdict:
        if updatemask:
            sdosubmpdict['mask'] = np.logical_and(
                sdosubmpdict['time'] >= trange.jd[0] + (Slider_trange.range[0] - 6.0) / 24. / 3600,
                sdosubmpdict['time'] <= trange.jd[0] + (Slider_trange.range[1] + 6.0) / 24. / 3600)
        maskidx = np.where(sdosubmpdict['mask'])[0]
        # print updatemask, sdosubmpdict['mask']
        sdofileinbound = [maskidx[0], maskidx[-1]]
        Slider_sdoidx.value = sdofileinbound[0] + 1
        # if Select_DiffImg.value != 'No diff images':
        sldertran = list(Slider_trange.range)
        sldertran[0] = trangesec[0] + sdofileinbound[0] * AIAcadence
        Slider_trange.range = tuple(sldertran)


def trange_change_handler(attr, old, new):
    trange_update()


def DropDn_slit_handler(attr, old, new):
    global slitfile
    if DropDn_slit.value == "Open":
        tkRoot = Tkinter.Tk()
        tkRoot.withdraw()  # Close the root window
        fin = tkFileDialog.askopenfilename(initialdir=database_dir, initialfile='cutslit-' + PlotID)
        LoadSlit(fin)
        slitfile = fin
    elif DropDn_slit.value == "Save As":
        tkRoot = Tkinter.Tk()
        tkRoot.withdraw()  # Close the root window
        fout = tkFileDialog.asksaveasfilename(initialdir=database_dir, initialfile='cutslit-' + PlotID)
        SaveSlit(fout)
        slitfile = fout
    elif DropDn_slit.value == "Load":
        LoadSlit(slitfile)
    elif DropDn_slit.value == "Save":
        SaveSlit(slitfile)
    else:
        pass


def exit_update():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


ports = []
'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
DButil.initconfig(suncasa_dir)
'''load config file'''
config_main = DButil.loadjsonfile(suncasa_dir + 'DataBrowser/config.json')
database_dir = config_main['datadir']['database']
database_dir = os.path.expandvars(database_dir) + '/aiaBrowserData/'
if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))

SDOdir = DButil.getSDOdir(config_main, database_dir, suncasa_dir)

if not os.path.exists(database_dir):
    os.system('mkdir {}'.format(database_dir))
infile = database_dir + 'MkPlot_args.json'
MkPlot_args_dict = DButil.loadjsonfile(infile)
try:
    imagetype = MkPlot_args_dict['imagetype']
except:
    imagetype = 'image'

trange = Time([MkPlot_args_dict['tst'], MkPlot_args_dict['ted']], format='iso', scale='utc')
PlotID = MkPlot_args_dict['PlotID']
sdofile = DButil.readsdofile(datadir=SDOdir, wavelength=MkPlot_args_dict['wavelength'], jdtime=trange.jd)
nsdofile = len(sdofile)
global sdofileidx, slitfile
if nsdofile == 0:
    raise SystemExit('No SDO file found at the select timestamp. Download the data with EvtBrowser first.')
sdofileidx = 0
sdofileidxstep = 1
smoth_factor0 = 0
smoth_factor1 = 0
slitfile = database_dir + 'cutslit-' + PlotID
cutslitplt = {}
sdosubmpdict = {}
AIAcadence = config_main['plot_config']['tab_MkPlot']['AIA_cadence'][MkPlot_args_dict['wavelength']]
datastep = 1
sdosubmp_quadselround = 0
sdo_RSPmap_quadselround = 0
sdomap = sunpy.map.Map(sdofile[sdofileidx])
sdomap = DButil.normalize_aiamap(sdomap)
MapRES = config_main['plot_config']['tab_MkPlot']['aia_RSPmap_RES']
dimensions = u.Quantity([MapRES, MapRES], u.pixel)
sdo_RSPmap = sdomap.resample(dimensions)
sdo_RSP_pfmap = PuffinMap(smap=sdo_RSPmap,
                          plot_height=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_hght'],
                          plot_width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
p_sdomap, r_sdomap = sdo_RSP_pfmap.PlotMap(DrawLimb=True, DrawGrid=True, grid_spacing=20 * u.deg,
                                           tools='crosshair,pan,wheel_zoom,reset,save')
mapx_sdo_RSPmap, mapy_sdo_RSPmap = sdo_RSP_pfmap.meshgrid(rescale=1.0)
mapx_sdo_RSPmap, mapy_sdo_RSPmap = mapx_sdo_RSPmap.value, mapy_sdo_RSPmap.value
ndy, ndx = mapx_sdo_RSPmap.shape
tapPointDF_sdo_RSPmap = pd.DataFrame({'xx': [], 'yy': []})  ## the selected point to fit in aia submap
dx = np.mean(np.diff(mapx_sdo_RSPmap[0, :]))
dy = np.mean(np.diff(mapy_sdo_RSPmap[:, 0]))
xLFE = mapx_sdo_RSPmap[0, :] - dx / 2.0
xRTE = np.append(xLFE[1:], xLFE[-1] + dx)
yBTE = mapy_sdo_RSPmap[:, 0] - dy / 2.0
yTPE = np.append(yBTE[1:], yBTE[-1] + dy)
SRC_sdo_RSPmap_quadx = ColumnDataSource(
    {'left': xLFE.ravel(), 'right': xRTE.ravel(), 'bottom': np.repeat(yBTE[0], ndx),
     'top': np.repeat(yBTE[-1] + dy, ndx)})
SRC_sdo_RSPmap_quady = ColumnDataSource(
    {'left': np.repeat(xLFE[0], ndy), 'right': np.repeat(xLFE[-1] + dx, ndy), 'bottom': yBTE.ravel(),
     'top': yTPE.ravel()})
r_sdo_RSPmap_quadx = p_sdomap.quad('left', 'right', 'top', 'bottom', source=SRC_sdo_RSPmap_quadx,
                                   fill_alpha=0.0, fill_color=None,
                                   line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                   selection_fill_color=None,
                                   nonselection_fill_alpha=0.5,nonselection_fill_color='black',
                                   selection_line_alpha=0.0, selection_line_color=None,
                                   nonselection_line_alpha=0.0)
r_sdo_RSPmap_quady = p_sdomap.quad('left', 'right', 'top', 'bottom', source=SRC_sdo_RSPmap_quady,
                                   fill_alpha=0.0, fill_color=None,
                                   line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                   selection_fill_color=None,
                                   nonselection_fill_alpha=0.5,nonselection_fill_color='black',
                                   selection_line_alpha=0.0, selection_line_color=None,
                                   nonselection_line_alpha=0.0)
p_sdomap.add_tools(TapTool(renderers=[r_sdo_RSPmap_quadx, r_sdo_RSPmap_quady]))

SRC_sdo_RSPmap_quadx.on_change('selected', sdosubmp_region_select)
SRC_sdo_RSPmap_quady.on_change('selected', sdosubmp_region_select)

x0 = sdo_RSP_pfmap.smap.center.x.value
y0 = sdo_RSP_pfmap.smap.center.y.value
lengthx = config_main['plot_config']['tab_MkPlot']['aia_submap_sz_max']
lengthy = config_main['plot_config']['tab_MkPlot']['aia_submap_sz_max']
x0, x1 = x0 - lengthx / 2, x0 + lengthx / 2
y0, y1 = y0 - lengthy / 2, y0 + lengthy / 2
sdosubmp = sdomap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]), u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
xaxis_sdosubmp, yaxis_sdosubmp = sdosubmp.data.shape
xaxis_sdosubmp_canvas, yaxis_sdosubmp_canvas = sdosubmp.data.shape
Canvas_scaleX, Canvas_scaleY = float(xaxis_sdosubmp) / xaxis_sdosubmp_canvas, float(
    yaxis_sdosubmp) / yaxis_sdosubmp_canvas
# Canvas_scale = np.sqrt(Canvas_scaleX ** 2 + Canvas_scaleY ** 2)
Canvas_scale = (Canvas_scaleX + Canvas_scaleY) / 2.0
SRC_sdo_RSPmap_patch = ColumnDataSource({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]})
r_sdo_RSPmap_square_patch = p_sdomap.patch('xx', 'yy', source=SRC_sdo_RSPmap_patch,
                                           fill_color=None, fill_alpha=0.0, line_color="white",
                                           line_alpha=1, line_width=1)
SRC_sdo_RSPmap_line = ColumnDataSource({'xs': [], 'ys': []})
r_sdo_RSPmap_line = p_sdomap.multi_line('xs', 'ys', source=SRC_sdo_RSPmap_line, line_color='cyan', line_width=1,
                                        alpha=0.8)

sdosubmp_pfmap = PuffinMap(smap=sdosubmp,
                           plot_height=config_main['plot_config']['tab_MkPlot']['aia_submap_hght'],
                           plot_width=config_main['plot_config']['tab_MkPlot']['aia_submap_wdth'])
p_sdosubmp, r_sdosubmp = sdosubmp_pfmap.PlotMap(DrawLimb=False, DrawGrid=False, ignore_coord=True,
                                                tools='crosshair,pan,wheel_zoom,reset,save', imagetype=imagetype)
mapx_sdosubmp, mapy_sdosubmp = sdosubmp_pfmap.meshgridpix(rescale=1.0)
mapx_sdosubmp, mapy_sdosubmp = mapx_sdosubmp.value, mapy_sdosubmp.value
ndy, ndx = mapx_sdosubmp.shape
tapPointDF_sdosubmp_canvas = pd.DataFrame({'xx': [], 'yy': []})  ## the selected point to fit in aia submap
dx = np.mean(np.diff(mapx_sdosubmp[0, :]))
dy = np.mean(np.diff(mapy_sdosubmp[:, 0]))
xLFE = mapx_sdosubmp[0, :] - dx / 2.0
xRTE = np.append(xLFE[1:], xLFE[-1] + dx)
yBTE = mapy_sdosubmp[:, 0] - dy / 2.0
yTPE = np.append(yBTE[1:], yBTE[-1] + dy)
SRC_sdosubmp_canvas_quadx = ColumnDataSource(
    {'left': xLFE.ravel(), 'right': xRTE.ravel(), 'bottom': np.repeat(yBTE[0], ndx),
     'top': np.repeat(yBTE[-1] + dy, ndx)})
SRC_sdosubmp_canvas_quady = ColumnDataSource(
    {'left': np.repeat(xLFE[0], ndy), 'right': np.repeat(xLFE[-1] + dx, ndy), 'bottom': yBTE.ravel(),
     'top': yTPE.ravel()})
r_sdosubmp_quadx = p_sdosubmp.quad('left', 'right', 'top', 'bottom', source=SRC_sdosubmp_canvas_quadx,
                                   fill_alpha=0.0, fill_color=None,
                                   line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                   selection_fill_color=None,
                                   nonselection_fill_alpha=0.0,
                                   selection_line_alpha=0.0, selection_line_color=None,
                                   nonselection_line_alpha=0.0)
r_sdosubmp_quady = p_sdosubmp.quad('left', 'right', 'top', 'bottom', source=SRC_sdosubmp_canvas_quady,
                                   fill_alpha=0.0, fill_color=None,
                                   line_color=None, line_alpha=0.0, selection_fill_alpha=0.0,
                                   selection_fill_color=None,
                                   nonselection_fill_alpha=0.0,
                                   selection_line_alpha=0.0, selection_line_color=None,
                                   nonselection_line_alpha=0.0)
p_sdosubmp.add_tools(TapTool(renderers=[r_sdosubmp_quadx, r_sdosubmp_quady]))
SRC_sdosubmp_canvas_quadx.on_change('selected', sdosubmp_quad_selection_change)
SRC_sdosubmp_canvas_quady.on_change('selected', sdosubmp_quad_selection_change)

SRC_sdosubmp_line = ColumnDataSource({'xx': [], 'yy': []})
SRC_sdosubmp_line0 = ColumnDataSource({'xx': [], 'yy': []})
SRC_sdosubmp_line1 = ColumnDataSource({'xx': [], 'yy': []})
r_sdosubmp_line = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line, line_color='white', line_width=2,
                                  alpha=0.7)
r_sdosubmp_line0 = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line0, line_color='white', line_width=1.5,
                                   alpha=0.7, line_dash='dashed')
r_sdosubmp_line1 = p_sdosubmp.line('xx', 'yy', source=SRC_sdosubmp_line1, line_color='white', line_width=1.5,
                                   alpha=0.7, line_dash='dashed')
r_sdosubmp_cross = p_sdosubmp.cross(x='xx', y='yy', size=15, line_width=2, alpha=0.7, line_color='firebrick',
                                    source=ColumnDataSource(pd.DataFrame({'xx': [], 'yy': []})))

TOOLS = "save"
p_lighcurve = figure(tools=TOOLS,
                     plot_width=config_main['plot_config']['tab_MkPlot']['light_curve_wdth'],
                     plot_height=config_main['plot_config']['tab_MkPlot']['light_curve_hght'],
                     toolbar_location="above")
p_lighcurve.axis.visible = True
p_lighcurve.title.text = "light Curve"
p_lighcurve.xaxis.axis_label = 'distance [arcsec]'
p_lighcurve.yaxis.axis_label = 'Count rate [DN/s]'
p_lighcurve.border_fill_alpha = 0.4
p_lighcurve.axis.major_tick_out = 0
p_lighcurve.axis.major_tick_in = 5
p_lighcurve.axis.minor_tick_out = 0
p_lighcurve.axis.minor_tick_in = 3
SRC_lightcurve = ColumnDataSource({'x': [], 'y': []})
r_lighcurve = p_lighcurve.line(x='x', y='y', alpha=0.8, line_width=1, line_color='black', source=SRC_lightcurve)

BUT_UndoDraw = Button(label='Undo', width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'],
                      button_type='warning')
BUT_UndoDraw.on_click(UndoDraw)
BUT_ClearDraw = Button(label='ClearDraw', width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'],
                       button_type='danger')
BUT_ClearDraw.on_click(ClearDraw)
BUT_loadchunk = Button(label='LoadChunk', width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'],
                       button_type='warning')
BUT_loadchunk.on_click(LoadChunk_handler)
BUT_Stackplt = Button(label='StackPlt', width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'],
                      button_type='success')
BUT_Stackplt.on_click(ViewStackplt)

ButtonsPlayCTRL = DButil.ButtonsPlayCTRL(plot_width=config_main['plot_config']['tab_MkPlot']['button_wdth_short'])
ButtonsPlayCTRL.buttons[0].on_click(ButtonFirst_handler)
ButtonsPlayCTRL.buttons[1].on_click(ButtonPrev_handler)
ButtonsPlayCTRL.buttons[2].on_click(ButtonPlay_handler)
ButtonsPlayCTRL.buttons[3].on_click(ButtonNext_handler)
ButtonsPlayCTRL.buttons[4].on_click(ButtonLast_handler)

Text_Cutwdth = TextInput(value='5.0', title="Cut Width (pix):")
Text_CutAng = TextInput(value='0.0', title="Cut Angle (deg):")
Text_SlitLgth = TextInput(value='{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4),
                          title="Slit Length (pix):")

Hideitemlist = ["Hide Slit&Points", "Hide Points only"]
Hideitem_checkbox = CheckboxGroup(labels=Hideitemlist, active=[])

ImgOption_checkbox = CheckboxButtonGroup(labels=['Track', 'Radial Filter'], active=[0])
ImgOption_checkbox.on_click(ImgOption_checkbox_handler)

fitmethoddict = {'0': "Polyfit", '1': "Spline", '2': "Param_Spline"}
fitMeth_radiogroup = RadioGroup(labels=["Polyfit", "Spline", "Parametric Spline"], active=0)
fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
ascending = True
clearpoint = False

Slider_smoothing_factor = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title='smoothing factor',
                                 width=config_main['plot_config']['tab_MkPlot']['slider_wdth'])
Slider_smoothing_factor.on_change('value', slider_smoothing_factor_update)

Hideitem_checkbox.on_click(HideItemUpdate)
fitMeth_radiogroup.on_click(fitMethUpdate)

Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b>
            first before closing the tab</p></b>""", width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
BUT_default_fitparam = Button(label='Default', width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'],
                              button_type='primary')
BUT_default_fitparam.on_click(default_fitparam)
BUT_exit = Button(label='Exit', width=config_main['plot_config']['tab_MkPlot']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit_update)

menu_slit = [("Open", "Open"), ("Save As", "Save As"), None, ("Load", "Load"), ("Save", "Save")]
DropDn_slit = Dropdown(label="Slit File", menu=menu_slit,
                       width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'])
DropDn_slit.on_change('value', DropDn_slit_handler)

trangesec = (0.0, float(int(np.diff(trange.jd)[0] * 24 * 3600)))
Slider_trange = RangeSlider(start=trangesec[0], end=trangesec[1], range=trangesec, step=AIAcadence,
                            title='Time range [seconds since {}]'.format(trange.iso[0]),
                            width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
Slider_trange.on_change('range', trange_change_handler)
sdofileinbound = [0, nsdofile - 1]

Slider_sdoidx = Slider(start=1, end=nsdofile, value=1, step=sdofileidxstep, title='frame',
                       width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
Slider_sdoidx.on_change('value', Slider_sdoidx_update)

Slider_sdoidxstep = Slider(start=1, end=nsdofile - 1, value=1, step=1, title='step',
                           width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
Slider_sdoidxstep.on_change('value', Slider_sdoidxstep_update)

aia_wv_list = ["1700", "1600", "304", "171", "193", "211", "335", "94", "131"]
Select_aia_wave = Select(title="AIA Wavelenght:", value=MkPlot_args_dict['wavelength'], options=aia_wv_list,
                         width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
Select_aia_wave.on_change('value', aiamap_wavelength_selection)

DiffImglabelsdict = config_main['plot_config']['tab_MkPlot']['imagetype']
DiffImglabels = []
DiffImglabelsdict_r = {}
for k, v in DiffImglabelsdict.items():
    DiffImglabels.append(v)
    DiffImglabelsdict_r[v] = k

Select_DiffImg = Select(title="Difference images:", value=DiffImglabelsdict[imagetype], options=DiffImglabels,
                        width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
Select_DiffImg.on_change('value', Select_DiffImg_update)

Slider_datadt = Slider(start=AIAcadence, end=int(nsdofile - 1) / 2 * AIAcadence, value=AIAcadence, step=AIAcadence,
                       title='dt [second]', width=config_main['plot_config']['tab_MkPlot']['slider_wdth'])
Slider_datadt.on_change('value', Slider_datadt_update)
try:
    LoadSlit(slitfile)
except:
    pass

SPCR_LFT_lbutt_play = Spacer(width=config_main['plot_config']['tab_MkPlot']['space_wdth80'])
lbutt_play = row(ButtonsPlayCTRL.buttons[0], ButtonsPlayCTRL.buttons[1], ButtonsPlayCTRL.buttons[2],
                 ButtonsPlayCTRL.buttons[3],
                 ButtonsPlayCTRL.buttons[4])
widgt0 = widgetbox(Slider_trange, Slider_sdoidx, Slider_sdoidxstep,
                   width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
loutc1 = column(p_sdomap, p_lighcurve, widgt0, lbutt_play)
loutc2 = p_sdosubmp
widgt1p1 = widgetbox(Select_aia_wave, Text_SlitLgth, Text_Cutwdth, Text_CutAng,
                     Slider_smoothing_factor,
                     Hideitem_checkbox, fitMeth_radiogroup, ImgOption_checkbox,
                     width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
tab1 = Panel(child=widgt1p1, title="cutslit")
widgt1p2 = widgetbox(Select_DiffImg, Slider_datadt,
                     width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
tab2 = Panel(child=widgt1p2, title="Diff")
tabs = Tabs(tabs=[tab1, tab2])
widgt2 = widgetbox(BUT_UndoDraw, DropDn_slit, BUT_loadchunk,
                   width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'])
widgt3 = widgetbox(BUT_ClearDraw, BUT_default_fitparam, BUT_Stackplt,
                   width=config_main['plot_config']['tab_MkPlot']['button_wdth_small'])
widgt4 = widgetbox(BUT_exit, Div_info, width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
loutc3 = column(tabs, row(widgt2, widgt3), widgt4)

lout = row(loutc1, loutc2, loutc3)

curdoc().add_root(lout)
curdoc().title = "MkPlot"
