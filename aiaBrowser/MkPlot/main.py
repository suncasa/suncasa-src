import glob
import json
import os
import astropy.units as u
import bokeh.palettes as bp
import pickle
import numpy as np
import pandas as pd
import sunpy.map
from astropy.time import Time
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (ColumnDataSource, Slider, Button, TextInput, CheckboxGroup, RadioGroup,
                          BoxSelectTool, TapTool, Div, Spacer, Range1d)
from bokeh.plotting import figure, curdoc
from sunpy.time import TimeRange
from suncasa.utils.puffin import PuffinMap
from suncasa.utils import DButil

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
        sdosubmp = sdosubmpdict['map'][fileidx]
    r_sdosubmp.data_source.data['image'] = [sdosubmp.data]
    p_sdosubmp.title.text = sdosubmp.name
    try:
        getimprofile(sdosubmp, cutslitplt)
    except:
        ValueError('make a cut slit first!')


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


def sdosubmp_region_select(attrname, old, new):
    global ImgDF0_sdo_RSPmap, sdosubmp, x0, x1, y0, y1
    global xaxis_sdosubmp, yaxis_sdosubmp
    global Canvas_scale, Canvas_scaleX, Canvas_scaleY
    r_sdo_RSPmap_square_selected = SRC_sdo_RSPmap_square.selected['1d']['indices']
    if r_sdo_RSPmap_square_selected:
        ImgDF_sdo_RSPmap = ImgDF0_sdo_RSPmap.iloc[r_sdo_RSPmap_square_selected, :]
        x0, x1 = ImgDF_sdo_RSPmap['xx'].min(), ImgDF_sdo_RSPmap['xx'].max()
        y0, y1 = ImgDF_sdo_RSPmap['yy'].min(), ImgDF_sdo_RSPmap['yy'].max()
        if x1 > x0 + mapx_sdo_RSPmap[0, 1] - mapx_sdo_RSPmap[0, 0] and y1 > y0 + mapy_sdo_RSPmap[0, 1] - \
                mapy_sdo_RSPmap[0, 0]:
            ## select at least 4 points
            patch_size = min(x1 - x0, y1 - y0, config_main['plot_config']['tab_MkPlot']['aia_submap_sz_max'])
            x1, y0 = x0 + patch_size, y1 - patch_size
            update_sdosubmp_region(x0, x1, y0, y1)
            Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)
            ClearDraw()


def ProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '=' * (length - filledLength)
    # return '%s |%s| %s%% %s' % (prefix, bar, percent, suffix)
    return '{} |{}| {}% {}'.format(prefix, bar, percent, suffix)


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
            ProgressBar(sidx + 1, nsdofile, suffix='Load', decimals=0, length=16, fill='#'))
    sdompdict = {'map': sdosubmplist, 'time': np.array(timestamps)}
    # Div_info.text = """<p><b>SDO submap chunk loaded.</b></p>"""
    return sdompdict


def LoadSubChunk(sdompdict):
    sdosubmplist = []
    for sidx, smap in enumerate(sdompdict['map']):
        sdosubmptmp = smap.submap(u.Quantity([x0 * u.arcsec, x1 * u.arcsec]),
                                  u.Quantity([y0 * u.arcsec, y1 * u.arcsec]))
        sdosubmplist.append(sdosubmptmp)
        Div_info.text = """<p>{}</p>""".format(
            ProgressBar(sidx + 1, nsdofile, suffix='Update', decimals=0, length=14, fill='#'))
    return {'map': sdosubmplist, 'time': sdompdict['time']}


def LoadChunk_handler():
    global sdosubmpdict
    if sdosubmpdict:
        if x0 >= sdosubmpdict['map'][0].xrange[0].value and x1 <= sdosubmpdict['map'][0].xrange[1].value and y0 >= \
                sdosubmpdict['map'][0].yrange[0].value and y1 <= sdosubmpdict['map'][0].yrange[1].value:
            sdosubmpdict = LoadSubChunk(sdosubmpdict)
        else:
            sdosubmpdict = LoadChunk()
    else:
        sdosubmpdict = LoadChunk()


def ButtonNext_handler():
    global sdofileidx
    sdofileidx += 1
    sdofileidx = sdofileidx % nsdofile
    update_sdosubmp_image(sdofileidx)


def ButtonPrev_handler():
    global sdofileidx
    sdofileidx -= 1
    sdofileidx = sdofileidx % nsdofile
    update_sdosubmp_image(sdofileidx)


def ButtonLast_handler():
    global sdofileidx
    sdofileidx = - 1
    sdofileidx = sdofileidx % nsdofile
    update_sdosubmp_image(sdofileidx)


def ButtonFirst_handler():
    global sdofileidx
    sdofileidx = 0
    sdofileidx = sdofileidx % nsdofile
    update_sdosubmp_image(sdofileidx)


def ButtonPlay_handler():
    global sdofileidx
    if ButtonsPlayCTRL.buttons[2].label == '>':
        ButtonsPlayCTRL.buttons[2].label = '||'
        curdoc().add_periodic_callback(sdosubmp_animate_update, 100)
    else:
        ButtonsPlayCTRL.buttons[2].label = '>'
        curdoc().remove_periodic_callback(sdosubmp_animate_update)


def sdosubmp_animate_update():
    if sdosubmpdict['map']:
        ButtonNext_handler()


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
        getimprofile(sdosubmp, cutslitplt)
        pushslit2plt(cutslitplt, clearpoint=clearpoint)


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
    Text_CutAng.value = '5.0'
    Text_SlitLgth.value = '{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4)


def LoadSlit():
    global cutslitplt, tapPointDF_sdosubmp_canvas, ascending
    global fitmethod, smoth_factor1, smoth_factor2
    global x0, x1, y0, y1
    fin = database_dir + 'cutslit-' + PlotID
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
        Div_info.text = """<p><b>cutslit settings</b> in {} loaded.</p>""".format(fin)
    else:
        Div_info.text = """<p>No cutslit file found!!</p>"""


def SaveSlit():
    if cutslitplt:
        fout = database_dir + 'cutslit-' + PlotID
        cutslitdict = {'tapPointDF_sdosubmp_canvas': tapPointDF_sdosubmp_canvas,
                       'ascending': ascending, 'fitmethod': fitmethod, 'smoth_factor1': smoth_factor1,
                       'smoth_factor2': smoth_factor1, 'smoth_factor': Slider_smoothing_factor.value,
                       'Cutwdth': Text_Cutwdth.value, 'CutAng': Text_CutAng.value,
                       'SlitLgth': Text_SlitLgth.value, 'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1}
        with open(fout, 'wb') as fp:
            pickle.dump(cutslitdict, fp)
        Div_info.text = """<p><b>cutslit settings</b> saved to {}.</p>""".format(fout)
    else:
        Div_info.text = """<p>make a slit first!!!</p>"""


def MkStackplt():
    global stackpltdict
    if sdosubmpdict:
        stackplt = []
        for sidx, smap in enumerate(sdosubmpdict['map']):
            intens = getimprofile(smap, cutslitplt, plot=False)
            stackplt.append(intens['y'])
        stackplt = np.vstack(stackplt)
        stackplt = stackplt.transpose()
        stackpltdict = {'zz': stackplt, 'x': sdosubmpdict['time'], 'y': intens['x'],
                        'wavelength': '{:.0f}'.format(smap.wavelength.value), 'observatory': smap.observatory,
                        'instrument': smap.instrument[0:3]}
        Div_info.text = """<p>click <b>Stackplt</b> to view the stack plot</p>"""
        return True
    else:
        Div_info.text = """<p>load the chunk first!!!</p>"""
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


def exit_update():
    Div_info.text = """<p><b>You may close the tab anytime you like.</b></p>"""
    raise SystemExit


ports = []
'''load config file'''
suncasa_dir = os.path.expandvars("${SUNCASA}") + '/'
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
with open(infile, 'rb') as fp:
    MkPlot_args_dict = json.load(fp)

trange = Time([MkPlot_args_dict['tst'], MkPlot_args_dict['ted']], format='iso', scale='utc')
PlotID = MkPlot_args_dict['PlotID']
sdofile = DButil.readsdofile(datadir=SDOdir, wavelength='171', jdtime=trange.jd)
nsdofile = len(sdofile)
global sdofileidx
if nsdofile == 0:
    raise SystemExit('No SDO file found at the select timestamp. Download the data with EvtBrowser first.')
sdofileidx = 0
cutslitplt = {}
sdosubmpdict = {}
sdomap = sunpy.map.Map(sdofile[sdofileidx])
sdomap = DButil.normalize_aiamap(sdomap)
MapRES = config_main['plot_config']['tab_MkPlot']['aia_RSPmap_RES']
dimensions = u.Quantity([MapRES, MapRES], u.pixel)
sdo_RSPmap = sdomap.resample(dimensions)
sdo_RSP_pfmap = PuffinMap(smap=sdo_RSPmap,
                          plot_height=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_hght'],
                          plot_width=config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'])
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
                                      size=max(config_main['plot_config']['tab_MkPlot']['aia_RSPmap_hght'] /
                                               config_main['plot_config']['tab_MkPlot']['aia_RSP_squareny'],
                                               config_main['plot_config']['tab_MkPlot']['aia_RSPmap_wdth'] /
                                               config_main['plot_config']['tab_MkPlot']['aia_RSP_squarenx']))

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
SRC_sdo_RSPmap_patch = ColumnDataSource(
    pd.DataFrame({'xx': [x0, x1, x1, x0], 'yy': [y0, y0, y1, y1]}))
r_sdo_RSPmap_square_patch = p_sdomap.patch('xx', 'yy', source=SRC_sdo_RSPmap_patch,
                                           fill_color=None, fill_alpha=0.5, line_color="white",
                                           line_alpha=1, line_width=1)
p_sdomap.add_tools(BoxSelectTool(renderers=[r_sdo_RSPmap_square]))

SRC_sdo_RSPmap_square.on_change('selected', sdosubmp_region_select)

sdosubmp_pfmap = PuffinMap(smap=sdosubmp,
                           plot_height=config_main['plot_config']['tab_MkPlot']['aia_submap_hght'],
                           plot_width=config_main['plot_config']['tab_MkPlot']['aia_submap_wdth'])
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
                                      size=max(config_main['plot_config']['tab_MkPlot']['aia_submap_hght'] /
                                               config_main['plot_config']['tab_MkPlot']['aia_RSP_squareny'],
                                               config_main['plot_config']['tab_MkPlot']['aia_submap_wdth'] /
                                               config_main['plot_config']['tab_MkPlot']['aia_RSP_squarenx']))
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
                     plot_width=config_main['plot_config']['tab_MkPlot']['time_plot_wdth'],
                     plot_height=config_main['plot_config']['tab_MkPlot']['time_plot_hght'],
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

BUT_UndoDraw = Button(label='Undo', width=config_main['plot_config']['tab_MkPlot']['button_wdth'],
                      button_type='warning')
BUT_UndoDraw.on_click(UndoDraw)
BUT_ClearDraw = Button(label='ClearDraw', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                       button_type='danger')
BUT_ClearDraw.on_click(ClearDraw)
BUT_loadchunk = Button(label='LoadChunk', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                       button_type='primary')
BUT_loadchunk.on_click(LoadChunk_handler)
BUT_Stackplt = Button(label='StackPlt', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                      button_type='primary')
BUT_Stackplt.on_click(ViewStackplt)

ButtonsPlayCTRL = DButil.ButtonsPlayCTRL(plot_width=config_main['plot_config']['tab_MkPlot']['button_wdth_short'])
ButtonsPlayCTRL.buttons[0].on_click(ButtonFirst_handler)
ButtonsPlayCTRL.buttons[1].on_click(ButtonPrev_handler)
ButtonsPlayCTRL.buttons[2].on_click(ButtonPlay_handler)
ButtonsPlayCTRL.buttons[3].on_click(ButtonNext_handler)
ButtonsPlayCTRL.buttons[4].on_click(ButtonLast_handler)

Text_Cutwdth = TextInput(value='5.0', title="Cut Width (pix):")
Text_CutAng = TextInput(value='5.0', title="Cut Angle (deg):")
Text_SlitLgth = TextInput(value='{:.0f}'.format(np.sqrt(xaxis_sdosubmp ** 2 + yaxis_sdosubmp ** 2) / 4),
                          title="Slit Length (pix):")

Hideitemlist = ["Hide Slit&Points", "Hide Points only"]
Hideitem_checkbox = CheckboxGroup(labels=Hideitemlist, active=[])
fitmethoddict = {'0': "Polyfit", '1': "Spline", '2': "Param_Spline"}
fitMeth_radiogroup = RadioGroup(labels=["Polyfit", "Spline", "Parametric Spline"], active=0)
fitmethod = fitmethoddict['{}'.format(fitMeth_radiogroup.active)]
ascending = True
clearpoint = False

SRC_sdosubmp_canvas_square.on_change('selected', sdosubmp_square_selection_change)
Slider_smoothing_factor = Slider(start=0.0, end=1.0, value=1.0, step=0.05, title='smoothing factor')
Slider_smoothing_factor.on_change('value', slider_smoothing_factor_update)

Hideitem_checkbox.on_click(HideItemUpdate)
fitMeth_radiogroup.on_click(fitMethUpdate)

Div_info = Div(text="""<p><b>Warning</b>: Click <b>Exit</b>
            first before closing the tab</p></b>""", width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
BUT_default_fitparam = Button(label='Default', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                              button_type='primary')
BUT_default_fitparam.on_click(default_fitparam)
BUT_loadslit = Button(label='LoadSlit', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                      button_type='success')
BUT_loadslit.on_click(LoadSlit)
BUT_saveslit = Button(label='SaveSlit', width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'],
                      button_type='success')
BUT_saveslit.on_click(SaveSlit)
BUT_exit = Button(label='Exit', width=config_main['plot_config']['tab_MkPlot']['button_wdth'], button_type='danger')
BUT_exit.on_click(exit_update)

SPCR_LFT_lbutt_play = Spacer(width=config_main['plot_config']['tab_MkPlot']['space_wdth80'])
# lbutt_play = row(BUT_first, BUT_prev, BUT_play, BUT_next, BUT_last)
# lbutt_play = buttons_play
lbutt_play = row(ButtonsPlayCTRL.buttons[0], ButtonsPlayCTRL.buttons[1], ButtonsPlayCTRL.buttons[2],
                 ButtonsPlayCTRL.buttons[3],
                 ButtonsPlayCTRL.buttons[4])
loutc1 = column(p_sdomap, p_lighcurve, lbutt_play)
loutc2 = p_sdosubmp
widgt1 = widgetbox(Text_SlitLgth, Text_Cutwdth, Text_CutAng, Slider_smoothing_factor, Hideitem_checkbox,
                   fitMeth_radiogroup, BUT_UndoDraw, width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
widgt2 = widgetbox(BUT_ClearDraw, BUT_loadslit, BUT_loadchunk,
                   width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'])
widgt3 = widgetbox(BUT_default_fitparam, BUT_saveslit, BUT_Stackplt,
                   width=config_main['plot_config']['tab_MkPlot']['button_wdth_half'])
widgt4 = widgetbox(BUT_exit, Div_info, width=config_main['plot_config']['tab_MkPlot']['button_wdth'])
loutc3 = column(widgt1, row(widgt2, widgt3), widgt4)

lout = row(loutc1, loutc2, loutc3)

# def timeout_callback():
#     print 'timeout'
#     raise SystemExit


curdoc().add_root(lout)
curdoc().title = "MkPlot"
