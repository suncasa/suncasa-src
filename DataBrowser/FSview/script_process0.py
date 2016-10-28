import os
import json
import numpy as np
import jdutil
from datetime import datetime

if os.path.exists('CASA_CLN_args.json'):
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))
    # mspath = '/srg/sjyu/20141101/'
    os.chdir(mspath)
    # msfile = 'sun_20141101_t191020-191040.50ms.cal.ms/'
    # ephemfile = 'horizons_sun_20141101.radecp'
    # msinfofile = 'sun_20141101_t191020-191040.50ms.cal.msinfo.npz'
    ms.open(vis)
    axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    freqInfo_ravel = freqInfo.ravel()
    timeInfo = axisInfo["axis_info"]["time_axis"]['MJDseconds']
    timran = ms.range(["time"])

    if not ('timerange' in locals()):
        timerange = [qa.time(qa.quantity(ll, 's'), prec=9)[0] for ll in timran['time']]
    ms_timeran = timerange
    ms_timeran_obj = [datetime.strptime(ll + '000', '%Y/%m/%d/%H:%M:%S.%f') for ll in ms_timeran.split('~')]
    ms_timeran_MJDseconds = [(jdutil.datetime_to_jd(ll) - 2400000.5) * 24. * 3600. for ll in ms_timeran_obj]

    for ll in ms_timeran_MJDseconds:
        if ll < timeInfo[0] or ll > timeInfo[-1]:
            raise ValueError('Selected time out of range!!!')
    timIdx0 = np.argmin(abs(timeInfo - ms_timeran_MJDseconds[0]))
    timIdx1 = np.argmin(abs(timeInfo - ms_timeran_MJDseconds[1]))

    if 'struct_id' in locals():
        structure_id = struct_id
    else:
        raise ValueError('define a struct_id!!!')

    if 'freqrange' in locals() and not ('spw' in locals()):
        freq0, freq1 = freqrange.split(' ')[0].split('~')
        freq0, freq1 = float(freq0), float(freq1)
        for ll in [freq0, freq1]:
            if not freqInfo_ravel[0] <= ll <= freqInfo_ravel[-1]:
                raise ValueError('Selected frequency out of range!!!')
        freqIdx0 = np.where(freqInfo == freq0)
        freqIdx1 = np.where(freqInfo == freq1)
        sz_freqInfo = freqInfo.shape
        ms_spw = ['{}'.format(ll) for ll in xrange(freqIdx0[0], freqIdx1[0] + 1)]
        if len(ms_spw) == 1:
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], freqIdx1[1][0])]
        else:
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1]-1)] + ['0~{}'.format(sz_freqInfo[1]-1) for ll in
                                                                          xrange(freqIdx0[0] + 1, freqIdx1[0])]
            ms_chan.append('0~{}'.format(freqIdx1[1][0]))
        spwchan = ','.join('%s:%s' % t for t in zip(ms_spw, ms_chan))
    elif 'spw' in locals():
        spwchan = spw
    else:
        spwchan = ''

    # tim0_char = jdutil.jd_to_datetime((timran['time'][0] / 3600. / 24. + 2400000.5) * 86400. / 3600. / 24.)
    # strdate = tim0_char.strftime('%Y-%m-%d')
    # strtmp = [
    #     jdutil.jd_to_datetime((ll / 3600. / 24. + 2400000.5) * 86400. / 3600. / 24.).strftime('%H%M%S') + '.{}'.format(
    #         round(tim0_char.microsecond / 1e3) * 1e3)[0:4] for ll in timran['time']]
    # timestr = 'T' + strtmp[0] + '-' + strtmp[1]
    # tofits = True

    print 'Script for calibrating 2014 Nov 1 data --- ' + structure_id
    print ''

    tb.open(vis, nomodify=False)
    colnames = tb.colnames()
    cols2rm = ["MODEL_DATA", 'CORRECTED_DATA']
    for l in range(len(cols2rm)):
        if cols2rm[l] in colnames:
            tb.removecols(cols2rm[l])
    tb.close()

    imgprefix = 'slfcal/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)
    imgprefix = 'slfcal/' + structure_id + '/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)


    t_int = int((timeInfo[1]-timeInfo[0])*1000)/1000.0
    # for tidx in xrange(timIdx0, timIdx1 + 1):
    #     t0 = timeInfo[tidx]
    #     t1 = timeInfo[tidx] + t_int
    #     t0str = qa.time(qa.quantity(t0, 's'), prec=9)[0]
    #     t1str = qa.time(qa.quantity(t1, 's'), prec=9)[0]
    #     timestr = t0str.translate(None, ':')
    #     timeran = t0str + '~' + t1str
    #     spwstr = 'spw{}'.format(ms_spw[spwcurr])

    '''first round of self calibration'''
    default('ptclean')
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))
    imagedir=imgprefix
    timerange=timeran
    spw = spwchan
    mask = ['slfcal/typeIII-prep/typeIII-prep_spw1.local.slfcal.mask']
    '''spw 0'''
    imsize=[128,128]
    # phasecenter = 'J2000 14h26m51.180 -14d29m30.983'
    # cell = ['5.0arcsec', '5.0arcsec']
    # '''spw 1'''
    # phasecenter = 'J2000 14h26m58.180 -14d30m35.983'
    # cell = ['2.0arcsec', '2.0arcsec']
    ptclean()
    # os.system('rm -fr ' + slfcal_img_local + '.flux')
    # os.system('rm -fr ' + slfcal_img_local + '.mask')
    # os.system('rm -fr ' + slfcal_img_local + '.model')
    # os.system('rm -fr ' + slfcal_img_local + '.psf')
    # os.system('rm -fr ' + slfcal_img_local + '.residual')
else:
    print 'CASA arguments config file not found!!'



# ia.open('slfcal/rU01-test2/rU01-test2_164619.150T.spw1.global.slfcal.image')
# ROI='centerbox[[14h26m59.250,-14d35m44.681], [384.0arcsec, 384.0arcsec]]'
# ia.fromimage(outfile='test', infile='slfcal/rU01-test2/rU01-test2_164619.150T.spw1.global.slfcal.image', region=ROI, overwrite=true)  
# ia.close()  


# phasecenter = me.direction('j2000','14h26m59.250','-14d35m44.681')
# ia.open('slfcal/rU01-test2/rU01-test2_164619.150T.spw1.global.slfcal.image')
# csys = ia.coordsys()
# cout=csys.convert(coordin=[qa.getvalue(phasecenter['m0'])[0],qa.getvalue(phasecenter['m1'])[0],0,0],absin=[T,T,T,T],  
#                   unitsin=["rad","rad","pix","pix"],  
#                   absout=[T,T,T,T],  
#                   unitsout=["pix","pix","","GHz"]) 
# ROI== rg.box(blc=[round(cout[0])-64, round(cout[1])-64, 0, 0], trc=[round(cout[0])+64, round(cout[1])+64, 0, 0])
# ia.fromimage(outfile='test', infile='slfcal/rU01-test2/rU01-test2_164619.150T.spw1.global.slfcal.image', region=ROI, overwrite=true) 


# ''' ----- step 4 ----- '''
# if tofits:
#     ephem=vla_prep.read_horizons(ephemfile=ephemfile)
#     reftime=[timeran]*2
#     helio=vla_prep.ephem_to_helio(msinfo=msinfofile,ephem=ephem,reftime=reftime)
#     imagenames=[slfcal_img1,slfcal_img2]
#     imagefile=[img+'.image' for img in imagenames]
#     vlafits=['./slfcal/fits/'+img.split('/')[-1]+'.image.fits' for img in imagenames]
#     vla_prep.imreg(imagefile=imagefile,fitsfile=vlafits,helio=helio,toTb=F,scl100=True)
