import os
import numpy as np
import glob

mspath = '/srg/sjyu/20141101/'
os.chdir(mspath)
msfile = 'SUN01_20141101T163940-164700.50ms.cal.ms'
ephemfile = 'horizons_sun_20141101.radecp'
msinfofile = 'SUN01_20141101.T163940-164700.50ms.cal.msinfo.npz'

##################  ----------- U04burst-prep -----------------#########
slfcalms_timeran = '16:46:12.375~16:46:32.475'
slfcalms_s = ['0', '1', '2', '3']
slfcalms_chan = ['0~63', '0~63', '0~63', '0~63']
structure_id = 'Uburst04-prep'
pol = 'LL'
refantenna = 'ea04'
antennas = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26'
spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))
##  ----------- spw 0 -----------------
struct_timerange = '2014/11/01/16:46:17.100~2014/11/01/16:46:17.600'
struct_freqrange = '1.080~1.116 GHz'
CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
phasecenter_local = 'J2000 14h26m59.300 -14d36m00.200'
phasecenter_global = 'J2000 14h26m22.7351 -14d29m29.801'
##  ----------- spw 1 -----------------
# struct_timerange = '2014/11/01/16:46:17.600~2014/11/01/16:46:18.400'
# struct_freqrange = '1.146~1.182 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# phasecenter_local = 'J2000 14h26m59.500 -14d36m01.100'
##  ----------- spw 2 -----------------
# struct_timerange = '2014/11/01/16:46:17.650~2014/11/01/16:46:18.250'
# struct_freqrange = '1.286~1.318 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# phasecenter_local = 'J2000 14h27m00.100 -14d35m24.200'
##  ----------- spw 3 -----------------
# struct_timerange = '2014/11/01/16:46:16.150~2014/11/01/16:46:19.650'
# struct_freqrange = '1.422~1.476 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# phasecenter_local = 'J2000 14h26m56.300 -14d35m11.000'


##################  ----------- typeIII01-prep -----------------#########
slfcalms_timeran = "16:42:48~16:42:54"
slfcalms_s = ['0', '1', '2', '3']
slfcalms_chan = ['0~63', '0~63', '0~63', '0~63']
structure_id = 'typeIII01-prep'
pol = 'LL'
refantenna = 'ea04'
antennas = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26'
spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))
##  ----------- spw 0 -----------------
struct_timerange = '2014/11/01/16:42:50.925~2014/11/01/16:42:53.825'
struct_freqrange = '1.016~1.102 GHz'
CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
phasecenter_local = 'J2000 14h26m59.300 -14d36m00.200'
phasecenter_global = 'J2000 14h26m22.7351 -14d29m29.801'

#
# ##################  ----------- typeIII01 -----------------#########
# slfcalms_timeran = '16:43:58.825~16:44:01.325'
# slfcalms_s = ['0', '1']
# slfcalms_chan = ['0~63', '0~63']
# structure_id = 'typeIII01-prep'
# pol = 'RR'
# spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))
# phasecenter_local = 'ljdfsldfj'
# # spwchan='0,1,2'
# ##  ----------- spw 0 -----------------
# struct_timerange = '2014/11/01/16:44:00.975~2014/11/01/16:44:01.125'
# struct_freqrange = '1.040~1.072 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# ##  ----------- spw 1 -----------------
# struct_timerange = '2014/11/01/16:44:00.475~2014/11/01/16:44:00.625'
# struct_freqrange = '1.160~1.208 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# ##  ----------- spw 2 -----------------
# # struct_timerange = '2014/11/01/16:46:17.650~2014/11/01/16:46:18.250'
# # struct_freqrange = '1.286~1.318 GHz'
# # CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# ##  ----------- spw 4 -----------------
# # struct_timerange = '2014/11/01/16:46:17.275~2014/11/01/16:46:19.450'
# # struct_freqrange = '1.426~1.472 GHz'
# # CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
#
#
# ##################  ----------- ZP01 -----------------#########
# slfcalms_timeran = '16:41:02~16:41:26'
# slfcalms_s = ['0', '1']
# slfcalms_chan = ['0~63', '0~63']
# structure_id = 'ZP01-prep'
# pol = 'RRLL'
# spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))
# phasecenter_local = 'J2000 14h26m43.082 -14d37m43.817'
# ##  ----------- spw 0 -----------------
# struct_timerange = '2014/11/01/16:41:11.125~2014/11/01/16:41:12.325'
# struct_freqrange = '1.044~1.066 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# ##  ----------- spw 1 -----------------
# # struct_timerange = '2014/11/01/16:41:12.575~2014/11/01/16:41:14.175'
# # struct_freqrange = '1.166~1.184 GHz'
# # CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)


# ##################  ----------- U02 -----------------#########
# slfcalms_timeran = '16:45:39~16:45:48'
# slfcalms_s = ['0', '1', '2']
# slfcalms_chan = ['0~63', '0~63', '0~63']
# structure_id = 'U02-prep'
# pol = 'LL'
# spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))
# phasecenter_local = 'J2000 14h26m58.891 -14d35m24.673'
# ##  ----------- spw 0 -----------------
# struct_timerange = '2014/11/01/16:45:41.925~2014/11/01/16:45:42.225'
# struct_freqrange = '1.042~1.060 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)
# ##  ----------- spw 1 -----------------
# struct_timerange = '2014/11/01/16:45:40.875~2014/11/01/16:45:41.175'
# struct_freqrange =  '1.168~1.196 GHz'
# CLNmask = 'slfcal/{}/region_spw0.rgn'.format(structure_id)

''' ----- step 1 ----- '''

dir_ms_split = './ms_split/'
if not os.path.exists(dir_ms_split):
    os.mkdir(dir_ms_split)

strtmp = [t.replace(':', '') for t in slfcalms_timeran.split('~')]
timestr = 'T' + strtmp[0] + '-' + strtmp[1]
slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.' + structure_id + '.50ms.slfcal.ms'
slfcal_dbkg_ms = dir_ms_split + 'SUN01_20141101T' + timestr + '.' + structure_id + '.50ms.slfcal.dbkg.ms'

debkg = False
prep = False
tofits = False

print 'Script for calibrating 2014 Nov 1 data --- ' + structure_id
print ''

if debkg:
    slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.' + structure_id + '.50ms.slfcal.dbkg.ms'
else:
    slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.' + structure_id + '.50ms.slfcal.ms'

if prep:
    if debkg:
        print '--de-background--'
        default('subvs')
        vis = msfile
        outputvis = slfcalms
        timerange = slfcalms_timeran
        spw = spwchan
        mode = 'linear'
        subtime1 = '16:46:11.525~16:46:11.825'
        subtime2 = '16:46:28.525~16:46:28.825'
        overwrite = True
        if os.path.exists(slfcalms):
            os.system('rm -fr %s*' % slfcalms)
        subvs()
    else:
        # Split the calibrated data
        print '--split--'
        default('split')
        vis = msfile
        datacolumn = 'data'
        field = 'SUN01'
        timerange = slfcalms_timeran
        spw = spwchan
        outputvis = slfcalms
        if os.path.exists(slfcalms):
            os.system('rm -fr %s*' % slfcalms)
        split()
        print 'split msfile to ' + slfcalms

ms.open(slfcalms)

axisInfo = ms.getdata(["axis_info"], ifraxis=True)
spwInfo = ms.getspectralwindowinfo()
freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
freqInfo_ravel = freqInfo.ravel()
timeInfo = axisInfo["axis_info"]["time_axis"]['MJDseconds']
timran = ms.range(["time"])
ms.close()

t_int = 0.05

tt = struct_timerange.split('~')
t0str = tt[0].split('/')[-1]
t1str = tt[1].split('/')[-1]
# timestr=qa.time(qa.quantity(tmid,'s'),prec=9)[0].translate(None, ':')
timestr = t0str.translate(None, ':')
timeran = t0str + '~' + t1str
f_int = 2.0  # MHz
# if 'freqrange' in locals():
freq0, freq1 = struct_freqrange.split(' ')[0].split('~')
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
    ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1] - 1)] + ['0~{}'.format(sz_freqInfo[1] - 1) for ll in
                                                                      xrange(freqIdx0[0] + 1, freqIdx1[0])]
    ms_chan.append('0~{}'.format(freqIdx1[1][0]))
spw = ','.join('%s:%s' % t for t in zip(ms_spw, ms_chan))
spwstr = spw
# spwstr = 'spw{}-{}'.format(slfcalms_s[0], slfcalms_s[-1])

calprefix = 'caltables/'
if not os.path.exists(calprefix):
    os.makedirs(calprefix)
imgprefix = 'slfcal/' + structure_id + '/'
if not os.path.exists(imgprefix):
    os.makedirs(imgprefix)

calprefix = calprefix + 'cal_SUN_' + timestr + 'T.spw' + ms_spw[0]
imgprefix = imgprefix + structure_id + '_'

slfcal_list = range(4)
# initial image and first slfcal table
slfcal_img_global = [imgprefix + timestr + 'T.' + spwstr.replace(':', '-') + '.global.slfcal{:d}'.format(ll) for ll in
                     slfcal_list]
slfcal_img_local = [imgprefix + timestr + 'T.' + spwstr.replace(':', '-') + '.local.slfcal{:d}'.format(ll) for ll in
                    slfcal_list]

slfcal_img_list = [imgprefix + timestr + 'T.' + spwstr.replace(':', '-') + '.slfcal{:d}'.format(ll) for ll in
                   slfcal_list]
slfcal_table_list = [calprefix + '.slfcal.G{:d}'.format(ll) for ll in slfcal_list]
uvrange_list = ['<5klambda', '<5klambda', '', '']
slfcalms_list = [slfcalms[:-3] + '{}.ms'.format(ll) for ll in slfcal_list]
# uvrange_list = ['<1klambda', '<3klambda', '<6klambda', '<9klambda', '','']




slfcal_iter = 0
clearcal(vis=slfcalms, spw=ms_spw[0])
tb.open(slfcalms, nomodify=False)
colnames = tb.colnames()
cols2rm = ["MODEL_DATA", 'CORRECTED_DATA']
for l in range(len(cols2rm)):
    if cols2rm[l] in colnames:
        tb.removecols(cols2rm[l])
tb.close()

for ll in slfcal_img_list:
    os.system('rm -rf {}*'.format(ll))
for ll in slfcal_table_list:
    os.system('rm -rf {}*'.format(ll))

## first round of self calibration
default('clean')
vis = slfcalms
imagename = slfcal_img_list[slfcal_iter]
spw = spwstr
timerange = timeran
mode = 'mfs'
imagermode = 'csclean'
weighting = 'briggs'
robust = 0.5
niter = 1000
antenna = antennas
uvrange = uvrange_list[slfcal_iter]
npercycle = 50
imsize = [512, 512]
cell = ['5arcsec', '5arcsec']
# imsize = [256, 256]
# cell = ['3arcsec', '3arcsec']
phasecenter = phasecenter_global
stokes = pol
# uvtaper = True
# outertaper = ['30.0arcsec']
interactive = True
usescratch = True
clean()
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')

##gain solution, phase only
gaincal(vis=slfcalms, refant=refantenna, caltable=slfcal_table_list[slfcal_iter], uvrange=uvrange_list[slfcal_iter],
        spw=ms_spw[0], selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p', minsnr=2,
        minblperant=2, antenna=antennas)
plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
        subplot=111, iteration='channel')
clearcal(vis=slfcalms, spw=ms_spw[0])
applycal(vis=slfcalms, gaintable=slfcal_table_list[slfcal_iter], spw=ms_spw[0], selectdata=True, \
         antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)

os.system('rm -fr ' + slfcalms_list[slfcal_iter + 1])
split(vis=slfcalms, outputvis=slfcalms_list[slfcal_iter + 1], datacolumn='corrected')

##second round of self calibration
slfcal_iter += 1
tget('clean')
vis = slfcalms_list[slfcal_iter]
imagename = slfcal_img_list[slfcal_iter]
# uvrange = uvrange_list[slfcal_iter]
clean()
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')

gaincal(vis=slfcalms_list[slfcal_iter], refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw=spwstr, \
        selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p',
        uvrange=uvrange_list[slfcal_iter], minsnr=2, minblperant=2, antenna=antennas)

plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', subplot=111)
clearcal(vis=slfcalms, spw=ms_spw[0])
applycal(vis=slfcalms, gaintable=slfcal_table_list[0:slfcal_iter + 1], spw=ms_spw[0], selectdata=True, \
         antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)

os.system('rm -fr ' + slfcalms_list[slfcal_iter + 1])
split(vis=slfcalms, outputvis=slfcalms_list[slfcal_iter + 1], spw=ms_spw[0], datacolumn='corrected')

##third round of self calibration
slfcal_iter += 1
tget('clean')
vis = slfcalms_list[slfcal_iter]
imagename = slfcal_img_list[slfcal_iter]
# uvrange = uvrange_list[slfcal_iter]
clean()
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')

gaincal(vis=slfcalms_list[slfcal_iter], refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw=spwstr, \
        selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p',
        uvrange=uvrange_list[slfcal_iter], minsnr=2, minblperant=2, antenna=antennas)

plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
        subplot=111)
clearcal(vis=slfcalms, spw=ms_spw[0])
applycal(vis=slfcalms, gaintable=slfcal_table_list[0:slfcal_iter + 1], spw=ms_spw[0], selectdata=True, \
         antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)

os.system('rm -fr ' + slfcalms_list[slfcal_iter + 1])
split(vis=slfcalms, outputvis=slfcalms_list[slfcal_iter + 1], spw=ms_spw[0], datacolumn='corrected')

##final image
slfcal_iter += 1
tget('clean')
vis = slfcalms_list[slfcal_iter]
imagename = slfcal_img_list[slfcal_iter]
clean()
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
# os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')

# split(vis=slfcalms_list[slfcal_iter], outputvis=slfcalms[:-3]+'.spw{}.ms'.format(ms_spw[0]), spw = ms_spw[0],datacolumn='data')

# os.system('cp -r {} caltables/{}.spw{}.G'.format(slfcal_table_list[slfcal_iter-1], structure_id, ms_spw[0]))

#
# for slfcal_iter in slfcal_list[1:]:
#     # slfcal_iter+=1
#     tget('clean')
#     vis = slfcalms
#     imagename = slfcal_img_list[slfcal_iter]
#     uvrange = uvrange_list[slfcal_iter]
#     # mask='rU01.mask1'
#     niter = 500
#     npercycle = 20
#     # interactive=False
#     clean()
#     os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
#     # os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
#     os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
#     os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
#     os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')
#
#     # gaincal(vis=slfcalms, refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw=spwstr, \
#     #         selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p', \
#     #         minsnr=3, gaintable=slfcal_table_list[0:slfcal_iter])
#     # clearcal(vis=slfcalms,spw=ms_spw[0])
#     # if os.path.exists(slfcal_table_list[slfcal_iter]):
#     #     plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
#     #             subplot=111)
#     #     applycal(vis=slfcalms, gaintable=slfcal_table_list[0:(slfcal_iter + 1)], spw=ms_spw[0], selectdata=True, \
#     #              antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
#     # else:
#     #     slfcal_iter -= 1
#     #     applycal(vis=slfcalms, gaintable=slfcal_table_list[0:(slfcal_iter + 1)], spw=ms_spw[0], selectdata=True, \
#     #              antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
#     #     break
#     gaincal(vis=slfcalms, refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw=spwstr, \
#             selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p',
#             uvrange=uvrange_list[slfcal_iter], \
#             minsnr=3)
#     clearcal(vis=slfcalms, spw=ms_spw[0])
#     if os.path.exists(slfcal_table_list[slfcal_iter]):
#         plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
#                 subplot=111)
#         applycal(vis=slfcalms, gaintable=slfcal_table_list[slfcal_iter], spw=ms_spw[0], selectdata=True, \
#                  antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
#     else:
#         slfcal_iter -= 1
#         applycal(vis=slfcalms, gaintable=slfcal_table_list[slfcal_iter], spw=ms_spw[0], selectdata=True, \
#                  antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
#         break
# print 'self calibration finished within {:d} iterations'.format(slfcal_iter + 1)
# # os.system('cp -r {} caltables/{}.spw{}.G'.format(slfcal_table_list[slfcal_iter], structure_id, ms_spw[0]))
#
# for slfcal_iter in slfcal_list:
#     try:
#         clearcal(vis=slfcalms, spw=ms_spw[0])
#         applycal(vis=slfcalms, gaintable=slfcal_table_list[slfcal_iter], spw=ms_spw[0], selectdata=True, \
#                  antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
#         # final imaging
#         tget('clean')
#         vis = slfcalms
#         imagename = slfcal_img_global[slfcal_iter]
#         uvrange = ''
#         # mask='rU01.mask2'
#         mask = CLNmask
#         niter = 200
#         npercycle = 10
#         imsize = [512, 512]
#         cell = ['5arcsec', '5arcsec']
#         phasecenter = phasecenter_global
#         uvtaper = True
#         outertaper = ['30.0arcsec']
#         interactive = False
#         usescratch = False
#         # outertaper = ['50arcsec']
#         clean()
#         os.system('rm -fr ' + slfcal_img_global[slfcal_iter] + '*.flux')
#         os.system('rm -fr ' + slfcal_img_global[slfcal_iter] + '*.mask')
#         os.system('rm -fr ' + slfcal_img_global[slfcal_iter] + '*.model')
#         os.system('rm -fr ' + slfcal_img_global[slfcal_iter] + '*.psf')
#         os.system('rm -fr ' + slfcal_img_global[slfcal_iter] + '*.residual')
#
#         tget('clean')
#         vis = slfcalms
#         imagename = slfcal_img_local[slfcal_iter]
#         uvrange = ''
#         mask = CLNmask
#         niter = 200
#         npercycle = 10
#         imsize = [64, 64]
#         cell = ['5.0arcsec', '5.0arcsec']
#         phasecenter = phasecenter_local
#         uvtaper = True
#         outertaper = ['30.0arcsec']
#         interactive = False
#         usescratch = False
#         clean()
#         os.system('rm -fr ' + slfcal_img_local[slfcal_iter] + '.flux')
#         os.system('rm -fr ' + slfcal_img_local[slfcal_iter] + '.mask')
#         os.system('rm -fr ' + slfcal_img_local[slfcal_iter] + '.model')
#         os.system('rm -fr ' + slfcal_img_local[slfcal_iter] + '.psf')
#         os.system('rm -fr ' + slfcal_img_local[slfcal_iter] + '.residual')
#     except:
#         pass
#
# for ll in slfcal_table_list:
#     os.system('rm -fr {}'.format(ll))
# for ll in slfcal_img_list:
#     os.system('rm -fr {}.*'.format(ll))

# clearcal(vis=slfcalms)
# caltabl = glob.glob('caltables/{}.spw?.slfcal.G?'.format(structure_id))
# for ll in caltabl:
#     stridx1 = ll.index('spw')
#     stridx2 = ll.index('.G')
#     applycal(vis=slfcalms, gaintable=ll, spw=ll[stridx1 + 3:stridx2], selectdata=True, \
#              antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)

clearcal(vis=slfcalms)
for ll in xrange(4):
    caltabl = glob.glob('caltables/cal_SUN*.spw{}.slfcal.G?'.format(ll))
    if len(caltabl) > 0:
        applycal(vis=slfcalms, gaintable=caltabl, selectdata=True, \
                 antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
    else:
        print 'caltable for spw {} not found'.format(ll)
