import os

mspath = '/srg/sjyu/20141101/'
os.chdir(mspath)
msfile = 'SUN01_20141101.T163940-164700.50ms.cal.ms'
ephemfile = 'horizons_sun_20141101.radecp'
msinfofile = 'SUN01_20141101.T163940-164700.50ms.cal.msinfo.npz'
# slfcalms_timeran='16:46:00~16:46:35'
slfcalms_timeran = '16:46:13.125~16:46:29.525'
slfcalms_s = ['0', '1', '2', '3']
slfcalms_chan = ['0~63', '0~63', '0~63', '0~63']
structure_id = 'U04-prep'
spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))

# specfile='/srg/sjyu/20141101/SUN01_20141101TT164613.125-164629.525.spw'+slfcalms_s[0]+'.50ms.slfcal.ms.spec.npz'
# specdata=np.load(specfile)
# spec=specdata['spec']
# npol=specdata['npol']
# nbl=specdata['nbl']
# ntim=specdata['ntim']
# nfreq=specdata['nfreq']
# tim=specdata['tim']
# freq=specdata['freq']
# bl=specdata['bl'].item()
# b=0
# pol='LL'
# sz_spec=spec.shape
# spec_plt=np.zeros((4,sz_spec[2],sz_spec[3]))
# spec_plt[0,:,:]=spec[1,b,:,:]
# # spec_plt[1,:,:]=spec[1,b,:,:]
# # spec_plt[2,:,:]=spec[1,b,:,:]
# spec_plt[3,:,:]=spec_med = medfilt2d(spec[1,b,:,:],kernel_size=[3,3])

# # threshold=140
# threshold=300
# idx_thresh1 = spec_plt[3,:,:] > threshold
# # idx_thresh2 = spec_plt[0,:,:] <= threshold
# spec_plt[1,:,:][idx_thresh1] = 400
# spec_plt[2,:,:][idx_thresh1] = spec_plt[0,:,:][idx_thresh1]
# idx_selec = np.where(spec_plt[3,:,:] > threshold)
# # idx_selec = np.where(spec_plt[0,:,:] > 0.8*np.max(spec_plt[0,:,:]))
# idx_tim=np.unique(idx_selec[1])
# idx_selec=[30,110]

''' ----- step 1 ----- '''

dir_ms_split = './ms_split/'
if not os.path.exists(dir_ms_split):
    os.mkdir(dir_ms_split)

strtmp = [t.replace(':', '') for t in slfcalms_timeran.split('~')]
timestr = 'T' + strtmp[0] + '-' + strtmp[1]
slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.50ms.slfcal.ms'
slfcal_dbkg_ms = dir_ms_split + 'SUN01_20141101T' + timestr + '.50ms.slfcal.dbkg.ms'

debkg = False
prep = False
tofits = True

print 'Script for calibrating 2014 Nov 1 data --- ' + structure_id
print ''

if debkg:
    slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.50ms.slfcal.dbkg.ms'
else:
    slfcalms = dir_ms_split + 'SUN01_20141101T' + timestr + '.50ms.slfcal.ms'

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

# tools_20141101_rU.slfcal(slfcalms=slfcalms, ephemfile=ephemfile, msinfofile=msinfofile,
#             timeran=timeran, s=s, 
#             chan=chan, structure_id=structure_id,
#             pol=pol,refantenna=refantenna,antennas=antennas,
#             prep=prep,debkg=debkg,tofits=tofits)
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------



pol = 'LL'
refantenna = 'ea04'
# antennas='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26'
antennas = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26'

# t_int=tim[idx_tim[-1]]-tim[idx_tim[len(idx_selec[0])/2]]
# ii=900#len(idx_selec[1])/2

t_int = 0.05
# t0=tim[idx_selec[1]]
# t1=t0+t_int
# tmid=(t0+t1)/2.0
# t0str=qa.time(qa.quantity(t0,'s'),prec=9)[0]
# t1str=qa.time(qa.quantity(t1,'s'),prec=9)[0]
t0str = '16:46:17.475'
t1str = '16:46:17.725'
# timestr=qa.time(qa.quantity(tmid,'s'),prec=9)[0].translate(None, ':')
timestr = t0str.translate(None, ':')
timeran = t0str + '~' + t1str
f_int = 2.0  # MHz
# f0=freq[idx_selec[0]]/1.e6
# f1=freq[idx_selec[0]]/1.e6+f_int
# fmid=(f0+f1)/2.0
# chn='{:d}~{:d}'.format(idx_selec[0],idx_selec[0]+1)
# freqrange='{0:5.1f}~{1:5.1f}MHz'.format(f0,f1)
# fmidstr='{:5.1f}MHz'.format(fmid)
spwstr = 'spw{}-{}'.format(slfcalms_s[0], slfcalms_s[-1])

calprefix = 'caltables/'
if not os.path.exists(calprefix):
    os.mkdir(calprefix)
imgprefix = 'slfcal/'
if not os.path.exists(imgprefix):
    os.mkdir(imgprefix)
imgprefix = 'slfcal/' + structure_id + '/'
if not os.path.exists(imgprefix):
    os.mkdir(imgprefix)

calprefix = calprefix + 'cal_SUN_' + timestr + 'T.' + spwstr
imgprefix = imgprefix + structure_id + '_'

# initial image and first slfcal table
slfcal_img_global = imgprefix + timestr + 'T.' + spwstr + '.global.slfcal'
slfcal_img_local = imgprefix + timestr + 'T.' + spwstr + '.local.slfcal'

slfcal_img_list = [imgprefix + timestr + 'T.' + spwstr + '.slfcal{:d}'.format(ll) for ll in range(5)]
slfcal_table_list = [calprefix + '.slfcal.G{:d}'.format(ll) for ll in range(5)]
uvrange_list = ['<3klambda', '<6klambda', '<9klambda', '', '']

slfcal_iter = 0
clearcal(vis=slfcalms)
tb.open(slfcalms, nomodify=False)
colnames = tb.colnames()
cols2rm = ["MODEL_DATA", 'CORRECTED_DATA']
for l in range(len(cols2rm)):
    if cols2rm[l] in colnames:
        tb.removecols(cols2rm[l])
tb.close()

# first round of self calibration
default('clean')
vis = slfcalms
imagename = slfcal_img_list[slfcal_iter]
spw = '0:30~40'
timerange = timeran
# mode='channel'
mode = 'mfs'
imagermode = 'csclean'
weighting = 'natural'
gain = 0.05
niter = 500
uvrange = uvrange_list[slfcal_iter]
npercycle = 50
# mask = ['circle [ [ 147pix, 181pix ], 28pix]']
# mask = 'box [ [ 119pix , 153pix] , [175pix, 211pix ] ]'
# mask='rU01.mask0'
imsize = [512, 512]
cell = ['5arcsec', '5arcsec']
phasecenter = 'J2000 14h26m22.7351 -14d29m29.801'
# imsize=[128,128]
# cell = ['3.0arcsec', '3.0arcsec']
# phasecenter = 'J2000 14h26m59.250 -14d35m44.681'
stokes = pol
uvtaper = False
# outertaper=['30.0arcsec']
interactive = True
# interactive=False
usescratch = True
clean()
# os.system('rm -fr '+slfcal_img_list[slfcal_iter]+'.flux')
# os.system('rm -fr '+slfcal_img_list[slfcal_iter]+'.mask')
# os.system('rm -fr '+slfcal_img_list[slfcal_iter]+'.model')
# os.system('rm -fr '+slfcal_img_list[slfcal_iter]+'.psf')
# os.system('rm -fr '+slfcal_img_list[slfcal_iter]+'.residual')




clearcal(vis=slfcalms)
# gain solution, phase only
gaincal(vis=slfcalms, refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw='1:32', \
        selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p', minsnr=3)
plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
        subplot=111, iteration='channel')
applycal(vis=slfcalms, gaintable=slfcal_table_list[0:(slfcal_iter + 1)], spw='0', selectdata=True, \
         antenna=antennas, interp='linear', flagbackup=False, applymode='calonly')

for slfcal_iter in range(1, 2):
    # slfcal_iter=1
    tget('clean')
    vis = slfcalms
    imagename = slfcal_img_list[slfcal_iter]
    uvrange = uvrange_list[slfcal_iter]
    # mask='rU01.mask1'
    niter = 500
    npercycle = 50
    # interactive=False
    clean()
    os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.flux')
    os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.mask')
    os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.model')
    os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.psf')
    os.system('rm -fr ' + slfcal_img_list[slfcal_iter] + '.residual')

    clearcal(vis=slfcalms)
    gaincal(vis=slfcalms, refant=refantenna, caltable=slfcal_table_list[slfcal_iter], spw='1:32', \
            selectdata=True, timerange=timeran, solint='inf', gaintype='G', calmode='p', \
            minsnr=3, gaintable=slfcal_table_list[0:slfcal_iter])
    if os.path.exists(slfcal_table_list[slfcal_iter]):
        plotcal(caltable=slfcal_table_list[slfcal_iter], antenna=antennas, xaxis='antenna', yaxis='phase', \
                subplot=111)
        applycal(vis=slfcalms, gaintable=slfcal_table_list[0:(slfcal_iter + 1)], spw='0', selectdata=True, \
                 antenna=antennas, interp='linear', flagbackup=False, applymode='calonly')
    else:
        slfcal_iter -= 1
        applycal(vis=slfcalms, gaintable=slfcal_table_list[0:(slfcal_iter + 1)], spw='0', selectdata=True, \
                 antenna=antennas, interp='linear', flagbackup=False, applymode='calonly')
        break
print 'self calibration finished within {:d} iterations'.format(slfcal_iter + 1)

# final imaging
tget('clean')
vis = slfcalms
imagename = slfcal_img_global
uvrange = ''
# mask='rU01.mask2'
niter = 2000
npercycle = 100
imsize = [512, 512]
cell = ['5arcsec', '5arcsec']
phasecenter = 'J2000 14h26m22.7351 -14d29m29.801'
outertaper = ['50.0arcsec']
interactive = False
usescratch = False
outertaper = ['50arcsec']
clean()
os.system('rm -fr ' + slfcal_img_global + '.flux')
os.system('rm -fr ' + slfcal_img_global + '.mask')
os.system('rm -fr ' + slfcal_img_global + '.model')
os.system('rm -fr ' + slfcal_img_global + '.psf')
os.system('rm -fr ' + slfcal_img_global + '.residual')

tget('clean')
vis = slfcalms
imagename = slfcal_img_local
uvrange = ''
# mask='rU01.mask2'
niter = 2000
npercycle = 100
imsize = [128, 128]
cell = ['3.0arcsec', '3.0arcsec']
phasecenter = 'J2000 14h26m59.250 -14d35m44.681'
outertaper = ['50arcsec']
interactive = False
usescratch = False
clean()
os.system('rm -fr ' + slfcal_img_local + '.flux')
os.system('rm -fr ' + slfcal_img_local + '.mask')
os.system('rm -fr ' + slfcal_img_local + '.model')
os.system('rm -fr ' + slfcal_img_local + '.psf')
os.system('rm -fr ' + slfcal_img_local + '.residual')


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
