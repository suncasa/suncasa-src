import os


d = {'a':'1', 'b':'2'}
for key,val in d.items():
    exec(key + '= {}'.format(val))


mspath = '/srg/sjyu/20141101/'
os.chdir(mspath)
msfile = 'sun_20141101_t191020-191040.50ms.cal.ms/'
ephemfile = 'horizons_sun_20141101.radecp'
msinfofile = 'sun_20141101_t191020-191040.50ms.cal.msinfo.npz'
# slfcalms_timeran='16:46:00~16:46:35'
slfcalms_timeran = '19:10:36.875~19:10:36.925'
slfcalms_s = ['0', '1']#, '2', '3']
slfcalms_chan = ['0~63', '0~63']#, '0~63', '0~63']
structure_id = 'typeIII-prep'
spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s, slfcalms_chan))


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
        field = 'SUN02'
        timerange = slfcalms_timeran
        spw = spwchan
        outputvis = slfcalms
        if os.path.exists(slfcalms):
            os.system('rm -fr %s*' % slfcalms)
        split()
        print 'split msfile to ' + slfcalms
