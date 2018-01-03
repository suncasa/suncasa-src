#!/common/casa/casa-release-5.0.0-218.el6/lib/casa/bin/casa
from datetime import datetime as dt
import glob
import os

tnow = dt.now()
# tnow = dt(2017, 12, 17)
yy = tnow.strftime("%Y")
ym = tnow.strftime("%Y%m")
ymd = tnow.strftime("%Y%m%d")
inpath = '/data1/eovsa/fits/UDB/{}/'.format(yy)
idbfiles = [os.path.basename(ll) for ll in glob.glob('{}UDB{}*'.format(inpath, ymd))]

outpath = '/data1/eovsa/fits/UDBms/{}/'.format(ym)
if not os.path.exists(outpath):
    os.makedirs(outpath)
msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB{}*.ms'.format(outpath, ymd))]

files2import = [inpath + ll for ll in list(set(idbfiles) - set(msfiles))]
if files2import:
    importeovsa(idbfiles=files2import, ncpu=1, timebin="0s", width=1,
                visprefix=outpath, nocreatms=False, doconcat=False, modelms="",doscaling=False)
else:
    print 'No new UDB files found. Quit.'
    # # add to crontab file
    # # cronjob to convert UDB data to CASA Measurement Sets every 10 minutes
    # */10 * * * * touch /data1/eovsa/fits/UDBms/LOG/UDB2MS$(date +\%Y\%m\%d).log;/bin/tcsh /home/user/sjyu/udb2ms.csh >> /data1/eovsa/fits/UDBms/LOG/UDB2MS$(date +\%Y\%m\%d).log 2>&1
