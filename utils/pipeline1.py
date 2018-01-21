#!/common/casa/casa-release-5.0.0-218.el6/lib/casa/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep
t = ep.Time('2017-07-10')
print t.iso
date = t.iso[:10]
print date
vis_corrected = ep.calib_pipeline(date, doimport=True)
