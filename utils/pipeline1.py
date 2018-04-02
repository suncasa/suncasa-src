#!/common/casa/casa-release-5.0.0-218.el6/lib/casa/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep

# Set to run 5 days earlier than the current date
mjdnow = ep.Time.now().mjd
t = ep.Time(mjdnow - 6, format='mjd')
print t.iso
date = t.iso[:10]
print date
vis_corrected = ep.calib_pipeline(date, doimport=True, synoptic=True)
