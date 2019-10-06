#!/common/casa/casa-release-5.4.1-31.el6/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep

# Set to run 5 days earlier than the current date
mjdnow = ep.Time.now().mjd
#t = ep.Time(mjdnow - 6, format='mjd')
# Uncomment below and set date to run for a given date
t = ep.Time('2019-09-30 20:00')
print(t.iso)
date = t.iso[:10]
print(date)
vis_corrected = ep.calib_pipeline(date, overwrite=True, doimport=True)
