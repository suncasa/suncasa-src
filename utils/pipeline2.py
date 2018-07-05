#!/common/casa/casa-release-5.0.0-218.el6/lib/casa/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep
from eovsapy.html_movie import html_movie

# Set to run 5 days earlier than current date
mjdnow = ep.Time.now().mjd
t = ep.Time(mjdnow - 6, format='mjd')
date = t.iso[:10]
ep.qlook_image_pipeline(date, synoptic=True)
html_movie(t, imgprefix='eovsa_qlimg_')
html_movie(t, imgprefix='eovsa_qlimg_', synoptic=True)
