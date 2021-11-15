#!/common/casa/casa-release-5.4.1-31.el6/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep
from eovsapy.html_movie import html_movie
import sys
syspath = sys.path
sys.path = [s for s in syspath if '.local' not in s]

# Set to run 5 days earlier than current date
mjdnow = ep.Time.now().mjd
t = ep.Time(mjdnow - 6, format='mjd')
# Uncomment below and set date to run for a given date
# t = ep.Time('2019-09-17 20:00')
date = t.iso[:10]
# ep.qlook_image_pipeline(date, ncpu=1, synoptic=True)
ep.qlook_image_pipeline(date, ncpu=1)
html_movie(t, imgprefix='eovsa_qlimg_')
# html_movie(t, imgprefix='eovsa_qlimg_', synoptic=True)
