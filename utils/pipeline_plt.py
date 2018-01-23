#!/common/casa/casa-release-5.0.0-218.el6/lib/casa/bin/casa
from suncasa.eovsa import eovsa_pipeline as ep
from eovsapy.html_movie import html_movie
t = ep.Time('2017-07-10')
date = t.iso[:10]
qlookfitsdir='/data1/eovsa/qlookfits/'
ep.imres=mk_qlook_image(date,imagedir=qlookfitsdir)
