## build suncasa tasks on inti
## export PATH="$PATH:/inti/software/casa-6.2.0-124/bin/"
## [ref]: https://casadocs.readthedocs.io/en/latest/api/casashell/buildmytasks.html#
## using buildmytasks requires a full installation of CASA. https://casadocs.readthedocs.io/en/stable/notebooks/usingcasa.html

import os
from glob import glob

xmlfiles = glob('ptclean6*.xml')

for xmlfile in xmlfiles:
    # os.system("buildmytasks --upgrade {}".format(xmlfile)) ## this only need to be run once.
    os.system("buildmytasks --module suncasatasks {}".format(xmlfile))
