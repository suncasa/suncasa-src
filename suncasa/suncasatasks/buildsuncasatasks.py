## build suncasa tasks on inti
## export PATH="$PATH:/inti/software/casa-6.2.0-124/bin/"

import os
from glob import glob

xmlfiles = glob('*.xml')

for xmlfile in xmlfiles:
    os.system("buildmytasks --upgrade {}".format(xmlfile))
    os.system("buildmytasks --module suncasatasks {}".format(xmlfile))
