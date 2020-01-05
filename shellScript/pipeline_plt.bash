#! /bin/bash -f

#set up python path variable, and path to Anaconda Python
export PYTHONPATH="/common/python/"
export PATH="/common/anaconda2/bin:$PATH"
source activate py3x

if ["$#" -ge 3]; then
    /common/anaconda2/envs/py3x/bin/python /common/python/suncasa/eovsa/eovsa_pltQlookImage.py "$@"
else
    /common/anaconda2/envs/py3x/bin/python /common/python/suncasa/eovsa/eovsa_pltQlookImage.py
fi