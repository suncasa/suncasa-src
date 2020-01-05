#! /bin/bash -f

source /home/user/.bashrc
source /home/user/.setenv_py3x

if($#argv >= 3) then
    /common/anaconda2/envs/py3x/bin/python /common/python/eovsa/eovsa_pltQlookImage.py $argv
else
    /common/anaconda2/envs/py3x/bin/python /common/python/eovsa/eovsa_pltQlookImage.py
endif