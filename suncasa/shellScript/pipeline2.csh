#! /bin/tcsh -f

source /home/user/.cshrc
cd /home/user/workdir
/common/casa/casa-release-5.4.1-31.el6/bin/casa --agg --nologfile -c /common/python/suncasa/utils/pipeline2.py