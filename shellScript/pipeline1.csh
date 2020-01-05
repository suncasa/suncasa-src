#! /bin/tcsh -f

source /home/user/.cshrc
# Do a daily update of the calibration calendar
/common/anaconda2/bin/python /common/python/current/cal_calendar.py `date +\%Y\ \%m`
# Now run the first pipeline casa task (importeovsa)
#/common/casa/casa-release-5.0.0-218.el6/bin/casa --nologfile -c /common/python/suncasa/utils/pipeline1.py
# /common/casa/casa-release-5.4.1-31.el6/bin/casa --nologfile --agg --nogui -c /common/python/suncasa/utils/pipeline1.py
if($#argv >= 3) then
    /common/casa/casa-release-5.4.1-31.el6/bin/casa --nologfile --agg --nogui -c /home/user/workdir/test_pipeline1.py $argv
else
    /common/casa/casa-release-5.4.1-31.el6/bin/casa --nologfile --agg --nogui -c /home/user/workdir/test_pipeline1.py
endif