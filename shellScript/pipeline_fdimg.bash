#! /bin/bash
/bin/tcsh /common/python/suncasa/shellScript/pipeline.csh -c 1 > /tmp/pipeline.log 2>&1
/bin/bash /common/python/suncasa/shellScript/pipeline_plt.bash > /tmp/pipeline_plt.log 2>&1
/bin/bash /common/python/suncasa/shellScript/pipeline_compress.bash -n 1 -O 1 > /tmp/pipeline_compress.log 2>&1

