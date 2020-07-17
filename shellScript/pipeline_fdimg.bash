#! /bin/bash
/bin/tcsh /common/python/suncasa/shellScript/pipeline1.csh --clearcache > /tmp/pipeline.log 2>&1
/bin/bash /common/python/suncasa/shellScript/pipeline_plt.bash > /tmp/pipeline_plt.log 2>&1
/bin/bash /common/python/suncasa/shellScript/pipeline_compress.bash > /tmp/pipeline_compress.log 2>&1

