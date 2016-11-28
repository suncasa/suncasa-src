#!/bin/tcsh -f
#
set port = `(netstat  -atn | awk '{printf "%s\n%s\n", $4, $4}' | grep -oE '[0-9]*$'; seq 32768 61000) | sort -n | uniq -u | head -n 1`

echo Starting Bokeh server on port $port
bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser --port $port