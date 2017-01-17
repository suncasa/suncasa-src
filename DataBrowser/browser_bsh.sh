#!/bin/bash -f
#
port=32768; while netstat -atn | grep -q $port; do port=$(expr $port + 1); done; echo $port
echo Starting Bokeh server on port $port
bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser --port $port
