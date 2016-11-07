import os
import json
import numpy as np
import jdutil
from datetime import datetime

if os.path.exists('CASA_CLN_args.json'):
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))

    if 'mspath' in locals():
        os.chdir(mspath)
    if 'vis' in locals():
        ms.open(vis)
    else:
        raise ValueError('define a vis file!!!')

    axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    freqInfo_ravel = freqInfo.ravel()
    timeInfo = axisInfo["axis_info"]["time_axis"]['MJDseconds']
    timran = ms.range(["time"])
    if 'struct_id' in locals():
        structure_id = struct_id
    else:
        raise ValueError('define a struct_id!!!')
    print 'Script for calibrating 2014 Nov 1 data --- ' + structure_id
    print ''

    tb.open(vis, nomodify=False)
    colnames = tb.colnames()
    cols2rm = ["MODEL_DATA", 'CORRECTED_DATA']
    for l in range(len(cols2rm)):
        if cols2rm[l] in colnames:
            tb.removecols(cols2rm[l])
    tb.close()

    imgprefix = 'slfcal/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)
    imgprefix = 'slfcal/' + structure_id + '/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)


    '''first round of self calibration'''
    default('ptclean')
    imagedir=imgprefix
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))
    if not ('timerange' in locals()):
        timerange = [qa.time(qa.quantity(ll, 's'), prec=9)[0] for ll in timran['time']]

    if 'freqrange' in locals() and not ('spw' in locals()):
        freq0, freq1 = freqrange.split(' ')[0].split('~')
        freq0, freq1 = float(freq0), float(freq1)
        for ll in [freq0, freq1]:
            if not freqInfo_ravel[0] <= ll <= freqInfo_ravel[-1]:
                raise ValueError('Selected frequency out of range!!!')
        freqIdx0 = np.where(freqInfo == freq0)
        freqIdx1 = np.where(freqInfo == freq1)
        sz_freqInfo = freqInfo.shape
        ms_spw = ['{}'.format(ll) for ll in xrange(freqIdx0[0], freqIdx1[0] + 1)]
        if len(ms_spw) == 1:
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], freqIdx1[1][0])]
        else:
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1]-1)] + ['0~{}'.format(sz_freqInfo[1]-1) for ll in
                                                                          xrange(freqIdx0[0] + 1, freqIdx1[0])]
            ms_chan.append('0~{}'.format(freqIdx1[1][0]))
        spw = ','.join('%s:%s' % t for t in zip(ms_spw, ms_chan))

    ptclean()
else:
    print 'CASA arguments config file not found!!'