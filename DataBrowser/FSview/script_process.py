import json
import numpy as np
import glob
import os
from  casac import *

database_dir = "${SUNCASADB}"
database_dir = os.path.expandvars(database_dir)+'/'
if os.path.exists('CASA_CLN_args.json'):
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))

    # if 'mspath' in locals():
    #     os.chdir(mspath)
    if 'vis' in locals():
        ms.open(vis)
        axisInfo = ms.getdata(["axis_info"], ifraxis=True)
        spwInfo = ms.getspectralwindowinfo()
        freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
        freqInfo_ravel = freqInfo.ravel()
        timeInfo = axisInfo["axis_info"]["time_axis"]['MJDseconds']
        mstimran = ms.range(["time"])
    else:
        raise ValueError('define a vis file!!!')

    if 'struct_id' in locals():
        structure_id = struct_id
    else:
        raise ValueError('define a struct_id!!!')
    print 'Script for calibrating --- {} in {}'.format(structure_id, event_id)
    print ''
    if not ('timerange' in locals()):
        timeran = [qa.time(qa.quantity(ll, 's'), prec=9)[0] for ll in mstimran['time']]
    else:
        timeran = timerange
    if not 'ncpu' in locals():
        ncpu = 10
    chunksize = ncpu
    (tstart, tend) = timeran.split('~')
    bt_s = qa.convert(qa.quantity(tstart, 's'), 's')['value']
    et_s = qa.convert(qa.quantity(tend, 's'), 's')['value']
    btidx = np.argmin(np.abs(timeInfo - bt_s))
    etidx = np.argmin(np.abs(timeInfo - et_s))
    dt = float('{:.3f}'.format(np.median(np.diff(timeInfo))))
    if not 'twidth' in locals():
        twidth = 1
    # timerans = []
    # if etidx <= btidx + twidth * chunksize:
    #     btstr = qa.time(qa.quantity(timeInfo[btidx] - dt / 2, 's'), prec=9, form='fits')[0]
    #     etstr = qa.time(qa.quantity(timeInfo[etidx] + dt / 2, 's'), prec=9, form='fits')[0]
    #     timerans.append('{}~{}'.format(btstr, etstr))
    # else:
    #     for iter in xrange(btidx, etidx + 1, twidth * chunksize):
    #         btstr = qa.time(qa.quantity(timeInfo[iter] - dt / 2, 's'), prec=9, form='fits')[0]
    #         if iter <= etidx - twidth * chunksize:
    #             etstr = qa.time(qa.quantity(timeInfo[iter + twidth * chunksize] - dt / 2, 's'), prec=9, form='fits')[0]
    #         else:
    #             etstr = qa.time(qa.quantity(et_s + dt / 2, 's'), prec=9, form='fits')[0]
    #         timerans.append('{}~{}'.format(btstr, etstr))

    imgprefix = 'slfcal/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)
    imgprefix = 'slfcal/' + structure_id + '/'
    if not os.path.exists(imgprefix):
        os.mkdir(imgprefix)
    '''set a loop to limit the number of timestamps in the time range'''
    # for TRang in timerans:
    # TRang=timerans[0]
    os.system('rm -rf {}'.format('cgrid_ft.im'))
    default('ptclean')
    with open('CASA_CLN_args.json', 'r') as fp:
        CASA_CLN_args = json.load(fp)
    for key, val in CASA_CLN_args.items():
        exec (key + '= {}'.format(val))
    timerange = timeran
    # width = 32
    if 'freqrange' in locals() and spw == '':
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
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1] - 1)] \
                      + ['0~{}'.format(sz_freqInfo[1] - 1) for ll in xrange(freqIdx0[0] + 1, freqIdx1[0])]
            ms_chan.append('0~{}'.format(freqIdx1[1][0]))
        spw = ','.join('{}:{}'.format(t[0], t[1]) for t in zip(ms_spw, ms_chan))
    # inp(ptclean)
    out = ptclean()

    imgdir = database_dir + event_id + '/' + struct_id + '/Synthesis_Image/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    imgdir = imgdir + 'local/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    print 'imgdir: {}'.format(imgdir)
    fitsfile = glob.glob('{}*.fits'.format(imageprefix))
    for fits in fitsfile:
        idxmms = fits.index('fits')
        mms = fits[idxmms - 4:idxmms - 1]
        fits1 = fits[0:idxmms - 4] + '{:03d}'.format(int(mms) + 25) + fits[idxmms - 1:]
        fits1 = fits1.split('/')[-1]
        # print imgdir+fits1
        os.system('mv {} {}'.format(fits, imgdir + fits1))

else:
    print 'CASA arguments config file not found!!'
