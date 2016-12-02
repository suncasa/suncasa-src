import json
import numpy as np
import glob
import os
from  casac import *
import pickle

database_dir = "${SUNCASADB}"
database_dir = os.path.expandvars(database_dir) + '/'
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
    print 'Script for clean --- {} in {}'.format(structure_id, event_id)
    print ''
    if not ('timerange' in locals()):
        timeran = [qa.time(qa.quantity(ll, 's'), prec=9)[0] for ll in mstimran['time']]
    else:
        timeran = timerange
    if not 'ncpu' in locals():
        ncpu = 10
    (tstart, tend) = timeran.split('~')
    bt_s = qa.convert(qa.quantity(tstart, 's'), 's')['value']
    et_s = qa.convert(qa.quantity(tend, 's'), 's')['value']
    btidx = np.argmin(np.abs(timeInfo - bt_s))
    etidx = np.argmin(np.abs(timeInfo - et_s))
    dt = float('{:.3f}'.format(np.median(np.diff(timeInfo))))
    if not 'twidth' in locals():
        twidth = 1
    # chunksize = ncpu
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
    with open(imgdir + 'CASA_CLN_out', 'w') as fp:
        pickle.dump(out, fp)
    imgdir = imgdir + 'local/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    print 'imgdir: {}'.format(imgdir)


    if not doreg:
        import suncasa.vla.vla_prep as vla_prep
        ms.open(vis)
        ms.selectinit()
        timfreq = ms.getdata(['time', 'axis_info'], ifraxis=True)
        tim = timfreq['time']
        dt = np.median(np.diff(tim))  # need to change to median of all time intervals
        freq = timfreq['axis_info']['freq_axis']['chan_freq'].flatten()
        ms.close()

        if twidth < 1 or twidth > len(tim):
            casalog.post('twidth not between 1 and # of time pixels in the dataset. Change to 1')
            twidth = 1

        # find out the start and end time index according to the parameter timerange
        # if not defined (empty string), use start and end from the entire time of the ms
        if not timerange:
            btidx = 1
            etidx = len(tim)-1
        else:
            try:
                (tstart, tend) = timerange.split('~')
                bt_s = qa.convert(qa.quantity(tstart, 's'), 's')['value']
                et_s = qa.convert(qa.quantity(tend, 's'), 's')['value']
                # only time is given but not date, add the date (at 0 UT) from the first record
                if bt_s < 86400. or et_s < 86400.:
                    bt_s += np.fix(qa.convert(qa.quantity(tim[0], 's'), 'd')['value']) * 86400.
                    et_s += np.fix(qa.convert(qa.quantity(tim[0], 's'), 'd')['value']) * 86400.
                btidx = np.argmin(np.abs(tim - bt_s))
                etidx = np.argmin(np.abs(tim - et_s))
                # make the indice back to those bracket by the timerange
                if tim[btidx] < bt_s:
                    btidx += 1
                if tim[etidx] > et_s:
                    etidx -= 1
                if etidx <= btidx:
                    print "ending time must be greater than starting time"
                    print "reinitiating to the entire time range"
                    btidx = 0
                    etidx = len(tim)
            except ValueError:
                print "keyword 'timerange' in wrong format"

        btstr = qa.time(qa.quantity(tim[btidx], 's'), prec=9)[0]
        etstr = qa.time(qa.quantity(tim[etidx], 's'), prec=9)[0]

        iterable = range(btidx, etidx + 1, twidth)
        print 'First time pixel: ' + btstr
        print 'Last time pixel: ' + etstr
        print str(len(iterable)) + ' images to clean...'
        timeranges = [
            qa.time(qa.quantity(tim[ll]-dt/2, 's'), prec=9)[0] + '~' + qa.time(qa.quantity(tim[ll + twidth]-dt/2, 's'), prec=9)[0] for
            ll in iterable]
        filestr = [qa.time(qa.quantity(tim[ll]-dt/2, 's'),form='fits', prec=9)[0].replace(':','').replace('-','') for ll in iterable]


        if not os.path.exists(database_dir + event_id + struct_id + 'Synthesis_Image'):
            os.system('mkdir {}'.format(database_dir + event_id + struct_id + 'Synthesis_Image/'))
            os.system('mkdir {}'.format(database_dir + event_id + struct_id + 'Synthesis_Image/local/'))
        # check if ephemfile and msinfofile exist
        if not ephemfile:
            print("ephemeris info does not exist!")
            raise ValueError
        for ll, timeran in enumerate(timeranges):
            ephem = vla_prep.read_horizons(ephemfile=ephemfile)
            reftime = [timeran]
            helio = vla_prep.ephem_to_helio(msinfo=msinfofile, ephem=ephem, reftime=reftime)
            imname = imgprefix+filestr[ll]
            imagefile = [imname + '.image']
            fitsfile = [imname + '.fits']
            try:
                vla_prep.imreg(imagefile=imagefile, fitsfile=fitsfile, helio=helio, toTb=False, scl100=True)
                # os.system('cp {} {}'.format(fitsfile[0], database_dir + event_id +'/'+ struct_id + '/Synthesis_Image/local/'))
                # print 'cp {} {}'.format(fitsfile[0], database_dir + event_id +'/'+ struct_id + '/Synthesis_Image/local/')
            except:
                '{} not found!'.format(imagefile)
    # else:
    fitsfile = glob.glob('{}*.fits'.format(imageprefix))
    for fits in fitsfile:
        idxmms = fits.index('fits')
        mms = fits[idxmms - 4:idxmms - 1]
        fits1 = fits[0:idxmms - 4] + '{:03d}'.format(int(mms) + int(dt/2*1000)) + fits[idxmms - 1:]
        fits1 = fits1.split('/')[-1]
        # print imgdir+fits1
        os.system('mv {} {}'.format(fits, imgdir + fits1))

else:
    print 'CASA arguments config file not found!!'
