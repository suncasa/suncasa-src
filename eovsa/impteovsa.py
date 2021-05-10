import numpy as np
import aipy
import time
import os
import scipy.constants as constants
from taskinit import smtool, me, casalog
from astropy.time import Time


def jd2mjds(tjd=None):
    tmjds = (tjd - 2400000.5) * 24. * 3600.
    return tmjds


def bl_list2(nant=16):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = -1.
    '''
    bl2ord = np.ones((nant, nant), dtype='int') * (-1)
    k = 0
    for i in range(nant):
        for j in range(i, nant):
            bl2ord[i, j] = k
            # bl2ord[j,i] = k
            k += 1
    return bl2ord


# def get_band_edge(nband=34):
#     # Input the frequencies from UV, returen the indices frequency edges of all bands
#     idx_start_freq = [0]
#     ntmp = 0
#     for i in range(1, nband + 1):
#         ntmp += len(chan_util_bc.start_freq(i))
#         idx_start_freq.append(ntmp)
#
#     return np.asarray(idx_start_freq)


def get_band(sfreq=None, sdf=None, date=None):
    # Input the frequencies from UV
    # return a dictionary contains the band information:
    # freq : center frequency of channels
    # df :  frequency resolution
    from operator import itemgetter
    from itertools import groupby
    # nband = 34
    bandlist = []
    if date.mjd > Time('2019-02-02 12:00:00').mjd:
        import eovsapy.chan_util_52 as chan_util
    else:
        import eovsapy.chan_util_bc as chan_util

    bands = chan_util.freq2bdname(sfreq)
    cidxs = range(len(sfreq))
    spwinfo = zip(bands,sfreq,sdf,cidxs)
    for k, g in groupby(sorted(spwinfo), key=itemgetter(0)):
        itm = map(itemgetter(1,2,3), g)
        freq=[]
        df =[]
        cidx = []
        for i in itm:
            freq.append(i[0])
            df.append(i[1])
            cidx.append(i[2])
        bandlist.append({'band':k,'freq':freq,'df':np.nanmean(df),'cidx':cidx})

    return bandlist


# def uv_hex_rm(uv=None):
#     # import re
#     uvs = {}
#     for ll in uv.vartable:
#         if type(uv[ll]) == str:
#             uvs[ll] = re.sub(r'[^\x20-\x7E].*', '', uv[ll])
#     return uvs


def creatms(idbfile, outpath, timebin=None, width=None):
    uv = aipy.miriad.UV(idbfile)
    uv.rewind()
    # if idbfile.split('/')[-1][0:3] == 'UDB':
    #     uv_str = uv_hex_rm(uv)
    # else:
    #     uv_str = uv

    # uv.select('antennae', 0, 1, include=True)
    # uv.select('polarization', -5, -5, include=True)
    times = []
    uv.rewind()
    for preamble, data in uv.all():
        uvw, t, (i, j) = preamble
        times.append(t)
    times = np.unique(times)

    uv.select('clear', -1, -1, include=True)
    times = jd2mjds(np.asarray(times))
    inttime = np.median((times - np.roll(times, 1))[1:])

    start_time = 0  # The start and stop times are referenced to ref_time_jd in second
    end_time = times[-1] - times[0] + inttime

    time0 = time.time()

    if 'antlist' in uv.vartable:
        ants = uv['antlist'].replace('\x00', '')
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    good_idx = np.where(uv['sfreq'] > 0)[0]

    ref_time_jd = uv['time']
    sfreq = uv['sfreq'][good_idx]
    sdf = uv['sdf'][good_idx]
    project = uv['proj'].replace('\x00', '')
    source_id = uv['source'].replace('\x00', '')
    chan_band = get_band(sfreq=sfreq, sdf=sdf, date=Time(ref_time_jd, format='jd'))
    msname = list(idbfile.split('/')[-1])
    msname = outpath + ''.join(msname) + '_tmp.ms'

    if os.path.exists(msname):
        os.system("rm -fr %s" % msname)

    """ Creates an empty measurement set using CASA simulate (sm) tool. """
    sm = smtool()
    sm.open(msname)

    enu = np.reshape(uv['antpos'], (16, 3)) * constants.speed_of_light / 1e9
    refpos_wgs84 = me.position('wgs84',
                               '-118.286952892965deg',
                               '37.2331698901026deg',
                               '1207.1339m')
    lon, lat, rad = [me.measure(refpos_wgs84, 'itrf')[x]['value'] for x in ['m0', 'm1', 'm2']]
    # 3x3 transform matrix. Each row is a normal vector, i.e. the rows are (dE,dN,dU)
    # ----------- local xyz ------------
    xform = np.array([
        [0, -np.sin(lat), np.cos(lat)],
        [1, 0, 0],
        [0, np.cos(lat), np.sin(lat)]])
    xyz = enu.dot(xform)  # + xyz0[np.newaxis,:]

    # ----------- global xyz ------------
    # xyz0 = rad*np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    # # 3x3 transform matrix. Each row is a normal vector, i.e. the rows are (dE,dN,dU)
    # xform = np.array([
    #     [-np.sin(lon),np.cos(lon),0],
    #     [-np.cos(lon)*np.sin(lat),-np.sin(lon)*np.sin(lat),np.cos(lat)],
    #     [np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
    # ])
    # xyz = xyz0[np.newaxis,:] + enu.dot(xform)

    dishdiam = np.asarray([2.1] * uv['nants'])
    dishdiam[-3:-1] = 27
    dishdiam[-1] = 2.1
    station = uv['telescop'].replace('\x00', '')
    mount = ['ALT-AZ'] * uv['nants']
    for l in [8, 9, 10, 12, 13, 14]:
        mount[l] = 'EQUATORIAL'
    sm.setconfig(telescopename=station,
                 x=np.asarray(xyz)[:, 0],
                 y=np.asarray(xyz)[:, 1],
                 z=np.asarray(xyz)[:, 2],
                 dishdiameter=dishdiam,
                 mount=mount,
                 antname=['eo' + "{0:02d}".format(l) for l in antlist],
                 padname=station,
                 coordsystem='local', referencelocation=refpos_wgs84)

    sm.setfield(sourcename=source_id,
                sourcedirection=me.direction('J2000',
                                             '{:22.19f}'.format(uv['obsra']) + 'rad',
                                             '{:22.19f}'.format(uv['obsdec']) + 'rad'))
    sm.setfeed(mode='perfect X Y')

    ref_time = me.epoch('tai',
                        '{:20.13f}'.format(ref_time_jd - 2400000.5) + 'd')

    sm.settimes(integrationtime='{:.3f}s'.format(inttime),
                usehourangle=False,
                referencetime=ref_time)

    for l, cband in enumerate(chan_band):
        nchannels = len(cband['freq'])
        stokes = 'XX YY XY YX'
        sm.setspwindow(spwname='band{:02d}'.format(cband['band']),
                       freq='{:22.19f}'.format(cband['freq'][0] - cband['df'] / 2.0) + 'GHz',
                       deltafreq='{:22.19f}'.format(cband['df']) + 'GHz',
                       freqresolution='{:22.19f}'.format(cband['df']) + 'GHz',
                       nchannels=nchannels,
                       stokes=stokes)

    for l, cband in enumerate(chan_band):
        print('sm-band{}'.format(cband['band']))
        sm.observe(source_id, 'band{:02d}'.format(cband['band']),
                   starttime=start_time, stoptime=end_time,
                   project=project,
                   state_obs_mode='')

    if sm.done():
        casalog.post('Empty MS {0} created in --- {1:10.2f} seconds ---'.format(msname, (time.time() - time0)))
    else:
        raise RuntimeError('Failed to create MS. Look at the log file. '
                           'Double check you settings.')

    sm.close()
    modelms = msname + '.MSmodel'
    os.system('mv {} {}'.format(msname, modelms))

    # if timebin != '0s' or width != 1:
    #     modelms = msname + '.MSmodel'
    #     split(vis=msname, outputvis=modelms, datacolumn='data', timebin=timebin, width=width)
    #     os.system('rm -rf {}'.format(msname))

    return modelms
