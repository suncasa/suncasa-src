import numpy as np
import matplotlib.pyplot as plt
import os
# from config import get_and_create_download_dir
import shutil
from astropy.io import fits
import urllib2
from split_cli import split_cli as split
from ptclean_cli import ptclean_cli as ptclean
from suncasa.utils import helioimage2fits as hf
import sunpy.map as smap
from sunpy.net import vso
from astropy import units as u
from astropy.time import Time
from taskinit import ms, tb, qa, iatool
from clean_cli import clean_cli as clean
from sunpy import lightcurve
from sunpy.time import TimeRange
from matplotlib.dates import DateFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord
from sunpy import lightcurve as lc
from sunpy.time import TimeRange, parse_time
import pickle
import datetime
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patches as patches
import glob
from suncasa.utils import DButil
import pdb


def mk_qlook_image(vis, ncpu=10, twidth=12, stokes='RR,LL', antenna='', imagedir=None,
                   spws=['0', '1', '2', '3', '4', '5', '6', '7'], toTb=True, overwrite=True, doslfcal=False,
                   phasecenter=''):
    vis = [vis]
    subdir = ['/']

    for idx, f in enumerate(vis):
        if f[-1] == '/':
            vis[idx] = f[:-1]
    if not stokes:
        stokes = 'RR'
    stokes = stokes.replace(',', '')

    if not imagedir:
        imagedir = './'
    imres = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': []}
    msfile = vis[0]
    ms.open(msfile)
    # axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    # freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    # freqInfo_ravel = freqInfo.ravel()
    ms.close()
    msfilebs = os.path.basename(msfile)
    imdir = imagedir + subdir[0]
    if not os.path.exists(imdir):
        os.makedirs(imdir)
    if doslfcal:
        slfcalms = './' + msfilebs + '.rr'
        split(msfile, outputvis=slfcalms, datacolumn='corrected', correlation='RR')
    for spw in spws:
        spwran = [s.zfill(2) for s in spw.split('~')]
        freqran = [
            (int(s) * spwInfo['0']['TotalWidth'] + spwInfo['0']['RefFreq'] + spwInfo['0']['TotalWidth'] / 2.0) / 1.0e9
            for s in spw.split('~')]
        cfreq = np.mean(freqran)
        bmsz = max(30. / cfreq, 30.)
        uvrange = '<3km'

        if cfreq < 10.:
            imsize = 512
            cell = ['5arcsec']
        else:
            imsize = 1024
            cell = ['2.5arcsec']
        if len(spwran) == 2:
            spwstr = spwran[0] + '~' + spwran[1]
        else:
            spwstr = spwran[0]

        restoringbeam = ['{0:.1f}arcsec'.format(bmsz)]
        imagesuffix = '.spw' + spwstr.replace('~', '-')
        # if cfreq > 10.:
        #     antenna = antenna + ';!0&1;!0&2'  # deselect the shortest baselines
        res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, uvrange=uvrange, spw=spw,
                      ncpu=ncpu, niter=1000, gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=stokes,
                      doreg=True, overwrite=overwrite, toTb=toTb, restoringbeam=restoringbeam, uvtaper=True,
                      outertaper=['30arcsec'], phasecenter=phasecenter)

        if res:
            imres['Succeeded'] += res['Succeeded']
            imres['BeginTime'] += res['BeginTime']
            imres['EndTime'] += res['EndTime']
            imres['ImageName'] += res['ImageName']
            imres['Spw'] += [spwstr] * len(res['ImageName'])
            imres['Vis'] += [msfile] * len(res['ImageName'])
        else:
            return None

    # save it for debugging purposes
    np.savez(os.path.join(imagedir, '{}.imres.npz'.format(os.path.basename(msfile))), imres=imres)

    return imres


def plt_qlook_image(imres, figdir=None, specdata=None, verbose=True, stokes='RR,LL',tojvscript=True):
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    from sunpy import sun
    import astropy.units as u
    if not figdir:
        figdir = './'

    polmap = {'RR': 0, 'LL': 1, 'I': 0, 'V': 1}
    pols = stokes.split(',')
    npols = len(pols)
    # SRL = set(['RR', 'LL'])
    # SXY = set(['XX', 'YY', 'XY', 'YX'])
    nspw = len(set(imres['Spw']))
    plttimes = list(set(imres['BeginTime']))
    ntime = len(plttimes)
    # sort the imres according to time
    images = np.array(imres['ImageName'])
    btimes = Time(imres['BeginTime'])
    etimes = Time(imres['EndTime'])
    spws = np.array(imres['Spw'])
    suc = np.array(imres['Succeeded'])
    inds = btimes.argsort()
    images_sort = images[inds].reshape(ntime, nspw)
    btimes_sort = btimes[inds].reshape(ntime, nspw)
    suc_sort = suc[inds].reshape(ntime, nspw)
    spws_sort = spws[inds].reshape(ntime, nspw)
    if verbose:
        print '{0:d} figures to plot'.format(ntime)
    plt.ioff()
    import matplotlib.gridspec as gridspec
    spec = specdata['spec']
    (npol, nbl, nfreq, ntim) = spec.shape
    tidx = range(ntim)
    fidx = range(nfreq)
    tim = specdata['tim']
    freq = specdata['freq']
    freqg = freq / 1e9
    pol = ''.join(pols)
    spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')
    timstrr = spec_tim.plot_date
    if npols == 1:
        if pol == 'RR':
            spec_plt = spec[0, 0, :, :]
        elif pol == 'LL':
            spec_plt = spec[1, 0, :, :]
        elif pol == 'I':
            spec_plt = (spec[0, 0, :, :] + spec[1, 0, :, :]) / 2.
        elif pol == 'V':
            spec_plt = (spec[0, 0, :, :] - spec[1, 0, :, :]) / 2.
        spec_plt = [spec_plt]
        print 'plot the dynamic spectrum in pol ' + pol  # ax1 = fig.add_subplot(211)

        hnspw = nspw / 2
        ncols = hnspw
        nrows = 2 + 2  # 1 image: 1x1, 1 dspec:2x4
        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(nrows, ncols)
        axs = [plt.subplot(gs[0, 0])]
        for ll in range(1, nspw):
            axs.append(plt.subplot(gs[ll / hnspw, ll % hnspw], sharex=axs[0], sharey=axs[0]))
        for ll in range(nspw):
            axs.append(plt.subplot(gs[ll / hnspw + 2, ll % hnspw], sharex=axs[0], sharey=axs[0]))
        axs_dspec = [plt.subplot(gs[2:, :])]
    elif npols == 2:
        R_plot = np.absolute(spec[0, 0, :, :])
        L_plot = np.absolute(spec[1, 0, :, :])
        if pol == 'RRLL':
            spec_plt = [R_plot, L_plot]
            polstr = ['RR', 'LL']
        if pol == 'IV':
            I_plot = (R_plot + L_plot) / 2.
            V_plot = (R_plot - L_plot) / 2.
            spec_plt = [I_plot, V_plot]
            polstr = ['I', 'V']
        print 'plot the dynamic spectrum in pol ' + pol

        hnspw = nspw / 2
        ncols = hnspw + 2  # 1 image: 1x1, 1 dspec:2x2
        nrows = 2 + 2
        fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(nrows, ncols)
        axs = [plt.subplot(gs[0, 0])]
        for ll in range(1, nspw):
            axs.append(plt.subplot(gs[ll / hnspw, ll % hnspw], sharex=axs[0], sharey=axs[0]))
        for ll in range(nspw):
            axs.append(plt.subplot(gs[ll / hnspw + 2, ll % hnspw], sharex=axs[0], sharey=axs[0]))
        axs_dspec = [plt.subplot(gs[:2, hnspw:])]
        axs_dspec.append(plt.subplot(gs[2:, hnspw:], sharex=axs_dspec[0], sharey=axs_dspec[0]))

    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    for i in range(ntime):
        plt.ioff()
        # plt.clf()
        for ax in axs:
            ax.cla()
        plttime = btimes_sort[i, 0]
        # tofd = plttime.mjd - np.fix(plttime.mjd)
        suci = suc_sort[i]
        # if tofd < 16. / 24. or sum(
        #         suci) < nspw - 2:  # if time of the day is before 16 UT (and 24 UT), skip plotting (because the old antennas are not tracking)
        #     continue
        # fig=plt.figure(figsize=(9,6))
        # fig.suptitle('EOVSA @ '+plttime.iso[:19])
        fig.text(0.01, 0.98, plttime.iso[:19], color='w', fontweight='bold', fontsize=12, ha='left', va='top')
        if verbose:
            print 'Plotting image at: ', plttime.iso

        if i == 0:
            dspecvspans = []
            for pol in range(npols):
                ax = axs_dspec[pol]
                ax.pcolormesh(timstrr, freqg, spec_plt[pol], cmap='jet')
                ax.xaxis_date()
                ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
                # plt.xticks(rotation=45)
                ax.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
                ax.set_ylim(freqg[fidx[0]], freqg[fidx[-1]])
                ax.set_xlabel('Time [UT]')
                ax.set_ylabel('Frequency [GHz]')
                ax.text(0.01, 0.98, pols[pol], color='w', transform=ax.transAxes, fontweight='bold', ha='left',
                        va='top')
                dspecvspans.append(ax.axvspan(btimes[i].plot_date, etimes[i].plot_date, color='w', alpha=0.5))
                ax_pos = ax.get_position().extents
                x0, y0, x1, y1 = ax_pos
                h, v = x1 - x0, y1 - y0
                x0_new = x0 + 0.15 * h
                y0_new = y0 + 0.15 * v
                x1_new = x1 - 0.05 * h
                y1_new = y1 - 0.05 * v
                ax.set_position(mpl.transforms.Bbox([[x0_new, y0_new], [x1_new, y1_new]]))
        else:
            for pol in range(npols):
                xy = dspecvspans[pol].get_xy()
                xy[:, 0][np.array([0, 1, 4])] = btimes[i].plot_date
                xy[:, 0][np.array([2, 3])] = etimes[i].plot_date
                dspecvspans[pol].set_xy(xy)

        for n in range(nspw):
            image = images_sort[i, n]
            # fig.add_subplot(nspw/3, 3, n+1)
            # fig.add_subplot(2, nspw / 2, n + 1)
            for pol in range(npols):
                if suci[n]:
                    try:
                        eomap = smap.Map(image)
                    except:
                        continue
                    sz = eomap.data.shape
                    # if eomap.meta['naxis4'] == 2 and len(SRL.intersection(Spols)) == 2:
                    #     pols = pols + ['I', 'V']
                    # if eomap.meta['naxis4'] == 4 and len(SXY.intersection(Spols)) == 4:
                    #     pols = pols + ['I', 'V']
                    # pdb.set_trace()
                    if len(sz) == 4:
                        eomap.data = eomap.data[max(polmap[pols[pol]], eomap.meta['naxis4'] - 1), 0, :, :].reshape(
                            (sz[2], sz[3]))
                    # resample the image for plotting
                    dim = u.Quantity([256, 256], u.pixel)
                    eomap = eomap.resample(dim)
                    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                    eomap.plot(axes=axs[n + nspw * pol])
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax = plt.gca()
                    ax.set_xlim([-1080, 1080])
                    ax.set_ylim([-1080, 1080])
                    spwran = spws_sort[i, n]
                    # freqran = [int(s) * 0.5 + 2.9 for s in spwran.split('~')]
                    # if len(freqran) == 1:
                    #     ax.text(0.98, 0.01, '{0:.1f} GHz'.format(freqran[0]), color='w',
                    #             transform=ax.transAxes, fontweight='bold', ha='right')
                    # else:
                    #     ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(freqran[0], freqran[1]), color='w',
                    #             transform=ax.transAxes, fontweight='bold', ha='right')
                    ax.text(0.98, 0.01, '{1}@{0:.3f} GHz'.format(eomap.meta['crval3'] / 1e9, pols[pol]), color='w',
                            transform=ax.transAxes, fontweight='bold', ha='right')
                    ax.set_title(' ')
                    # ax.set_title('spw '+spws_sort[i,n])
                    # ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                else:
                    # make an empty map
                    data = np.zeros((512, 512))
                    header = {"DATE-OBS": plttime.isot, "EXPTIME": 0., "CDELT1": 5., "NAXIS1": 512, "CRVAL1": 0.,
                              "CRPIX1": 257, "CUNIT1": "arcsec", "CTYPE1": "HPLN-TAN", "CDELT2": 5., "NAXIS2": 512,
                              "CRVAL2": 0., "CRPIX2": 257, "CUNIT2": "arcsec", "CTYPE2": "HPLT-TAN",
                              "HGLT_OBS": sun.heliographic_solar_center(plttime)[1].value, "HGLN_OBS": 0.,
                              "RSUN_OBS": sun.solar_semidiameter_angular_size(plttime).value,
                              "RSUN_REF": sun.constants.radius.value,
                              "DSUN_OBS": sun.sunearth_distance(plttime).to(u.meter).value, }
                    eomap = smap.Map(data, header)
                    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                    eomap.plot(axes=axs[n + nspw * pol])
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax = plt.gca()
                    ax.set_xlim([-1080, 1080])
                    ax.set_ylim([-1080, 1080])
                    # ax.set_title('spw '+spwran+'( )'))
                    spwran = spws_sort[i, n]
                    freqran = [int(s) * 0.5 + 2.9 for s in spwran.split('~')]
                    spwran = spws_sort[i, n]
                    # ax.set_title('{0:.1f} - {1:.1f} GHz'.format(freqran[0],freqran[1]))
                    # ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(freqran[0], freqran[1]), color='w',
                    #         transform=ax.transAxes, fontweight='bold', ha='right')
                    ax.text(0.98, 0.01, '{1}@{0:.3f} GHz'.format(eomap.meta['crval3'] / 1e9, pols[pol]), color='w',
                            transform=ax.transAxes, fontweight='bold', ha='right')
                    ax.set_title(' ')

                    # ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
        figname = eomap.meta['telescop']+'_qlimg_' + plttime.isot.replace(':', '').replace('-', '')[:19] + '.png'
        fig_tdt = plttime.to_datetime()
        # fig_subdir = fig_tdt.strftime("%Y/%m/%d/")
        figdir_ = figdir  # + fig_subdir
        if not os.path.exists(figdir_):
            os.makedirs(figdir_)
        if verbose:
            print 'Saving plot to :' + os.path.join(figdir_, figname)
        plt.savefig(os.path.join(figdir_, figname))
    plt.close(fig)
    DButil.img2html_movie(figdir_)


def dspec_external(vis, workdir='./', specfile=None):
    dspecscript = os.path.join(workdir, 'dspec.py')
    if not specfile:
        specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
    os.system('rm -rf {}'.format(dspecscript))
    fi = open(dspecscript, 'wb')
    fi.write('from suncasa.utils import dspec2 as ds \n')
    fi.write('specdata = ds.get_dspec("{0}", specfile="{1}", domedian=True, verbose=True, savespec=True) \n'.format(vis,
                                                                                                                    specfile))
    fi.close()
    os.system('casa --nologger -c {}'.format(dspecscript))


def svplot(vis, timerange=None, freqrange='', workdir='./', specfile=None, stokes='RR,LL', dmin=None, dmax=None,
           goestime=None, reftime=None, fov=None, aiawave=171, imagefile=None, savefig=False, aiafits=None,
           changeheader=True, redoclean=False, fitsfile=None, mkmovie=False, overwrite=True, twidth=1):
    '''
    Required inputs:
            vis: calibrated CASA measurement set
            timerange: timerange for clean. Standard CASA time selection format
            freqrange: frequency selection for clean. example: '1.25~1.35 GHz' or '1250~1350 MHz'
    Optional inputs:
            specfile: supply dynamic spectrum save file (from suncasa.utils.dspec2.get_dspec()). Otherwise
                      generate on the fly
            pol: pol of the dynamic spectrum,can be 'RR','LL','I','V','IV','RRLL', default is 'RRLL'
            dmin,dmax: color bar parameter
            goestime: goes plot time, example ['2016/02/18 18:00:00','2016/02/18 23:00:00']
            rhessisav: rhessi savefile
            reftime: reftime for the image
            fov: field of view in aia image, in unit of arcsec, example:[[-400,-200],[100,300]]
            aiawave: wave length of aia file in a
            imagefile: cleaned image file
            fitsfile: exist vla fitsfile
            savefig: whether to save the figure
    Example:

    '''
    # first part of dyn spec

    pol = ''.join(stokes.split(','))
    if pol != 'RR' and pol != 'LL' and pol != 'I' and pol != 'V' and pol != 'RRLL' and pol != 'IV':
        print 'wrong pol(only LL,RR,RRLL,I,V and IV)'
        return 0

    if vis[-1] == '/':
        vis = vis[:-1]
    if not os.path.exists(vis):
        print 'input measurement not exist'
        return -1
    if aiafits is None:
        aiafits = ''
    # split the data
    # generating dynamic spectrum
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    if specfile:
        try:
            specdata = np.load(specfile)
        except:
            print('Provided dynamic spectrum file not numpy npz. Generating one from the visibility data')
            specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
            dspec_external(vis, workdir=workdir, specfile=specfile)
            specdata = np.load(specfile)  # specdata = ds.get_dspec(vis, domedian=True, verbose=True)
    else:
        print('Dynamic spectrum file not provided; Generating one from the visibility data')
        # specdata = ds.get_dspec(vis, domedian=True, verbose=True)
        specfile = os.path.join(workdir, os.path.basename(vis) + '.dspec.npz')
        dspec_external(vis, workdir=workdir, specfile=specfile)
        specdata = np.load(specfile)

    tb.open(vis)
    starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
    tb.close()
    datstr = starttim.iso[:10]

    try:
        (tstart, tend) = timerange.split('~')
        if tstart[2] == ':':
            starttim1 = Time(datstr + 'T' + tstart)
            endtim1 = Time(datstr + 'T' + tend)
            timerange = '{0}/{1}~{0}/{2}'.format(datstr.replace('-', '/'), tstart, tend)
        else:
            starttim1 = Time(qa.quantity(tstart, 'd')['value'], format='mjd')
            endtim1 = Time(qa.quantity(tend, 'd')['value'], format='mjd')
    except ValueError:
        print "keyword 'timerange' in wrong format"
    # I do not really like the following section. All we need here is a timerange and frequency range, which does
    # not require splitting and gridding the data at all.  --Bin
    if vis.endswith('/'):
        vis = vis[:-1]
    visname = os.path.basename(vis)
    ret1 = starttim1.plot_date
    ret2 = endtim1.plot_date

    ms.open(vis)
    axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    freqInfo_ravel = freqInfo.ravel()
    ms.close()
    if freqrange != '':
        unit = freqrange.split(' ')[1]
        freq0, freq1 = freqrange.split(' ')[0].split('~')
        if unit == 'MHz':
            freq0, freq1 = float(freq0) / 1e3, float(freq1) / 1e3
        elif unit == 'GHz':
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
            ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1] - 1)] + ['0~{}'.format(sz_freqInfo[1] - 1) for ll
                                                                              in xrange(freqIdx0[0] + 1, freqIdx1[0])]
            ms_chan.append('0~{}'.format(freqIdx1[1][0]))
        spw = ','.join('{}:{}'.format(t[0], t[1]) for t in zip(ms_spw, ms_chan))
        req1, req2 = freq0, freq1
    else:
        spw = ''
        req1 = freqInfo_ravel[0]
        req2 = freqInfo_ravel[-1]

    midtime_mjd = (starttim1.mjd + endtim1.mjd) / 2.

    if mkmovie:
        plt.ioff()
        # fig = plt.figure(figsize=(12, 7.5), dpi=100)
        if fitsfile:
            pass
        else:
            if not imagefile:
                # from ptclean_cli import ptclean_cli as ptclean
                eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
                phasecenter = 'J2000 ' + str(eph['ra'][0])[:15] + 'rad ' + str(eph['dec'][0])[:15] + 'rad'
                print 'use phasecenter: ' + phasecenter
                qlookfigdir = os.path.join(workdir, 'qlookfits/')
                imresfile = os.path.join(qlookfigdir, '{}.imres.npz'.format(os.path.basename(vis)))
                if overwrite:
                    imres = mk_qlook_image(vis, twidth=twidth, ncpu=8, imagedir=qlookfigdir, phasecenter=phasecenter,
                                           stokes=stokes)
                else:
                    if os.path.exists(imresfile):
                        imres = np.load(imresfile)
                        imres = imres['imres'].item()
                    else:
                        print('Image results file not found; Creating new images.')
                        imres = mk_qlook_image(vis, twidth=twidth, ncpu=8, imagedir=qlookfigdir,
                                               phasecenter=phasecenter, stokes=stokes)
                if not os.path.exists(qlookfigdir):
                    os.makedirs(qlookfigdir)
                plt_qlook_image(imres, figdir=qlookfigdir, specdata=specdata, verbose=True, stokes=stokes)

    else:
        spec = specdata['spec']
        (npol, nbl, nfreq, ntim) = spec.shape
        tidx = range(ntim)
        fidx = range(nfreq)
        tim = specdata['tim']
        freq = specdata['freq']
        freqg = freq / 1e9
        spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')
        timstrr = spec_tim.plot_date
        plt.ion()
        fig = plt.figure(figsize=(10, 7.5), dpi=100)
        if pol != 'RRLL' and pol != 'IV':
            if pol == 'RR':
                spec_plt = spec[0, 0, :, :]
            elif pol == 'LL':
                spec_plt = spec[1, 0, :, :]
            elif pol == 'I':
                spec_plt = (spec[0, 0, :, :] + spec[1, 0, :, :]) / 2.
            elif pol == 'V':
                spec_plt = (spec[0, 0, :, :] - spec[1, 0, :, :]) / 2.

            print 'plot the dynamic spectrum in pol ' + pol
            ax1 = fig.add_subplot(221)
            ax1.pcolormesh(timstrr, freqg, spec_plt, cmap='jet')
            # f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
            ax1.add_patch(patches.Rectangle((ret1, req1), ret2 - ret1, req2 - req1, alpha=0.4))
            ax1.xaxis_date()
            ax1.xaxis.set_major_formatter(DateFormatter("%H:%M"))
            ax1.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
            ax1.set_ylim(freqg[fidx[0]], freqg[fidx[-1]])

            ax1.set_ylabel('Frequency (GHz)', fontsize=10)
            ax1.set_title('Dynamic spectrum @ pol ' + pol, fontsize=12)
            for tick in ax1.get_xticklabels():
                # tick.set_fontname('Comic Sans MS')
                tick.set_fontsize(8)
            for tick in ax1.get_yticklabels():
                tick.set_fontsize(8)
            ax1.set_autoscale_on(False)
        else:
            R_plot = np.absolute(spec[0, 0, :, :])
            L_plot = np.absolute(spec[1, 0, :, :])
            if pol == 'RRLL':
                spec_plt_1 = R_plot
                spec_plt_2 = L_plot
                polstr = ['RR', 'LL']
            if pol == 'IV':
                I_plot = (R_plot + L_plot) / 2.
                V_plot = (R_plot - L_plot) / 2.
                spec_plt_1 = I_plot
                spec_plt_2 = V_plot
                polstr = ['I', 'V']

            print 'plot the dynamic spectrum in pol ' + pol
            ax1 = fig.add_subplot(321)
            ax1.pcolormesh(timstrr, freqg, spec_plt_1, cmap='jet', vmin=dmin, vmax=dmax)
            ax1.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
            ax1.xaxis_date()
            ax1.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
            ax1.set_ylim(freqg[fidx[0]], freqg[fidx[-1]])
            ax1.set_ylabel('Frequency (GHz)', fontsize=10)
            ax1.set_title('Dynamic spectrum @pol ' + polstr[0] + ' in upper and ' + polstr[1] + ' in bottom',
                          fontsize=12)
            ax1.set_autoscale_on(False)
            # f1.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
            ax1.add_patch(patches.Rectangle((ret1, req1), ret2 - ret1, req2 - req1, alpha=0.4))
            for tick in ax1.get_xticklabels():
                tick.set_fontsize(8)
            for tick in ax1.get_yticklabels():
                tick.set_fontsize(8)
            ax2 = fig.add_subplot(323)
            ax2.pcolormesh(timstrr, freqg, spec_plt_2, cmap='jet', vmin=dmin, vmax=dmax)
            ax2.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
            ax2.xaxis_date()
            ax2.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
            ax2.set_ylim(freqg[fidx[0]], freqg[fidx[-1]])
            ax2.set_ylabel('Frequency (GHz)', fontsize=10)
            for tick in ax2.get_xticklabels():
                tick.set_fontsize(8)
            for tick in ax2.get_yticklabels():
                tick.set_fontsize(8)
            ax2.set_autoscale_on(False)
            # f2.add_patch(patches.Rectangle((ret1, req1),ret2-ret1,req2-req1,fill=False))
            ax2.add_patch(patches.Rectangle((ret1, req1), ret2 - ret1, req2 - req1, alpha=0.4))

        # second part of Goes plot
        print 'plot the Goes soft X-ray curve and derivate'
        if goestime:
            goesplottim = TimeRange(goestime[0], goestime[1])
        else:
            datstrg = datstr.replace('-', '/')
            goestime0 = datstrg + ' ' + qa.time(qa.quantity(tim[0] - 1800, 's'), form='clean', prec=9)[0]
            goestime1 = datstrg + ' ' + qa.time(qa.quantity(tim[tidx[-1] - 1] + 1800, 's'), form='clean', prec=9)[0]
            goesplottim = TimeRange(goestime0, goestime1)

        if pol == 'RRLL' or pol == 'IV':
            ax3 = fig.add_subplot(427)
        else:
            ax3 = fig.add_subplot(223)

        goesscript = os.path.join(workdir, 'goes.py')
        goesdatafile = os.path.join(workdir, 'goes.dat')
        os.system('rm -rf {}'.format(goesscript))
        fi = open(goesscript, 'wb')
        fi.write('import os \n')
        fi.write('from sunpy.time import TimeRange \n')
        fi.write('from sunpy import lightcurve as lc \n')
        fi.write('import pickle \n')
        fi.write('goes = lc.GOESLightCurve.create(goesplottim) \n')
        fi.write('fi2 = open(goesdatafile, "wb") \n')
        fi.write('pickle.dump(goes, fi2) \n')
        fi.write('fi2.close()')
        fi.close()

        try:
            execfile(goesscript)

            fi1 = file(goesdatafile, 'rb')
            goest = pickle.load(fi1)
            fi1.close()

            dates = mpl.dates.date2num(parse_time(goest.data.index))
            goesdif = np.diff(goest.data['xrsb'])
            gmax = np.nanmax(goesdif)
            gmin = np.nanmin(goesdif)
            ra = gmax - gmin
            db = 3e-4 / ra
            goesdifp = goesdif * db + gmin + 1e-4
            ax3.plot_date(dates, goest.data['xrsb'], '-', label='1.0--8.0 $\AA$', color='red', lw=2)
            ax3.plot_date(dates[0:-1], goesdifp, '-', label='derivate', color='blue', lw=0.4)

            ax3.set_yscale("log")
            ax3.set_ylim(1e-7, 1e-3)
            ax3.set_title('Goes Soft X-ray and derivate')
            ax3.set_ylabel('Watts m$^{-2}$')
            ax3.set_xlabel(datetime.datetime.isoformat(goest.data.index[0])[0:10])
            ax3.axvspan(dates[899], dates[dates.size - 899], alpha=0.2)

            ax2 = ax3.twinx()
            ax2.set_yscale("log")
            ax2.set_ylim(1e-7, 1e-3)
            ax2.set_yticks((1e-7, 1e-6, 1e-5, 1e-4, 1e-3))
            ax2.set_yticklabels(('B', 'C', 'M', 'X'))

            # ax3 = f3.twinx()
            # # ax3.set_yscale("linear")
            # ax3.plot_date(dates[0:-1], goesdifp, '-', label='derivate', color='blue', lw=0.4)
            # ax3.set_visible(False)

            ax3.yaxis.grid(True, 'major')
            ax3.xaxis.grid(False, 'major')
            ax3.legend(prop={'size': 6})

            formatter = mpl.dates.DateFormatter('%H:%M')
            ax3.xaxis.set_major_formatter(formatter)

            ax3.fmt_xdata = mpl.dates.DateFormatter('%H:%M')
        except:
            print 'error in plotting GOES. Continue without GOES'

        # third part
        # start to download the fits files
        if not aiafits:
            newlist = []
            items = glob.glob('*.fits')
            for names in items:
                str1 = starttim1.iso[:4] + '_' + starttim1.iso[5:7] + '_' + starttim1.iso[8:10] + 't' + starttim1.iso[
                                                                                                        11:13] + '_' + starttim1.iso[
                                                                                                                       14:16]
                str2 = str(aiawave)
                if names.endswith(".fits"):
                    if names.find(str1) != -1 and names.find(str2) != -1:
                        newlist.append(names)
                newlist.append('0')
            if os.path.exists(newlist[0]):
                aiafits = newlist[0]
            else:
                print 'downloading the aiafits file'
                client = vso.VSOClient()
                wave1 = aiawave - 3
                wave2 = aiawave + 3
                t1 = Time(starttim1.mjd - 0.02 / 24., format='mjd')
                t2 = Time(endtim1.mjd + 0.02 / 24., format='mjd')
                qr = client.query(vso.attrs.Time(t1.iso, t2.iso), vso.attrs.Instrument('aia'),
                                  vso.attrs.Wave(wave1 * u.AA, wave2 * u.AA))
                res = client.get(qr, path='{file}')
        # Here something is needed to check whether it has finished downloading the fits files or not

        if not aiafits:
            newlist = []
            items = glob.glob('*.fits')
            for nm in items:
                str1 = starttim1.iso[:4] + '_' + starttim1.iso[5:7] + '_' + starttim1.iso[8:10] + 't' + starttim1.iso[
                                                                                                        11:13] + '_' + starttim1.iso[
                                                                                                                       14:16]
                str2 = str(aiawave)
                if nm.find(str1) != -1 and nm.find(str2) != -1:
                    newlist.append(nm)
            if newlist:
                aiafits = newlist[0]
                print 'AIA fits ' + aiafits + ' selected'
            else:
                print 'no AIA fits files found. Proceed without AIA'

        try:
            aiamap = smap.Map(aiafits)
        except:
            print 'error in reading aiafits. Proceed without AIA'

        ax4 = fig.add_subplot(222)
        ax5 = fig.add_subplot(224)

        if fitsfile:
            pass
        else:
            if not imagefile:
                eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
                phasecenter = 'J2000 ' + str(eph['ra'][0])[:15] + 'rad ' + str(eph['dec'][0])[:15] + 'rad'
                imagename = os.path.join(workdir, visname)
                if os.path.exists(imagename + '.image') or os.path.exists(imagename + '.flux'):
                    os.system('rm -rf ' + imagename + '*')
                print 'do clean for ' + timerange + ' in spw ' + spw
                print 'use phasecenter: ' + phasecenter
                clean(vis=vis, imagename=imagename, selectdata=True, spw=spw, timerange=timerange, niter=500,
                      interactive=False, npercycle=50, imsize=[512, 512], cell=['5.0arcsec'], phasecenter=phasecenter)
                os.system('rm -rf ' + imagename + '.psf')
                os.system('rm -rf ' + imagename + '.flux')
                os.system('rm -rf ' + imagename + '.model')
                os.system('rm -rf ' + imagename + '.mask')
                os.system('rm -rf ' + imagename + '.residual')
                imagefile = imagename + '.image'
            fitsfile = imagefile + '.fits'
            hf.imreg(vis=vis, ephem=eph, imagefile=imagefile, timerange=timerange, reftime=reftime, fitsfile=fitsfile,
                     verbose=True)
            if changeheader:
                data, header = fits.getdata(fitsfile, header=True)
                header['date-obs'] = starttim1.iso[:10] + 'T' + starttim1.iso[12:]
                fits.writeto(fitsfile, data, header, clobber=True)
        print 'fits file ' + fitsfile + ' selected'
        ax4.cla()
        ax5.cla()

        vlafits = fitsfile
        vlamap = smap.Map(vlafits)
        vlamap.data = vlamap.data.reshape(vlamap.meta['naxis1'], vlamap.meta['naxis2'])
        # vlamap = vlamap.submap([-1200, 1200] * u.arcsec, [-1200, 1200] * u.arcsec)
        XX, YY = np.meshgrid(np.arange(vlamap.data.shape[1]), np.arange(vlamap.data.shape[0]))
        vlamapx, vlamapy = vlamap.pixel_to_data(XX * u.pix, YY * u.pix)

        if fov:
            fov = fov
        else:
            raw, col = vlamap.data.shape
            positon = np.nanargmax(vlamap.data)
            m, n = divmod(positon, col)
            length = 200 * u.arcsec
            x0 = vlamap.xrange[0] + vlamap.scale[1] * (n + 0.5) * u.pix
            y0 = vlamap.yrange[0] + vlamap.scale[0] * (m + 0.5) * u.pix
            x1 = x0 - length
            x2 = x0 + length
            y1 = y0 - length
            y2 = y0 + length
            fov = [[x1.value, x2.value], [y1.value, y2.value]]

        vmax, vmin = np.nanmax(vlamap.data), np.nanmin(vlamap.data)
        if aiafits:
            aiamap.plot(axes=ax4, title='overview image')
            ax4.contour(vlamapx.value, vlamapy.value, vlamap.data,
                        levels=np.linspace(0.5, 0.9, 3) * np.nanmax(vlamap.data), cmap=cm.jet)
            aiamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
        else:
            vlamap.plot(axes=ax4, title='overview image', cmap=cm.jet, vmax=vmax, vmin=(vmax + vmin) / 2.0)
            vlamap.draw_limb()
            vlamap.draw_grid()
            ax4.contour(vlamapx.value, vlamapy.value, vlamap.data,
                        levels=np.linspace(0.5, 0.9, 3) * np.nanmax(vlamap.data), cmap=cm.gray)
            vlamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
        ax4.set_xlim(-1200, 1200)
        ax4.set_ylim(-1200, 1200)

        # cmap = smap.CompositeMap(aiamap)
        # cmap.add_map(vlamap, levels=np.array([0.5, 0.7, 0.9]) * np.nanmax(vlamap.data))
        # cmap.set_colors(1, cm=cm.rainbow)
        # f4 = cmap.plot(title='overview image')
        # f4[2].set_autoscale_on(False)  gca

        subvlamap = vlamap.submap(fov[0] * u.arcsec, fov[1] * u.arcsec)
        XX, YY = np.meshgrid(np.arange(subvlamap.data.shape[1]), np.arange(subvlamap.data.shape[0]))
        subvlamapx, subvlamapy = subvlamap.pixel_to_data(XX * u.pix, YY * u.pix)
        if aiafits:
            subaiamap = aiamap.submap(fov[0] * u.arcsec, fov[1] * u.arcsec)
            subaiamap.plot(axes=ax5, title='zoomed view')
            ax5.contour(subvlamapx.value, subvlamapy.value, subvlamap.data,
                        levels=np.linspace(0.5, 0.9, 3) * np.nanmax(subvlamap.data), cmap=cm.jet)
            subaiamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
            ax5.text(0.01, 0.02,
                     'AIA {} @ '.format(aiamap.wavelength.value) + aiamap.date.strftime('%Y %h %d %H:%M:%S'),
                     verticalalignment='bottom', horizontalalignment='left', transform=ax5.transAxes, color='w',
                     fontsize=8)
        else:
            subvlamap.plot(axes=ax5, title='zoomed view', cmap=cm.jet, vmax=vmax, vmin=(vmax + vmin) / 2.0)
            subvlamap.draw_limb()
            subvlamap.draw_grid()
            ax5.contour(subvlamapx.value, subvlamapy.value, subvlamap.data,
                        levels=np.linspace(0.5, 0.9, 3) * np.nanmax(subvlamap.data), cmap=cm.gray)
            subvlamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
        ax5.set_xlim(fov[0])
        ax5.set_ylim(fov[1])
        ax5.text(0.01, 0.06,
                 'VLA @ {} GHz '.format((req1 + req2) / 2.0) + vlamap.date.strftime('%Y %h %d %H:%M:%S.%f')[:-3],
                 verticalalignment='bottom', horizontalalignment='left', transform=ax5.transAxes, color='w', fontsize=8)

        # scmap = smap.CompositeMap(subaiamap)
        # scmap.add_map(subvlamap, levels=np.array([0.5, 0.7, 0.9]) * np.nanmax(subvlamap.data))
        # scmap.set_colors(1, cm=cm.rainbow)
        # f5 = scmap.plot(title='zoomed view')
        fig.show()

        # os.system('rm -rf '+vis_sp)
        os.system('rm -rf {}'.format(goesscript))
        os.system('rm -rf {}'.format(goesdatafile))
