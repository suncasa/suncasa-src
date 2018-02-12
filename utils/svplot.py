import numpy as np
import matplotlib.pyplot as plt
import os, sys
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
from astropy.io import fits
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
from matplotlib import gridspec
import glob
from suncasa.utils import DButil
import copy
import pdb


def uniq(lst):
    last = object()
    nlst = []
    for item in lst:
        if item == last:
            continue
        nlst.append(item)
        last = item
    return nlst


def mk_qlook_image(vis, ncpu=10, twidth=12, stokes='I,V', antenna='', imagedir=None, spws=[], toTb=True, overwrite=True, doslfcal=False,
                   phasecenter='', c_external=True):
    vis = [vis]
    subdir = ['/']

    for idx, f in enumerate(vis):
        if f[-1] == '/':
            vis[idx] = f[:-1]

    if not imagedir:
        imagedir = './'
    imres = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': [], 'Freq': []}
    msfile = vis[0]
    ms.open(msfile)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    # axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    # freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    # freqInfo_ravel = freqInfo.ravel()
    ms.close()

    if not spws:
        if observatory == 'EVLA':
            spws = ['0', '1', '2', '3', '4', '5', '6', '7']
        if observatory == 'EOVSA':
            spws = ['1~5', '6~10', '11~15', '16~25']
    if observatory == 'EOVSA':
        print 'Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be IV'
        stokes = 'I,V'

    msfilebs = os.path.basename(msfile)
    imdir = imagedir + subdir[0]
    if not os.path.exists(imdir):
        os.makedirs(imdir)
    if doslfcal:
        slfcalms = './' + msfilebs + '.rr'
        split(msfile, outputvis=slfcalms, datacolumn='corrected', correlation='RR')
    for spw in spws:
        spwran = [s.zfill(2) for s in spw.split('~')]
        # freqran = [
        #     (int(s) * spwInfo['0']['TotalWidth'] + spwInfo['0']['RefFreq'] + spwInfo['0']['TotalWidth'] / 2.0) / 1.0e9
        #     for s in spw.split('~')]
        spw_ = spw.split('~')
        if len(spw_) == 2:
            freqran = [(spwInfo['{}'.format(s)]['RefFreq'] + spwInfo['{}'.format(s)]['TotalWidth'] / 2.0) / 1.0e9 for s in spw.split('~')]
        elif len(spw_) == 1:
            s = spw_[0]
            freqran = np.array([0, spwInfo['{}'.format(s)]['TotalWidth']]) + spwInfo['{}'.format(s)]['RefFreq']
            freqran = freqran / 1.0e9
            freqran = list(freqran)
        else:
            raise ValueError("Keyword 'spw' in wrong format")

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
        sto = stokes.replace(',', '')
        if c_external:
            cleanscript = os.path.join(imdir, 'ptclean_external.py')
            resfile = os.path.join(imdir, os.path.basename(msfile) + '.res.npz')
            os.system('rm -rf {}'.format(cleanscript))
            inpdict = {'vis': msfile, 'imageprefix': imdir, 'imagesuffix': imagesuffix, 'twidth': twidth, 'uvrange': uvrange, 'spw': spw,
                       'ncpu': ncpu, 'niter': 1000, 'gain': 0.05, 'antenna': antenna, 'imsize': imsize, 'cell': cell, 'stokes': sto, 'doreg': True,
                       'overwrite': overwrite, 'toTb': toTb, 'restoringbeam': restoringbeam, 'uvtaper': True, 'outertaper': ['30arcsec'],
                       'phasecenter': phasecenter}
            for key, val in inpdict.items():
                if type(val) is str:
                    inpdict[key] = '"{}"'.format(val)
            fi = open(cleanscript, 'wb')
            fi.write('from ptclean_cli import ptclean_cli as ptclean \n')
            fi.write('import numpy as np \n')
            fi.write(
                'res = ptclean(vis={i[vis]},imageprefix={i[imageprefix]},imagesuffix={i[imagesuffix]},twidth={i[twidth]},uvrange={i[uvrange]},spw={i[spw]},ncpu={i[ncpu]},niter={i[niter]},gain={i[gain]},antenna={i[antenna]},imsize={i[imsize]},cell={i[cell]},stokes={i[stokes]},doreg={i[doreg]},overwrite={i[overwrite]},toTb={i[toTb]},restoringbeam={i[restoringbeam]},uvtaper={i[uvtaper]},outertaper={i[outertaper]},phasecenter={i[phasecenter]}) \n'.format(
                    i=inpdict))
            fi.write('np.savez("{}",res=res) \n'.format(resfile))
            fi.close()

            os.system('casa --nologger -c {}'.format(cleanscript))
            res = np.load(resfile)
            res = res['res'].item()
        else:
            res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, uvrange=uvrange, spw=spw, ncpu=ncpu, niter=1000,
                          gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=sto, doreg=True, overwrite=overwrite, toTb=toTb,
                          restoringbeam=restoringbeam, uvtaper=True, outertaper=['30arcsec'], phasecenter=phasecenter)

        if res:
            imres['Succeeded'] += res['Succeeded']
            imres['BeginTime'] += res['BeginTime']
            imres['EndTime'] += res['EndTime']
            imres['ImageName'] += res['ImageName']
            imres['Spw'] += [spwstr] * len(res['ImageName'])
            imres['Vis'] += [msfile] * len(res['ImageName'])
            imres['Freq'] += [freqran] * len(res['ImageName'])
        else:
            return None

    # save it for debugging purposes
    np.savez(os.path.join(imagedir, '{}.imres.npz'.format(os.path.basename(msfile))), imres=imres)

    return imres


def plt_qlook_image(imres, figdir=None, specdata=None, verbose=True, stokes='I,V', fov=None):
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    from sunpy import sun
    import astropy.units as u
    if not figdir:
        figdir = './'

    observatory = 'EOVSA'
    polmap = {'RR': 0, 'LL': 1, 'I': 0, 'V': 1}
    pols = stokes.split(',')
    npols = len(pols)
    # SRL = set(['RR', 'LL'])
    # SXY = set(['XX', 'YY', 'XY', 'YX'])
    Spw = sorted(list(set(imres['Spw'])))
    nspw = len(Spw)
    # Freq = set(imres['Freq']) ## list is an unhashable type
    Freq = sorted(uniq(imres['Freq']))

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
    freqghz = freq / 1e9
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
        cmaps = ['jet']
    elif npols == 2:
        R_plot = np.absolute(spec[0, 0, :, :])
        L_plot = np.absolute(spec[1, 0, :, :])
        if pol == 'RRLL':
            spec_plt = [R_plot, L_plot]
            polstr = ['RR', 'LL']
            cmaps = ['jet'] * 2
        if pol == 'IV':
            I_plot = (R_plot + L_plot) / 2.
            V_plot = (R_plot - L_plot) / 2.
            spec_plt = [I_plot, V_plot]
            polstr = ['I', 'V']
            cmaps = ['jet', 'RdBu']
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
    timetext = fig.text(0.01, 0.98, '', color='w', fontweight='bold', fontsize=12, ha='left', va='top')
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
        timetext.set_text(plttime.iso[:19])
        if verbose:
            print 'Plotting image at: ', plttime.iso

        if i == 0:
            dspecvspans = []
            for pol in range(npols):
                ax = axs_dspec[pol]
                ax.pcolormesh(timstrr, freqghz, spec_plt[pol], cmap=cmaps[pol])
                ax.xaxis_date()
                ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
                # plt.xticks(rotation=45)
                ax.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
                ax.set_xlabel('Time [UT]')
                ax.set_ylabel('Frequency [GHz]')
                for idx, freq in enumerate(Freq):
                    ax.axhspan(freq[0], freq[1], linestyle='dotted', edgecolor='w', alpha=0.7, facecolor='none')
                    xtext, ytext = ax.transAxes.inverted().transform(ax.transData.transform([timstrr[tidx[0]], np.mean(freq)]))
                    ax.text(xtext + 0.01, ytext, 'spw ' + Spw[idx], color='w', transform=ax.transAxes, fontweight='bold', ha='left', va='center',
                            fontsize=8, alpha=0.5)
                ax.text(0.01, 0.98, 'Stokes ' + pols[pol], color='w', transform=ax.transAxes, fontweight='bold', ha='left', va='top')
                dspecvspans.append(ax.axvspan(btimes[i].plot_date, etimes[i].plot_date, color='w', alpha=0.4))
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
                    if len(sz) == 4:
                        eomap.data = eomap.data[min(polmap[pols[pol]], eomap.meta['naxis4'] - 1), 0, :, :].reshape((sz[2], sz[3]))
                    # resample the image for plotting
                    if fov:
                        fov = [np.array(ll) for ll in fov]
                        pad = max(np.diff(fov[0])[0], np.diff(fov[1])[0])
                        eomap = eomap.submap((fov[0] + np.array([-1.0, 1.0]) * pad) * u.arcsec, (fov[1] + np.array([-1.0, 1.0]) * pad) * u.arcsec)
                    else:
                        dim = u.Quantity([256, 256], u.pixel)
                        eomap = eomap.resample(dim)
                    eomap.plot_settings['cmap'] = plt.get_cmap(cmaps[pol])
                    eomap.plot(axes=axs[n + nspw * pol])
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax = plt.gca()
                    ax.set_autoscale_on(False)
                    if fov:
                        # pass
                        ax.set_xlim(fov[0])
                        ax.set_ylim(fov[1])
                    else:
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
                    ax.text(0.98, 0.01, 'Stokes {1} @ {0:.3f} GHz'.format(eomap.meta['crval3'] / 1e9, pols[pol]), color='w', transform=ax.transAxes,
                            fontweight='bold', ha='right')
                    ax.set_title(' ')
                    # ax.set_title('spw '+spws_sort[i,n])
                    # ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                else:
                    # make an empty map
                    data = np.zeros((512, 512))
                    header = {"DATE-OBS": plttime.isot, "EXPTIME": 0., "CDELT1": 5., "NAXIS1": 512, "CRVAL1": 0., "CRPIX1": 257, "CUNIT1": "arcsec",
                              "CTYPE1": "HPLN-TAN", "CDELT2": 5., "NAXIS2": 512, "CRVAL2": 0., "CRPIX2": 257, "CUNIT2": "arcsec",
                              "CTYPE2": "HPLT-TAN", "HGLT_OBS": sun.heliographic_solar_center(plttime)[1].value, "HGLN_OBS": 0.,
                              "RSUN_OBS": sun.solar_semidiameter_angular_size(plttime).value, "RSUN_REF": sun.constants.radius.value,
                              "DSUN_OBS": sun.sunearth_distance(plttime).to(u.meter).value, }
                    eomap = smap.Map(data, header)
                    # resample the image for plotting
                    if fov:
                        fov = [np.array(ll) for ll in fov]
                        pad = max(np.diff(fov[0])[0], np.diff(fov[1])[0])
                        try:
                            eomap = eomap.submap((fov[0] + np.array([-1.0, 1.0]) * pad) * u.arcsec, (fov[1] + np.array([-1.0, 1.0]) * pad) * u.arcsec)
                        except:
                            x0, x1 = fov[0] + np.array([-1.0, 1.0]) * pad
                            y0, y1 = fov[1] + np.array([-1.0, 1.0]) * pad
                            bl = SkyCoord(x0 * u.arcsec, y0 * u.arcsec, frame=eomap.coordinate_frame)
                            tr = SkyCoord(x1 * u.arcsec, y1 * u.arcsec, frame=eomap.coordinate_frame)
                            eomap = eomap.submap(bl, tr)
                    else:
                        dim = u.Quantity([256, 256], u.pixel)
                        eomap = eomap.resample(dim)
                    eomap.plot_settings['cmap'] = plt.get_cmap(cmaps[pol])
                    eomap.plot(axes=axs[n + nspw * pol])
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax = plt.gca()
                    ax.set_autoscale_on(False)
                    if fov:
                        # pass
                        ax.set_xlim(fov[0])
                        ax.set_ylim(fov[1])
                    else:
                        ax.set_xlim([-1080, 1080])
                        ax.set_ylim([-1080, 1080])
                    # ax.set_title('spw '+spwran+'( )'))
                    spwran = spws_sort[i, n]
                    freqran = [int(s) * 0.5 + 2.9 for s in spwran.split('~')]
                    spwran = spws_sort[i, n]
                    # ax.set_title('{0:.1f} - {1:.1f} GHz'.format(freqran[0],freqran[1]))
                    # ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(freqran[0], freqran[1]), color='w',
                    #         transform=ax.transAxes, fontweight='bold', ha='right')
                    ax.text(0.98, 0.01, 'Stokes {1} @ {0:.3f} GHz'.format(0., pols[pol]), color='w', transform=ax.transAxes, fontweight='bold',
                            ha='right')
                    ax.set_title(' ')

                    # ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
        figname = observatory + '_qlimg_' + plttime.isot.replace(':', '').replace('-', '')[:19] + '.png'
        fig_tdt = plttime.to_datetime()
        # fig_subdir = fig_tdt.strftime("%Y/%m/%d/")
        figdir_ = figdir  # + fig_subdir
        if not os.path.exists(figdir_):
            os.makedirs(figdir_)
        if verbose:
            print 'Saving plot to: ' + os.path.join(figdir_, figname)
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
    fi.write('specdata = ds.get_dspec("{0}", specfile="{1}", domedian=True, verbose=True, savespec=True) \n'.format(vis, specfile))
    fi.close()
    os.system('casa --nologger -c {}'.format(dspecscript))


def svplot(vis, timerange=None, spw='', workdir='./', specfile=None, stokes='RR,LL', dmin=None, dmax=None, goestime=None, reftime=None, fov=None,
           usephacenter=True, imagefile=None, fitsfile=None, plotaia=True, aiawave=171, aiafits=None, savefig=False, mkmovie=False, overwrite=True,
           ncpu=10, twidth=1, verbose=True):
    '''
    Required inputs:
            vis: calibrated CASA measurement set
    Important optional inputs:
            timerange: timerange for clean. Standard CASA time selection format. 
                       If not provided, use the entire range (*BE CAREFUL, COULD BE VERY SLOW*)
            spw: spectral window selection following the CASA syntax. 
                 Examples: spw='1:2~60' (spw id 1, channel range 2-60); spw='*:1.2~1.3GHz' (selects all channels within 1.2-1.3 GHz; note the *) 
            specfile: supply dynamic spectrum save file (from suncasa.utils.dspec2.get_dspec()). Otherwise
                      generate a median dynamic spectrum on the fly
    Optional inputs:
            bl: baseline to generate dynamic spectrum
            uvrange: uvrange to select baselines for generating dynamic spectrum 
            stokes: polarization of the clean image, can be 'RR,LL' or 'I,V'
            dmin,dmax: color bar parameter
            goestime: goes plot time, example ['2016/02/18 18:00:00','2016/02/18 23:00:00']
            rhessisav: rhessi savefile
            reftime: reftime for the image
            fov: field of view in aia image, in the format of [[xcenter,ycenter],[xlength,ylength]], in unit of arcsec, example:[[-400,-200],[100,300]]
            aiawave: wave length of aia file in a
            imagefile: if imagefile provided, use it. Otherwise do clean and generate a new one.
            fitsfile: if fitsfile provided, use it. Otherwise generate a new one
            savefig: whether to save the figure
    Example:

    '''

    stokes_allowed = ['RR,LL', 'I,V', 'RRLL', 'IV']
    if not stokes in stokes_allowed:
        print 'wrong stokes parameter ' + str(stokes) + '. Allowed values are ' + ', '.join(stokes_allowed)
        return -1
    if stokes == 'RRLL':
        stokes = 'RR,LL'
    if stokes == 'IV':
        stokes = 'I,V'

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
    endtim = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
    tb.close()
    datstr = starttim.iso[:10]

    if timerange is None or timerange == '':
        starttim1 = starttim
        endtim1 = endtim
        timerange = '{0}~{1}'.format(starttim.iso.replace('-', '/').replace(' ', '/'), endtim.iso.replace('-', '/').replace(' ', '/'))
    else:
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
    midtime_mjd = (starttim1.mjd + endtim1.mjd) / 2.

    if vis.endswith('/'):
        vis = vis[:-1]
    visname = os.path.basename(vis)
    bt = starttim1.plot_date
    et = endtim1.plot_date

    # find out min and max frequency for plotting in dynamic spectrum
    ms.open(vis)
    metadata = ms.metadata()
    observatory = metadata.observatorynames()[0]
    spwInfo = ms.getspectralwindowinfo()
    nspw = len(spwInfo)
    if not spw:
        spw = '0~' + str(nspw - 1)
    staql = {'timerange': timerange, 'spw': spw}
    if ms.msselect(staql, onlyparse=True):
        ndx = ms.msselectedindices()
        chan_sel = ndx['channel']
        nspw = chan_sel.shape[0]
        bspw = chan_sel[0, 0]
        bchan = chan_sel[0, 1]
        espw = chan_sel[-1, 0]
        echan = chan_sel[-1, 2]
        bfreq = spwInfo[str(bspw)]['Chan1Freq'] + spwInfo[str(bspw)]['ChanWidth'] * bchan
        efreq = spwInfo[str(espw)]['Chan1Freq'] + spwInfo[str(espw)]['ChanWidth'] * echan
        bfreqghz = bfreq / 1e9
        efreqghz = efreq / 1e9
        if verbose:
            print 'selected timerange {}'.format(timerange)
            print 'selected frequency range {0:6.3f} to {1:6.3f} GHz'.format(bfreqghz, efreqghz)
    else:
        print "spw or timerange selection failed. Aborting..."
        ms.close()
        return -1
    ms.close()

    if observatory == 'EOVSA':
        print 'Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be IV'
        stokes = 'I,V'

    if mkmovie:
        plt.ioff()
        # fig = plt.figure(figsize=(12, 7.5), dpi=100)
        if fitsfile:
            pass
        else:
            if not imagefile:
                # from ptclean_cli import ptclean_cli as ptclean
                eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
                if observatory == 'EOVSA' or (not usephacenter):
                    phasecenter = ''
                else:
                    phasecenter = 'J2000 ' + str(eph['ra'][0])[:15] + 'rad ' + str(eph['dec'][0])[:15] + 'rad'
                print 'use phasecenter: ' + phasecenter
                qlookfitsdir = os.path.join(workdir, 'qlookfits/')
                qlookfigdir = os.path.join(workdir, 'qlookimgs/')
                imresfile = os.path.join(qlookfitsdir, '{}.imres.npz'.format(os.path.basename(vis)))
                if overwrite:
                    imres = mk_qlook_image(vis, twidth=twidth, ncpu=ncpu, imagedir=qlookfitsdir, phasecenter=phasecenter, stokes=stokes,
                                           c_external=True)
                else:
                    if os.path.exists(imresfile):
                        imres = np.load(imresfile)
                        imres = imres['imres'].item()
                    else:
                        print('Image results file not found; Creating new images.')
                        imres = mk_qlook_image(vis, twidth=twidth, ncpu=ncpu, imagedir=qlookfigdir, phasecenter=phasecenter, stokes=stokes,
                                               c_external=True)
                if not os.path.exists(qlookfigdir):
                    os.makedirs(qlookfigdir)
                plt_qlook_image(imres, figdir=qlookfigdir, specdata=specdata, verbose=True, stokes=stokes, fov=fov)

    else:
        spec = specdata['spec']
        (npol, nbl, nfreq, ntim) = spec.shape
        tidx = range(ntim)
        fidx = range(nfreq)
        tim = specdata['tim']
        freq = specdata['freq']
        freqghz = freq / 1e9
        spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')
        timstrr = spec_tim.plot_date
        plt.ion()
        fig = plt.figure(figsize=(12, 7), dpi=100)
        gs1 = gridspec.GridSpec(3, 1)
        gs1.update(left=0.08, right=0.32, wspace=0.05)
        gs2 = gridspec.GridSpec(2, 2)
        gs2.update(left=0.38, right=0.98, hspace=0.02, wspace=0.02)

        spec_1 = np.absolute(spec[0, 0, :, :])
        spec_2 = np.absolute(spec[1, 0, :, :])
        if observatory == 'EVLA':
            # circular feeds
            polstr = ['RR', 'LL']
        if observatory == 'EOVSA' or observatory == 'ALMA':
            # linear feeds
            polstr = ['XX', 'YY']

        print 'plot the dynamic spectrum in pol ' + ' & '.join(polstr)
        ax1 = plt.subplot(gs1[0])
        ax1.pcolormesh(timstrr, freqghz, spec_1, cmap='jet', vmin=dmin, vmax=dmax)
        ax1.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
        ax1.xaxis_date()
        ax1.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
        # ax1.set_xticklabels(['']*10)
        ax1.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
        ax1.set_ylabel('Frequency (GHz)', fontsize=10)
        ax1.set_title(observatory + ' ' + datstr + ' ' + polstr[0] + ' & ' + polstr[1], fontsize=12)
        ax1.set_autoscale_on(False)
        ax1.add_patch(patches.Rectangle((bt, bfreqghz), et - bt, efreqghz - bfreqghz, ec='w', fill=False))
        ax1.plot([(bt + et) / 2.], [(bfreqghz + efreqghz) / 2.], '*w', ms=12)
        for tick in ax1.get_xticklabels():
            tick.set_fontsize(8)
        for tick in ax1.get_yticklabels():
            tick.set_fontsize(8)
        ax2 = plt.subplot(gs1[1])
        ax2.pcolormesh(timstrr, freqghz, spec_2, cmap='jet', vmin=dmin, vmax=dmax)
        ax2.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
        ax2.xaxis_date()
        ax2.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
        ax2.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
        ax2.set_ylabel('Frequency (GHz)', fontsize=10)
        for tick in ax2.get_xticklabels():
            tick.set_fontsize(8)
        for tick in ax2.get_yticklabels():
            tick.set_fontsize(8)
        ax2.set_autoscale_on(False)
        ax2.add_patch(patches.Rectangle((bt, bfreqghz), et - bt, efreqghz - bfreqghz, ec='w', fill=False))
        ax2.plot([(bt + et) / 2.], [(bfreqghz + efreqghz) / 2.], '*w', ms=12)

        # Second part: GOES plot
        if goestime:
            btgoes = goestime[0]
            etgoes = goestime[1]
        else:
            datstrg = datstr.replace('-', '/')
            btgoes = datstrg + ' ' + qa.time(qa.quantity(tim[0] - 1800, 's'), form='clean', prec=9)[0]
            etgoes = datstrg + ' ' + qa.time(qa.quantity(tim[tidx[-1] - 1] + 1800, 's'), form='clean', prec=9)[0]
        if verbose:
            print 'Acquire GOES soft X-ray data in from ' + btgoes + ' to ' + etgoes

        ax3 = plt.subplot(gs1[2])

        try:
            from sunpy import lightcurve as lc
            from sunpy.time import TimeRange
            goest = lc.GOESLightCurve.create(TimeRange(btgoes, etgoes))
        except:
            goesscript = os.path.join(workdir, 'goes.py')
            goesdatafile = os.path.join(workdir, 'goes.dat')
            os.system('rm -rf {}'.format(goesscript))
            fi = open(goesscript, 'wb')
            fi.write('import os \n')
            fi.write('from sunpy.time import TimeRange \n')
            fi.write('from sunpy import lightcurve as lc \n')
            fi.write('import pickle \n')
            fi.write('goesplottim = TimeRange("{0}", "{1}") \n'.format(btgoes, etgoes))
            fi.write('goes = lc.GOESLightCurve.create(goesplottim) \n')
            fi.write('fi2 = open("{}", "wb") \n'.format(goesdatafile))
            fi.write('pickle.dump(goes, fi2) \n')
            fi.write('fi2.close()')
            fi.close()

            try:
                os.system('python {}'.format(goesscript))
            except NameError:
                print "Bad input names"
            except ValueError:
                print "Bad input values"
            except:
                print "Unexpected error:", sys.exc_info()[0]
                print "Error in generating GOES light curves. Proceed without GOES..."

            if os.path.exists(goesdatafile):
                fi1 = file(goesdatafile, 'rb')
                goest = pickle.load(fi1)
                fi1.close()

        try:
            dates = mpl.dates.date2num(parse_time(goest.data.index))
            goesdif = np.diff(goest.data['xrsb'])
            gmax = np.nanmax(goesdif)
            gmin = np.nanmin(goesdif)
            ran = gmax - gmin
            db = 2.8 / ran
            goesdifp = goesdif * db + gmin + (-6)
            ax3.plot_date(dates, np.log10(goest.data['xrsb']), '-', label='1.0--8.0 $\AA$', color='red', lw=2)
            ax3.plot_date(dates[0:-1], goesdifp, '-', label='derivate', color='blue', lw=0.4)

            ax3.set_ylim([-7, -3])
            ax3.set_yticks([-7, -6, -5, -4, -3])
            ax3.set_yticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$'])
            ax3.set_title('Goes Soft X-ray', fontsize=12)
            ax3.set_ylabel('Watts m$^{-2}$')
            ax3.set_xlabel(datetime.datetime.isoformat(goest.data.index[0])[0:10])
            ax3.axvspan(dates[899], dates[dates.size - 899], alpha=0.2)

            ax2 = ax3.twinx()
            # ax2.set_yscale("log")
            ax2.set_ylim([-7, -3])
            ax2.set_yticks([-7, -6, -5, -4, -3])
            ax2.set_yticklabels(['B', 'C', 'M', 'X', ''])

            ax3.yaxis.grid(True, 'major')
            ax3.xaxis.grid(False, 'major')
            ax3.legend(prop={'size': 6})

            formatter = mpl.dates.DateFormatter('%H:%M')
            ax3.xaxis.set_major_formatter(formatter)

            ax3.fmt_xdata = mpl.dates.DateFormatter('%H:%M')
        except:
            print 'Error in downloading GOES soft X-ray data. Proceeding with out soft X-ray plot.'

        # third part
        # start to download the fits files
        if plotaia:
            if not aiafits:
                newlist = []
                items = glob.glob('*.fits')
                for names in items:
                    str1 = starttim1.iso[:4] + '_' + starttim1.iso[5:7] + '_' + starttim1.iso[8:10] + 't' + starttim1.iso[
                                                                                                            11:13] + '_' + starttim1.iso[14:16]
                    str2 = str(aiawave)
                    if names.endswith(".fits"):
                        if names.find(str1) != -1 and names.find(str2) != -1:
                            newlist.append(names)
                    newlist.append('0')
                if newlist and os.path.exists(newlist[0]):
                    aiafits = newlist[0]
                else:
                    print 'downloading the aiafits file'
                    client = vso.VSOClient()
                    wave1 = aiawave - 3
                    wave2 = aiawave + 3
                    t1 = Time(starttim1.mjd - 0.02 / 24., format='mjd')
                    t2 = Time(endtim1.mjd + 0.02 / 24., format='mjd')
                    qr = client.query(vso.attrs.Time(t1.iso, t2.iso), vso.attrs.Instrument('aia'), vso.attrs.Wave(wave1 * u.AA, wave2 * u.AA))
                    res = client.get(qr, path='{file}')

            # Here something is needed to check whether it has finished downloading the fits files or not

            if not aiafits:
                newlist = []
                items = glob.glob('*.fits')
                for nm in items:
                    str1 = starttim1.iso[:4] + '_' + starttim1.iso[5:7] + '_' + starttim1.iso[8:10] + 't' + starttim1.iso[
                                                                                                            11:13] + '_' + starttim1.iso[14:16]
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

        # RCP or I
        ax4 = plt.subplot(gs2[0, 0])
        ax5 = plt.subplot(gs2[1, 0])
        # LCP or V
        ax6 = plt.subplot(gs2[0, 1])
        ax7 = plt.subplot(gs2[1, 1])

        if fitsfile:
            pass
        else:
            if not imagefile:
                eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
                if observatory == 'EOVSA' or (not usephacenter):
                    print 'This is EOVSA data'
                    phasecenter = ''
                    if stokes == 'RRLL' or stokes == 'RR,LL':
                        print 'Provide stokes: ' + str(stokes) + '. However EOVSA has linear feeds. Force stokes to be IV'
                        stokes = 'I,V'
                else:
                    phasecenter = 'J2000 ' + str(eph['ra'][0])[:15] + 'rad ' + str(eph['dec'][0])[:15] + 'rad'
                imagename = os.path.join(workdir, visname + '.outim')
                if os.path.exists(imagename + '.image') or os.path.exists(imagename + '.flux'):
                    os.system('rm -rf ' + imagename + '.*')
                sto = stokes.replace(',', '')
                print 'do clean for ' + timerange + ' in spw ' + spw + ' stokes ' + sto
                print 'use phasecenter: ' + phasecenter
                clean(vis=vis, imagename=imagename, selectdata=True, spw=spw, timerange=timerange, stokes=sto, niter=500, interactive=False,
                      npercycle=50, imsize=[512, 512], cell=['5.0arcsec'], phasecenter=phasecenter)
                os.system('rm -rf ' + imagename + '.psf')
                os.system('rm -rf ' + imagename + '.flux')
                os.system('rm -rf ' + imagename + '.model')
                os.system('rm -rf ' + imagename + '.mask')
                os.system('rm -rf ' + imagename + '.residual')
                imagefile = imagename + '.image'
            fitsfile = imagefile + '.fits'
            hf.imreg(vis=vis, ephem=eph, imagefile=imagefile, timerange=timerange, reftime=reftime, fitsfile=fitsfile, verbose=True, overwrite=True)
        print 'fits file ' + fitsfile + ' selected'
        ax4.cla()
        ax5.cla()
        ax6.cla()
        ax7.cla()

        rfits = fitsfile
        try:
            hdulist = fits.open(rfits)
            hdu = hdulist[0]
            (npol, nf, nx, ny) = hdu.data.shape
            rmap = smap.Map(hdu.data[0, 0, :, :], hdu.header)
        except:
            print 'radio fits file not recognized by sunpy.map. Aborting...'
            return -1
        if npol > 1:
            rmap1 = smap.Map(hdu.data[0, 0, :, :], hdu.header)
            rmap2 = smap.Map(hdu.data[1, 0, :, :], hdu.header)

        XX, YY = np.meshgrid(np.arange(rmap.data.shape[1]), np.arange(rmap.data.shape[0]))
        try:
            rmapx, rmapy = rmap.pixel_to_data(XX * u.pix, YY * u.pix)
        except:
            rmapxy = rmap.pixel_to_data(XX * u.pix, YY * u.pix)
            rmapx = rmapxy.Tx
            rmapy = rmapxy.Ty

        if fov:
            xc, yc = fov[0]
            xlen, ylen = fov[1]
            fov = [[xc - xlen / 2.0, yc - ylen / 2.0], [xc + xlen / 2.0, yc + ylen / 2.0]]
        else:
            row, col = rmap1.data.shape
            positon = np.nanargmax(rmap1.data)
            m, n = divmod(positon, col)
            length = 200 * u.arcsec
            x0 = rmap1.xrange[0] + rmap1.scale[1] * (n + 0.5) * u.pix
            y0 = rmap1.yrange[0] + rmap1.scale[0] * (m + 0.5) * u.pix
            x1 = x0 - length
            x2 = x0 + length
            y1 = y0 - length
            y2 = y0 + length
            fov = [[x1.value, x2.value], [y1.value, y2.value]]

        clevels1 = np.linspace(0.2, 0.9, 5)
        if stokes.split(',')[1] == 'V':
            clevels2 = np.array([-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8])
        else:
            clevels2 = np.linspace(0.2, 0.9, 5)
        if 'aiamap' in vars():
            aiamap.plot_settings['cmap'] = plt.get_cmap('binary')
            if rmap:
                title = 'AIA {0:.0f} + {1} {2:6.3f} GHz'.format(aiamap.wavelength.value, observatory, (bfreqghz + efreqghz) / 2.0)
            else:
                title = 'AIA {0:.0f}'.format(aiamap.wavelength.value)
            aiamap.plot(axes=ax4)
            ax4.set_title(title + ' ' + stokes.split(',')[0], fontsize=12)
            aiamap.draw_limb()
            aiamap.draw_grid()
            aiamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
            aiamap.plot(axes=ax6)
            ax6.set_title(title + ' ' + stokes.split(',')[1], fontsize=12)
            aiamap.draw_limb()
            aiamap.draw_grid()
            aiamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
            if rmap:
                ax4.contour(rmapx.value, rmapy.value, rmap1.data, levels=clevels1 * np.nanmax(rmap1.data), cmap=cm.jet)
                ax6.contour(rmapx.value, rmapy.value, rmap2.data, levels=clevels2 * np.nanmax(rmap2.data), cmap=cm.jet)
            ax4.text(0.02, 0.02, 'AIA {0:.0f} '.format(aiamap.wavelength.value) + aiamap.date.strftime('%H:%M:%S'), verticalalignment='bottom',
                     horizontalalignment='left', transform=ax4.transAxes, color='k', fontsize=10)
            ax6.text(0.02, 0.02, 'AIA {0:.0f} '.format(aiamap.wavelength.value) + aiamap.date.strftime('%H:%M:%S'), verticalalignment='bottom',
                     horizontalalignment='left', transform=ax6.transAxes, color='k', fontsize=10)
        else:
            title = '{0} {1:6.3f} GHz'.format(observatory, (bfreqghz + efreqghz) / 2.0)
            rmap1.plot(axes=ax4, cmap=cm.jet)
            ax4.set_title(title + ' ' + stokes.split(',')[0], fontsize=12)
            rmap1.draw_limb()
            rmap1.draw_grid()
            rmap1.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
            rmap2.plot(axes=ax6, cmap=cm.jet)
            ax6.set_title(title + ' ' + stokes.split(',')[1], fontsize=12)
            rmap2.draw_limb()
            rmap2.draw_grid()
            # ax4.contour(rmapx.value, rmapy.value, rmap1.data, levels=np.linspace(0.2, 0.9, 5) * np.nanmax(rmap1.data),
            #            cmap=cm.gray)
            # ax6.contour(rmapx.value, rmapy.value, rmap2.data, levels=np.linspace(0.2, 0.9, 5) * np.nanmax(rmap2.data),
            #            cmap=cm.gray)
            rmap2.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec,
                                 400 * u.arcsec)  # ax4.contour(rmapx.value, rmapy.value, rmap1.data,  #            levels=clevels1 * np.nanmax(rmap1.data), cmap=cm.gray)  # ax6.contour(rmapx.value, rmapy.value, rmap2.data,  #            levels=clevels2 * np.nanmax(rmap2.data), cmap=cm.gray)
        ax4.set_xlim(-1200, 1200)
        ax4.set_ylim(-1200, 1200)
        ax6.set_xlim(-1200, 1200)
        ax6.set_ylim(-1200, 1200)

        try:
            subrmap1 = rmap1.submap(fov[0] * u.arcsec, fov[1] * u.arcsec)
            subrmap2 = rmap2.submap(fov[0] * u.arcsec, fov[1] * u.arcsec)
        except:
            bl = SkyCoord(fov[0][0] * u.arcsec, fov[1][0] * u.arcsec, frame=rmap1.coordinate_frame)
            tr = SkyCoord(fov[0][1] * u.arcsec, fov[1][1] * u.arcsec, frame=rmap1.coordinate_frame)
            subrmap1 = rmap1.submap(bl, tr)
            subrmap2 = rmap2.submap(bl, tr)

        XX, YY = np.meshgrid(np.arange(subrmap1.data.shape[1]), np.arange(subrmap1.data.shape[0]))
        try:
            subrmapx, subrmapy = subrmap1.pixel_to_data(XX * u.pix, YY * u.pix)
        except:
            subrmapxy = subrmap1.pixel_to_data(XX * u.pix, YY * u.pix)
            subrmapx = subrmapxy.Tx
            subrmapy = subrmapxy.Ty

        if 'aiamap' in vars():
            try:
                subaiamap = aiamap.submap(fov[0] * u.arcsec, fov[1] * u.arcsec)
            except:
                bl = SkyCoord(fov[0][0] * u.arcsec, fov[1][0] * u.arcsec, frame=aiamap.coordinate_frame)
                tr = SkyCoord(fov[0][1] * u.arcsec, fov[1][1] * u.arcsec, frame=aiamap.coordinate_frame)
                subaiamap = aiamap.submap(bl, tr)

            subaiamap.plot(axes=ax5, title='')
            subaiamap.draw_limb()
            subaiamap.draw_grid()
            subaiamap.plot(axes=ax7, title='')
            subaiamap.draw_limb()
            subaiamap.draw_grid()
            ax5.contour(subrmapx.value, subrmapy.value, subrmap1.data, levels=clevels1 * np.nanmax(subrmap1.data), cmap=cm.jet)
            ax7.contour(subrmapx.value, subrmapy.value, subrmap2.data, levels=clevels2 * np.nanmax(subrmap2.data),
                        cmap=cm.jet)  # subaiamap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
        else:
            subrmap1.plot(axes=ax5, cmap=cm.jet, title='')
            subrmap1.draw_limb()
            subrmap1.draw_grid()
            subrmap2.plot(axes=ax7, cmap=cm.jet, title='')
            subrmap2.draw_limb()
            subrmap2.draw_grid()  # ax5.contour(subrmapx.value, subrmapy.value, subrmap1.data,  #            levels=clevels1 * np.nanmax(subrmap1.data), cmap=cm.gray)  # ax7.contour(subrmapx.value, subrmapy.value, subrmap2.data,  #            levels=clevels2 * np.nanmax(subrmap2.data), cmap=cm.gray)  # subrmap1.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)  # subrmap2.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, 400 * u.arcsec, 400 * u.arcsec)
        ax5.set_xlim(fov[0])
        ax5.set_ylim(fov[1])
        ax5.text(0.02, 0.02, observatory + ' ' + rmap.date.strftime('%H:%M:%S.%f')[:-3], verticalalignment='bottom', horizontalalignment='left',
                 transform=ax5.transAxes, color='k', fontsize=10)
        ax7.set_xlim(fov[0])
        ax7.set_ylim(fov[1])
        ax7.text(0.02, 0.02, observatory + ' ' + rmap.date.strftime('%H:%M:%S.%f')[:-3], verticalalignment='bottom', horizontalalignment='left',
                 transform=ax7.transAxes, color='k', fontsize=10)

        fig.show()

        os.system('rm -rf {}'.format(goesscript))  # os.system('rm -rf {}'.format(goesdatafile))
