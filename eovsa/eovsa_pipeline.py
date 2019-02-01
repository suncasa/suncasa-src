from suncasa.eovsa import eovsa_prep as ep
from ptclean_cli import ptclean_cli as ptclean
from eovsapy.util import Time
from importeovsa_cli import importeovsa_cli as importeovsa
from suncasa.tasks import task_importeovsa as timporteovsa
from split_cli import split_cli as split
from clean_cli import clean_cli as clean
from delmod_cli import delmod_cli as delmod
from clearcal_cli import clearcal_cli as clearcal
from gaincal_cli import gaincal_cli as gaincal
from applycal_cli import applycal_cli as applycal
from suncasa.tasks import task_calibeovsa as calibeovsa
from eovsapy import refcal_anal as ra
from taskinit import ms
import astropy.units as u
import os
import numpy as np
import pdb

udbmsscldir = os.getenv('EOVSAUDBMSSCL')
udbmsdir = os.getenv('EOVSAUDBMS')
udbdir = os.getenv('EOVSAUDB')

if not udbmsdir:
    print 'Environmental variable for EOVSA udbms path not defined'
    print 'Use default path on pipeline'
    udbmsdir = '/data1/eovsa/fits/UDBms/'
if not udbmsscldir:
    print 'Environmental variable for scaled EOVSA udbms path not defined'
    print 'Use default path on pipeline'
    udbmsscldir = '/data1/eovsa/fits/UDBms_scl/'
if not udbdir:
    print 'Environmental variable for EOVSA udb path not defined'
    print 'Use default path on pipeline'
    udbdir = '/data1/eovsa/fits/UDB/'


def trange2ms(trange=None, doimport=False, verbose=False, doscaling=False):
    '''This finds all solar UDBms files within a timerange; If the UDBms file does not exist 
       in EOVSAUDBMSSCL, create one by calling importeovsa
       Required inputs:
       trange - can be 1) a single string or Time() object in UTC: use the entire day, e.g., '2017-08-01' or Time('2017-08-01')
                          if just a date, find all scans withing the same date in local time. 
                          if a complete time stamp, find the local date first (which may be different from that provided, 
                            and return all scans within that day
                       2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                       3) None -- use current date Time.now()
       doimport - Boolean. If true, call importeovsa to import UDB files that are missing from 
                  those found in the directory specified in EOVSAUDBMSSCL. Otherwise, return
                  a list of ms files it has found.
       doscaling - Boolean. If true, scale cross-correlation amplitudes by using auto-correlations
       verbose - Boolean. If true, return more information
    '''
    import glob
    import pytz
    from datetime import datetime
    if trange is None:
        trange = Time.now()
    if type(trange) == list or type(trange) == str:
        try:
            trange = Time(trange)
        except:
            print('trange format not recognised. Abort....')
            return None
    # if type(trange) == Time:
    try:
        # if single Time object, the following line would report an error
        nt = len(trange)
        if len(trange) > 1:
            # more than one value
            trange = Time([trange[0], trange[-1]])
            tdatetime = trange[0].to_datetime()
        else:
            # single value in a list
            if trange[0].mjd == np.fix(trange[0].mjd):
                # if only date is given, move the time from 00 to 12 UT
                trange[0] = Time(trange[0].mjd + 0.5, format='mjd')

            tdatetime = trange[0].to_datetime()
            dhr = trange[0].LocalTime.utcoffset().total_seconds() / 60 / 60 / 24
            btime = Time(np.fix(trange[0].mjd + dhr) - dhr, format='mjd')
            etime = Time(btime.mjd + 1, format='mjd')
            trange = Time([btime, etime])
    except:
        # the case of a single Time object
        if trange.mjd == np.fix(trange.mjd):
            # if only date is given, move the time from 00 to 12 UT
            trange = Time(trange.mjd + 0.5, format='mjd')

        tdatetime = trange.to_datetime()
        dhr = trange.LocalTime.utcoffset().total_seconds() / 60 / 60 / 24
        btime = Time(np.fix(trange.mjd + dhr) - dhr, format='mjd')
        etime = Time(btime.mjd + 1, format='mjd')
        trange = Time([btime, etime])

    print 'Selected timerange in UTC: ', trange.iso

    sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
    udbfilelist = sclist['scanlist']
    udbfilelist = [os.path.basename(ll) for ll in udbfilelist]
    if doscaling:
        udbmspath = udbmsscldir
    else:
        udbmspath = udbmsdir
    inpath = '{}{}/'.format(udbdir, tdatetime.strftime("%Y"))
    outpath = '{}{}/'.format(udbmspath, tdatetime.strftime("%Y%m"))
    if not os.path.exists(outpath):
        if verbose:
            print outpath + ' does not exist. Making a new directory.'
        os.makedirs(outpath)
        msfiles = []
    else:
        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]

    msfile_synoptic = os.path.join(outpath, 'UDB' + tdatetime.strftime("%Y%m%d") + '.ms')
    if os.path.exists(msfile_synoptic):
        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': [], 'ms': [msfile_synoptic],
                'tstlist': sclist['tstlist'], 'tedlist': sclist['tedlist']}
    else:
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))
        if filelist and doimport:
            import multiprocessing as mprocs
            #ncpu = mprocs.cpu_count()
            #if ncpu > 10:
            #    ncpu = 10
            #if ncpu > len(filelist):
            #    ncpu = len(filelist)
            ncpu = 1
            timporteovsa.importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1, visprefix=outpath, nocreatms=False,
                                     doconcat=False, modelms="", doscaling=doscaling, keep_nsclms=False, udb_corr=True)

        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))

        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': filelist,
                'ms': [outpath + ll + '.ms' for ll in sorted(list(msfiles))], 'tstlist': sclist['tstlist'], 'tedlist': sclist['tedlist']}


def calib_pipeline(trange, doimport=False, synoptic=False):
    ''' 
       trange: can be 1) a single Time() object: use the entire day
                      2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                      3) a single or a list of UDBms file(s)
                      4) None -- use current date Time.now()
    '''
    if type(trange) == Time:
        mslist = trange2ms(trange=trange, doimport=doimport)
        invis = mslist['ms']
        tsts = [l.to_datetime() for l in mslist['tstlist']]
        subdir = [tst.strftime("%Y/%m/%d/") for tst in tsts]
    if type(trange) == str:
        try:
            date = Time(trange)
            mslist = trange2ms(trange=trange, doimport=doimport)
            invis = mslist['ms']
        except:
            invis = [trange]
        subdir = ['/']

    for idx, f in enumerate(invis):
        if f[-1] == '/':
            invis[idx] = f[:-1]

    if synoptic:
        vis = calibeovsa.calibeovsa(invis, caltype=['refpha', 'phacal'], interp='nearest', doflag=True, flagant='13~15', doimage=False, doconcat=True,
                                    msoutdir=os.path.dirname(invis[0]), concatvis=os.path.basename(invis[0])[:11] + '.ms', keep_orig_ms=False)
    else:
        vis = calibeovsa.calibeovsa(invis, caltype=['refpha', 'phacal'], interp='nearest', doflag=True, flagant='13~15', doimage=False, doconcat=True,
                                    msoutdir=os.path.dirname(invis[0]), keep_orig_ms=False)
    return vis


def mk_qlook_image(trange, doimport=False, docalib=False, ncpu=10, twidth=12, stokes=None, antenna='0~12',
                   #imagedir=None, spws=['1~3','4~6','7~9','10~13','14~18','19~28'],verbose=False):
                   imagedir=None, spws=['1~5', '6~10', '11~15', '16~25'], toTb=True, overwrite=True, doslfcal=False, verbose=False):
    '''
       trange: can be 1) a single Time() object: use the entire day
                      2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                      3) a single or a list of UDBms file(s)
                      4) None -- use current date Time.now()
    '''
    antenna0 = antenna
    if type(trange) == Time:
        mslist = trange2ms(trange=trange, doimport=doimport)
        vis = mslist['ms']
        tsts = [l.to_datetime() for l in mslist['tstlist']]
        subdir = [tst.strftime("%Y/%m/%d/") for tst in tsts]
    if type(trange) == str:
        try:
            date = Time(trange)
            mslist = trange2ms(trange=trange, doimport=doimport)
            vis = mslist['ms']
        except:
            vis = [trange]
        subdir = ['/']

    for idx, f in enumerate(vis):
        if f[-1] == '/':
            vis[idx] = f[:-1]
    if not stokes:
        stokes = 'XX'

    if not imagedir:
        imagedir = './'
    imres = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': [],
             'Synoptic': {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': []}}
    for n, msfile in enumerate(vis):
        msfilebs = os.path.basename(msfile)
        imdir = imagedir + subdir[n]
        if not os.path.exists(imdir):
            os.makedirs(imdir)
        if doslfcal:
            slfcalms = './' + msfilebs + '.xx'
            split(msfile, outputvis=slfcalms, datacolumn='corrected', correlation='XX')
        for spw in spws:
            antenna = antenna0
            spwran = [s.zfill(2) for s in spw.split('~')]
            freqran = [int(s) * 0.5 + 2.9 for s in spw.split('~')]
            cfreq = np.mean(freqran)
            bmsz = max(150. / cfreq, 20.)
            uvrange = '<10klambda'
            if doslfcal:
                slfcal_img = './' + msfilebs + '.slf.spw' + spw.replace('~', '-') + '.slfimg'
                slfcal_tb = './' + msfilebs + '.slf.spw' + spw.replace('~', '-') + '.slftb'
                try:
                    clean(vis=slfcalms, antenna=antenna, imagename=slfcal_img, spw=spw, mode='mfs', timerange='', imagermode='csclean',
                          psfmode='clark', imsize=[512, 512], cell=['5arcsec'], niter=100, gain=0.05, stokes='I', weighting='natural',
                          restoringbeam=[str(bmsz) + 'arcsec'], pbcor=False, interactive=False, usescratch=True)
                except:
                    print 'error in cleaning spw: ' + spw
                    break
                gaincal(vis=slfcalms, refant='0', antenna=antenna, caltable=slfcal_tb, spw=spw, uvrange='', gaintable=[], selectdata=True,
                        timerange='', solint='600s', gaintype='G', calmode='p', combine='', minblperant=3, minsnr=2, append=False)
                if not os.path.exists(slfcal_tb):
                    print 'No solution found in spw: ' + spw
                    break
                else:
                    clearcal(slfcalms)
                    delmod(slfcalms)
                    applycal(vis=slfcalms, gaintable=[slfcal_tb], spw=spw, selectdata=True, antenna=antenna, interp='nearest', flagbackup=False,
                             applymode='calonly', calwt=False)
                    msfile = slfcalms

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
            if cfreq > 10.:
                antenna = antenna + ';!0&1;!0&2'  #deselect the shortest baselines

            res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, uvrange=uvrange, spw=spw, ncpu=ncpu, niter=1000,
                          gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=stokes, doreg=True, usephacenter=False, overwrite=overwrite,
                          toTb=toTb, restoringbeam=restoringbeam, uvtaper=True, outertaper=['30arcsec'])

            if res:
                imres['Succeeded'] += res['Succeeded']
                imres['BeginTime'] += res['BeginTime']
                imres['EndTime'] += res['EndTime']
                imres['ImageName'] += res['ImageName']
                imres['Spw'] += [spwstr] * len(res['ImageName'])
                imres['Vis'] += [msfile] * len(res['ImageName'])
            else:
                continue

    if len(vis) == 1:
        # produce the band-by-band whole-day images
        ms.open(msfile)
        ms.selectinit()
        timfreq = ms.getdata(['time', 'axis_info'], ifraxis=True)
        tim = timfreq['time']
        ms.close()

        imdir = imagedir + subdir[0]
        if not os.path.exists(imdir):
            os.makedirs(imdir)
        for spw in spws:
            spwran = [s.zfill(2) for s in spw.split('~')]
            freqran = [int(s) * 0.5 + 2.9 for s in spw.split('~')]
            cfreq = np.mean(freqran)
            bmsz = max(150. / cfreq, 20.)
            uvrange = '<10klambda'
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
            imagesuffix = '.synoptic.spw' + spwstr.replace('~', '-')
            # if cfreq > 10.:
            #     antenna = antenna + ';!0&1;!0&2'  #deselect the shortest baselines

            res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=len(tim), uvrange=uvrange, spw=spw, ncpu=1, niter=1000,
                          gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=stokes, doreg=True, usephacenter=False, overwrite=overwrite,
                          toTb=toTb, restoringbeam=restoringbeam, uvtaper=True, outertaper=['30arcsec'])
            if res:
                imres['Synoptic']['Succeeded'] += res['Succeeded']
                imres['Synoptic']['BeginTime'] += res['BeginTime']
                imres['Synoptic']['EndTime'] += res['EndTime']
                imres['Synoptic']['ImageName'] += res['ImageName']
                imres['Synoptic']['Spw'] += [spwstr] * len(res['ImageName'])
                imres['Synoptic']['Vis'] += [msfile] * len(res['ImageName'])
            else:
                continue

    #save it for debugging purposes
    np.savez('imres.npz', imres=imres)

    return imres


def plt_qlook_image(imres, figdir=None, verbose=True, synoptic=False):
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    from sunpy import sun
    from matplotlib import colors
    import astropy.units as u
    if not figdir:
        figdir = './'
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
    fig = plt.figure(figsize=(8, 8))
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    for i in range(ntime):
        plt.ioff()
        plt.clf()
        plttime = btimes_sort[i, 0]
        tofd = plttime.mjd - np.fix(plttime.mjd)
        suci = suc_sort[i]
        if not synoptic:
            if tofd < 16. / 24. or sum(
                    suci) < nspw - 2:  # if time of the day is before 16 UT (and 24 UT), skip plotting (because the old antennas are not tracking)
                continue
        #fig=plt.figure(figsize=(9,6))
        #fig.suptitle('EOVSA @ '+plttime.iso[:19])
        if synoptic:
            fig.text(0.01, 0.98, plttime.iso[:10], color='w', fontweight='bold', fontsize=12, ha='left')
        else:
            fig.text(0.01, 0.98, plttime.iso[:19], color='w', fontweight='bold', fontsize=12, ha='left')
        if verbose:
            print 'Plotting image at: ', plttime.iso
        for n in range(nspw):
            plt.ioff()
            image = images_sort[i, n]
            #fig.add_subplot(nspw/3, 3, n+1)
            fig.add_subplot(nspw / 2, 2, n + 1)
            if suci[n]:
                try:
                    eomap = smap.Map(image)
                except:
                    continue
                sz = eomap.data.shape
                if len(sz) == 4:
                    eomap.data = eomap.data.reshape((sz[2], sz[3]))
                eomap.data[np.isnan(eomap.data)] = 0.0
                #resample the image for plotting
                dim = u.Quantity([256, 256], u.pixel)
                eomap = eomap.resample(dim)
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot_settings['norm'] = colors.Normalize(vmin=-1e5, vmax=1e6)
                eomap.plot()
                if not synoptic:
                    eomap.draw_limb()
                eomap.draw_grid()
                ax = plt.gca()
                ax.set_xlim([-1080, 1080])
                ax.set_ylim([-1080, 1080])
                spwran = spws_sort[i, n]
                freqran = [int(s) * 0.5 + 2.9 for s in spwran.split('~')]
                ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(freqran[0], freqran[1]), color='w', transform=ax.transAxes, fontweight='bold',
                        ha='right')
                ax.set_title(' ')
                #ax.set_title('spw '+spws_sort[i,n])
                #ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticklabels([''])
                ax.set_yticklabels([''])
            else:
                #make an empty map
                data = np.zeros((512, 512))
                header = {"DATE-OBS": plttime.isot, "EXPTIME": 0., "CDELT1": 5., "NAXIS1": 512, "CRVAL1": 0., "CRPIX1": 257, "CUNIT1": "arcsec",
                          "CTYPE1": "HPLN-TAN", "CDELT2": 5., "NAXIS2": 512, "CRVAL2": 0., "CRPIX2": 257, "CUNIT2": "arcsec", "CTYPE2": "HPLT-TAN",
                          "HGLT_OBS": sun.heliographic_solar_center(plttime)[1].value, "HGLN_OBS": 0.,
                          "RSUN_OBS": sun.solar_semidiameter_angular_size(plttime).value, "RSUN_REF": sun.constants.radius.value,
                          "DSUN_OBS": sun.sunearth_distance(plttime).to(u.meter).value, }
                eomap = smap.Map(data, header)
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot_settings['norm'] = colors.Normalize(vmin=-1e5, vmax=1e6)
                eomap.plot()
                if not synoptic:
                    eomap.draw_limb()
                eomap.draw_grid()
                ax = plt.gca()
                ax.set_xlim([-1080, 1080])
                ax.set_ylim([-1080, 1080])
                #ax.set_title('spw '+spwran+'( )'))
                spwran = spws_sort[i, n]
                freqran = [int(s) * 0.5 + 2.9 for s in spwran.split('~')]
                spwran = spws_sort[i, n]
                #ax.set_title('{0:.1f} - {1:.1f} GHz'.format(freqran[0],freqran[1]))
                ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(freqran[0], freqran[1]), color='w', transform=ax.transAxes, fontweight='bold',
                        ha='right')
                ax.set_title(' ')

                #ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticklabels([''])
                ax.set_yticklabels([''])
        fig_tdt = plttime.to_datetime()
        if synoptic:
            fig_subdir = fig_tdt.strftime("%Y/")
            figname = 'eovsa_qlimg_' + plttime.iso[:10].replace('-', '') + '.png'
        else:
            fig_subdir = fig_tdt.strftime("%Y/%m/%d/")
            figname = 'eovsa_qlimg_' + plttime.isot.replace(':', '').replace('-', '')[:15] + '.png'
        figdir_ = figdir + fig_subdir
        if not os.path.exists(figdir_):
            os.makedirs(figdir_)
        if verbose:
            print 'Saving plot to :' + figdir_ + figname
        plt.savefig(figdir_ + figname)
    plt.close(fig)


def qlook_image_pipeline(date, twidth=10, ncpu=15, doimport=False, docalib=False, synoptic=False):
    ''' date: date string or Time object. e.g., '2017-07-15' or Time('2017-07-15')
    '''
    import pytz
    from datetime import datetime
    if date is None:
        date = Time.now()
    try:
        date = Time(date)
    except:
        print('date format not recognised. Abort....')
        return None

    qlookfitsdir = os.getenv('EOVSAQLOOKFITS')
    qlookfigdir = os.getenv('EOVSAQLOOKFIG')
    synopticfigdir = os.getenv('EOVSASYNOPTICFIG')
    if not qlookfitsdir:
        qlookfitsdir = '/data1/eovsa/fits/qlook_10m/'
    if not qlookfigdir:
        qlookfigdir = '/common/webplots/qlookimg_10m/'
    if not synopticfigdir:
        synopticfigdir = '/common/webplots/SynopticImg/'

    imagedir = qlookfitsdir
    if synoptic:
        vis_synoptic = os.path.join(udbmsdir, date.datetime.strftime("%Y%m"), 'UDB' + date.datetime.strftime("%Y%m%d") + '.ms')
        if os.path.exists(vis_synoptic):
            date = vis_synoptic
        else:
            print('Whole-day ms file {} not existed. About..... Use pipeline1.py to make one.'.format(vis_synoptic))
            return None
    if docalib:
        vis = calib_pipeline(date, doimport=doimport, synoptic=synoptic)

    imres = mk_qlook_image(date, twidth=twidth, ncpu=ncpu, doimport=doimport, docalib=docalib, imagedir=imagedir, verbose=True)
    figdir = qlookfigdir
    plt_qlook_image(imres, figdir=figdir, verbose=True)
    if imres['Synoptic']['Succeeded']:
        figdir = synopticfigdir
        plt_qlook_image(imres['Synoptic'], figdir=figdir, verbose=True, synoptic=True)
