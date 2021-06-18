from ptclean3_cli import ptclean3_cli as ptclean3
from eovsapy.util import Time
from suncasa.tasks import task_importeovsa as timporteovsa
from split_cli import split_cli as split
from clean_cli import clean_cli as clean
from delmod_cli import delmod_cli as delmod
from clearcal_cli import clearcal_cli as clearcal
from gaincal_cli import gaincal_cli as gaincal
from applycal_cli import applycal_cli as applycal
from suncasa.tasks import task_calibeovsa as calibeovsa
from eovsapy import refcal_anal as ra
from taskinit import ms, tb
import os
import numpy as np
from suncasa.eovsa import eovsa_diskmodel as ed
from suncasa.utils import mstools as mstl

udbmsdir = os.getenv('EOVSAUDBMS')
udbmsscldir = os.getenv('EOVSAUDBMSSCL')
udbmsslfcaleddir = os.getenv('EOVSAUDBMSSLFCALED')
udbdir = os.getenv('EOVSAUDB')
caltbdir = os.getenv('EOVSACAL')
slfcaltbdir = os.getenv('EOVSASLFCAL')
qlookfitsdir = os.getenv('EOVSAQLOOKFITS')
qlookfigdir = os.getenv('EOVSAQLOOKFIG')
synopticfigdir = os.getenv('EOVSASYNOPTICFIG')

if not udbmsdir:
    print('Environmental variable for EOVSA udbms path not defined')
    print('Use default path on pipeline')
    udbmsdir = '/data1/eovsa/fits/UDBms/'
if not udbmsscldir:
    print('Environmental variable for scaled EOVSA udbms path not defined')
    print('Use default path on pipeline')
    udbmsscldir = '/data1/eovsa/fits/UDBms_scl/'
if not udbmsslfcaleddir:
    print('Environmental variable for EOVSA udbms path not defined')
    print('Use default path on pipeline')
    udbmsslfcaleddir = '/data1/eovsa/fits/UDBms_slfcaled/'
if not udbdir:
    print('Environmental variable for EOVSA udb path not defined')
    print('Use default path on pipeline')
    udbdir = '/data1/eovsa/fits/UDB/'
# check if the calibration table directory is defined

if not qlookfitsdir:
    qlookfitsdir = '/data1/eovsa/fits/synoptic/'
    if not os.path.exists(qlookfitsdir): os.makedirs(qlookfitsdir)
if not qlookfigdir:
    qlookfigdir = '/common/webplots/qlookimg_10m/'
    if not os.path.exists(qlookfigdir): os.makedirs(qlookfigdir)
if not synopticfigdir:
    synopticfigdir = '/common/webplots/SynopticImg/'
    if not os.path.exists(synopticfigdir): os.makedirs(synopticfigdir)

if not caltbdir:
    print('Task calibeovsa')
    caltbdir = '/data1/eovsa/caltable/'
    print('Environmental variable for EOVSA calibration table path not defined')
    print('Use default path on pipeline ' + caltbdir)

if not slfcaltbdir:
    print('Task calibeovsa')
    slfcaltbdir = '/data1/eovsa/slfcaltable/'
    print('Environmental variable for EOVSA disk calibration table path not defined')
    print('Use default path on pipeline ' + slfcaltbdir)


def getspwfreq(vis):
    '''

    :param vis:
    :return: mid frequencies in GHz of each spw in the vis
    '''
    tb.open(vis + '/SPECTRAL_WINDOW')
    reffreqs = tb.getcol('REF_FREQUENCY')
    bdwds = tb.getcol('TOTAL_BANDWIDTH')
    cfreqs = reffreqs + bdwds / 2.
    tb.close()
    cfreqs = cfreqs / 1.0e9
    return cfreqs


def trange2ms(trange=None, doimport=False, verbose=False, doscaling=False, overwrite=True):
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

    print('Selected timerange in UTC: ', trange.iso)

    if doscaling:
        udbmspath = udbmsscldir
    else:
        udbmspath = udbmsdir
    inpath = '{}{}/'.format(udbdir, tdatetime.strftime("%Y"))
    outpath = '{}{}/'.format(udbmspath, tdatetime.strftime("%Y%m"))
    if not os.path.exists(outpath):
        if verbose:
            print(outpath + ' does not exist. Making a new directory.')
        os.makedirs(outpath)
        msfiles = []
    else:
        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]

    msfile_synoptic = os.path.join(outpath, 'UDB' + tdatetime.strftime("%Y%m%d") + '.ms')
    if overwrite and doimport:
        if os.path.exists(msfile_synoptic):
            os.system('rm -rf {}'.format(msfile_synoptic))

    sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
    udbfilelist = sclist['scanlist']
    udbfilelist = [os.path.basename(ll) for ll in udbfilelist]

    if os.path.exists(msfile_synoptic):
        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': [],
                'ms': [msfile_synoptic],
                'tstlist': sclist['tstlist'], 'tedlist': sclist['tedlist']}
    else:
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))
        if filelist and doimport:
            import multiprocessing as mprocs
            # ncpu = mprocs.cpu_count()
            # if ncpu > 10:
            #    ncpu = 10
            # if ncpu > len(filelist):
            #    ncpu = len(filelist)
            ncpu = 1
            timporteovsa.importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1,
                                     visprefix=outpath, nocreatms=False,
                                     doconcat=False, modelms="", doscaling=doscaling, keep_nsclms=False, udb_corr=True)

        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))

        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': filelist,
                'ms': [outpath + ll + '.ms' for ll in sorted(list(msfiles))], 'tstlist': sclist['tstlist'],
                'tedlist': sclist['tedlist']}


def calib_pipeline(trange, workdir=None, doimport=False, overwrite=False, clearcache=False, verbose=False, pols='XX'):
    ''' 
       trange: can be 1) a single Time() object: use the entire day
                      2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                      3) a single or a list of UDBms file(s)
                      4) None -- use current date Time.now()
    '''

    if workdir is None:
        workdir = '/data1/workdir'
    os.chdir(workdir)
    if type(trange) == Time:
        mslist = trange2ms(trange=trange, doimport=False)
        invis = mslist['ms']
    if type(trange) == str:
        try:
            mslist = trange2ms(trange=trange, doimport=False)
            invis = mslist['ms']
        except:
            invis = [trange]

    for idx, f in enumerate(invis):
        if f[-1] == '/':
            invis[idx] = f[:-1]

    if overwrite or (invis == []):
        if type(trange) == Time:
            mslist = trange2ms(trange=trange, doimport=doimport, overwrite=overwrite)
            invis = mslist['ms']
        if type(trange) == str:
            try:
                mslist = trange2ms(trange=trange, doimport=doimport, overwrite=overwrite)
                invis = mslist['ms']
            except:
                invis = [trange]

        for idx, f in enumerate(invis):
            if f[-1] == '/':
                invis[idx] = f[:-1]

        outputvis = os.path.join(os.path.dirname(invis[0]), os.path.basename(invis[0])[:11] + '.ms')
        vis = calibeovsa.calibeovsa(invis, caltype=['refpha', 'phacal'], caltbdir=caltbdir, interp='nearest',
                                    doflag=True,
                                    flagant='13~15',
                                    doimage=False, doconcat=True,
                                    concatvis=outputvis, keep_orig_ms=False)
    else:
        vis = invis[0]

    udbmspath = udbmsslfcaleddir
    tdate = mstl.get_trange(vis)[0]
    outpath = os.path.join(udbmspath, tdate.datetime.strftime('%Y%m')) + '/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    imgoutdir = os.path.join(qlookfitsdir, tdate.datetime.strftime("%Y/%m/%d/"))
    if not os.path.exists(imgoutdir):
        os.makedirs(imgoutdir)
    figoutdir = os.path.join(synopticfigdir, tdate.datetime.strftime("%Y/"))
    if not os.path.exists(figoutdir):
        os.makedirs(figoutdir)
    if verbose:
        print('input of pipeline_run:')
        print({'vis': vis,
               'outputvis': outpath + os.path.basename(invis[0])[:11] + '.ms',
               'workdir': workdir,
               'slfcaltbdir': os.path.join(slfcaltbdir, tdate.datetime.strftime('%Y%m')) + '/',
               'imgoutdir': imgoutdir,
               'figoutdir': figoutdir})
    vis = ed.pipeline_run(vis, outputvis=outpath + os.path.basename(invis[0])[:11] + '.ms',
                          workdir=workdir,
                          slfcaltbdir=os.path.join(slfcaltbdir, tdate.datetime.strftime('%Y%m')) + '/',
                          imgoutdir=imgoutdir, figoutdir=figoutdir, clearcache=clearcache, pols=pols)
    return vis


def mk_qlook_image(trange, doimport=False, docalib=False, ncpu=10, twidth=12, stokes=None, antenna='0~12',
                   lowcutoff_freq=3.7, imagedir=None, spws=['1~5', '6~10', '11~15', '16~25'], toTb=True, overwrite=True,
                   doslfcal=False, verbose=False):
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
    if type(trange) == str:
        try:
            date = Time(trange)
            mslist = trange2ms(trange=trange, doimport=doimport)
            vis = mslist['ms']
            tsts = [l.to_datetime() for l in mslist['tstlist']]
        except:
            vis = [trange]
            tsts = []
            for v in vis:
                tb.open(v + '/OBSERVATION')
                tsts.append(Time(tb.getcell('TIME_RANGE')[0] / 24 / 3600, format='mjd').datetime)
                tb.close()
    subdir = [tst.strftime("%Y/%m/%d/") for tst in tsts]

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
        cfreqs = getspwfreq(msfile)
        for spw in spws:
            antenna = antenna0
            if spw == '':
                continue
            spwran = [s.zfill(2) for s in spw.split('~')]
            freqran = [cfreqs[int(s)] for s in spw.split('~')]
            cfreq = np.mean(freqran)
            bmsz = max(150. / cfreq, 20.)
            uvrange = '<10klambda'
            if doslfcal:
                slfcal_img = './' + msfilebs + '.slf.spw' + spw.replace('~', '-') + '.slfimg'
                slfcal_tb = './' + msfilebs + '.slf.spw' + spw.replace('~', '-') + '.slftb'
                try:
                    clean(vis=slfcalms, antenna=antenna, imagename=slfcal_img, spw=spw, mode='mfs', timerange='',
                          imagermode='csclean',
                          psfmode='clark', imsize=[512, 512], cell=['5arcsec'], niter=100, gain=0.05, stokes='I',
                          weighting='natural',
                          restoringbeam=[str(bmsz) + 'arcsec'], pbcor=False, interactive=False, usescratch=True)
                except:
                    print('error in cleaning spw: ' + spw)
                    break
                gaincal(vis=slfcalms, refant='0', antenna=antenna, caltable=slfcal_tb, spw=spw, uvrange='',
                        gaintable=[], selectdata=True,
                        timerange='', solint='600s', gaintype='G', calmode='p', combine='', minblperant=3, minsnr=2,
                        append=False)
                if not os.path.exists(slfcal_tb):
                    print('No solution found in spw: ' + spw)
                    break
                else:
                    clearcal(slfcalms)
                    delmod(slfcalms)
                    applycal(vis=slfcalms, gaintable=[slfcal_tb], spw=spw, selectdata=True, antenna=antenna,
                             interp='nearest', flagbackup=False,
                             applymode='calonly', calwt=False)
                    msfile = slfcalms

            imsize = 512
            cell = ['5arcsec']
            if len(spwran) == 2:
                spwstr = spwran[0] + '~' + spwran[1]
            else:
                spwstr = spwran[0]

            restoringbeam = ['{0:.1f}arcsec'.format(bmsz)]
            imagesuffix = '.spw' + spwstr.replace('~', '-')
            if cfreq > 10.:
                antenna = antenna + ';!0&1;!0&2'  # deselect the shortest baselines
            # else:
            #     antenna = antenna + ';!0&1'  # deselect the shortest baselines

            res = ptclean3(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, uvrange=uvrange,
                           spw=spw, ncpu=ncpu, niter=1000,
                           gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=stokes, doreg=True,
                           usephacenter=False, overwrite=overwrite,
                           toTb=toTb, restoringbeam=restoringbeam, specmode="mfs", deconvolver="hogbom",
                           datacolumn='data', pbcor=True)

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

        cfreqs = getspwfreq(msfile)
        imdir = imagedir + subdir[0]
        if not os.path.exists(imdir):
            os.makedirs(imdir)
        for spw in spws:
            antenna = antenna0
            if spw == '':
                spw = '{:d}~{:d}'.format(next(x[0] for x in enumerate(cfreqs) if x[1] > lowcutoff_freq),
                                         len(cfreqs) - 1)
            spwran = [s.zfill(2) for s in spw.split('~')]
            freqran = [cfreqs[int(s)] for s in spw.split('~')]
            cfreq = np.mean(freqran)
            bmsz = max(150. / cfreq, 20.)
            uvrange = ''
            imsize = 512
            cell = ['5arcsec']
            if len(spwran) == 2:
                spwstr = spwran[0] + '~' + spwran[1]
            else:
                spwstr = spwran[0]

            restoringbeam = ['{0:.1f}arcsec'.format(bmsz)]
            imagesuffix = '.synoptic.spw' + spwstr.replace('~', '-')
            antenna = antenna + ';!0&1'  # deselect the shortest baselines

            res = ptclean3(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=len(tim), uvrange=uvrange,
                           spw=spw, ncpu=1, niter=0,
                           gain=0.05, antenna=antenna, imsize=imsize, cell=cell, stokes=stokes, doreg=True,
                           usephacenter=False, overwrite=overwrite,
                           toTb=toTb, restoringbeam=restoringbeam, specmode="mfs", deconvolver="hogbom",
                           datacolumn='data', pbcor=True)
            if res:
                imres['Synoptic']['Succeeded'] += res['Succeeded']
                imres['Synoptic']['BeginTime'] += res['BeginTime']
                imres['Synoptic']['EndTime'] += res['EndTime']
                imres['Synoptic']['ImageName'] += res['ImageName']
                imres['Synoptic']['Spw'] += [spwstr] * len(res['ImageName'])
                imres['Synoptic']['Vis'] += [msfile] * len(res['ImageName'])
            else:
                continue

    # save it for debugging purposes
    np.savez('imres.npz', imres=imres)

    return imres


def plt_qlook_image(imres, figdir=None, verbose=True, synoptic=False):
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    from sunpy import sun
    from matplotlib import colors
    import astropy.units as u
    from suncasa.utils import plot_mapX as pmX
    # from matplotlib import gridspec as gridspec

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
    if verbose:
        print('{0:d} figures to plot'.format(ntime))
    plt.ioff()
    fig = plt.figure(figsize=(8, 8))

    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    axs = []
    ims = []
    pltst = 0
    for i in range(ntime):
        plt.ioff()
        plttime = btimes_sort[i, 0]
        tofd = plttime.mjd - np.fix(plttime.mjd)
        suci = suc_sort[i]
        if not synoptic:
            if tofd < 16. / 24. or sum(
                    suci) < nspw - 2:  # if time of the day is before 16 UT (and 24 UT), skip plotting (because the old antennas are not tracking)
                continue
            else:
                if pltst == 0:
                    i0 = i
                    pltst = 1
        else:
            if pltst == 0:
                i0 = i
                pltst = 1
        if i == i0:
            if synoptic:
                timetext = fig.text(0.01, 0.98, plttime.iso[:10], color='w', fontweight='bold', fontsize=12, ha='left')
            else:
                timetext = fig.text(0.01, 0.98, plttime.iso[:19], color='w', fontweight='bold', fontsize=12, ha='left')
        else:
            if synoptic:
                timetext.set_text(plttime.iso[:10])
            else:
                timetext.set_text(plttime.iso[:19])
        if verbose:
            print('Plotting image at: ', plttime.iso)
        for n in range(nspw):
            plt.ioff()
            if i == i0:
                if nspw == 1:
                    ax = fig.add_subplot(111)
                else:
                    ax = fig.add_subplot(nspw / 2, 2, n + 1)
                axs.append(ax)
            else:
                ax = axs[n]
            image = images_sort[i, n]
            if suci[n] or os.path.exists(image):
                try:
                    eomap = smap.Map(image)
                except:
                    continue
                data = eomap.data
                sz = data.shape
                if len(sz) == 4:
                    data = data.reshape((sz[2], sz[3]))
                data[np.isnan(data)] = 0.0
                # add a basin flux to the image to avoid negative values
                data = data + 0.8e5
                data[data < 0] = 0.0
                data = np.sqrt(data)
                eomap = smap.Map(data, eomap.meta)
                # resample the image for plotting
                dim = u.Quantity([256, 256], u.pixel)
                eomap = eomap.resample(dim)
            else:
                # make an empty map
                data = np.zeros((256, 256))
                header = {"DATE-OBS": plttime.isot, "EXPTIME": 0., "CDELT1": 10., "NAXIS1": 256, "CRVAL1": 0.,
                          "CRPIX1": 128.5, "CUNIT1": "arcsec",
                          "CTYPE1": "HPLN-TAN", "CDELT2": 10., "NAXIS2": 256, "CRVAL2": 0., "CRPIX2": 128.5,
                          "CUNIT2": "arcsec", "CTYPE2": "HPLT-TAN",
                          "HGLT_OBS": sun.heliographic_solar_center(plttime)[1].value, "HGLN_OBS": 0.,
                          "RSUN_OBS": sun.solar_semidiameter_angular_size(plttime).value,
                          "RSUN_REF": sun.constants.radius.value,
                          "DSUN_OBS": sun.sunearth_distance(plttime).to(u.meter).value, }
                eomap = smap.Map(data, header)
            if i == i0:
                eomap_ = pmX.Sunmap(eomap)
                # im = eomap_.imshow(axes=ax, cmap='jet', norm=colors.LogNorm(vmin=0.1, vmax=1e8))
                im = eomap_.imshow(axes=ax, cmap='jet', norm=colors.Normalize(vmin=150, vmax=700))
                ims.append(im)
                if not synoptic:
                    eomap_.draw_limb(axes=ax)
                eomap_.draw_grid(axes=ax)
                ax.set_xlim([-1080, 1080])
                ax.set_ylim([-1080, 1080])
                try:
                    cfreq = eomap.meta['crval3'] / 1.0e9
                    bdwid = eomap.meta['cdelt3'] / 1.0e9
                    ax.text(0.98, 0.01, '{0:.1f} - {1:.1f} GHz'.format(cfreq - bdwid / 2.0, cfreq + bdwid / 2.0),
                            color='w', transform=ax.transAxes, fontweight='bold', ha='right')
                except:
                    pass
                ax.set_title(' ')
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticklabels([''])
                ax.set_yticklabels([''])
            else:
                ims[n].set_data(eomap.data)

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
            print('Saving plot to :' + figdir_ + figname)

        plt.savefig(figdir_ + figname)
    plt.close(fig)


def qlook_image_pipeline(date, twidth=10, ncpu=15, doimport=False, docalib=False, synoptic=False, overwrite=True):
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

    if date.mjd > Time('2019-02-02 12:00:00').mjd:
        ## the last '' window is for fullBD synthesis image. Now obsolete.
        # spws = ['6~10', '11~20', '21~30', '31~43', '']
        # spws = ['6~10', '11~20', '21~30', '31~43']
        spws = ['0~1', '2~5', '6~10', '11~20', '21~30', '31~49']
    else:
        ## the last '' window is for fullBD synthesis image. Now obsolete.
        # spws = ['1~5', '6~10', '11~15', '16~25', '']
        spws = ['1~3', '4~6', '7~10', '10~14', '15~20', '21~30']

    if docalib:
        vis = calib_pipeline(date, doimport=doimport, synoptic=synoptic)

    imagedir = qlookfitsdir
    if synoptic:
        vis_synoptic = os.path.join(udbmsdir, date.datetime.strftime("%Y%m"),
                                    'UDB' + date.datetime.strftime("%Y%m%d") + '.ms')
        if os.path.exists(vis_synoptic):
            date = vis_synoptic
        else:
            print('Whole-day ms file {} not existed. About..... Use pipeline to make one.'.format(vis_synoptic))
            return None

    imres = mk_qlook_image(date, twidth=twidth, ncpu=ncpu, doimport=doimport, docalib=docalib, imagedir=imagedir,
                           spws=spws, verbose=True, overwrite=overwrite)

    figdir = qlookfigdir
    plt_qlook_image(imres, figdir=figdir, verbose=True)
    # if imres['Synoptic']['Succeeded']:
    #     figdir = synopticfigdir
    #     imres_bds = {}
    #     # imres_allbd = {}
    #     for k, v in imres['Synoptic'].items():
    #         imres_bds[k] = v  # [:4]
    #         # imres_allbd[k] = v[4:]
    #
    #     plt_qlook_image(imres_bds, figdir=figdir, verbose=True, synoptic=True)
    #     # plt_qlook_image(imres_allbd, figdir=figdir + 'FullBD/', verbose=True, synoptic=True)


def pipeline(year=None, month=None, day=None, ndays=1, clearcache=True, overwrite=True, doimport=True, pols='XX'):
    workdir = '/data1/workdir/'
    os.chdir(workdir)
    # Set to run 5 days earlier than the current date
    if year is None:
        mjdnow = Time.now().mjd
        t = Time(mjdnow - 2, format='mjd')
    else:
        # Uncomment below and set date to run for a given date
        t = Time('{}-{:02d}-{:02d} 20:00'.format(year, month, day))
    for d in range(ndays):
        t1 = Time(t.mjd - d, format='mjd')
        datestr = t1.iso[:10]
        subdir = os.path.join(workdir, t1.datetime.strftime('%Y%m%d/'))
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        else:
            os.system('rm -rf {}/*'.format(subdir))
        vis_corrected = calib_pipeline(datestr, overwrite=overwrite, doimport=doimport,
                                       workdir=subdir, clearcache=False, pols=pols)
        if clearcache:
            os.chdir(workdir)
            os.system('rm -rf {}'.format(subdir))


if __name__ == '__main__':
    '''
    Name: 
    eovsa_pipeline --- main pipeline for importing and calibrating EOVSA visibility data.
    
    Synopsis:
    eovsa_pipeline.py [options]... [DATE_IN_YY_MM_DD]
    
    Description:
    Import and calibrate EOVSA visibility data of the date specified
    by DATE_IN_YY_MM_DD (or from ndays before the DATE_IN_YY_MM_DD if option --ndays/-n is provided).
    If DATE_IN_YY_MM_DD is omitted, it will be set to 2 days before now by default. 
    The are no mandatory arguments in this command.
    
    -c, --clearcache
            Remove temporary files
            
    -n, --ndays
            Processing the date spanning from DATE_IN_YY_MM_DD-ndays to DATE_IN_YY_MM_DD. Default is 30.
            
    -o, --overwrite
            If True, overwrite imported and calibrated ms. Reprocess the date from scratch.
            Syntax: True, False, T, F, 1, 0
            
    -i, --doimport
            If False, skip the import step. overwrite the calibrated ms. Reprocess the date from the imported ms.
            Syntax: True, False, T, F, 1, 0
                      
    
    Example: 
    eovsa_pipeline.py -c True -n 1 -o True -i True 2020 06 10
    '''
    import sys
    import numpy as np
    import getopt

    year = None
    month = None
    day = None
    ndays = 1
    clearcache = True
    overwrite = True
    doimport = True
    pols = 'XX'

    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(argv, "c:n:o:i:p:", ['clearcache=', 'ndays=', 'overwrite=', 'doimport=', 'pols='])
        print(opts, args)
        for opt, arg in opts:
            print(opt, arg, type(arg))
            if opt in ['-c', '--clearcache']:
                if arg in ['True', 'T', '1']:
                    clearcache = True
                elif arg in ['False', 'F', '0']:
                    clearcache = False
                else:
                    clearcache = np.bool(arg)
            elif opt in ('-n', '--ndays'):
                ndays = np.int(arg)
            elif opt in ('-o', '--overwrite'):
                if arg in ['True', 'T', '1']:
                    overwrite = True
                elif arg in ['False', 'F', '0']:
                    overwrite = False
                else:
                    overwrite = np.bool(arg)
            elif opt in ('-i', '--doimport'):
                if arg in ['True', 'T', '1']:
                    doimport = True
                elif arg in ['False', 'F', '0']:
                    doimport = False
                else:
                    doimport = np.bool(arg)
            elif opt in ('-p', '--pols'):
                if arg in ['XX', 'XXYY']:
                    pols = arg
        nargs = len(args)
        if nargs == 3:
            year = np.int(args[0])
            month = np.int(args[1])
            day = np.int(args[2])
        else:
            year = None
            month = None
            day = None
    except getopt.GetoptError as err:
        print(err)
        print('Error interpreting command line argument')
        year = None
        month = None
        day = None
        ndays = 1
        clearcache = True
        overwrite = True
        doimport = True
        pols = 'XX'

    print("Running pipeline_plt for date {}-{}-{}.".format(year, month, day))
    kargs = {'ndays': ndays,
             'clearcache': clearcache,
             'overwrite': overwrite,
             'doimport': doimport,
             'pols': pols}
    for k, v in kargs.items():
        print(k, v)
    pipeline(year, month, day, ndays=ndays, clearcache=clearcache, overwrite=overwrite, doimport=doimport, pols=pols)
