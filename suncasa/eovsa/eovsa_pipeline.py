import argparse
from datetime import datetime, timedelta
# from astropy.time import Time
from suncasa.suncasatasks import ptclean6 as ptclean
from suncasa.suncasatasks import calibeovsa
from suncasa.suncasatasks import importeovsa

import re
import sys
import numpy as np
from eovsapy.dump_tsys import findfiles
from eovsapy.util import Time
import os
from suncasa.eovsa import eovsa_diskmodel as ed
from suncasa.utils import mstools as mstl

from suncasa.casa_compat import import_casatasks, import_casatools
from suncasa.eovsa.update_log import  EOVSA15_UPGRADE_DATE,DCM_IF_FILTER_UPGRADE_DATE

tasks = import_casatasks('split', 'tclean', 'gencal', 'clearcal', 'applycal', 'gaincal',
                         'delmod')
split = tasks.get('split')
tclean = tasks.get('tclean')
gencal = tasks.get('gencal')
clearcal = tasks.get('clearcal')
applycal = tasks.get('applycal')
gaincal = tasks.get('gaincal')
delmod = tasks.get('delmod')

tools = import_casatools(['qatool', 'iatool', 'mstool', 'tbtool'])
qatool = tools['qatool']
iatool = tools['iatool']
mstool = tools['mstool']
tbtool = tools['tbtool']
ms = mstool()
tb = tbtool()


def get_tdate_from_basename(vis):
    # Define the regular expression pattern
    pattern = r'UDB(\d{8})\.ms(\.tar\.gz)?'

    # Extract the basename from the vis path
    basename = os.path.basename(vis)

    # Search for the pattern in the basename
    match = re.search(pattern, basename)

    if match:
        # Extract the date string
        date_str = match.group(1)

        # Convert the date string to a datetime object
        tdate = datetime.strptime(date_str, '%Y%m%d').replace(hour=20, minute=0, second=0)

        return tdate
    else:
        raise ValueError("The basename does not match the expected format.")


import socket
import os

hostname = socket.gethostname()
is_on_server = hostname in ['pipeline', 'inti.hpcnet.campus.njit.edu']
is_on_inti = hostname == 'inti.hpcnet.campus.njit.edu'
if is_on_server:
    base_dir = "/inti/data/pipeline_mirror" if is_on_inti else ""
else:
    base_dir = './'


class Path_config:
    def __init__(self, base_dir=base_dir):
        self.paths = {}

        self.base_dir = base_dir
        # Setting paths
        self.udbmsdir = self._get_env_var('EOVSAUDBMS', f'{base_dir}/data1/eovsa/fits/UDBms/')
        self.udbmsscldir = self._get_env_var('EOVSAUDBMSSCL', f'{base_dir}/data1/eovsa/fits/UDBms_scl/')
        self.udbmsslfcaleddir = self._get_env_var('EOVSAUDBMSSLFCALED', f'{base_dir}/data1/eovsa/fits/UDBms_slfcaled/')
        self.udbdir = self._get_env_var('EOVSAUDB', f'{base_dir}/data1/eovsa/fits/UDB/')
        self.caltbdir = self._get_env_var('EOVSACAL', f'{base_dir}/data1/eovsa/caltable/')
        self.slfcaltbdir = self._get_env_var('EOVSASLFCAL', f'{base_dir}/data1/eovsa/slfcaltable/')
        self.qlookfitsdir = self._get_env_var('EOVSAQLOOKFITS', f'{base_dir}/data1/eovsa/fits/synoptic/')
        self.qlookfigdir = self._get_env_var('EOVSAQLOOKFIG', f'{base_dir}/common/webplots/qlookimg_10m/')
        self.synopticfigdir = self._get_env_var('EOVSASYNOPTICFIG', f'{base_dir}/common/webplots/SynopticImg/')
        self.workdir_default = self._get_env_var('EOVSAWORKDIR', f'{base_dir}/data1/workdir/')

        # Print a summary of paths
        self._print_summary()

    def _get_env_var(self, env_var, default_path):
        path = os.getenv(env_var) or default_path
        if not os.path.exists(path):
            # if not is_on_server:
            #     path = os.path.basename(default_path.rstrip('/')) + '/'
            #     if not os.path.exists(path):
            #         os.makedirs(path)
            # else:
            os.makedirs(path)
        self.paths[env_var] = path
        return path

    def _print_summary(self):
        print("Paths Configuration Summary:")
        for env_var, path in self.paths.items():
            print(f"  {env_var}: {path}")


# Usage
pathconfig = Path_config()
print(pathconfig.udbmsdir)  # Accessing the directory path

udbmsdir = pathconfig.udbmsdir
udbmsscldir = pathconfig.udbmsscldir
udbmsslfcaleddir = pathconfig.udbmsslfcaleddir
udbdir = pathconfig.udbdir
caltbdir = pathconfig.caltbdir
slfcaltbdir = pathconfig.slfcaltbdir
qlookfitsdir = pathconfig.qlookfitsdir
qlookfigdir = pathconfig.qlookfigdir
synopticfigdir = pathconfig.synopticfigdir
workdir_default = pathconfig.workdir_default


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

    # # in case of a single Time object was passed, adjusting the time range to start at the local beginning of the day.
    # # Initially, 'trange' is set to start at 20:00:00 UTC of the given day.
    # # After adjustment, 'trange' spans from 08:00:00 UTC of the same day (start of the local day)
    # # to 08:00:00 UTC of the following day, covering the entire local day.
    # # Example of change:
    # # Before: trange = "yyyy-mm-dd 20:00:00.000" (single starting point)
    # # After:  trange = ["yyyy-mm-dd 08:00:00.000", "yyyy-mm-dd+1 08:00:00.000"] (full day range)
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

    print('Selected idb files in the  time range (UTC): ', trange.iso)

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
        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms*'.format(outpath)) if
                   ll.endswith('.ms') or ll.endswith('.ms.tar.gz')]

    msfile_synoptic = os.path.join(outpath, 'UDB' + tdatetime.strftime("%Y%m%d") + '.ms')

    if os.path.exists(msfile_synoptic):
        if overwrite and doimport:
            os.system(f'rm -rf {msfile_synoptic}*')

    # sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
    sclist = findfiles(trange, projid='NormalObserving', srcid='Sun')
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
            # import multiprocessing as mprocs
            # ncpu = mprocs.cpu_count()
            # if ncpu > 10:
            #    ncpu = 10
            # if ncpu > len(filelist):
            #    ncpu = len(filelist)
            ncpu = 1
            importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1,
                        visprefix=outpath, nocreatms=False,
                        doconcat=False, modelms="", doscaling=doscaling, keep_nsclms=False, udb_corr=True)

        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms*'.format(outpath)) if
                   ll.endswith('.ms') or ll.endswith('.ms.tar.gz')]
        udbfilelist_set = set(udbfilelist)
        msfiles = udbfilelist_set.intersection(msfiles)
        filelist = udbfilelist_set - msfiles
        filelist = sorted(list(filelist))

        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': filelist,
                'ms': [outpath + ll + '.ms' for ll in sorted(list(msfiles))], 'tstlist': sclist['tstlist'],
                'tedlist': sclist['tedlist']}


def calib_pipeline(trange, workdir=None, doimport=False, overwrite=False, clearcache=False, verbose=False, pols='XX',
                   version='v3.0', ncpu='auto', caltype=['refpha', 'phacal'], interp='nearest'):
    ''' 
       trange: can be 1) a single Time() object: use the entire day
                      2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                      3) a single or a list of UDBms file(s)
                      4) None -- use current date Time.now()
    '''

    if workdir is None:
        workdir = workdir_default
    os.chdir(workdir)
    if isinstance(trange, Time):
        mslist = trange2ms(trange=trange, doimport=False)
        invis = mslist['ms']
    if type(trange) == str:
        try:
            mslist = trange2ms(trange=trange, doimport=False)
            invis = mslist['ms']
        except:
            invis = [trange]

    for idx, f in enumerate(invis):
        invis[idx] = f.rstrip('/')



    if overwrite or (invis == []):
        if isinstance(trange, Time):
            mslist = trange2ms(trange=trange, doimport=doimport, overwrite=overwrite)
            invis = mslist['ms']
        if type(trange) == str:
            try:
                mslist = trange2ms(trange=trange, doimport=doimport, overwrite=overwrite)
                invis = mslist['ms']
            except:
                invis = [trange]

        for idx, f in enumerate(invis):
            invis[idx] = f.rstrip('/')


        outputvis = os.path.join(os.path.dirname(invis[0]), os.path.basename(invis[0])[:11] + '.ms')
        tdate = get_tdate_from_basename(outputvis)
        flagant = '13~15' if Time(tdate).mjd >= EOVSA15_UPGRADE_DATE.mjd else '15'
        vis = calibeovsa(invis, caltype=caltype, caltbdir=caltbdir, interp=interp,
                         doflag=True,
                         flagant=flagant,
                         doimage=False, doconcat=True,
                         concatvis=outputvis, keep_orig_ms=False)
    else:
        vis = invis[0]

    udbmspath = udbmsslfcaleddir
    # tdate = mstl.get_trange(vis)[0]
    tdate = get_tdate_from_basename(vis)

    outpath = os.path.join(udbmspath, tdate.strftime('%Y%m')) + '/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    imgoutdir = os.path.join(qlookfitsdir, tdate.strftime("%Y/%m/%d/"))
    if not os.path.exists(imgoutdir):
        os.makedirs(imgoutdir)
    figoutdir = os.path.join(synopticfigdir, tdate.strftime("%Y/"))
    if not os.path.exists(figoutdir):
        os.makedirs(figoutdir)

    if version == 'v1.0':
        output_file_path = outpath + os.path.basename(invis[0])[:11] + '.ms'
    else:
        output_file_path = outpath + os.path.basename(invis[0])[:11] + f'.{version}.ms'
    slfcaltbdir_path = os.path.join(slfcaltbdir, tdate.strftime('%Y%m')) + '/'

    if verbose:
        print('input of pipeline_run:')
        print({'vis': vis,
               'outputvis': output_file_path,
               'workdir': workdir,
               'slfcaltbdir': slfcaltbdir_path,
               'imgoutdir': imgoutdir,
               'figoutdir': figoutdir,
               'overwrite': overwrite,
               'clearcache': clearcache,
               'pols': pols, 'ncpu': ncpu})
    if version == 'v1.0':
        vis = ed.pipeline_run(vis, outputvis=output_file_path,
                              workdir=workdir,
                              slfcaltbdir=slfcaltbdir_path,
                              imgoutdir=imgoutdir, figoutdir=figoutdir, clearcache=clearcache, pols=pols)
    elif version == 'v2.0':
        from suncasa.eovsa import eovsa_synoptic_imaging_pipeline as esip
        vis = esip.pipeline_run(vis, outputvis=output_file_path,
                                workdir=workdir,
                                slfcaltbdir=slfcaltbdir_path,
                                imgoutdir=imgoutdir, figoutdir=figoutdir, clearcache=clearcache, pols=pols, ncpu=ncpu,
                                overwrite=overwrite)
    elif version == 'v3.0':
        from suncasa.eovsa import eovsa_synoptic_imaging_pipeline_wsclean as esip
        vis = esip.pipeline_run(vis, outputvis=output_file_path,
                                workdir=workdir,
                                slfcaltbdir=slfcaltbdir_path,
                                imgoutdir=imgoutdir, pols=pols, overwrite=overwrite)
        if clearcache:
            os.system(f'rm -rf {workdir}/*')
    else:
        print(f'Version {version} is not supported. Valid versions are v1.0, v2.0, and v3.0. Use the default version 1.0.')
        vis = ed.pipeline_run(vis, outputvis=output_file_path,
                              workdir=workdir,
                              slfcaltbdir=slfcaltbdir_path,
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
    if isinstance(trange, Time):
        mslist = trange2ms(trange=trange, doimport=doimport)
        vis = mslist['ms']
        tsts = [l.to_datetime() for l in mslist['tstlist']]
    if isinstance(trange, str):
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
                    tclean(vis=slfcalms, antenna=antenna, imagename=slfcal_img, spw=spw, mode='mfs', timerange='',
                           deconvolver='hogbom',
                           imsize=[512, 512], cell=['5arcsec'], niter=100, gain=0.05, stokes='I',
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

            res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, uvrange=uvrange,
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
        ms.done()

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

            res = ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=len(tim), uvrange=uvrange,
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

    if date.mjd >= DCM_IF_FILTER_UPGRADE_DATE.mjd:
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


def pipeline(year=None, month=None, day=None, ndays=1, clearcache=True, overwrite=False, doimport=True, pols='XX',
             version='v1.0', ncpu='auto', debugging=False, caltype=['refpha', 'phacal'], interp='nearest'):
    """
    Main pipeline for importing and calibrating EOVSA visibility data.

    Name:
        eovsa_pipeline --- main pipeline for importing and calibrating EOVSA visibility data.

    Synopsis:
        eovsa_pipeline.py [options]... [DATE_IN_YY_MM_DD]

    Description:
        Import and calibrate EOVSA visibility data of the date specified
        by DATE_IN_YY_MM_DD (or from ndays before the DATE_IN_YY_MM_DD if option --ndays/-n is provided).
        If DATE_IN_YY_MM_DD is omitted, it will be set to 2 days before now by default.
        There are no mandatory arguments in this command.

    :param year: The year for which data should be processed, defaults to None.
    :type year: int, optional
    :param month: The month for which data should be processed, defaults to None.
    :type month: int, optional
    :param day: The day for which data should be processed, defaults to None.
    :type day: int, optional
    :param ndays: Number of days before the specified date to include in the processing, defaults to 1.
    :type ndays: int, optional
    :param clearcache: Whether to clear cache after processing, defaults to True.
    :type clearcache: bool, optional
    :param overwrite: Whether to overwrite existing files, defaults to True.
    :type overwrite: bool, optional
    :param doimport: Whether to perform the import step, defaults to True.
    :type doimport: bool, optional
    :param pols: Polarizations to process, can be 'XX' or 'XXYY', defaults to 'XX'.
    :type pols: str, optional
    :param version: Version of the pipeline to use, choices are 'v1.0' or 'v2.0', defaults to 'v1.0'.
    :type version: str, optional
    :param ncpu: Number of CPUs to use for processing, defaults to 'auto'.
    :type ncpu: str, optional
    :param debugging: Whether to run the pipeline in debugging mode, defaults to False.
    :type debugging: bool, optional
    :param caltype: Calibration types to use, defaults to ['refpha','phacal'].
    :type caltype: list, optional
    :param interp: Interpolation method to use for calibration tables, defaults to 'nearest'. Options are 'nearest', 'linear'
    :type interp: str, optional

    :raises ValueError: Raises an exception if the date parameters are out of the valid Gregorian calendar range.

    Example:
    --------
    To process data for November 24th, 2021 using version 2.0 of the pipeline, with all options enabled:

    >>> python eovsa_pipeline.py --date 2021-11-24T20:00 --clearcache --overwrite --doimport --pols XX --version v2.0 --ndays 2

    If you want to see the help message, you can run:

    >>> python eovsa_pipeline.py -h
    """
    workdir = workdir_default
    os.chdir(workdir)
    if year is None:
        # Default behavior: Process data from one day prior to the current date.
        # Calculate the Modified Julian Date (MJD) for yesterday.
        mjdnow = Time.now().mjd - 1
        # Convert MJD to a datetime object and format it to start processing at 20:00 UT.
        t = Time(Time(mjdnow, format='mjd').to_datetime().strftime('%Y-%m-%dT20:00'))
    else:
        t = Time('{}-{:02d}-{:02d} 20:00'.format(year, month, day))
    for d in range(ndays):
        t1 = Time(t.mjd - d, format='mjd')
        datestr = t1.iso[:10]
        subdir = os.path.join(workdir, t1.datetime.strftime('%Y%m%d/'))
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        else:
            if overwrite:
                os.system('rm -rf {}/*'.format(subdir))
        # ##debug
        # vis_corrected = calib_pipeline(datestr, overwrite=overwrite, doimport=doimport,
        #                                workdir=subdir, clearcache=False, pols=pols)

        if debugging:
            vis_corrected = calib_pipeline(t1, overwrite=overwrite, doimport=doimport,
                                           workdir=subdir, clearcache=False, pols=pols, version=version, ncpu=ncpu,
                                           caltype=caltype, interp=interp)
        else:
            try:
                vis_corrected = calib_pipeline(t1, overwrite=overwrite, doimport=doimport,
                                               workdir=subdir, clearcache=False, pols=pols, version=version, ncpu=ncpu,
                                               caltype=caltype, interp=interp)
            except Exception as e:
                print(f'error in processing {datestr}. Error message: {e}')
        if clearcache:
            os.chdir(workdir)
            os.system('rm -rf {}'.format(subdir))


if __name__ == '__main__':
    # Define the parser
    parser = argparse.ArgumentParser(description='EOVSA Pipeline for importing and calibrating visibility data.')
    # Default date is set to one day before the current date, formatted as YYYY-MM-DDT20:00
    default_date = (datetime.now() - timedelta(days=1)).strftime('%Y-%m-%dT20:00')
    parser.add_argument('--date', type=str, default=default_date,
                        help='Date to process in YYYY-MM-DDT20:00 format, defaults to 20:00 UT of one day before the current date.')
    parser.add_argument('--clearcache', action='store_true', default=False,
                        help='Remove temporary files after processing')
    parser.add_argument('--ndays', type=int, default=1,
                        help='Process data from DATE_IN_YY_MM_DD-ndays to DATE_IN_YY_MM_DD, default is 1.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite existing processed data')
    parser.add_argument('--doimport', action='store_true', default=True, help='Perform import step before processing')
    parser.add_argument('--pols', type=str, default='XX', choices=['XX', 'XXYY'], help='Polarizations to process')
    parser.add_argument('--ncpu', type=str, default='auto', help='Number of CPUs to use for processing')
    parser.add_argument('--version', type=str, default='v3.0', choices=['v1.0', 'v2.0', 'v3.0'],
                        help='Version of the EOVSA pipeline to use')
    parser.add_argument('--debugging', action='store_true', default=False, help='Run the pipeline in debugging mode')
    parser.add_argument('--caltype', type=str, nargs='+', default=['refpha', 'phacal'],
                        help='Calibration types to use, defaults to refpha and phacal')
    parser.add_argument('--interp', type=str, default='nearest', choices=['nearest', 'linear', 'auto'],
                        help='Interpolation method for calibration tables. Options: "nearest", "linear", "auto". '
                             'If "auto" is selected, calibeovsa automatically chooses the method: '
                             'if the time difference between the observation and the phacal calibrations is less than 1 hour, it uses "nearest"; '
                             'otherwise, it uses "linear".')

    # Parse the arguments
    args = parser.parse_args()

    # Convert --date argument to an astropy Time object to ensure consistency in time handling
    t = Time(args.date)

    # Extract year, month, day from the --date argument
    year, month, day = t.datetime.year, t.datetime.month, t.datetime.day

    # Run the main pipeline function
    pipeline(year, month, day, args.ndays, args.clearcache, args.overwrite, args.doimport, args.pols,
             args.version, args.ncpu, args.debugging, args.caltype, args.interp)
