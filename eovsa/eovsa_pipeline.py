from suncasa.eovsa import eovsa_prep as ep
from ptclean_cli import ptclean_cli as ptclean
from eovsapy.util import Time
from importeovsa_cli import importeovsa_cli as importeovsa
from suncasa.tasks import task_calibeovsa as calibeovsa
from eovsapy import refcal_anal as ra
import os
import numpy as np

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

def trange2ms(trange=None, doimport=False, verbose=False, doscaling=True):
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
    if type(trange) == Time:
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

    print 'Selected timerange in UTC: ',trange.iso

    sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
    udbfilelist = sclist['scanlist']
    udbfilelist = [os.path.basename(ll) for ll in udbfilelist]
    if doscaling:
        udbmspath = udbmsscldir
    else:
        udbmspath = udbmsdir
    outpath = '{}{}/'.format(udbmspath, tdatetime.strftime("%Y%m"))
    if not os.path.exists(outpath):
        if verbose:
            print outpath+' does not exist. Making a new directory.'
        os.makedirs(outpath)
        msfiles = []
    else:
        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]
    udbfilelist_set = set(udbfilelist)
    msfiles = udbfilelist_set.intersection(msfiles)
    filelist = udbfilelist_set - msfiles
    filelist = sorted(list(filelist))

    inpath = '{}{}/'.format(udbdir, tdatetime.strftime("%Y"))
    if filelist and doimport:
        import multiprocessing as mprocs
        ncpu = mprocs.cpu_count()
        if ncpu > 10:
            ncpu = 10
        if ncpu > len(filelist):
            ncpu = len(filelist)
        importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1,
                    visprefix=outpath, nocreatms=False, doconcat=False, modelms="", doscaling=doscaling,
                    keep_nsclms=False)

    msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]
    udbfilelist_set = set(udbfilelist)
    msfiles = udbfilelist_set.intersection(msfiles)
    filelist = udbfilelist_set - msfiles
    filelist = sorted(list(filelist))

    if verbose:
        return {'mspath': outpath, 'udbpath': inpath, 'udbfile': sorted(udbfilelist), 'udb2ms': filelist,
                'ms': [outpath + ll + '.ms' for ll in sorted(list(msfiles))], 
                'tstlist':sclist['tstlist'], 'tedlist':sclist['tedlist']}
    else:
        return {'ms': [outpath + ll + '.ms' for ll in sorted(list(msfiles))], 
                'tstlist':sclist['tstlist'], 'tedlist':sclist['tedlist']}

def mk_qlook_image(trange, doimport=False, docalib=False, ncpu=10, twidth=10, stokes=None, antenna='0~12', 
        imagedir=None, spws=['1~5','6~10','11~15','16~20','21~25','26~30']):

        
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
        invis = [trange]
        subdir = ['/']

    for idx, f in enumerate(invis):
        if f[-1] == '/':
            invis[idx] = f[:-1]

    if docalib:
        vis=calibeovsa.calibeovsa(invis, caltype=None, interp='nearest', docalib=docalib, 
                       doflag=True, flagant='13~15', doimage=False, doconcat=False)
    else:
        vis=invis 

    if not stokes:
        stokes = 'XX'
     
    if not imagedir:
        imagedir='./'
    imres = {'Succeeded': [], 'BeginTime': [], 'EndTime': [], 'ImageName': [], 'Spw': [], 'Vis': []}
    for n, msfile in enumerate(vis):
        imdir = imagedir + subdir[n]
        if not os.path.exists(imdir):
            os.makedirs(imdir)
        for spw in spws:
            spwran = [s.zfill(2) for s in spw.split('~')]
            if len(spwran) == 2:
                spwstr = spwran[0]+'~'+spwran[1]
            else:
                spwstr = spwran[0]
            imagesuffix='.spw'+spwstr.replace('~','-')
            res=ptclean(vis=msfile, imageprefix=imdir, imagesuffix=imagesuffix, twidth=twidth, 
                        spw=spw, ncpu=ncpu, niter=200, imsize=[512], cell=['5arcsec'], 
                        stokes=stokes, doreg=True, usephacenter=False, restoringbeam=['25arcsec'])

            if res:
                imres['Succeeded'] += res['Succeeded']
                imres['BeginTime'] += res['BeginTime']
                imres['EndTime'] += res['EndTime']
                imres['ImageName'] += res['ImageName']
                imres['Spw'] += [spwstr]*len(res['ImageName'])
                imres['Vis'] += [msfile]*len(res['ImageName'])
            else:
                return None

    return imres

def plt_qlook_image(imres,figdir=None,verbose=True): 
    from matplotlib import pyplot as plt
    from sunpy import map as smap
    from sunpy import sun
    import astropy.units as u
    if not figdir:
        figdir='./'
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
    images_sort = images[inds].reshape(ntime,nspw)
    btimes_sort = btimes[inds].reshape(ntime,nspw)
    suc_sort = suc[inds].reshape(ntime,nspw)
    spws_sort = spws[inds].reshape(ntime,nspw)
    if verbose:
        print '{0:d} figures to plot'.format(ntime)
    for i in range(ntime): 
    #for i in range(1): 
        plt.ioff()
        fig=plt.figure(figsize=(11,6))
        plttime=btimes_sort[i,0]
        fig.suptitle('EOVSA @ '+plttime.iso[:19])
        if verbose:
            print 'Plotting image at: ', plttime.iso
            suci = suc_sort[i]
        for n in range(nspw):
            plt.ioff()
            image = images_sort[i, n]
            fig.add_subplot(nspw/3, 3, n+1)
            if suci[n]:
                try:
                    eomap = smap.Map(image)
                except:
                    continue
                sz = eomap.data.shape
                if len(sz) == 4:
                    eomap.data = eomap.data.reshape((sz[2], sz[3]))
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax=plt.gca()
                ax.set_xlim([-1050,1050])
                ax.set_ylim([-1050,1050])
                ax.set_title('spw '+spws_sort[i,n])
                #ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                if n != nspw-3:
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    ax.set_xticklabels([''])
                    ax.set_yticklabels([''])
            else:
                #make an empty map
                data = np.zeros((512,512))
                header = {"DATE-OBS": plttime.isot, "EXPTIME": 0.,
                        "CDELT1": 5., "NAXIS1": 512, "CRVAL1": 0., "CRPIX1": 257, "CUNIT1": "arcsec", "CTYPE1": "HPLN-TAN",
                        "CDELT2": 5., "NAXIS2": 512, "CRVAL2": 0., "CRPIX2": 257, "CUNIT2": "arcsec", "CTYPE2": "HPLT-TAN",
                        "HGLT_OBS": sun.heliographic_solar_center(plttime)[1].value, 
                        "HGLN_OBS": 0., 
                        "RSUN_OBS": sun.solar_semidiameter_angular_size(plttime).value,
                        "RSUN_REF": sun.constants.radius.value,
                        "DSUN_OBS": sun.sunearth_distance(plttime).to(u.meter).value,
                          }
                eomap = smap.Map(data, header)
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax=plt.gca()
                ax.set_xlim([-1050,1050])
                ax.set_ylim([-1050,1050])
                ax.set_title('spw '+spws_sort[i,n])
                #ax.text(0.01,0.02, plttime.isot,transform=ax.transAxes,color='white')
                if n != 6:
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    ax.set_xticklabels([''])
                    ax.set_yticklabels([''])
        figname='eovsa_qlimg_'+plttime.isot.replace(':','').replace('-','')[:15]+'.png'
        fig_tdt= plttime.to_datetime()
        fig_subdir = fig_tdt.strftime("%Y/%m/%d/")  
        figdir_ = figdir + fig_subdir
        if not os.path.exists(figdir_):
            os.makedirs(figdir_)
        if verbose:
            print 'Saving plot to :'+figdir_+figname
        plt.savefig(figdir_+figname)
        plt.close(fig)

def qlook_image_pipeline(date, twidth=10, ncpu=15, doimport=False, docalib=False):
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
    if not qlookfitsdir:
        qlookfitsdir='/data1/eovsa/qlookfits/'
    if not qlookfigdir:
        qlookfigdir='/common/webplots/qlookimg_10m/'

    imagedir=qlookfitsdir
    imres=mk_qlook_image(date, twidth=ncpu, ncpu=ncpu, doimport=doimport, docalib=docalib, imagedir=imagedir)
    figdir=qlookfigdir
    plt_qlook_image(imres,figdir=figdir,verbose=True)
    




    
