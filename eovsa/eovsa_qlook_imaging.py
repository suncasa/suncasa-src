from suncasa.eovsa import eovsa_prep as ep
from ptclean_cli import ptclean_cli as ptclean
from eovsapy.util import Time
from sunpy import map as smap
from importeovsa_cli import importeovsa_cli as importeovsa
from eovsapy import refcal_anal as ra
import os

udbmsscldir = os.getenv('EOVSAUDBMSSCL')
udbdir = os.getenv('EOVSAUDB')

if not udbmsscldir:
    print 'Environmental variable for EOVSA udbms path not defined'
    print 'Use default path on pipeline'
    udbmsscldir = '/data1/eovsa/fits/UDBms_scl/'
if not udbdir:
    print 'Environmental variable for EOVSA udb path not defined'
    print 'Use default path on pipeline'
    udbdir = '/data1/eovsa/fits/UDB/'

def trange2ms(trange=None, verbose=False, doscaling=True):
    '''This finds all solar UDBms files within a timerange; If the UDBms file does not exist 
       in EOVSAUDBMSSCL, create one by calling importeovsa
       Required inputs:
       trange - can be 1) a single Time() object: use the entire day
                       2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                       4) a list of UDBms files
                       3) None -- use current date Time.now()
    '''
    import glob
    import pytz
    if trange is None:
        trange = Time.now()
    if type(trange) == list:
        try:
            trange = Time(trange)
        except:
            print('trange format not recognised. Abort....')
            return None
    local_tz = pytz.timezone('America/Los_Angeles')
    try:
        if len(trange) > 1:
            trange = Time([trange[0], trange[-1]])
            tdatetime = trange[0].to_datetime()
        else:
            tdatetime = trange[0].to_datetime()
            btime = Time(local_tz.localize(tdatetime, is_dst=None).astimezone(pytz.utc))
            etime = Time(btime.mjd + 1.0, format='mjd')
            trange = Time([btime, etime])
    except:
        tdatetime = trange.to_datetime()
        btime = Time(local_tz.localize(tdatetime, is_dst=None).astimezone(pytz.utc))
        etime = Time(btime.mjd + 1.0, format='mjd')
        trange = Time([btime, etime])

    sclist = ra.findfiles(trange, projid='NormalObserving', srcid='Sun')
    udbfilelist = sclist['scanlist']
    udbfilelist = [os.path.basename(ll) for ll in udbfilelist]
    outpath = '{}{}/'.format(udbmsscldir, tdatetime.strftime("%Y%m"))
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        msfiles = []
    else:
        msfiles = [os.path.basename(ll).split('.')[0] for ll in glob.glob('{}UDB*.ms'.format(outpath))]
    udbfilelist_set = set(udbfilelist)
    msfiles = udbfilelist_set.intersection(msfiles)
    filelist = udbfilelist_set - msfiles
    filelist = sorted(list(filelist))

    if filelist:
        import multiprocessing as mp
        ncpu = mp.cpu_count()
        if ncpu > 10:
            ncpu = 10
        if ncpu > len(filelist):
            ncpu = len(filelist)
        inpath = '{}{}/'.format(udbdir, tdatetime.strftime("%Y"))
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
                'ms': [ll + '.ms' for ll in sorted(list(msfiles))]}
    else:
        return [outpath + ll + '.ms' for ll in sorted(list(msfiles))]

def mk_qlook_image(trange, docalib=False, twidth=30, stokes=None, antenna='0~12', bds=['1~3'], savejpg=True, imagedir=None):
        
    ''' 
       trange: can be 1) a single Time() object: use the entire day
                   2) a range of Time(), e.g., Time(['2017-08-01 00:00','2017-08-01 23:00'])
                   3) a single or a list of UDBms file(s)
                   4) None -- use current date Time.now()
    '''
    if type(trange) == Time:
        invis = trange2ms(trange=trange)
    if type(invis) == str:
        invis = [trange]

    for idx, f in enumerate(invis):
        if f[-1] == '/':
            invis[idx] = f[:-1]

    if docalib:
        vis=ocalibeovsa(invis, caltype=None, interp='nearest', docalib=docalib, 
                   doflag=True, flagant='13~15', qlookimage=False, doconcat=False)
    else:
        vis=invis 

    if not stokes:
        stokes = 'XX'
     
    if not imagedir:
        imagedir='./'
    outimgs=[]
    for msfile in vis:
        for bd in bds:
            imageprefix=imagedir+'bd'+str(bd).zfill(2)+'_'
            outimg=ptclean(vis=msfile, imageprefix=imageprefix, twidth=twidth, spw=bd, ncpu=4, niter=200, imsize=[512], cell=['5arcsec'], 
                    stokes=stokes, doreg=True, usephacenter=False)
            outimgs.append(outimg)
    return outimgs

    
