import os,sys
import numpy as np
import pickle
import warnings
from suncasa.dspec import dspec
from suncasa.utils import qlookplot
from suncasa.utils import helioimage2fits as hf
from suncasa.utils import fitsutils
import sunpy
from sunpy import map as smap
from matplotlib import gridspec as gridspec
from matplotlib import pyplot as plt
import math
import astropy.units as u
from astropy.coordinates import SkyCoord

warnings.filterwarnings("ignore")

'''
Example script for self-calibrating EOVSA flare data. 
- Input of the script is a CASA visibility measurement set which has prior flux and gain calibration
    The example can be downloaded from the following link:
        https://drive.google.com/file/d/1i46fvsCS4B7DfoTwwL1JQkBiB5MTVVVE/view?usp=share_link
- Output of the script is a self-calibrated CASA visibility measurement set, self-calibration tables, and intermediate
    products such as images after each self-calibration rounds 
- The script involves interactive mask drawing with CASA's imview, and hence need to be run under 
    a **monolithic** CASA 6.5 environment (both Python 3.6 and Python 3.8 are tested). 
    A modular CASA can be used if interactive mask drawing is not required (e.g., you already have the masks)
- Also needed is a Python 3.6/3.8 installation of suncasa, astropy, and sunpy, which can be side loaded by 
    including your Python search path, which can be found by doing the following under your Python environment:
    > import sys
    > sys.path()
  In Anaconda Python, the search path to the installed packages is usually something like
    ~/anaconda3/envs/suncasa/lib/python3.6/site-packages/
- In CASA, append the search path as follows
    > import sys
    > sys.path.append("YOUR_PYTHON_3_SEARCH_PATH")
- Run the script in CASA
    > execfile('slfcal_flare_example.py')

History:
    2019-May-16 Bin Chen (bin.chen@njit.edu)
       Initial creation of an example script as part of the EOVSA tutorial at
       the RHESSI XVIII Workshop (http://rhessi18.umn.edu/). 
       The original script used CASA 5.4 with Python 2.7. Link to the script:
            https://github.com/binchensun/eovsa-tutorial/blob/master/rhessi18/slfcal_example.py
       
    2023-Jan-08 Bin Chen (bin.chen@njit.edu)
       Major updates for compatibility with CASA 6.5 / Python 3, along with many feature enhancements.
       This script is used for the flare self-calibration tutorial at the 2023 EOVSA/GX/FASR Workshop. 
       Link to this script on Github:
            https://github.com/binchensun/suncasa/blob/master/examples/eovsa_flare_slfcal_example.py
                         
'''


# ################ The following are function blocks needed to run the script. Change at your own risk. ###########
###################################################################################################################
def draw_masks(slfcalms, spwrans_mask=None, maskdir=None, maskfile='masks.p',
               restoringbeam=['6arcsec'], timerange=None, phasecenter='', npix=512, gain=0.05, uvrange=''):
    """Function for drawing masks interactively"""
    antennas = ''
    pol = 'I'
    imgprefix = maskdir + 'slf_t0'

    # step 1: set up the clean masks
    img_init = imgprefix + '_mask_'
    os.system('rm -rf ' + img_init + '*')

    masks = []
    imnames = []
    for spwran in spwrans_mask:
        imname = img_init + spwran.replace('~', '-')
        try:
            tclean(vis=slfcalms, antenna=antennas, imagename=imname, spw=spwran,
                   uvrange=uvrange, specmode='mfs', timerange=timerange, imsize=[npix],
                   cell=['1arcsec'], niter=1000, gain=gain, stokes=stokes,
                   restoringbeam=restoringbeam, phasecenter=phasecenter, weighting='briggs', robust=1.0,
                   interactive=True, datacolumn='corrected', pbcor=True)
            imnames.append(imname + '.image')
            masks.append(imname + '.mask')
            clnjunks = ['.flux', '.model', '.psf', '.residual']
            for clnjunk in clnjunks:
                if os.path.exists(imname + clnjunk):
                    os.system('rm -rf ' + imname + clnjunk)
        except:
            print('error in cleaning spw: ' + spwran)

    pickle.dump(masks, open(maskfile, 'wb'))
    return maskfile


def make_images(slfcalms, msinfo=None, spws=None, spwrans_mask=None, usemask=True, maskfile=None,
                imgprefix='./slf', round=0, trange='', uvrange='', antennas='',
                npix=512, cell=['1arcsec'], niter=1000, gain=0.05, robust=1.0, sbeam_1GHz=40.,
                sclfactor=1, xycen=None, fov=[256, 256], stokes='I', nrow=3,
                usemsphacenter=True, observatory='EOVSA'):
    """Function for making images based on CASA tclean and registration with suncasa"""
    if not msinfo:
        msinfo = hf.read_msinfo(slfcalms)

    phasecenter, midtime = hf.calc_phasecenter_from_solxy(slfcalms, timerange=trange, xycen=xycen,
                                                          usemsphacenter=usemsphacenter, observatory=observatory)
    msmd.open(slfcalms)
    nspw = msmd.nspw()
    cfreqsghz = [msmd.meanfreq(s, unit='GHz') for s in range(nspw)]
    msmd.done()
    # Make new images after each round
    nspw = len(spws)
    ncol = math.ceil(nspw / nrow)
    xran = [xycen[0] - fov[0] / 2., xycen[0] + fov[0] / 2.]
    yran = [xycen[1] - fov[1] / 2., xycen[1] + fov[1] / 2.]
    if xran or yran:
        xyasp = (xran[1] - xran[0]) / (yran[1] - yran[0])
    ysize = 8.
    fig = plt.figure(figsize=(ncol / nrow * ysize * xyasp, ysize))
    gs = gridspec.GridSpec(nrow, ncol)
    fitsfiles = []
    for s, sp in enumerate(spws):
        print('cleaning spw: ' + sp)
        cfreqghz = cfreqsghz[int(sp)]
        # setting restoring beam size (not very useful for selfcal anyway, but just to see the results)
        bm = max(sbeam_1GHz / cfreqghz, 3.)
        print('Beam size {0:.1f}" at ref freq {1:.2f} GHz'.format(bm, cfreqghz))
        slfcal_img = imgprefix + '.spw' + sp.zfill(2) + '.slfcal' + str(round)
        # only the first round uses nearby spws for getting initial model
        spwran = sp
        if usemask:
            try:
                masks = pickle.load(open(maskfile, 'rb'))
            except:
                print('Mask file {} does not exist. Do not use any mask'.format(maskfile))
                mask = ''
            if len(masks) > 0:
                foundmask = False
                for m, spwran_mask in enumerate(spwrans_mask):
                    sbg, sed = spwran_mask.split('~')
                    if (int(sp) >= int(sbg)) and (int(sp) <= int(sed)):
                        mask = masks[m]
                        print('using mask {0:s}'.format(os.path.basename(mask)))
                        foundmask = True
                if not foundmask:
                    print('mask not defined for {}. Do use any masks.'.format(sp))
                    mask = ''
            else:
                mask = ''
        else:
            mask = ''
        try:
            tclean(vis=slfcalms, antenna=antennas, imagename=slfcal_img, uvrange=uvrange, spw=spwran,
                   specmode='mfs', timerange='', imsize=[npix], cell=cell, niter=niter,
                   gain=gain, stokes=stokes, weighting='briggs', robust=robust, phasecenter=phasecenter, mask=mask,
                   restoringbeam=[str(bm) + 'arcsec'], pbcor=False, interactive=False)
            if os.path.exists(slfcal_img + '.image'):
                fitsfile = slfcal_img + '.fits'
                hf.imreg(vis=slfcalms, msinfo=msinfo, imagefile=slfcal_img + '.image', fitsfile=fitsfile,
                         timerange=trange, usephacenter=True, geocentric=True, sclfactor=sclfactor, toTb=True,
                         verbose=False,
                         overwrite=True)
            clnjunks = ['.mask', '.flux', '.model', '.psf', '.residual', '.image', '.pb', '.image.pbcor', '.sumwt']
            for clnjunk in clnjunks:
                if os.path.exists(slfcal_img + clnjunk):
                    os.system('rm -rf ' + slfcal_img + clnjunk)
            eomap = smap.Map(fitsfile)
            ax = fig.add_subplot(gs[s], projection=eomap)
            top_right = SkyCoord(xran[1] * u.arcsec, yran[1] * u.arcsec, frame=eomap.coordinate_frame)
            bottom_left = SkyCoord(xran[0] * u.arcsec, yran[0] * u.arcsec, frame=eomap.coordinate_frame)
            eomap_sub = eomap.submap(bottom_left, top_right=top_right)
            eomap_sub.plot_settings['cmap'] = plt.get_cmap('jet')
            eomap_sub.plot(axes=ax)
            lon = ax.coords[0]
            lat = ax.coords[1]
            ax.set_title(' ')
            ax.text(0.96, 0.96, 'SPW {0:s} at {1:.1f} GHz'.format(sp, cfreqghz), transform=ax.transAxes,
                    ha='right', va='top', color='w', fontweight='bold')
            plt.subplots_adjust(left=0.1, bottom=0.08, right=0.96, top=0.94, wspace=0, hspace=0)
            if s / ncol < nrow - 1:
                ax.set_xlabel(' ')
                lon.set_ticklabel_visible(False)
            else:
                ax.set_xlabel('Solar X (arcsec)')
            if s % ncol != 0:
                ax.set_ylabel(' ')
                lat.set_ticklabel_visible(False)
            else:
                ax.set_ylabel('Solar Y (arcsec)')
            fitsfiles.append(fitsfile)
        except:
            print('error in cleaning spw: ' + sp)
            continue

    figname = imgprefix + '_t0_n{:d}.png'.format(round)
    plt.subplots_adjust(left=0.1, bottom=0.08, right=0.96, top=0.94, wspace=0, hspace=0)
    plt.suptitle("Images after round {0:d}".format(round))
    plt.savefig(figname)
    plt.close()
    return figname, fitsfiles


def plt_slftable(slfcalms, slftbs, caltypes, docombine=False, markersize=2.0):
    """
    Function for plotting self-calibration tables

    :param slfcalms: reference ms file for the spw and antenna information
    :param slftbs: a list of self-calibration gain tables
    :param caltypes: a list of self-calibration types ('a', 'p', or 'ap')
    :param docombine: if True, combine all the calibration tables to make a combined phase plot
    :param markersize: size of the plot symbols
    :return:
        1: succeeded. -1: failed
    """
    if type(slftbs) == str:
        slftbs = [slftbs]
    if type(caltypes) == str:
        caltypes = [caltypes]

    if len(slftbs) != len(caltypes):
        print('The length of slftbs must equal to caltypes. Abort...')
        return -1

    msmd.open(slfcalms)
    nspw_ms = msmd.nspw()
    cfreqsghz = [msmd.meanfreq(s, unit='GHz') for s in range(nspw_ms)]
    nant_ms = msmd.nantennas()
    msmd.done()

    nslftb = len(slftbs)
    # now read the self-calibration tables
    gains = []
    mgains = []
    flags = []
    errs = []
    snrs = []

    nrow = 4
    for n, slftb in enumerate(slftbs):
        caltype = caltypes[n]
        tb.open(slftb)
        spws = np.unique(tb.getcol('SPECTRAL_WINDOW_ID'))
        nspw = len(spws)
        ants = np.unique(tb.getcol('ANTENNA1'))
        nant = len(ants)
        ncol = math.ceil(nant / nrow)
        cparm = tb.getcol('CPARAM')
        npol, nval, ndata = cparm.shape
        gain = cparm.reshape(npol, nspw, nant)
        flag = tb.getcol('FLAG').reshape(npol, nspw, nant)
        mgain = np.ma.masked_array(gain, mask=flag)
        err = tb.getcol('PARAMERR').reshape(npol, nspw, nant)
        snr = tb.getcol('SNR').reshape(npol, nspw, nant)
        tb.close()
        if n == 0:
            gain_tot = np.copy(gain)
            flag_tot = np.copy(flag)
        else:
            gain_tot *= gain
            flag_tot *= flag

        if 'p' in caltype:
            fig = plt.figure(figsize=(ncol / nrow * 10., 8.))
            gs = gridspec.GridSpec(nrow, ncol, left=0.08, right=0.98, wspace=0.3, hspace=0.5)
            for i in range(nant):
                ax = fig.add_subplot(gs[i])
                ax.plot(spws, np.angle(mgain[0, :, i], deg=True), 'o', ms=markersize, color='b')
                ax.plot(spws, np.angle(mgain[1, :, i], deg=True), 'o', ms=markersize, color='r')
                if i / ncol >= nrow - 1:
                    ax.set_xlabel('SPW ID')
                if i % ncol == 0:
                    ax.set_ylabel('Phase (Deg)')
                ax.set_xlim([0, nspw_ms - 1])
                ax.set_ylim([-180, 180])
                ax.set_title('Ant {0}'.format(i))
            plt.suptitle("Phase solutions in gain table {}".format(os.path.basename(slftb)))
            plt.savefig(slftb + '.pha.png')
            plt.close()

        if 'a' in caltype:
            fig = plt.figure(figsize=(ncol / nrow * 10., 8.))
            gs = gridspec.GridSpec(nrow, ncol, left=0.08, right=0.98, wspace=0.3, hspace=0.5)
            for i in range(nant):
                ax = fig.add_subplot(gs[i])
                ax.plot(spws, np.abs(mgain[0, :, i]), 'o', ms=markersize, color='b')
                ax.plot(spws, np.abs(mgain[1, :, i]), 'o', ms=markersize, color='r')
                if i / ncol >= nrow - 1:
                    ax.set_xlabel('SPW ID')
                if i % ncol == 0:
                    ax.set_ylabel('Amplitude')
                ax.set_xlim([0, nspw_ms - 1])
                ax.set_ylim([0.1, 3])
                ax.set_title('Ant {0}'.format(i))
            plt.suptitle("Amplitude solutions in gain table {}".format(os.path.basename(slftb)))
            plt.savefig(slftb + '.amp.png')
            plt.close()

        gains.append(gain)
        mgains.append(mgain)
        flags.append(flag)
        errs.append(err)
        snrs.append(snr)

    if nslftb > 1 and docombine:
        mgain_tot = np.ma.masked_array(gain_tot, mask=flag_tot)
        fig = plt.figure(figsize=(ncol / nrow * 10., 8.))
        gs = gridspec.GridSpec(nrow, ncol, left=0.08, right=0.98, wspace=0.3, hspace=0.5)
        for i in range(nant):
            ax = fig.add_subplot(gs[i])
            ax.plot(spws, np.angle(mgain_tot[0, :, i], deg=True), 'o', ms=markersize, color='b')
            ax.plot(spws, np.angle(mgain_tot[1, :, i], deg=True), 'o', ms=markersize, color='r')
            if i / ncol >= nrow - 1:
                ax.set_xlabel('SPW ID')
            if i % ncol == 0:
                ax.set_ylabel('Phase (Deg)')
            ax.set_xlim([0, nspw_ms - 1])
            ax.set_ylim([-180, 180])
            ax.set_title('Ant {0}'.format(i))
        plt.suptitle("Combined Phases")
        plt.savefig(slftbs[0] + '-' + os.path.basename(slftbs[-1]) + '.pha.comb.png')
        plt.close()

        fig = plt.figure(figsize=(ncol / nrow * 10., 8.))
        gs = gridspec.GridSpec(nrow, ncol, left=0.08, right=0.98, wspace=0.3, hspace=0.5)
        for i in range(nant):
            ax = fig.add_subplot(gs[i])
            ax.plot(spws, np.abs(mgain_tot[0, :, i]), 'o', ms=markersize, color='b')
            ax.plot(spws, np.abs(mgain_tot[1, :, i]), 'o', ms=markersize, color='r')
            if i / ncol >= nrow - 1:
                ax.set_xlabel('SPW ID')
            if i % ncol == 0:
                ax.set_ylabel('Amplitude')
            ax.set_xlim([0, nspw_ms - 1])
            ax.set_ylim([0.1, 3])
            ax.set_title('Ant {0}'.format(i))
        plt.suptitle("Combined Amplitudes")
        plt.savefig(slftbs[0] + '-' + os.path.basename(slftbs[-1]) + '.amp.comb.png')
        plt.close()

    return 1


# ####################### End of function blocks ##############################################
################################################################################################


# ############################### Step 0: Initial Definitions #####################################
# ========== declare input visibility in CASA measurement set format ==========
ms_in = 'IDB20170821201800-202300.4s.corrected.ms'

# define spectral windows used for self-calibration
spws_slf = ['3', '6', '9', '12', '15', '18', '21', '24', '27']
# Uncomment the following to do all of them (spw 0 does not have prior calibration)
# spws_slf = [str(s+1) for s in range(30)]

# ============ Prior definitions for EOVSA data ==============
# define polarization and antennas to self-calibration (no need to change for current EOVSA)
stokes = 'XX'  # polarization to be used for self-calibration
antennas = '0~12'  # antennas numbers used

# =========== TASK HANDLERS =============
dofindflare = 0  # try to find flare time and phase center?

# ============ declaring the working directories (current directory assumed. change as necessary) ============
workdir = os.getcwd() + '/'  # main working directory. Using current directory in this example
slfcaldir = workdir + 'slfcal/'  # place to put all selfcalibration products
imagedir = slfcaldir + 'images_slfcal/'  # place to put all selfcalibration images
maskdir = slfcaldir + 'masks/'  # place to put clean masks
imagedir_slfcaled = slfcaldir + 'images_slfcaled/'  # place to put final self-calibrated images
caltbdir = slfcaldir + 'caltbs/'  # place to put calibration tables
# make these directories if they do not already exist
dirs = [workdir, slfcaldir, imagedir, maskdir, imagedir_slfcaled, caltbdir]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

# ################################ Step 1: Examine data for flare timing and location ##################
if dofindflare:
    # First, make a cross-power dynamic spectrum by doing median over some appropriate uv range
    d = dspec.Dspec(ms_in, uvrange='0.2~0.8km', usetbtool=True, domedian=True)
    d.plot()
    d.tofits(ms_in + '.dspec.fits')

    # From the previous step, take an appropriate time range for doing self-calibration.
    # Here I select an 8-s piece around the flare peak.
    trange = '2017/08/21/20:21:10~2017/08/21/20:21:18'
    # Do a full-disk imaging to find a new phase center and an appropriate image size.
    qlookplot.qlookplot(ms_in, timerange=trange, spw='3GHz~10GHz', imsize=[512], cell=['5.0arcsec'],
                        niter=100, uvrange='>500lambda', plotaia=False, stokes=stokes)

# Take an appropriate time range for doing self-calibration. Here I select 10-s piece around the flare peak.
trange = '2017/08/21/20:21:10~2017/08/21/20:21:18'
# Using the image it produced, define a new phase center of the image
xycen = [375., 45.]  # in solar coordinates
# Also define the field of view in arcsec
fov = [256., 256.]
# Convert to RA & DEC coordinates (will be used in TCLEAN)
phasecenter, midtime = hf.calc_phasecenter_from_solxy(ms_in, timerange=trange, xycen=xycen,
                                                      usemsphacenter=False)


# ############################ Step 2: Split a short period for self-calibration #################
msinfo = hf.read_msinfo(ms_in, verbose=True)
# output self-calibrated visibility
ms_slfcaled = workdir + os.path.basename(ms_in).replace('corrected', 'slfcaled')

# The following does the split
slfcalms = slfcaldir + 'slfcalms.slfcal'
slfcaledms = slfcaldir + 'slfcalms.slfcaled'
if not os.path.exists(slfcalms):
    print('splitting a small piece for ', trange)
    split(vis=ms_in, outputvis=slfcalms, datacolumn='data', timerange=trange, correlation=stokes, antenna=antennas)


# ################################# Step 3: Main self-calibration loop #############################
# More control parameters for self-calibration
# define ranges of spectral windows you wish to draw masks. Need to change on case-by-case basis
spwran_masks = ['1~5', '6~12', '13~20', '21~30']
half_wid = 2  # how many nearby spws (+- to the spw of choice) to generate model image for the first round of selfcal

# Optional. No need to change in usual cases
maxnround = 10  # maximum rounds of self-calibration allowed
npix = 512
uvrange = '>500lambda'
robust = 1.0
sbeam_1GHz = 80.
refantenna = '0'
sclfactor = 1.  # this applies to VLA images only (with 20 dB attenuators)

# Optional. Consider only this range of spws for imaging/generating models
spw_min = 1  # first spw for consideration
spw_max = 30  # last spw for consideration

os.system('rm -rf ' + imagedir + '*')
os.system('rm -rf ' + caltbdir + '*')
# first step: make a mock caltable for the entire database
print('Processing ' + trange)
msmd.open(slfcalms)
nspw = msmd.nspw()
cfreqsghz = [msmd.meanfreq(s, unit='GHz') for s in range(nspw)]
msmd.done()
spws_ms = [str(s) for s in np.arange(nspw)]

calprefix = caltbdir + 'slf'
imgprefix = imagedir + 'slf'
# starting beam size at 1.0 GHz in arcsec
strtmp = [m.replace(':', '') for m in trange.split('~')]
timestr = 't' + strtmp[0] + '-' + strtmp[1]

# make initial images without self-calibration
n = 0
print('==========Doing initial imaging==========')
clearcal(slfcalms)
prompt3 = input('Wish to draw masks? ')
maskfile = slfcaldir + 'masks.p'
if prompt3.lower() == 'y':
    maskfile = draw_masks(slfcalms, spwrans_mask=spwran_masks, maskdir=maskdir,
                          maskfile=maskfile, restoringbeam=['20arcsec'], timerange=trange,
                          phasecenter=phasecenter, npix=npix, gain=0.05, uvrange=uvrange)
    masks = pickle.load(open(maskfile, 'rb'))
if prompt3.lower() == 'n':
    print('Chose not to draw new masks. Looking for existing ones.')
    if os.path.exists(slfcaldir + 'masks.p'):
        masks = pickle.load(open(slfcaldir + 'masks.p', 'rb'))
    if not os.path.exists(slfcaldir + 'masks.p'):
        print('masks do not exist. Use default mask')
        masks = []

figname, fitsfiles = make_images(slfcalms, msinfo=msinfo, spws=spws_slf, spwrans_mask=spwran_masks,
                                 usemask=True, maskfile=maskfile, imgprefix=imgprefix, round=n,
                                 trange=trange, xycen=xycen, fov=fov,
                                 uvrange=uvrange, antennas='', npix=npix, cell=['1arcsec'],
                                 niter=500, gain=0.05, robust=1.0, sbeam_1GHz=sbeam_1GHz, nrow=3,
                                 usemsphacenter=False, observatory='EOVSA',sclfactor=sclfactor)
print('Resulting image is {}'.format(figname))
slftbs = []
caltypes = []
niters = []

# Question on whether phase or amplitude selfcalibration
prompt2 = input('Phase only (p), amplitude only (a), or amplitude & phase (ap)? ')
while prompt2.lower() != 'p' and prompt2.lower() != 'a' and prompt2.lower() != 'ap':
    print('Your input is {}. Your input has to be "p", "a", or "ap"'.format(prompt2))
    prompt2 = input('Phase only (p), amplitude only (a), or amplitude & phase (ap)? ')
if prompt2.lower() == 'p':
    caltype = 'p'
    print('Doing phase-only calibration')
if prompt2.lower() == 'a':
    caltype = 'a'
    print('Doing amplitude-only calibration')
if prompt2.lower() == 'ap':
    caltype = 'ap'
    print('Doing amplitude & phase calibration')
caltypes.append(caltype)
# Question on how many iterations
while True:
    try:
        niter = int(input('How many iterations (e.g., 100)? '))
    except ValueError:
        print('Your input has to be an integer number. Try again.')
    else:
        break
print('Using {0:d} iterations for this self-calibration round'.format(niter))
niters.append(niter)

while n < maxnround:
    # Name of the self-calibration gain table, round is 1 based
    slfcal_tb_g = calprefix + '.G' + str(n + 1)

    solv_success = False
    for s, sp in enumerate(spws_slf):
        print('Self-calibrating spw: ' + sp)
        cfreqghz = cfreqsghz[int(sp)]
        # setting restoring beam size (not very useful for selfcal anyway, but just to see the results)
        bm = max(sbeam_1GHz / cfreqghz, 4.)
        print('Beam size {0:.1f}" at ref freq {1:.2f} GHz'.format(bm, cfreqghz))
        slfcal_img = imgprefix + '.spw' + sp.zfill(2) + '.slfcal' + str(n) + '.tmp'
        # only the first and second round uses nearby spws for getting initial model
        if n < 1:
            if not spw_min:
                spw_min = int(spws_ms[0])
            if not spw_max:
                spw_min = int(spws_ms[-1])
            spbg = int(max(int(sp) - half_wid, spw_min))
            sped = int(min(int(sp) + half_wid, spw_max))
            spwran = str(spbg) + '~' + str(sped)
            print('using spw {0:s} as model'.format(spwran))
        else:
            spwran = sp
            print('using spw {0:s} as model'.format(spwran))
        foundmask = False
        if len(masks) > 0:
            for m, spwran_mask in enumerate(spwran_masks):
                sbg, sed = spwran_mask.split('~')
                if (int(sp) >= int(sbg)) and (int(sp) <= int(sed)):
                    mask = masks[m]
                    print('using mask {0:s}'.format(os.path.basename(mask)))
                    foundmask = True
            if not foundmask:
                print('mask not defined for {}. Do use any masks.'.format(sp))
                mask = ''
        else:
            mask = ''
        try:
            tclean(vis=slfcalms, antenna=antennas, imagename=slfcal_img, uvrange=uvrange, spw=spwran,
                   specmode='mfs', timerange=trange, imsize=[npix], cell=['1arcsec'], niter=niter,
                   gain=0.05, stokes=stokes, weighting='briggs', robust=robust, phasecenter=phasecenter, mask=mask,
                   restoringbeam=[str(bm) + 'arcsec'], pbcor=False, interactive=False, savemodel='modelcolumn')

            clnjunks = ['.mask', '.flux', '.model', '.psf', '.residual', '.image', '.pb', '.image.pbcor', '.sumwt']
            for clnjunk in clnjunks:
                if os.path.exists(slfcal_img + clnjunk):
                    os.system('rm -rf ' + slfcal_img + clnjunk)
        except:
            print('error in cleaning spw: ' + sp)
            print('using nearby spws for initial model')
            sp_e = min(int(sp) + 2, int(spws_ms[-1]))
            sp_i = max(int(sp) - 2, int(spws_ms[0]))
            sp_ = str(sp_i) + '~' + str(sp_e)
            try:
                tclean(vis=slfcalms, antenna=antennas, imagename=slfcal_img, uvrange=uvrange, spw=sp_,
                       specmode='mfs', timerange=trange, imsize=[npix], cell=['1arcsec'], niter=niter,
                       gain=0.05, stokes=stokes, weighting='briggs', robust=robust, phasecenter=phasecenter,
                       mask=mask, restoringbeam=[str(bm) + 'arcsec'],
                       pbcor=False, interactive=False, savemodel='modelcolumn')
                print('using spw {0:s} as model'.format(sp_))
            except:
                print('still not successful. abort...')
                solv_success = False
                continue

        # Check what the calibration type is. If it involves amplitude, need to use the same uvrange at used in
        # generating the model image.
        if caltype == 'p':
            uvrange_for_gaincal = ''
        else:
            uvrange_for_gaincal = uvrange

        if s == 0 or not solv_success:
            append = False
        else:
            append = True
        try:
            # Tested (in CASA 6.5) that flags in gaintables being applied on the fly do not affect the solution
            # (equivalent to applymode='calonly')
            gaincal(vis=slfcalms, refant=refantenna, antenna=antennas, caltable=slfcal_tb_g, spw=sp,
                    uvrange=uvrange_for_gaincal, gaintable=slftbs, selectdata=True, timerange=trange,
                    solint='inf', gaintype='G', calmode=caltype, combine='', minblperant=3, minsnr=2,
                    append=append)
        except:
            print('No solution found in spw: ' + sp)
            solv_success = False
            continue

        if os.path.exists(slfcal_tb_g):
            solv_success = True
        else:
            solv_success = False

    if os.path.exists(slfcal_tb_g):
        slftbs.append(slfcal_tb_g)
        clearcal(slfcalms)
        delmod(slfcalms)
        applycal(vis=slfcalms, gaintable=slftbs, spw=','.join(spws_slf), selectdata=True,
                 antenna=antennas, interp='nearest', flagbackup=False, applymode='calonly', calwt=False)
        res = plt_slftable(slfcalms, slfcal_tb_g, caltype, docombine=False)

    if solv_success:
        # make images after this self-calibration round
        figname, fitsfiles = make_images(slfcalms, msinfo=msinfo, spws=spws_slf, spwrans_mask=spwran_masks,
                                         usemask=False, imgprefix=imgprefix, round=n + 1,
                                         trange=trange, xycen=xycen, fov=fov, uvrange=uvrange, antennas='',
                                         npix=npix,
                                         cell=['1arcsec'], niter=500, gain=0.1, robust=1.0, sbeam_1GHz=sbeam_1GHz,
                                         nrow=3, usemsphacenter=False, observatory='EOVSA', sclfactor=sclfactor)
        print('Resulting image is {}'.format(figname))
        print('Check the resulting image {0} and calibration table plot in {1}.'.format(figname,
                                                                                        os.path.dirname(
                                                                                            slfcal_tb_g)))

    prompt = input('Continue selfcal? (y or n) ')
    while prompt.lower() != 'y' and prompt.lower() != 'n':
        print('Your input is {}. Your input has to be "y" or "n"'.format(prompt))
        prompt = input('Continue selfcal? (y or n) ')
    if prompt.lower() == 'n':
        if os.path.exists(slfcaledms):
            os.system('rm -rf ' + slfcaledms)
        split(slfcalms, slfcaledms, datacolumn='corrected')
        print('Final self-calibrated ms is {0:s}'.format(slfcaledms))
        success = plt_slftable(slfcalms, slftbs, caltypes, docombine=True)
        outfits = '{0:s}/slfcaled_image_comb.fits'.format(imagedir_slfcaled)
        print('Final self-calibrated image is {}'.format(outfits))
        if len(fitsfiles) > 1:
            fitsutils.fits_wrap_spwX(fitsfiles, outfits)
        break
    if prompt.lower() == 'y':
        n = n + 1
        if n < maxnround:
            print('Continue to the No. {0:d} round'.format(n + 1))
            # Question on whether phase or amplitude selfcalibration
            prompt2 = input('Phase only (p), amplitude only (a), or amplitude & phase (ap)? ')
            while prompt2.lower() != 'p' and prompt2.lower() != 'a' and prompt2.lower() != 'ap':
                print('Your input is {}. Your input has to be "p", "a", or "ap"'.format(prompt2))
                prompt2 = input('Phase only (p), amplitude only (a), or amplitude & phase (ap)? ')
            if prompt2.lower() == 'p':
                caltype = 'p'
                print('Doing phase-only calibration')
            if prompt2.lower() == 'a':
                caltype = 'a'
                print('Doing amplitude-only calibration')
            if prompt2.lower() == 'ap':
                caltype = 'ap'
                print('Doing amplitude & phase calibration')
            caltypes.append(caltype)
            # Question on making new masks
            prompt3 = input('Wish to update masks? (y or n) ')
            while prompt3.lower() != 'y' and prompt3.lower() != 'n':
                print('Your input is {}. Your input has to be "y" or "n"'.format(prompt3))
                prompt3 = input('Wish to update masks? (y or n) ')
            if prompt3.lower() == 'y':
                maskfile = draw_masks(slfcalms, spwrans_mask=spwran_masks, maskdir=maskdir,
                                      maskfile=maskfile, restoringbeam=['20arcsec'], timerange=trange,
                                      phasecenter=phasecenter, npix=npix, gain=0.05, uvrange=uvrange)
                masks = pickle.load(open(maskfile, 'rb'))
            else:
                print('Chose not to draw new masks. Looking for existing ones.')
                if os.path.exists(slfcaldir + 'masks.p'):
                    masks = pickle.load(open(slfcaldir + 'masks.p', 'rb'))
                if not os.path.exists(slfcaldir + 'masks.p'):
                    print('masks do not exist. Use default mask')
                    masks = []
            # Question on how many iterations
            while True:
                try:
                    niter = int(input('How many iterations (e.g., 100)? '))
                except ValueError:
                    print('Your input has to be an integer number. Try again.')
                else:
                    break
            print('Using {0:d} iterations for this self-calibration round'.format(niter))
            niters.append(niter)
        else:
            print('Hitting the max rounds of self-calibration of {0:d}'.format(n))
            print('Splitting final calibrated ms is {0:s}'.format(slfcaledms))
            if os.path.exists(slfcaledms):
                os.system('rm -rf ' + slfcaledms)
            split(slfcalms, slfcaledms, datacolumn='corrected')
            success = plt_slftable(slfcalms, slftbs, caltypes, docombine=True)
            outfits = '{0:s}/slfcaled_image_comb.fits'.format(imagedir_slfcaled)
            print('Final self-calibrated image is {}'.format(outfits))
            if len(fitsfiles) > 1:
                fitsutils.fits_wrap_spwX(fitsfiles, outfits)

# ####################### Step 4: Apply self-calibration tables to the original visibility dataset ####################
prompt4 = input('Do you want to apply the self-calibration tables to the original ms? (y or n) ')
while prompt4.lower() != 'y' and prompt4.lower() != 'n':
    print('Your input is {}. Your input has to be "y" or "n"'.format(prompt4))
    prompt4 = input('Do you want to apply the self-calibration tables to the original ms? (y or n) ')
if prompt4.lower() == 'y':
    clearcal(ms_in)
    applycal(vis=ms_in, gaintable=slftbs, spw=','.join(spws_slf), selectdata=True,
             antenna=antennas, interp='nearest', flagbackup=False, applymode='calonly', calwt=False)
    if os.path.exists(ms_slfcaled):
        print('Found an existing self-calibrated dataset.')
        prompt5 = input('Wish to overwrite this dataset? (y or n) ')
        while prompt5.lower() != 'y' and prompt5.lower() != 'n':
            print('Your input is {}. Your input has to be "y" or "n"'.format(prompt5))
            prompt5 = input('Wish to overwrite this dataset? (y or n) ')
        if prompt5.lower() == 'y':
            os.system('rm -rf ' + ms_slfcaled)
            ms_out = ms_slfcaled
            print('Splitting self-calibrated dataset as {}'.format(ms_out))
        else:
            ms_out = ms_slfcaled + '.copy'
    else:
        ms_out = ms_slfcale
    print('Splitting self-calibrated dataset as {}'.format(ms_out))
    split(ms_in, ms_out, spw=','.join(spws_slf), correlation=stokes, antenna=antennas, datacolumn='corrected')
else:
    print('Chose not to apply self-calibration tables. You can apply them later using applycal().')
    print('Self-calibration tables are:', slftbs)


