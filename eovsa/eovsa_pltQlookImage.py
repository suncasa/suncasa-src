import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sunpy import map as smap
import sunpy.cm.cm as cm_smap
from suncasa.utils import plot_mapX as pmX
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar as colorbar
import matplotlib.patches as patches
from datetime import timedelta
from datetime import datetime
from glob import glob
import numpy as np
from astropy.time import Time
import urllib.request
import socket

socket.setdefaulttimeout(180)

imgfitsdir = '/data1/eovsa/fits/synoptic/'
imgfitstmpdir = '/data1/workdir/fitstmp/'
pltfigdir = '/common/webplots/SynopticImg/eovsamedia/eovsa-browser/'


def clearImage():
    for (dirpath, dirnames, filenames) in os.walk(pltfigdir):
        for filename in filenames:
            for k in ['0094', '0193', '0335', '4500', '0171', '0304', '0131', '1700', '0211', '1600', '_HMIcont',
                      '_HMImag']:
                # for k in ['_Halph_fr']:
                if k in filename:
                    print(os.path.join(dirpath, filename))
                    os.system('rm -rf ' + os.path.join(dirpath, filename))


def pltEmptyImage2(dpis_dict={'t': 32.0}):
    imgoutdir = './nodata/'

    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    rect_bkg = patches.Rectangle((-1227, -1227), 1227 * 2, 1227 * 2, linewidth=0, edgecolor='none', facecolor='k',
                                 alpha=0.9)
    rect_bar = patches.Rectangle((-1227, -300), 1227 * 2, 300 * 2, linewidth=0, edgecolor='none', facecolor='w',
                                 alpha=0.5)

    ax.add_patch(rect_bkg)
    ax.add_patch(rect_bar)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.text(0.5, 0.5, 'No Data',
            transform=ax.transAxes, color='w', ha='center', va='center', fontsize=120)
    ax.set_xlim(-1227, 1227)
    ax.set_ylim(-1227, 1227)

    for l, dpi in dpis_dict.items():
        figname = 'nodata.jpg'
        fig.savefig(figname, dpi=np.int(dpi), quality=85)
    return


def pltEmptyImage(datestr, spws, vmaxs, vmins, dpis_dict={'t': 32.0}):
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datestrdir
    imgoutdir = './nodata/'

    cmap = cm_smap.get_cmap('sdoaia304')

    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    rect = patches.Rectangle((-1227, -300), 1227 * 2, 300 * 2, linewidth=0, edgecolor='none', facecolor='k', alpha=0.5)

    for s, sp in enumerate(spws):
        ax.cla()
        spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
        eofile = imgindir + 'eovsa_{}.spw{}.tb.disk.fits'.format(dateobj.strftime('%Y%m%d'), spwstr)
        if not os.path.exists(eofile): continue
        if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
        eomap = smap.Map(eofile)
        norm = colors.Normalize(vmin=vmins[s], vmax=vmaxs[s])
        eomap_ = pmX.Sunmap(eomap)
        eomap_.imshow(axes=ax, cmap=cmap, norm=norm, alpha=0.75)
        eomap_.draw_limb(axes=ax, lw=0.5, alpha=0.5)
        eomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(0.5, 0.5, 'No Data',
                transform=ax.transAxes, color='w', ha='center', va='center', fontsize=120)
        ax.add_patch(rect)
        ax.set_xlim(-1227, 1227)
        ax.set_ylim(-1227, 1227)

        for l, dpi in dpis_dict.items():
            figname = os.path.join(imgoutdir, '{}_eovsa_bd{:02d}.jpg'.format(l, s + 1))
            fig.savefig(figname, dpi=np.int(dpi), quality=85)
    return


def pltEovsaQlookImage(datestr, spws, vmaxs, vmins, dpis_dict, fig=None, ax=None, overwrite=False, verbose=False):
    from astropy.visualization.stretch import AsinhStretch
    from astropy.visualization import ImageNormalize
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datestrdir
    imgoutdir = pltfigdir + datestrdir

    cmap = cm_smap.get_cmap('sdoaia304')
    cmap.set_bad(color='k')

    if fig is None or ax is None:
        mkfig = True
    else:
        mkfig = False

    if mkfig:
        fig, ax = plt.subplots(figsize=(8, 8))
        fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    if verbose: print('Processing EOVSA images for date {}'.format(dateobj.strftime('%Y-%m-%d')))
    for s, sp in enumerate(spws):
        fexists = []
        for l, dpi in dpis_dict.items():
            figname = os.path.join(imgoutdir, '{}_eovsa_bd{:02d}.jpg'.format(l, s + 1))
            fexists.append(os.path.exists(figname))

        if overwrite or (False in fexists):
            ax.cla()
            spwstr = '-'.join(['{:02d}'.format(int(sp_)) for sp_ in sp.split('~')])
            eofile = imgindir + 'eovsa_{}.spw{}.tb.disk.fits'.format(dateobj.strftime('%Y%m%d'), spwstr)
            if not os.path.exists(eofile): continue
            if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
            try:
                eomap = smap.Map(eofile)
                stretch = AsinhStretch(a=0.15)
                norm = ImageNormalize(vmin=vmins[s], vmax=vmaxs[s], stretch=stretch)
                # norm = colors.Normalize(vmin=vmins[s], vmax=vmaxs[s])
                eomap_ = pmX.Sunmap(eomap)
                eomap_.imshow(axes=ax, cmap=cmap, norm=norm)
                eomap_.draw_limb(axes=ax, lw=0.5, alpha=0.5)
                eomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.5)
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.text(0.02, 0.02,
                        'EOVSA {:.1f} GHz  {}'.format(eomap.meta['CRVAL3'] / 1e9, eomap.date.strftime('%d-%b-%Y 20:00 UT')),
                        transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=9)
                ax.text(0.98, 0.02, 'Max Tb {:.0f} K'.format(np.nanmax(eomap.data)),
                        transform=ax.transAxes, color='w', ha='right', va='bottom', fontsize=9)
                ax.set_xlim(-1227, 1227)
                ax.set_ylim(-1227, 1227)

                for l, dpi in dpis_dict.items():
                    figname = os.path.join(imgoutdir, '{}_eovsa_bd{:02d}.jpg'.format(l, s + 1))
                    fig.savefig(figname, dpi=np.int(dpi), quality=85)
            except Exception as err:
                print('Fail to plot {}'.format(eofile))
                print(err)
    if mkfig:
        pass
    else:
        plt.close(fig)
    return


def pltSdoQlookImage(datestr, dpis_dict, fig=None, ax=None, overwrite=False, verbose=False, clearcache=False):
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = os.path.join(imgfitstmpdir, datestr)
    imgoutdir = pltfigdir + datestrdir
    if not os.path.exists(imgindir):
        os.makedirs(imgindir)

    aiaDataSource = {"0094": 8,
                     "0193": 11,
                     "0335": 14,
                     # "4500": 17,
                     "0171": 10,
                     "0304": 13,
                     "0131": 9,
                     "1700": 16,
                     "0211": 12,
                     # "1600": 15,
                     "_HMIcont": 18,
                     "_HMImag": 19}

    if fig is None or ax is None:
        mkfig = True
    else:
        mkfig = False

    if mkfig:
        fig, ax = plt.subplots(figsize=(8, 8))
        fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    if verbose: print('Processing SDO images for date {}'.format(dateobj.strftime('%Y-%m-%d')))
    for key, sourceid in aiaDataSource.items():
        fexists = []
        for l, dpi in dpis_dict.items():
            figname = os.path.join(imgoutdir, '{}{}.jpg'.format(l, key))
            fexists.append(os.path.exists(figname))

        if overwrite or (False in fexists):
            sdourl = 'https://api.helioviewer.org/v2/getJP2Image/?date={}T20:00:00Z&sourceId={}'.format(datestr,
                                                                                                        sourceid)
            sdofile = os.path.join(imgindir, key + '.jp2')
            if overwrite and os.path.exists(sdofile):
                os.system('rm -rf {}'.format(sdofile))
            if not os.path.exists(sdofile):
                try:
                    urllib.request.urlretrieve(sdourl, sdofile)
                except:
                    print('The connection with {} has timed out. Skipped!'.format(sdourl))
            ax.cla()

            if not os.path.exists(sdofile): continue
            if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
            try:
                sdomap = smap.Map(sdofile)
                norm = colors.Normalize()
                sdomap_ = pmX.Sunmap(sdomap)
                if "HMI" in key:
                    cmap = plt.get_cmap('gray')
                else:
                    cmap = cm_smap.get_cmap('sdoaia' + key.lstrip('0'))
                sdomap_.imshow(axes=ax, cmap=cmap, norm=norm)
                sdomap_.draw_limb(axes=ax, lw=0.5, alpha=0.5)
                sdomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.5)
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.text(0.02, 0.02,
                        '{}/{} {}  {}'.format(sdomap.observatory, sdomap.instrument.split(' ')[0], sdomap.measurement,
                                              sdomap.date.strftime('%d-%b-%Y %H:%M UT')),
                        transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=9)
                ax.set_xlim(-1227, 1227)
                ax.set_ylim(-1227, 1227)

                for l, dpi in dpis_dict.items():
                    figname = os.path.join(imgoutdir, '{}{}.jpg'.format(l, key))
                    fig.savefig(figname, dpi=np.int(dpi), quality=85)
            except Exception as err:
                print('Fail to plot {}'.format(sdofile))
                print(err)
    if clearcache:
        os.system('rm -rf ' + imgindir)

    if mkfig:
        pass
    else:
        plt.close(fig)
    return


def pltBbsoQlookImage(datestr, dpis_dict, fig=None, ax=None, overwrite=False, verbose=False, clearcache=False):
    from astropy.io import fits
    from html.parser import HTMLParser
    class MyHTMLParser(HTMLParser):
        def __init__(self, prefix='bbso_halph_fr_', suffix='.fts'):
            HTMLParser.__init__(self)
            self.prefix = prefix
            self.suffix = suffix

        def handle_starttag(self, tag, attrs):
            if tag != 'a':
                return
            for name, value in attrs:
                if name == "href":
                    if value.startswith(self.prefix) and value.endswith(self.suffix):
                        self.links.append(value)

    def extract(url, prefix='bbso_halph_fr_', suffix='.fts'):
        import urllib.request
        with urllib.request.urlopen(url) as response:
            f = response.read()

        parser = MyHTMLParser(prefix, suffix)
        parser.links = []
        parser.feed(str(f))
        return parser.links

    bbsodir = 'http://www.bbso.njit.edu/pub/archive/'
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = os.path.join(imgfitstmpdir, datestr)
    imgoutdir = pltfigdir + datestrdir
    if not os.path.exists(imgindir):
        os.makedirs(imgindir)

    bbsoDataSource = {"_Halph_fr": ["bbso_halph_fr_", ".fts"]}

    if fig is None or ax is None:
        mkfig = True
    else:
        mkfig = False

    if mkfig:
        fig, ax = plt.subplots(figsize=(8, 8))
        fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    if verbose: print('Processing BBSO images for date {}'.format(dateobj.strftime('%Y-%m-%d')))
    for key, sourceid in bbsoDataSource.items():
        fexists = []
        for l, dpi in dpis_dict.items():
            figname = os.path.join(imgoutdir, '{}{}.jpg'.format(l, key))
            fexists.append(os.path.exists(figname))

        if overwrite or (False in fexists):
            bbsosite = os.path.join(bbsodir, datestrdir)
            filelist = extract(bbsosite, sourceid[0], sourceid[1])
            if filelist:
                tfilelist = Time(
                    [datetime.strptime(tf.replace(sourceid[0], '').replace(sourceid[1], ''), "%Y%m%d_%H%M%S") for tf in
                     filelist])
                bbsourl = os.path.join(bbsosite, filelist[
                    np.nanargmin(np.abs(np.array(tfilelist.mjd - (Time(dateobj).mjd + 20. / 24.))))])

                bbsofile = os.path.join(imgindir, key + '.fits')
                if not os.path.exists(bbsofile):
                    try:
                        urllib.request.urlretrieve(bbsourl, bbsofile)
                    except:
                        print('The connection with {} has timed out. Skipped!'.format(bbsourl))
                ax.cla()
                if not os.path.exists(bbsofile): continue
                if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
                try:
                    hdu = fits.open(bbsofile)[0]
                    header = hdu.header
                    header['WAVELNTH'] = 6562.8
                    header['WAVEUNIT'] = 'angstrom'
                    header['WAVE_STR'] = 'Halph'
                    header['CTYPE1'] = 'HPLN-TAN'
                    header['CUNIT1'] = 'arcsec'
                    header['CTYPE2'] = 'HPLT-TAN'
                    header['CUNIT2'] = 'arcsec'
                    header['DATE-OBS'] = header['DATE_OBS']
                    for k in ['CONTRAST', 'WAVE ERR']:
                        try:
                            header.remove(k)
                        except:
                            pass

                    bbsomap = smap.Map(hdu.data, header)
                    med = np.nanmean(bbsomap.data)
                    norm = colors.Normalize(vmin=med - 1500, vmax=med + 1500)
                    bbsomap_ = pmX.Sunmap(bbsomap)
                    cmap = cm_smap.get_cmap('sdoaia304')
                    bbsomap_.imshow(axes=ax, cmap=cmap, norm=norm)
                    bbsomap_.draw_limb(axes=ax, lw=0.5, alpha=0.5)
                    bbsomap_.draw_grid(axes=ax, grid_spacing=10. * u.deg, lw=0.5)
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                    ax.text(0.02, 0.02,
                            '{}  {}'.format(bbsomap.instrument, bbsomap.date.strftime('%d-%b-%Y %H:%M UT')),
                            transform=ax.transAxes, color='w', ha='left', va='bottom', fontsize=9)
                    ax.set_xlim(-1227, 1227)
                    ax.set_ylim(-1227, 1227)
                    ax.set_facecolor('k')

                    for l, dpi in dpis_dict.items():
                        figname = os.path.join(imgoutdir, '{}{}.jpg'.format(l, key))
                        fig.savefig(figname, dpi=np.int(dpi), quality=85)
                except Exception as err:
                    print('Fail to plot {}'.format(bbsofile))
                    print(err)
    if clearcache:
        os.system('rm -rf ' + imgindir)

    if mkfig:
        pass
    else:
        plt.close(fig)
    return


def main(year=None, month=None, day=None, ndays=1, clearcache=False, ovwrite_eovsa=False, ovwrite_sdo=False,
         ovwrite_bbso=False, show_warning=False):
    # tst = datetime.strptime("2017-04-01", "%Y-%m-%d")
    # ted = datetime.strptime("2019-12-31", "%Y-%m-%d")
    if not show_warning:
        import warnings
        warnings.filterwarnings("ignore")
    if year:
        ted = datetime(year, month, day)
    else:
        ted = datetime.now() - timedelta(days=2)
    tst = Time(np.fix(Time(ted).mjd) - ndays, format='mjd').datetime
    tsep = datetime.strptime('2019-02-22', "%Y-%m-%d")

    # vmaxs = [22.0e4, 8.0e4, 5.4e4, 3.5e4, 2.3e4, 1.8e4, 1.5e4]
    # vmins = [-9.0e3, -5.5e3, -3.4e3, -2.5e3, -2.5e3, -2.5e3, -2.5e3]
    vmaxs = [70.0e4, 30e4, 18e4, 13e4, 8e4, 6e4, 6e4]
    vmins = [-18.0e3, -8e3, -4.8e3, -3.4e3, -2.1e3, -1.6e3, -1.6e3]

    dpis = np.array([256, 512, 1024]) / 8
    dpis_dict_eo = {'t': dpis[0], 'l': dpis[1], 'f': dpis[2]}
    dpis = np.array([256, 512, 2048]) / 8
    dpis_dict_sdo = {'t': dpis[0], 'l': dpis[1], 'f': dpis[2]}
    dpis = np.array([256, 512, 2048]) / 8
    dpis_dict_bbso = {'t': dpis[0], 'l': dpis[1], 'f': dpis[2]}

    plt.ioff()
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)

    dateobs = tst
    while dateobs < ted:
        dateobs = dateobs + timedelta(days=1)

        # dateobs = datetime.strptime("2019-12-21", "%Y-%m-%d")

        if dateobs > tsep:
            spws = ['0~1', '2~5', '6~10', '11~20', '21~30', '31~43', '44~49']
        else:
            spws = ['1~3', '4~9', '10~16', '17~24', '25~30']

        datestr = dateobs.strftime("%Y-%m-%d")
        pltEovsaQlookImage(datestr, spws, vmaxs, vmins, dpis_dict_eo, fig, ax, overwrite=ovwrite_eovsa, verbose=True)
        pltSdoQlookImage(datestr, dpis_dict_sdo, fig, ax, overwrite=ovwrite_sdo, verbose=True, clearcache=clearcache)
        pltBbsoQlookImage(datestr, dpis_dict_bbso, fig, ax, overwrite=ovwrite_bbso, verbose=True,
                          clearcache=clearcache)


if __name__ == '__main__':
    '''
    Name: 
    eovsa_pltQlookImage --- pipeline for plotting EOVSA daily full-disk images at multi frequencies.

    Synopsis:
    eovsa_pltQlookImage.py [options]... [DATE_IN_YY_MM_DD]

    Description:
    Plot EOVSA daily full-disk images at multi frequencies of the date specified
    by DATE_IN_YY_MM_DD (or from ndays before the DATE_IN_YY_MM_DD if option --ndays/-n is provided).
    If DATE_IN_YY_MM_DD is omitted, it will be set to 2 days before now by default. 
    The are no mandatory arguments in this command.

    -c, --clearcache
            Remove temporary files

    -n, --ndays
            Processing the date spanning from DATE_IN_YY_MM_DD-ndays to DATE_IN_YY_MM_DD. Default is 30

    -e, --ovwrite_eovsa
            If True, overwrite eovsa images.
            Syntax: True, False, T, F, 1, 0

    -s, --ovwrite_sdo
            If True, overwrite sdo images.
            Syntax: True, False, T, F, 1, 0
            
    -b, --ovwrite_bbso
            If True, overwrite bbso images.
            Syntax: True, False, T, F, 1, 0                                


    Example: 
    eovsa_pltQlookImage.py -c True -n 20 -e True -s False -b False 2020 06 10
    '''

    import sys
    import numpy as np
    import getopt

    year = None
    month = None
    day = None
    ndays = 1
    clearcache = True
    ovwrite_eovsa = False
    ovwrite_sdo = False
    ovwrite_bbso = False
    show_warning = False
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(argv, "c:n:e:s:b:w:",
                                   ['clearcache=', 'ndays=', 'ovwrite_eovsa=', 'ovwrite_sdo=', 'ovwrite_bbso=',
                                    'show_warning='])
        print(opts, args)
        for opt, arg in opts:
            if opt in ['-c', '--clearcache']:
                if arg in ['True', 'T', '1']:
                    clearcache = True
                elif arg in ['False', 'F', '0']:
                    clearcache = False
                else:
                    clearcache = np.bool(arg)
            elif opt in ('-n', '--ndays'):
                ndays = np.int(arg)
            elif opt in ('-e', '--ovwrite_eovsa'):
                if arg in ['True', 'T', '1']:
                    ovwrite_eovsa = True
                elif arg in ['False', 'F', '0']:
                    ovwrite_eovsa = False
                else:
                    ovwrite_eovsa = np.bool(arg)
            elif opt in ('-s', '--ovwrite_sdo'):
                if arg in ['True', 'T', '1']:
                    ovwrite_sdo = True
                elif arg in ['False', 'F', '0']:
                    ovwrite_sdo = False
                else:
                    ovwrite_sdo = np.bool(arg)
            elif opt in ('-s', '--ovwrite_bbso'):
                if arg in ['True', 'T', '1']:
                    ovwrite_bbso = True
                elif arg in ['False', 'F', '0']:
                    ovwrite_bbso = False
                else:
                    ovwrite_bbso = np.bool(arg)
            elif opt in ('-w', '--show_warning'):
                if arg in ['True', 'T', '1']:
                    show_warning = True
                elif arg in ['False', 'F', '0']:
                    show_warning = False
                else:
                    show_warning = np.bool(arg)
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
        ovwrite_eovsa = False
        ovwrite_sdo = False
        ovwrite_bbso = False
        show_warning = False

    print("Running pipeline_plt for date {}-{}-{}.".format(year, month, day))
    kargs = {'ndays': ndays,
             'clearcache': clearcache,
             'ovwrite_eovsa': ovwrite_eovsa,
             'ovwrite_sdo': ovwrite_sdo,
             'ovwrite_bbso': ovwrite_bbso,
             'show_warning': show_warning}
    for k, v in kargs.items():
        print(k, v)

    main(year=year, month=month, day=day, ndays=ndays,
         clearcache=clearcache,
         ovwrite_eovsa=ovwrite_eovsa,
         ovwrite_sdo=ovwrite_sdo,
         ovwrite_bbso=ovwrite_bbso,
         show_warning=show_warning)
