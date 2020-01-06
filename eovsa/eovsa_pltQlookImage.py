import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sunpy import map as smap
import sunpy.cm.cm as cm
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

imgfitsdir = '/data1/eovsa/fits/qlook_10m/'
imgfitstmpdir = '/data1/workdir/fitstmp/'
pltfigdir = '/common/webplots/SynopticImg/eovsamedia/eovsa-browser/'


def clearImage():
    for (dirpath, dirnames, filenames) in os.walk(pltfigdir):
        for filename in filenames:
            if '_eovsa_bd01' in filename:
                print(os.path.join(dirpath, filename))
                # os.system('rm -rf '+os.path.join(dirpath,filename))


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
    datastrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datastrdir
    imgoutdir = './nodata/'

    cmap = cm.get_cmap('sdoaia304')

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
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datastrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datastrdir
    imgoutdir = pltfigdir + datastrdir

    cmap = cm.get_cmap('sdoaia304')

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
            eomap = smap.Map(eofile)
            norm = colors.Normalize(vmin=vmins[s], vmax=vmaxs[s])
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
    if mkfig:
        pass
    else:
        plt.close(fig)
    return


def pltSdoQlookImage(datestr, dpis_dict, fig=None, ax=None, overwrite=False, verbose=False, clearcache=False):
    plt.ioff()
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datastrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = os.path.join(imgfitstmpdir, datestr)
    imgoutdir = pltfigdir + datastrdir
    if not os.path.exists(imgindir):
        os.makedirs(imgindir)

    aiaDataSource = {"0094": 8,
                     "0193": 11,
                     "0335": 14,
                     "4500": 17,
                     "0171": 10,
                     "0304": 13,
                     "0131": 9,
                     "1700": 16,
                     "0211": 12,
                     "1600": 15,
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
            if not os.path.exists(sdofile):
                urllib.request.urlretrieve(sdourl, sdofile)
            ax.cla()

            if not os.path.exists(sdofile): continue
            if not os.path.exists(imgoutdir): os.makedirs(imgoutdir)
            sdomap = smap.Map(sdofile)
            norm = colors.Normalize()
            sdomap_ = pmX.Sunmap(sdomap)
            if "HMI" in key:
                cmap = cm.get_cmap('gray')
            else:
                cmap = cm.get_cmap('sdoaia' + key.lstrip('0'))
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
    if clearcache:
        os.rmdir(imgindir)

    if mkfig:
        pass
    else:
        plt.close(fig)
    return


def main(year=None, month=None, day=None, clearcache=False):
    # tst = datetime.strptime("2017-04-01", "%Y-%m-%d")
    # ted = datetime.strptime("2019-12-31", "%Y-%m-%d")
    if year:
        ted = datetime(year, month, day)
    else:
        ted = datetime.now()
    tst = Time(np.fix(Time(ted).mjd) - 30, format='mjd').datetime
    tsep = datetime.strptime('2019-02-22', "%Y-%m-%d")

    vmaxs = [22.0e4, 8.0e4, 5.4e4, 3.5e4, 2.3e4, 1.8e4, 1.5e4]
    vmins = [-9.0e3, -5.5e3, -3.4e3, -2.5e3, -2.5e3, -2.5e3, -2.5e3]

    dpis = np.array([256, 512, 1024]) / 8
    dpis_dict = {'t': dpis[0], 'l': dpis[1], 'f': dpis[2]}

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
        pltEovsaQlookImage(datestr, spws, vmaxs, vmins, dpis_dict, fig, ax, overwrite=False, verbose=True)
        pltSdoQlookImage(datestr, dpis_dict, fig, ax, overwrite=False, verbose=True, clearcache=clearcache)


if __name__ == '__main__':
    import sys
    import numpy as np

    # import subprocess
    # shell = subprocess.check_output('echo $0', shell=True).decode().replace('\n', '').split('/')[-1]
    # print("shell " + shell + " is using")

    print(sys.argv)
    try:
        argv = sys.argv[1:]
        if '--clearcache' in argv:
            clearcache = True
            argv.remove('--clearcache')  # Allows --clearcache to be either before or after date items
        else:
            clearcache = False
        year = np.int(argv[0])
        month = np.int(argv[1])
        day = np.int(argv[2])
    except:
        print('Error interpreting command line argument')
        year = None
        month = None
        day = None
        clearcache = True
    print("Running pipeline_plt for date {}-{}-{}. clearcache {}".format(year, month, day, clearcache))
    main(year, month, day, clearcache)
