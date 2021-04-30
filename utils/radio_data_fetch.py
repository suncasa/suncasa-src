import pandas as pd
# import urllib2
import urllib.request
import datetime as dt
from astropy.time import Time
import gzip
from tqdm import tqdm
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pdb


def get_rstn_data(time, outdir='./RSTN/', ylim=[0, 800]):
    rstn_csv_dir = outdir + 'csv/'
    if not os.path.exists(os.path.join(rstn_csv_dir, 'png/')):
        os.makedirs(os.path.join(rstn_csv_dir, 'png/'))
    if not os.path.exists(os.path.join(outdir, 'data/')):
        os.makedirs(os.path.join(outdir, 'data/'))
    OBSCID = {'k7o': 'sagamore-hill', 'phf': 'palehua', 'apl': 'learmonth', 'lis': 'san-vito'}

    date = time.datetime
    for obscid, obsc in OBSCID.items():
        gzfilenameout = '{}.{}'.format(date.strftime('%Y-%m-%d'), obscid.upper())
        gzfile = os.path.join(outdir, 'data/', gzfilenameout + '.gz')
        if not os.path.exists(gzfile):
            try:
                gzfilename = '{}.{}'.format(date.strftime('%d%b%y').upper(), obscid.upper())
                gzfile_url = "ftp://ftp.ngdc.noaa.gov/STP/space-weather/solar-data/solar-features/solar-radio/rstn-1-second/{p[obsc]}/{p[yy]:04d}/{p[mm]:02d}/{p[gzfile]}.gz".format(
                    p={'gzfile': gzfilename, 'yy': date.year, 'mm': date.month, 'obsc': obsc})
                # gzfiledata = urllib2.urlopen(gzfile_url)
                with urllib.request.urlopen(gzfile_url) as response, open(gzfile, 'wb') as out_file:
                    gzfiledata = response.read()  # a `bytes` object
                    out_file.write(gzfiledata)
                print('Loading {}'.format(gzfile_url))
            except:
                try:
                    gzfilename = '{}.{}'.format(date.strftime('%d%b%y').lower(), obscid.lower())
                    gzfile_url = "ftp://ftp.ngdc.noaa.gov/STP/space-weather/solar-data/solar-features/solar-radio/rstn-1-second/{p[obsc]}/{p[yy]:04d}/{p[mm]:02d}/{p[gzfile]}.gz".format(
                        p={'gzfile': gzfilename, 'yy': date.year, 'mm': date.month, 'obsc': obsc})
                    # gzfiledata = urllib2.urlopen(gzfile_url)
                    with urllib.request.urlopen(gzfile_url) as response, open(gzfile, 'wb') as out_file:
                        gzfiledata = response.read()  # a `bytes` object
                        out_file.write(gzfiledata)
                    print('Loading {}'.format(gzfile_url))
                except:
                    print('Fail to load {}'.format(gzfile_url))
                    continue
            # with open(gzfile, 'wb') as f:
            #     f.write(gzfiledata.read())
        with gzip.open(gzfile, 'rb') as f:
            # lines = f.readlines()
            lines = f.read().splitlines()
        YYYY = lines[0][4:8]
        MM = lines[0][8:10]
        DD = lines[0][10:12]
        HH = lines[0][12:14]
        mm = lines[0][14:16]
        SS = lines[0][16:18]
        tst = Time(dt.datetime(int(YYYY), int(MM), int(DD), int(HH), int(mm), int(SS)))
        YYYY = lines[-1][4:8]
        MM = lines[-1][8:10]
        DD = lines[-1][10:12]
        HH = lines[-1][12:14]
        mm = lines[-1][14:16]
        SS = lines[-1][16:18]
        ted = Time(dt.datetime(int(YYYY), int(MM), int(DD), int(HH), int(mm), int(SS)))
        # pdb.set_trace()
        if tst.mjd < time.mjd < ted.mjd:
            rstn_flux = {'time': [], 'OBSC': [], 'FFF245': [], 'FFF410': [], 'FFF610': [], 'FF1415': [], 'FF2695': [],
                         'FF4995': [], 'FF8800': [],
                         'F15400': []}
            gzfile_csv = gzfilenameout + '.csv'
            for idx, ll in enumerate(lines):
                line = str(ll, 'utf-8')
                OBSC = line[0:4]
                YYYY = line[4:8]
                MM = line[8:10]
                DD = line[10:12]
                HH = line[12:14]
                mm = line[14:16]
                SS = line[16:18]
                ltime = Time(dt.datetime(int(YYYY), int(MM), int(DD), int(HH), int(mm), int(SS)))
                if time.mjd - 0.5 / 24 <= ltime.mjd <= time.mjd + 1.0 / 24:
                    ll_list = []
                    for l in line.split(' '):
                        if l:
                            ll_list.append(l)
                    # FFF245 = ll[18:24]
                    # FFF410 = ll[24:30]
                    # FFF610 = ll[30:36]
                    # FF1415 = ll[36:42]
                    # FF2695 = ll[42:48]
                    # FF4995 = ll[48:54]
                    # FF8800 = ll[54:60]
                    # F15400 = ll[60:66]
                    FFF245 = ll_list[1 + 0]
                    FFF410 = ll_list[1 + 1]
                    FFF610 = ll_list[1 + 2]
                    FF1415 = ll_list[1 + 3]
                    FF2695 = ll_list[1 + 4]
                    FF4995 = ll_list[1 + 5]
                    FF8800 = ll_list[1 + 6]
                    F15400 = ll_list[1 + 7]
                    rstn_flux['time'].append(ltime.iso)
                    rstn_flux['OBSC'].append(OBSC)
                    try:
                        rstn_flux['FFF245'].append(float(FFF245))
                    except:
                        rstn_flux['FFF245'].append(np.nan)
                    try:
                        rstn_flux['FFF410'].append(float(FFF410))
                    except:
                        rstn_flux['FFF410'].append(np.nan)
                    try:
                        rstn_flux['FFF610'].append(float(FFF610))
                    except:
                        rstn_flux['FFF610'].append(np.nan)
                    try:
                        rstn_flux['FF1415'].append(float(FF1415))
                    except:
                        rstn_flux['FF1415'].append(np.nan)
                    try:
                        rstn_flux['FF2695'].append(float(FF2695))
                    except:
                        rstn_flux['FF2695'].append(np.nan)
                    try:
                        rstn_flux['FF4995'].append(float(FF4995))
                    except:
                        rstn_flux['FF4995'].append(np.nan)
                    try:
                        rstn_flux['FF8800'].append(float(FF8800))
                    except:
                        rstn_flux['FF8800'].append(np.nan)
                    try:
                        rstn_flux['F15400'].append(float(F15400))
                    except:
                        rstn_flux['F15400'].append(np.nan)
            rstn_flux_df = pd.DataFrame(rstn_flux)
            rstn_flux_df.time = pd.to_datetime(rstn_flux_df['time'])


            try:
                ax = rstn_flux_df.plot(x='time')
                ax.get_xaxis().set_major_formatter(mpl.dates.DateFormatter('%H:%M:%S'))
                ax.set_xlabel(rstn_flux_df.time[0].strftime('%Y-%m-%d'))
                ax.set_ylabel('solar radio flux [sfu]')
                ax.set_ylim(ylim)
                ax.set_xlim(time.plot_date - 0.5 / 24, time.plot_date + 1.0 / 24)
                ax.axvline(time.plot_date, ls=':')
                fig = plt.gcf()
                figfile = os.path.join(rstn_csv_dir, 'png/', gzfilenameout + '.png')
                fig.savefig(figfile)
                plt.close(fig)
            except:
                print('Fail to plot')
            # with open(rstn_df_dir + gzfile_csv, 'wb') as f:
            #     pickle.dump(rstn_flux_df, f)
            rstn_flux_df.to_csv(os.path.join(rstn_csv_dir, gzfile_csv))
            print('RSTN {} data saved to {}'.format(obsc, os.path.join(rstn_csv_dir, gzfile_csv)))
        else:
            print('Flare time {} is out the time range ({}~{}) in the downloaded file: {}.'.format(time.iso, tst.iso,
                                                                                                   ted.iso, gzfile))
    return
