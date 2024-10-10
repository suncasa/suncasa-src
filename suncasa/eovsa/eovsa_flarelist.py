##ipython eovsa_flarelist.py
##=============
import os
import sys
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from scipy.signal import find_peaks
from datetime import datetime, timedelta
from astropy.time import Time
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# Global settings and initializations
# EO_WIKI_URL = "http://www.ovsa.njit.edu/wiki/index.php/Recent_Flare_List_(2021-)"
EO_WIKI_URLs = [
    "https://www.ovsa.njit.edu/wiki/index.php/2017",
    "https://www.ovsa.njit.edu/wiki/index.php/2018",
    "https://www.ovsa.njit.edu/wiki/index.php/2019",
    "https://www.ovsa.njit.edu/wiki/index.php/2020",
    "https://www.ovsa.njit.edu/wiki/index.php/2021",
    "https://www.ovsa.njit.edu/wiki/index.php/2022",
    "https://www.ovsa.njit.edu/wiki/index.php/2023",
    "https://www.ovsa.njit.edu/wiki/index.php/2024",
    "https://www.ovsa.njit.edu/wiki/index.php/2025",
    "https://www.ovsa.njit.edu/wiki/index.php/2026"
]


def fetch_flare_data_from_wiki(eo_wiki_urls, given_date_strp, outcsvfile):
    """
    Fetches flare data from the given wiki URL and saves it to a CSV file in the specified directory.

    Parameters:
    - eo_wiki_url: URL of the EO wiki page to fetch data from.
    - given_date_strp: A tuple of datetime objects specifying the start and end dates to filter the data.
    - outcsvfile: path to the resulting CSV file.

    Returns:
    - None
    """
    date_data = []
    time_ut_data = []
    flare_class_data = []
    flare_ID_data = []
    depec_file = []

    for eo_wiki_url in eo_wiki_urls:
        response = requests.get(eo_wiki_url)

        # Check if the request was successful
        if response.status_code == 200:
            # Parse the HTML content of the page using BeautifulSoup
            soup = BeautifulSoup(response.text, "html.parser")
            tables = soup.find_all("table", {"class": "wikitable"})

            for table in tables:
                for row in table.find_all("tr"):
                    cells = row.find_all("td")

                    if len(cells) >= 3:
                        date = cells[0].text.strip()
                        time_ut = cells[1].text.strip()
                        flare_class = cells[2].text.strip()
                        datetime_strp = datetime.strptime(date + ' ' + time_ut, '%Y-%m-%d %H:%M')

                        if given_date_strp[0] <= datetime_strp <= given_date_strp[1]:
                            # get Date, Time (UT) and Flare Class
                            date_data.append(date)
                            time_ut_data.append(f"{time_ut}:00")
                            flare_class_data.append(flare_class)
                            # get Flare_ID
                            eotime = f"{date.replace('/', '-')}T{time_ut}:00"
                            flare_ID_data.append(eotime.replace('-', '').replace('T', '').replace(':', ''))
                            print("Fetched: ", eotime)

                            depec_file_tmp = ''
                            for cell in cells:
                                link_cell = cell.find('a', class_='external text', href=True, rel='nofollow')
                                if link_cell:
                                    url = link_cell['href']
                                    if url.endswith(".dat") or url.endswith(".fits"):
                                        depec_file_tmp = url.split('/')[-1]
                                        break
                            # get depec_file
                            depec_file.append(str(depec_file_tmp))
        else:
            print("Failed to retrieve the webpage. Status code:", response.status_code)

    # reformate the date and time above
    data = {
        "ID": np.arange(len(date_data)) + 1,
        "Flare_ID": flare_ID_data,
        "Date": date_data,
        "Time (UT)": time_ut_data,
        "Flare Class": flare_class_data,
        "depec_file": depec_file
    }

    df = pd.DataFrame(data)
    df.to_csv(outcsvfile, index=False)
    print(f"Date and Time (UT) data saved to {outcsvfile}")


def rd_datfile(file):
    ''' Read EOVSA binary spectrogram file and return a dictionary with times
        in Julian Date, frequencies in GHz, and cross-power data in sfu.

        Return Keys:
          'time'     Numpy array of nt times in JD format
          'fghz'     Numpy array of nf frequencies in GHz
          'data'     Numpy array of size [nf, nt] containing cross-power data

        Returns empty dictionary ({}) if file size is not compatible with inferred dimensions
    '''
    import struct
    import numpy as np
    def dims(file):
        # Determine time and frequency dimensions (assumes the file has fewer than 10000 times)
        f = open(file, 'rb')
        tmp = f.read(83608)  # max 10000 times and 451 frequencies
        f.close()
        nbytes = len(tmp)
        tdat = np.array(struct.unpack(str(int(nbytes / 8)) + 'd', tmp[:nbytes]))
        nt = np.where(tdat < 2400000.)[0]
        nf = np.where(np.logical_or(tdat[nt[0]:] > 18, tdat[nt[0]:] < 1))[0]
        return nt[0], nf[0]

    nt, nf = dims(file)
    f = open(file, 'rb')
    tmp = f.read(nt * 8)
    times = struct.unpack(str(nt) + 'd', tmp)
    tmp = f.read(nf * 8)
    fghz = struct.unpack(str(nf) + 'd', tmp)
    tmp = f.read()
    f.close()
    if len(tmp) != nf * nt * 4:
        print('File size is incorrect for nt=', nt, 'and nf=', nf)
        return {}
    data = np.array(struct.unpack(str(nt * nf) + 'f', tmp)).reshape(nf, nt)
    return {'time': times, 'fghz': fghz, 'data': data}


def moving_average(data, window_size):
    # Create a convolution kernel for the moving average
    kernel = np.ones(window_size) / window_size
    return np.convolve(data, kernel, mode='valid')


def get_eoflarelist(timerange=None, work_dir='./', file_out='EOVSA_flare_list_from_wiki_sub.csv'):
    '''
    Get eovsa flarelist from given timerange.
    On pipeline or ovsa server.

    Parameters:
    timerange:

    Example:
    final_csv_file = get_eoflarelist(timerange=["2024-01-01 20:00:00", "2024-01-01 21:00:00"])

    '''
    # Convert timerange strings to datetime objects for internal use

    timerange_strp = [datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S') for date_str in timerange]
    print(f"Fetching data for time range {timerange[0]} to {timerange[1]}")

    ##=============Step 1: Fetching data from wiki=============
    print("##=============Step 1: capture the radio peak times from flare list wiki webpage")

    init_csv = os.path.join(work_dir, "get_time_from_wiki_given_date.csv")
    fetch_flare_data_from_wiki(EO_WIKI_URLs, timerange_strp, init_csv)

    ##=============Step 2: read the spectrum data=============
    print("##=============Step 2: read the spectrum data")

    spec_data_dir = "/common/webplots/events/"  # YYYY/
    print("Spec data in ", spec_data_dir)

    df = pd.read_csv(init_csv)
    flare_id = df['Flare_ID']
    depec_file = df['depec_file']

    files_wiki = [spec_data_dir + str(flare_id[i])[0:4] + "/" + str(file_name) for i, file_name in
                  enumerate(depec_file)]

    ##=============
    tpk_spec_wiki = [] 
    tst_mad_spec_wiki, ted_mad_spec_wiki = [], []
    tst_thrd_spec_wiki, ted_thrd_spec_wiki = [], []
    depec_file = []

    ##=============
    for ww, file_wiki in enumerate(files_wiki):  # len(files_wiki)

        filename1 = os.path.basename(file_wiki)
        depec_file.append(filename1)
        print("Reading spec data: ", filename1)

        try:
            if filename1.split('.')[-1] == 'dat':
                filename = filename1.split('.dat')[0]
                data1 = rd_datfile(file_wiki)
                spec = np.array(data1['data'])
                fghz = np.array(data1['fghz'])
                time1 = data1['time']
            if filename1.split('.')[-1] == 'fits':
                filename = filename1.split('.fits')[0]
                eospecfits = fits.open(file_wiki)
                spec = eospecfits[0].data  # [freq, time]
                fghz = np.array(eospecfits[1].data['FGHZ'])  # in GHz
                time1 = np.array(eospecfits[2].data['TIME'])  # in jd format

        except Exception as e:
            temp = datetime.strptime(str(flare_id[ww]), "%Y%m%d%H%M%S")
            temp_st = (temp - timedelta(minutes=2)).strftime("%Y-%m-%d %H:%M:%S")
            temp_ed = (temp + timedelta(minutes=2)).strftime("%Y-%m-%d %H:%M:%S")

            tpk_spec_wiki.append(temp.strftime("%Y-%m-%d %H:%M:%S"))
            tst_mad_spec_wiki.append(temp_st)
            ted_mad_spec_wiki.append(temp_ed)
            tst_thrd_spec_wiki.append(temp_st)
            ted_thrd_spec_wiki.append(temp_ed)

            print('no data found for:', file_wiki)
            print('st/ed time will be tpk \u00B1 2mins: ', temp_st, temp_ed)
            continue

        time_spec = Time(time1, format='jd')
        time_str = time_spec.isot

        ##=============try MAD method
        tpk_tot_ind = []

        tst_tot_ind = []
        ted_tot_ind = []
        mad_threshd = []

        tst_tot_ind_thrd = []
        ted_tot_ind_thrd = []

        spec = np.nan_to_num(spec, nan=0.0)
        spec[spec < 0] = 0.01
        spec_median = np.median(spec, axis=1, keepdims=True)
        spec_abs_deviations = np.abs(spec - spec_median)
        spec_mad = np.median(spec_abs_deviations, axis=1, keepdims=True)
        mad_threshd = 3.0 * spec_mad

        outliers = spec_abs_deviations > mad_threshd

        for ff in range(len(fghz)):
            good_channel = False
            flux_array = spec[ff, :]

            ##=============try MAD method
            outlier_ind = np.where(outliers[ff])[0]
            if len(outlier_ind) > 0:
                tst_tot_ind.append(outlier_ind[0])
                ted_tot_ind.append(outlier_ind[-1])

            ##=============try to set threshold
            window_size = 10
            y = moving_average(flux_array, window_size) + 0.001  ##flux_array
            peaks, _ = find_peaks(y, height=1.)
            noise_thrd_st = np.mean(y[0:5])
            noise_thrd_ed = np.mean(y[-5:])

            if noise_thrd_st == 0:
                noise_thrd_st = 0.005 * np.max(y)
            if noise_thrd_ed == 0:
                noise_thrd_ed = 0.01 * np.max(y)
            if np.max(y) / noise_thrd_ed > 100:
                noise_thrd_ed = 0.02 * np.max(y)
                # noise_thrd_st = np.max([np.mean(y[0:10]),0.01*np.max(y)])
            # noise_thrd_ed = np.max([np.mean(y[-10:])*np.median(y[peaks]),0.01*np.max(y)])
            # if ff == 5:
            #     print(np.max(y), np.median(y[peaks]), np.mean(y[0:10]), np.mean(y[-10:]))

            for ind in range(len(y) - 5):
                if y[ind] < y[ind + 1] < y[ind + 2] < y[ind + 3] < y[ind + 4] < y[ind + 5]:
                    if y[ind + 5] >= 2 * y[ind]:
                        if all(y[i] > noise_thrd_st for i in range(ind, ind + 6)):
                            tst_tot_ind_thrd.append(ind)
                            break
            ind_tmp = np.argmax(flux_array) - 30
            for ind in range(len(y) - 5):
                if y[ind] > y[ind + 1] > y[ind + 2] > y[ind + 3] > y[ind + 4] > y[ind + 5]:
                    if y[ind + 3] <= 2 * y[ind]:
                        if all(abs(y[i]) > noise_thrd_ed for i in range(ind, ind + 6)):
                            ind_tmp = ind + 5
            ted_tot_ind_thrd.append(ind_tmp)

            ##=============tpeak
            tpk_tot_ind.append(np.argmax(flux_array))

        time_st_mad = time_spec[int(np.round(np.median(np.array(tst_tot_ind))))]
        time_ed_mad = time_spec[int(np.round(np.median(np.array(ted_tot_ind))))]

        time_st_thrd = time_spec[int(np.round(np.median(np.array(tst_tot_ind_thrd))))]
        time_ed_thrd = time_spec[int(np.round(np.median(np.array(ted_tot_ind_thrd))))]

        time_st = time_st_thrd
        time_ed = time_ed_thrd

        time_pk = time_spec[int(np.median(np.array(tpk_tot_ind)))]

        time_pk_obj = Time(time_pk, format='jd')

        tpk_spec_wiki.append(time_pk_obj.strftime('%Y-%m-%d %H:%M:%S'))

        tst_mad_spec_wiki.append(Time(time_st_mad, format='jd').strftime('%Y-%m-%d %H:%M:%S'))
        ted_mad_spec_wiki.append(Time(time_ed_mad, format='jd').strftime('%Y-%m-%d %H:%M:%S'))

        tst_thrd_spec_wiki.append(Time(time_st_thrd, format='jd').strftime('%Y-%m-%d %H:%M:%S'))
        ted_thrd_spec_wiki.append(Time(time_ed_thrd, format='jd').strftime('%Y-%m-%d %H:%M:%S'))


    #apply the thrd method
    EO_tpeak = tpk_spec_wiki
    EO_tstart_thrd = tst_thrd_spec_wiki
    EO_tend_thrd = ted_thrd_spec_wiki

    data_csv = {
        "ID": np.arange(len(flare_id)) + 1,
        "Flare_ID": flare_id,
        "Date": df['Date'],
        "Time (UT)": df['Time (UT)'],
        "Flare Class": df['Flare Class'],
        'EO_tstart': EO_tstart_thrd,
        'EO_tpeak': EO_tpeak,
        'EO_tend': EO_tend_thrd,
        'depec_file': df['depec_file']
    }

    df = pd.DataFrame(data_csv)
    final_csv_file = os.path.join(work_dir, file_out)
    df.to_csv(final_csv_file, index=False)
    print(f"Success: {len(flare_id)} flares saved to {file_out}")

    return final_csv_file


def run_eoflare_pipeline(flarelist_csv=None, to_web=True):
    ##=========================reading flarelist from .csv=========================
    import pandas as pd
    import os
    import shutil

    df = pd.read_csv(flarelist_csv)
    flare_id_tot = df['Flare_ID']
    EO_tstart_tot = df['EO_tstart']
    EO_tend_tot = df['EO_tend']

    for i,j in enumerate(flare_id_tot):
        print(f"Flare_ID, st, ed time: {flare_id_tot[i]}, {EO_tstart_tot[i]}, {EO_tend_tot[i]}")

    ##=========================running flare pipeline=========================
    ans = 'Y'
    # ans = input('Do you want to continue to run flare pipeline for Flare_IDs above? (say no if you want to stop) [y/n]?')

    root_dir = os.getcwd()

    if ans.upper() == 'Y':
        from suncasa.eovsa import eovsa_flare_pipeline
        from eovsapy.util import Time

        for ff, flare_id in enumerate(flare_id_tot):#flare_id_tot, len(flare_id_tot)
            os.chdir(root_dir)
            try:
                EO_tstart = EO_tstart_tot[ff]
                EO_tend = EO_tend_tot[ff]
                flare_id = str(flare_id)
                print(f"To run flare pipeline for Flare_ID: {flare_id}")

                work_dir = os.path.join(root_dir, flare_id)
                os.makedirs(work_dir, exist_ok=True)
                os.chdir(work_dir)

                trange_str = [EO_tstart, EO_tend]
                trange = Time(trange_str)
                fp = eovsa_flare_pipeline.FlareSelfCalib(vis=trange)
                fp.imaging_start = trange_str[0]
                fp.imaging_end = trange_str[1]
                fp.slfcal_pipeline(doselfcal=True, doimaging=True)
                if to_web:
                    fp.rename_move_files(flare_id=flare_id, fitsdir_web_tp='/data1/eovsa/fits/flares/', 
                                        movdir_web_tp='/common/webplots/SynopticImg/eovsamedia/eovsa-browser/',
                                        msdir_web_tp='/data1/eovsa/fits/flares/', 
                                        dorename_fits=True, domove_fits=True, 
                                        dorename_mov=True, domove_mov=True, domove_ms=True, dormworkdir=True, docopy=True)
                print(f"Success for Flare_ID: {flare_id}")

            except Exception as e:
                print(f"Error at ff = {ff} for Flare_id = {flare_id}: {e}")
                continue
        os.chdir(root_dir)
        print('Done!')

# if __name__ == "__main__":
#     main()
