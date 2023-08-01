from astropy.time import Time
import os
import requests
from bs4 import BeautifulSoup

'''
t=Time(['2023-05-09T18:35:19','2023-05-09T18:46:00'])
files=eovsa_filedownloader(t,outpath='/home/surajit')
print (files)
'''

def get_times_from_web(link):
    pieces=link.split('/')
    if pieces[-1]=='':
        ymd=pieces[-2]
        sep=1
    else:
        ymd=pieces[-1]
        sep=0
    year=ymd[:4] 
    page=requests.get(link).text
    soup = BeautifulSoup(page, 'html.parser')
    if sep==0:
        times=[node.get('href')[3:] for node in soup.find_all('a') if "IDB"+year in node.get('href')]
    else:
        times=[node.get('href')[3:-1] for node in soup.find_all('a') if "IDB"+year in node.get('href')]
    return times

def eovsa_filedownloader(trange, outpath='./'):
    #Given a timerange, this routine will take all relevant IDBfiles from
    #  that time range, put them in a list, and return that list.
    

        trange.format='datetime'
        year_month_date1=trange[0].strftime("%Y%m%d")
        year_month_date2=trange[1].strftime("%Y%m%d")

        if year_month_date1>='20170400' and year_month_date1<='20210400':
            link='https://research.ssl.berkeley.edu/data/eovsa/'
        elif year_month_date1>='20210300' and year_month_date1<='20221231':
            link='http://ovsa.njit.edu/IDB2/'
        else:
            link='http://www.ovsa.njit.edu/fits/IDB/'
            
            
        if year_month_date1==year_month_date2:
            ymd=[year_month_date1]
        else:
            diff=trange[1]-trange[0]
            diff_seconds=abs(diff.value)
            if diff_seconds>1:
                raise RuntimeError("The supplied two times differby more than 1 day!!!")
            else:
                ymd=[year_month_date1,year_month_date2]

        for year_month_date in ymd:
            #str1='curl -s '+link+year_month_date+'/ | grep href=\\\"IDB | awk \'{print $5}\' | cut -d\'\"\' -f2 | cut -d\'I\' -d\'B\' -f2 | cut -d\'/\' -f1'
            #print (str1)
            #filelist=os.popen(str1).read()
            #times1=filelist.split('\n')[:-1]
            times1=get_times_from_web(link+year_month_date+'/')
            times=Time([t1[:4]+"-"+t1[4:6]+"-"+t1[6:8]+"T"+t1[8:10]+":"+t1[10:12]+":"+t1[12:] for t1 in times1])
            mint=trange[0]
            maxt=trange[1]
            down_times=[]
            for t in times:
                if t>=mint and t<=maxt:
                    down_times.append(t)
        down_times=Time(down_times)
        down_times.format='datetime'
        files=[]
        os.chdir(outpath)
        for t in down_times:
            filename='IDB'+t.strftime("%Y%m%d%H%M%S")
            str1="wget -r -np -e robots=\"off\" -R \"index.html\" "+link+\
                    year_month_date+"/"+filename+"/"
            os.system(str1)
            os.system("mv "+link.split('//')[1]+year_month_date+'/'+filename+" ./")
            os.system("rm -rf "+link.split('//')[1].split('/')[0])
            files.append(filename)
        files.sort()
        return files


