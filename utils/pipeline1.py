#!/common/casa/casa-release-5.4.1-31.el6/bin/casa

def pipeline1(year=None, month=None, day=None, clearcache=True):
    from suncasa.eovsa import eovsa_pipeline as ep
    import os
    from astropy.time import Time

    os.chdir('/data1/workdir/')
    # Set to run 5 days earlier than the current date
    if year is None:
        mjdnow = Time.now().mjd
        t = Time(mjdnow - 2, format='mjd')
    else:
        # Uncomment below and set date to run for a given date
        t = Time('{}-{:02d}-{:02d} 20:00'.format(year, month, day))
    print(t.iso)
    date = t.iso[:10]

    workdir = t.datetime.strftime('%Y%m%d/')
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    vis_corrected = ep.calib_pipeline(date, overwrite=True, doimport=True, workdir=workdir)
    if clearcache:
        os.system('rm -rf ' + workdir)


if __name__ == '__main__':
    import sys
    import numpy as np

    print(sys.argv)
    try:
        argv = sys.argv[3:]
        year = np.int(argv[0])
        month = np.int(argv[1])
        day = np.int(argv[2])
        if '--clearcache' in argv:
            clearcache = True
        else:
            clearcache = False
    except:
        print('Error interpreting command line argument')
        year = None
        month = None
        day = None
        clearcache = True
    pipeline1(year, month, day, clearcache=clearcache)
