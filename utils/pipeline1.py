#!/common/casa/casa-release-5.4.1-31.el6/bin/casa

def pipeline1(year=None, month=None, day=None, clearcache=True,overwrite=True, doimport=True):
    from suncasa.eovsa import eovsa_pipeline as ep
    import os
    from eovsapy.util import Time

    workdir = '/data1/workdir/'
    os.chdir(workdir)
    # Set to run 5 days earlier than the current date
    if year is None:
        mjdnow = Time.now().mjd
        t = Time(mjdnow - 2, format='mjd')
    else:
        # Uncomment below and set date to run for a given date
        t = Time('{}-{:02d}-{:02d} 20:00'.format(year, month, day))
    print(t.iso)
    datestr = t.iso[:10]

    subdir = t.datetime.strftime('%Y%m%d/')
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    vis_corrected = ep.calib_pipeline(datestr, overwrite=overwrite, doimport=doimport,
                                      workdir=os.path.join(workdir, subdir))
    if clearcache:
        os.chdir(workdir)
        os.system('rm -rf ' + subdir)


if __name__ == '__main__':
    import sys
    import numpy as np

    print(sys.argv)
    try:
        argv = sys.argv[3:]
        if '--clearcache' in argv:
            clearcache = True
            argv.remove('--clearcache')  # Allows --clearcache to be either before or after date items
        else:
            clearcache = False
        try:
            year = np.int(argv[0])
            month = np.int(argv[1])
            day = np.int(argv[2])
        except:
            print('Error parsing date string in the command line argument')
            year = None
            month = None
            day = None
    except:
        print('Error interpreting command line argument')
        year = None
        month = None
        day = None
        clearcache = True
    print("Running pipeline1 for date {}-{}-{}. clearcache {}".format(year, month, day, clearcache))
    pipeline1(year, month, day, clearcache=clearcache)
