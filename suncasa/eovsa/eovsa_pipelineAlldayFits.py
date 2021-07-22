if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")

    import pipeline_cal as pc
    import eovsa_fits as ef
    import glob
    from astropy.time import Time
    import sys, os

    print(sys.argv)
    try:
        argv = sys.argv[1:]
        if '--clearcache' in argv:
            clearcache = True
            argv.remove('--clearcache')   # Allows --clearcache to be either before or after date items
        else:
            clearcache = False
        year = argv[0]
        month = argv[1]
        day = argv[2]
        t = Time(year+'-'+month+'-'+day+' 20:00:00')
    except:
        print('Error interpreting command line arguments--will analyze data from yesterday.')
        # No arguments (or no arguments given), so default to yesterday's data to analyze
        mjdnow = Time.now().mjd
        t = Time(mjdnow-1,format='mjd')
        clearcache = True
    # Reread year month day to preserve leading 0, and make time an array
    year, month, day = t.iso.split('-')
    day = day.split(' ')[0]
    t = Time([t.iso])
    # Change to standard working directory and delete any existing IDB files there
    datstr = t[0].iso[:10].replace('-','')+'/'
    outpath = '/data1/dgary/HSO/'+datstr
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    fitsoutpath = '/data1/eovsa/fits/synoptic/'
    os.chdir(outpath)
    os.system('rm -rf IDB*')
    # Run first (and lengthy!) task to create corrected IDB files for the entire day
    pc.allday_udb_corr(t, outpath=outpath)
    # Process the entire day's IDB files to create fits files
    pc.allday_process(path=outpath)
    print(outpath,year,month,day,'Finding files:',outpath+year+'/'+month+'/'+day+'/*_TP_*.fts')
    files = glob.glob(outpath+year+'/'+month+'/'+day+'/*_TP_*.fts')
    files.sort()
    print(len(files),'files found')
    spec = ef.eovsa_combinefits(files, freqgaps=True, outpath=fitsoutpath, ac_corr=True, savfig=True)
    print(outpath,year,month,day,'Finding files:',outpath+year+'/'+month+'/'+day+'/*_XP_*.fts')
    files = glob.glob(outpath+year+'/'+month+'/'+day+'/*_XP_*.fts')
    files.sort()
    print(len(files),'files found')
    spec = ef.eovsa_combinefits(files, freqgaps=True, outpath=fitsoutpath, ac_corr=True, savfig=True)
    if clearcache:
        os.chdir('..')
        os.system('rm -rf '+datstr)
