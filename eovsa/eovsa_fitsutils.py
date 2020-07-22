from astropy.io import fits
import os
from datetime import timedelta
from datetime import datetime
from glob import glob
import numpy as np
from astropy.time import Time
from suncasa.utils import fitsutils as fu

imgfitsdir = '/data1/eovsa/fits/synoptic/'
# imgfitsbkdir = '/data1/workdir/synoptic_bk/'

# imgfitsdir = '/data1/workdir/synoptic_bk/'
imgfitsbkdir = '/data1/workdir/synoptic_newbk/'


#
# def write_compress_image_fits(fname, data, header, mask=None, **kwargs):
#     """
#     Take a data header pair and write a compressed FITS file.
#
#     Parameters
#     ----------
#     fname : `str`
#         File name, with extension.
#     data : `numpy.ndarray`
#         n-dimensional data array.
#     header : `dict`
#         A header dictionary.
#     compression_type: `str`, optional
#         Compression algorithm: one of 'RICE_1', 'RICE_ONE', 'PLIO_1', 'GZIP_1', 'GZIP_2', 'HCOMPRESS_1'
#     hcomp_scale: `float`, optional
#         HCOMPRESS scale parameter
#     """
#     if kwargs is {}:
#         kwargs.update({'compression_type': 'RICE_1', 'quantize_level': 4.0})
#     if isinstance(fname, str):
#         fname = os.path.expanduser(fname)
#
#     hdunew = fits.CompImageHDU(data=data, header=header, **kwargs)
#     if mask is None:
#         hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew])
#     else:
#         hdumask = fits.CompImageHDU(data=mask.astype(np.uint8), **kwargs)
#         hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, hdumask])
#     hdulnew.writeto(fname, output_verify='fix')

def rewriteImageFits(datestr, verbose=False, writejp2=False):
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datestrdir
    imgbkdir = imgfitsbkdir + datestrdir
    if not os.path.exists(imgbkdir):
        os.makedirs(imgbkdir)

    if verbose: print('Processing EOVSA image fits files for date {}'.format(dateobj.strftime('%Y-%m-%d')))
    files = glob(os.path.join(imgindir, '*.tb.*fits'))
    for fl in files:
        if not os.path.exists(fl): continue
        filein = os.path.join(imgbkdir, os.path.basename(fl))
        if not os.path.exists(filein):
            os.system('mv {} {}'.format(fl, filein))
        hdul = fits.open(filein)
        for hdu in hdul:
            if hdu.header['NAXIS'] ==0:
                continue
            else:
                break
        data = np.squeeze(hdu.data).copy()
        # if verbose: print('Processing {}'.format(fl))
        if not os.path.exists(fl):
            data[np.isnan(data)] = 0.0
            fu.write_compress_image_fits(fl, data, hdu.header, compression_type='RICE_1', quantize_level=4.0)
        fj2name = fl.replace('.fits', '.jp2')
        if writejp2:
            if not os.path.exists(fj2name):
                data = np.squeeze(hdu.data).copy()
                fu.write_j2000_image(fj2name, data[::-1,:], hdu.header)
    return


def main(year=None, month=None, day=None, ndays=1):
    # tst = datetime.strptime("2017-04-01", "%Y-%m-%d")
    # ted = datetime.strptime("2019-12-31", "%Y-%m-%d")
    if year:
        ted = datetime(year, month, day)
    else:
        ted = datetime.now() - timedelta(days=2)
    tst = Time(np.fix(Time(ted).mjd) - ndays+1, format='mjd').datetime
    print("Running pipeline_fitsutils for date from {} to {}".format(tst.strftime("%Y-%m-%d"),
                                                                 ted.strftime("%Y-%m-%d")))
    dateobs = tst
    while dateobs <= ted:
        datestr = dateobs.strftime("%Y-%m-%d")
        rewriteImageFits(datestr, verbose=True, writejp2=True)
        dateobs = dateobs + timedelta(days=1)


if __name__ == '__main__':
    '''
    usage: eovsa_fitsutils -n 2 2020 06 10
    '''
    import sys
    import numpy as np
    import getopt
    from datetime import datetime, timedelta

    # import subprocess
    # shell = subprocess.check_output('echo $0', shell=True).decode().replace('\n', '').split('/')[-1]
    # print("shell " + shell + " is using")

    print(sys.argv)
    year = None
    month = None
    day = None
    ndays = 1
    clearcache = True
    opts = []
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(argv, "c:n:", ['clearcache=', 'ndays='])
        print(opts, args)
        for opt, arg in opts:
            if opt in ['-c', '--clearcache']:
                if arg is 'True':
                    clearcache = True
                elif arg is 'False':
                    clearcache = False
                else:
                    clearcache = np.bool(arg)
            elif opt in ('-n', '--ndays'):
                ndays = np.int(arg)
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
        opts = []



    main(year, month, day, ndays)
