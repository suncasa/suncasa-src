from astropy.io import fits
import os
from datetime import timedelta
from datetime import datetime
from glob import glob
import numpy as np
from astropy.time import Time
from suncasa.io import ndfits

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


def rewriteImageFits(datestr, verbose=False, writejp2=False, overwritejp2=False, overwritefits=False):
    dateobj = datetime.strptime(datestr, "%Y-%m-%d")
    datestrdir = dateobj.strftime("%Y/%m/%d/")
    imgindir = imgfitsdir + datestrdir
    imgbkdir = imgfitsbkdir + datestrdir
    if not os.path.exists(imgbkdir):
        os.makedirs(imgbkdir)

    if verbose: print('Processing EOVSA image fits files for date {}'.format(dateobj.strftime('%Y-%m-%d')))
    files = glob(os.path.join(imgindir, '*.tb.*fits'))
    files = sorted(files)
    for fl in files:
        if not os.path.exists(fl): continue
        hdul = fits.open(fl)
        if len(hdul) > 1:
            hdul.close()
            continue
        else:
            hdul.close()
        filein = os.path.join(imgbkdir, os.path.basename(fl))
        if not os.path.exists(filein):
            os.system('mv {} {}'.format(fl, filein))
        hdul = fits.open(filein)
        for hdu in hdul:
            if hdu.header['NAXIS'] == 0:
                continue
            else:
                break
        data = np.squeeze(hdu.data).copy()
        # if verbose: print('Processing {}'.format(fl))
        if overwritefits:
            if os.path.exists(fl):
                os.system('rm -f {}'.format(fl))
        if not os.path.exists(fl):
            data[np.isnan(data)] = 0.0
            ndfits.write(fl, data, hdu.header, compression_type='RICE_1', quantize_level=4.0)

        fj2name = fl.replace('.fits', '.jp2')
        if writejp2:
            if overwritejp2:
                if os.path.exists(fj2name):
                    os.system('rm -f {}'.format(fj2name))
            if not os.path.exists(fj2name):
                data = np.squeeze(hdu.data).copy()
                data[np.isnan(data)] = 0.0
                ndfits.write_j2000_image(fj2name, data[::-1, :], hdu.header)
    return

def main(dateobj=None, ndays=1, overwritejp2=False, overwritefits=False):
    """
    Main pipeline for creating compressed FITS and JP2 files of EOVSA daily full-disk images.

    :param dateobj: The starting datetime for processing. If None, defaults to two days before now.
    :type dateobj: datetime, optional
    :param ndays: Number of days to process (spanning from dateobj - ndays + 1 to dateobj), defaults to 1.
    :type ndays: int, optional
    :param overwritejp2: If True, overwrite existing EOVSA JP2 files, defaults to False.
    :type overwritejp2: bool, optional
    :param overwritefits: If True, overwrite existing EOVSA FITS files, defaults to False.
    :type overwritefits: bool, optional
    :raises Exception: If an error occurs during processing.
    :return: None
    :rtype: None
    """
    from datetime import timedelta
    import numpy as np
    from astropy.time import Time

    # Use dateobj if provided; otherwise, default to two days before now.
    if dateobj is None:
        ted = datetime.now() - timedelta(days=2)
    else:
        ted = dateobj

    # Compute the start date (tst) for processing based on ndays.
    tst = Time(np.fix(Time(ted).mjd) - ndays + 1, format='mjd').datetime
    print("Running pipeline_fitsutils for date from {} to {}.".format(
        tst.strftime("%Y-%m-%d"), ted.strftime("%Y-%m-%d")))
    dateobs = tst
    while dateobs <= ted:
        datestr = dateobs.strftime("%Y-%m-%d")
        rewriteImageFits(datestr, verbose=True, writejp2=True,
                          overwritejp2=overwritejp2, overwritefits=overwritefits)
        dateobs = dateobs + timedelta(days=1)


if __name__ == '__main__':
    import argparse
    from datetime import datetime, timedelta
    from astropy.time import Time

    parser = argparse.ArgumentParser(
        description='Pipeline for creating compressed FITS and JP2 files of EOVSA daily full-disk images.'
    )
    # Default date is set to two days before the current date at 20:00 UT (YYYY-MM-DDT20:00).
    default_date = (datetime.now() - timedelta(days=2)).strftime('%Y-%m-%dT20:00')
    parser.add_argument(
        '--date', type=str, default=default_date,
        help='Date to process in YYYY-MM-DDT20:00 format, defaults to 20:00 UT two days before the current date.'
    )
    parser.add_argument(
        '--ndays', type=int, default=1,
        help='Process data spanning from DATE minus ndays to DATE (default: 1 day).'
    )
    parser.add_argument(
        '--overwritejp2', action='store_true',
        help='Overwrite existing EOVSA JP2 files.'
    )
    parser.add_argument(
        '--overwritefits', action='store_true',
        help='Overwrite existing EOVSA FITS files.'
    )
    # Optional positional date arguments: year month day (overrides --date if provided)
    parser.add_argument(
        'date_args', type=int, nargs='*',
        help='Optional date arguments: year month day. If provided, overrides --date.'
    )

    args = parser.parse_args()

    # Determine the processing date.
    if len(args.date_args) == 3:
        year, month, day = args.date_args
        dateobj = datetime(year, month, day, 20)  # Use 20:00 UT for the specified date.
    else:
        dateobj = Time(args.date).datetime

    print(f"Running eovsa_fitsutils for date {dateobj.strftime('%Y-%m-%d')}.")
    print("Arguments:")
    print(f"  ndays: {args.ndays}")
    print(f"  overwritejp2: {args.overwritejp2}")
    print(f"  overwritefits: {args.overwritefits}")

    # Call the main function with the parsed datetime object.
    main(dateobj, args.ndays, overwritejp2=args.overwritejp2, overwritefits=args.overwritefits)
