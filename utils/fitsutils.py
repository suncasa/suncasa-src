import os
import numpy as np
from astropy.io import fits

def write_compress_image_fits(fname, data, header, mask=None, **kwargs):
    """
    Take a data header pair and write a compressed FITS file.

    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        n-dimensional data array.
    header : `dict`
        A header dictionary.
    compression_type: `str`, optional
        Compression algorithm: one of 'RICE_1', 'RICE_ONE', 'PLIO_1', 'GZIP_1', 'GZIP_2', 'HCOMPRESS_1'
    hcomp_scale: `float`, optional
        HCOMPRESS scale parameter
    """
    if kwargs is {}:
        kwargs.update({'compression_type': 'RICE_1', 'quantize_level': 4.0})
    if isinstance(fname, str):
        fname = os.path.expanduser(fname)

    hdunew = fits.CompImageHDU(data=data, header=header, **kwargs)
    if mask is None:
        hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew])
    else:
        hdumask = fits.CompImageHDU(data=mask.astype(np.uint8), **kwargs)
        hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, hdumask])
    hdulnew.writeto(fname, output_verify='fix')


def fits_wrap_spwX(fitsfiles, outfitsfile='output.fits'):
    fitsfiles = sorted(fitsfiles)
    nband = len(fitsfiles)
    fits_exist = []
    idx_fits_exist = []
    for sidx, fitsf in enumerate(fitsfiles):
        if os.path.exists(fitsf):
            fits_exist.append(fitsf)
            idx_fits_exist.append(sidx)
    if len(fits_exist) == 0: raise ValueError('None of the input fitsfiles exists!')
    os.system('cp {} {}'.format(fits_exist[0], outfitsfile))
    hdu0 = fits.open(outfitsfile, mode='update')
    header = hdu0[0].header
    npol, nbd, ny, nx = int(header['NAXIS4']), nband, int(header['NAXIS2']), int(header['NAXIS1'])
    data = np.zeros((npol, nbd, ny, nx))
    cfreqs = []
    for sidx, fitsf in enumerate(fits_exist):
        hdu = fits.open(fitsf)
        cfreqs.append(hdu[0].header['CRVAL3'])
        for pidx in range(npol):
            if len(hdu[0].data.shape) == 2:
                data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data
            else:
                data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data[pidx, 0, :, :]
    df = np.nanmean(np.diff(cfreqs)/np.diff(idx_fits_exist)) ## in case some of the band is missing
    header['cdelt3'] = df
    header['NAXIS3'] = nband
    header['NAXIS'] = 4
    header['CRVAL3'] = header['CRVAL3'] - df * idx_fits_exist[0]
    if os.path.exists(outfitsfile):
        os.system('rm -rf {}'.format(outfitsfile))
    fits.writeto(outfitsfile, data, header)
    print('wrapped fits writed to ' + outfitsfile)
    return


def fits_wrap_spw(fitsfiles, outfitsfile='output.fits', df=0.5e9, nband=30):
    fitsfiles = sorted(fitsfiles)
    for sidx, fitsf in enumerate(fitsfiles):
        if not os.path.exists(fitsf):
            raise ValueError('The {} does not exist!'.format(fitsf))
    os.system('cp {} {}'.format(fitsfiles[0], outfitsfile))
    hdu0 = fits.open(outfitsfile, mode='update')
    header = hdu0[0].header
    freqref = header['crval3']
    header['cdelt3'] = df
    header['NAXIS3'] = nband
    header['NAXIS'] = 4
    # if 'NAXIS4' not in header.keys():
    #     header['NAXIS4'] = npol
    npol, nbd, ny, nx = int(header['NAXIS4']), int(header['NAXIS3']), int(header['NAXIS2']), int(header['NAXIS1'])
    data = np.zeros((npol, nbd, ny, nx))
    for sidx, fitsf in enumerate(fitsfiles):
        hdu = fits.open(fitsf)
        for pidx in range(npol):
            if nband == 30:
                bdidx = int(np.round((hdu[0].header['crval3'] - freqref) / df))
            else:
                bdidx = sidx
            if len(hdu[0].data.shape) == 2:
                data[pidx, bdidx, :, :] = hdu[0].data
            else:
                data[pidx, bdidx, :, :] = hdu[0].data[pidx, 0, :, :]
    if os.path.exists(outfitsfile):
        os.system('rm -rf {}'.format(outfitsfile))
    fits.writeto(outfitsfile, data, header)
    print('wrapped fits writed to ' + outfitsfile)
    return
