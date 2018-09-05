import os
import numpy as np
import astropy.units as u
from astropy.io import fits


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
        hdu = fits.open(fitsf, mode='update')
        for pidx in range(npol):
            if nband == 30:
                bdidx = int(np.round((hdu[0].header['crval3'] - freqref) / df))
            else:
                bdidx = sidx
            if len(hdu[0].data.shape) == 2:
                data[pidx, bdidx, :, :] = hdu[0].data
            else:
                data[pidx, bdidx, :, :] = hdu[0].data[0, 0, :, :]
    if os.path.exists(outfitsfile):
        os.system('rm -rf {}'.format(outfitsfile))
    fits.writeto(outfitsfile, data, header)
    print('wrapped fits writed to ' + outfitsfile)
    return
