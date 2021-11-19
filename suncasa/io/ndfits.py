import os
import numpy as np
from astropy.io import fits
from sunpy import map as smap
import warnings
warnings.simplefilter("ignore")

stokesval = {'1': 'I', '2': 'Q', '3': 'U', '4': 'V', '-1': 'RR', '-2': 'LL', '-3': 'RL', '-4': 'LR', '-5': 'XX',
             '-6': 'YY', '-7': 'XY', '-8': 'YX'}


def headerfix(header, PC_coor=True):
    '''
	this code fixes the header problem of fits out from CASA 5.4+ which leads to a streched solar image.
    Setting PC_coor equal to True will reset the rotation matrix.
    '''

    keys2remove = []
    for k in header:
        if k.upper().startswith('PC'):
            if not k.upper().startswith('PC0'):
                pcidxs = k.upper().replace('PC', '')
                hd_ = 'PC0' + pcidxs
                keys2remove.append(k)
                if PC_coor:
                    pcidx0, pcidx1 = pcidxs.split('_')
                    if pcidx0 == pcidx1:
                        header[hd_] = 1.0
                    else:
                        header[hd_] = 0.0
                else:
                    header[hd_] = header[k]
    for k in keys2remove:
        header.remove(k)
    return header


def headersqueeze(header, data):
    '''
        Only 1D, 2D, or 3D images are currently supported by
        the astropy fits compression.
        This code remove single-dimensional entries from the n-dimensional
        image data array and update fits header
    '''
    dshape = data.shape
    ndim = data.ndim
    ## nonsdim: non-single-dimensional entries
    nonsdim = ndim - np.count_nonzero(np.array(dshape) == 1)
    if nonsdim > 3:
        return None, None
    else:
        keys2chng = ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CUNIT']  # ,'PC01_', 'PC02_', 'PC03_', 'PC04_']
        idx_nonsdim = 0
        for idx, dim in enumerate(dshape[::-1]):
            # if dim>1: continue
            if dim > 1:
                idx_nonsdim = idx_nonsdim + 1
            for k in keys2chng:
                k_ = '{}{}'.format(k, idx + 1)
                v = header[k_]
                header.remove(k_)
                if dim > 1:
                    k_new = '{}{}'.format(k, idx_nonsdim)
                    header[k_new] = v
                else:
                    if k is 'CTYPE' and v.startswith('STOKES'):
                        header['STOKES'] = header['CRVAL{}'.format(idx + 1)]

        idx_nonsdim1 = 0

        for idx1, dim1 in enumerate(dshape[::-1]):
            if dim1 > 1:
                idx_nonsdim1 = idx_nonsdim1 + 1
            idx_nonsdim2 = 0
            for idx2, dim2 in enumerate(dshape[::-1]):
                if dim2 > 1:
                    idx_nonsdim2 = idx_nonsdim2 + 1
                k_ = 'PC{:02d}_{:d}'.format(idx1 + 1, idx2 + 1)
                v = header[k_]
                header.remove(k_)
                if dim1 > 1 and dim2 > 1:
                    k_new = 'PC{:02d}_{:d}'.format(idx_nonsdim1, idx_nonsdim2)
                    header[k_new] = v

        header['NAXIS'] = nonsdim
        data = np.squeeze(data)
        return header, data


def get_bdinfo(freq, bw):
    """
    get band information from center frequencies and band widths.

    Parameters
    ----------
    freq : array_like
        an array of the center frequencies of all frequency bands in Hz
    bw: array_like
        an array of the band widths of all frequency bands in Hz

    Returns
    -------
    fbounds : `dict`
        A dict of band information
    """
    fghz = freq / 1.e9
    bwghz = bw / 1.e9
    bounds_lo = fghz - bwghz / 2.0
    bounds_hi = fghz + bwghz / 2.0
    bounds_all = np.hstack([bounds_lo, bounds_hi[-1]])
    fbounds = {'cfreqs': fghz, 'cfreqs_all': fghz, 'bounds_lo': bounds_lo,
               'bounds_hi': bounds_hi, 'bounds_all': bounds_all}
    return fbounds


def read(filepath, hdus=None, memmap=None, verbose=False, **kwargs):
    """
    Read a fits file.

    Parameters
    ----------
    filepath : `str`
        The fits file to be read.
    hdus: `int` or iterable
        The HDU indexes to read from the file.
    verbose: `bool`
        if verbose

    Returns
    -------
    pairs : `list`
        A list of (data, header) tuples

    Notes
    -----
    This routine reads all the HDU's in a fits file and returns a list of the
    data and a FileHeader instance for each one.

    Also all comments in the original file are concatenated into a single
    "comment" key in the returned FileHeader.
    """
    with fits.open(filepath, ignore_blank=True, memmap=memmap) as hdulist:
        if hdus is not None:
            if isinstance(hdus, int):
                hdulist = hdulist[hdus]
            elif isinstance(hdus, collections.Iterable):
                hdulist = [hdulist[i] for i in hdus]

        hdulist = fits.hdu.HDUList(hdulist)
        for h in hdulist:
            h.verify('silentfix+warn')

        meta = {}
        for i, hdu in enumerate(hdulist):
            try:
                ndim = hdu.data.ndim
                header = hdu.header
                slc = [slice(None)] * ndim
                freq_axis = None
                pol_axis = None
                _ = None
                npol = 1
                nfreq_fits = 1
                for idx in range(ndim):
                    v = header['CTYPE{}'.format(idx + 1)]
                    if v.startswith('FREQ'):
                        freq_axis = ndim - (idx + 1)
                        nfreq_fits = header['NAXIS{}'.format(idx + 1)]
                    if v.startswith('STOKES'):
                        pol_axis = ndim - (idx + 1)
                        npol = header['NAXIS{}'.format(idx + 1)]
                if freq_axis is not None:
                    slc[freq_axis] = slice(0, 1)
                    meta['ref_cfreqs'] = (hdu.header['CRVAL{}'.format(ndim - freq_axis)] + hdu.header[
                        'CDELT{}'.format(ndim - freq_axis)] * np.arange(
                        hdu.header['NAXIS{}'.format(ndim - freq_axis)]))
                    meta['ref_freqdelts'] = np.ones_like(meta['ref_cfreqs']) \
                                            * hdu.header['CDELT{}'.format(ndim - freq_axis)]
                else:
                    meta['ref_cfreqs'] = np.array([hdu.header['RESTFRQ']])

                if pol_axis is not None:
                    slc[pol_axis] = slice(0, 1)
                    meta['pol_idxs'] = (hdu.header['CRVAL{}'.format(ndim - pol_axis)] + hdu.header[
                        'CDELT{}'.format(ndim - pol_axis)] * np.arange(
                        hdu.header['NAXIS{}'.format(ndim - pol_axis)]))
                    meta['pol_names'] = [stokesval['{0:d}'.format(int(p))] for p in meta['pol_idxs']]
                if _ is not None:
                    slc[_] = slice(0, 1)
                hdu.data[np.isnan(hdu.data)] = 0.0
                rmap = smap.Map(np.squeeze(hdu.data[slc]), hdu.header)
                data = hdu.data.copy()
                meta['header'] = hdu.header.copy()
                meta['refmap'] = rmap  # this is a sunpy map of the first slice
                meta['naxis'] = ndim
                meta['hgln_axis'] = ndim - 1  # solar X
                meta['nx'] = header['NAXIS1']
                meta['hglt_axis'] = ndim - 2  # solar Y
                meta['ny'] = header['NAXIS2']
                meta['freq_axis'] = freq_axis
                meta['nfreq'] = nfreq_fits
                meta['pol_axis'] = pol_axis
                meta['npol'] = npol
                break
            except Exception as e:
                if verbose:
                    if hasattr(e, 'message'):
                        print(e.message)
                    else:
                        print(e)
                    print('skipped HDU {}'.format(i))
                meta, data = {}, None

        # Check if an additional frequency axis exists. If so, they should be in the last hdu.
        if hasattr(hdulist[-1].data, 'cfreqs'):
            if verbose:
                print('FITS file contains an additional frequency axis. '
                      'Update the frequency information in cfreqs and cdelts. '
                      'Ignore the original version with equal spacings.')
            meta['ref_cfreqs'] = np.array(hdulist[-1].data['cfreqs'])
            meta['ref_freqdelts'] = np.array(hdulist[-1].data['cdelts'])
        else:
            if verbose:
                print('FITS file does not have an additional frequency axis. '
                      'Use the original version with equal spacings.')
        hdulist.close()

        return meta, data


def write(fname, data, header, mask=None, fix_invalid=True, filled_value=0.0, **kwargs):
    """
    Take a data header pair and write a compressed FITS file.
    Caveat: only 1D, 2D, or 3D images are currently supported by Astropy fits compression.
    To be compressed, the image data array (n-dimensional) must have
    at least n-3 single-dimensional entries.

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

    dshape = data.shape
    dim = data.ndim
    if dim - np.count_nonzero(np.array(dshape) == 1) > 3:
        return 0
    else:
        if fix_invalid:
            data[np.isnan(data)] = filled_value
        if kwargs is {}:
            kwargs.update({'compression_type': 'RICE_1', 'quantize_level': 4.0})
        if isinstance(fname, str):
            fname = os.path.expanduser(fname)

        header, data = headersqueeze(header, data)
        hdunew = fits.CompImageHDU(data=data, header=header, **kwargs)
        if mask is None:
            hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew])
        else:
            hdumask = fits.CompImageHDU(data=mask.astype(np.uint8), **kwargs)
            hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, hdumask])
        hdulnew.writeto(fname, output_verify='fix')
        return 1


def wrap(fitsfiles, outfitsfile='output.fits', docompress=False, mask=None, fix_invalid=True, filled_value=0.0,
         **kwargs):
    if len(fitsfiles)<=1:
        print('There is only one files in the fits file list. wrap is aborted!')
        return 0
    else:
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
        cdelts = []
        cfreqs = []
        for sidx, fitsf in enumerate(fits_exist):
            hdu = fits.open(fitsf)
            cdelts.append(hdu[0].header['CDELT3'])
            cfreqs.append(hdu[0].header['CRVAL3'])
            for pidx in range(npol):
                if len(hdu[0].data.shape) == 2:
                    data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data
                else:
                    data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data[pidx, 0, :, :]
        cfreqs = np.array(cfreqs)
        cdelts = np.array(cdelts)
        df = np.nanmean(np.diff(cfreqs) / np.diff(idx_fits_exist))  ## in case some of the band is missing
        header['cdelt3'] = df
        header['NAXIS3'] = nband
        header['NAXIS'] = 4
        header['CRVAL3'] = header['CRVAL3'] - df * idx_fits_exist[0]
        if os.path.exists(outfitsfile):
            os.system('rm -rf {}'.format(outfitsfile))

        col1 = fits.Column(name='cfreqs', format='E', array=cfreqs)
        col2 = fits.Column(name='cdelts', format='E', array=cdelts)
        tbhdu = fits.BinTableHDU.from_columns([col1, col2])
        if docompress:
            dshape = data.shape
            dim = data.ndim
            if dim - np.count_nonzero(np.array(dshape) == 1) > 3:
                pass
            else:
                if fix_invalid:
                    data[np.isnan(data)] = filled_value
                if kwargs is {}:
                    kwargs.update({'compression_type': 'RICE_1', 'quantize_level': 4.0})
                if isinstance(outfitsfile, str):
                    outfitsfile = os.path.expanduser(outfitsfile)

                header, data = headersqueeze(header, data)
                hdunew = fits.CompImageHDU(data=data, header=header, **kwargs)
                if mask is None:
                    hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, tbhdu])
                else:
                    hdumask = fits.CompImageHDU(data=mask.astype(np.uint8), **kwargs)
                    hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, tbhdu, hdumask])
                hdulnew.writeto(outfitsfile, output_verify='fix')
                return 1

        hdulnew = fits.HDUList([fits.PrimaryHDU(data=data, header=header), tbhdu])
        hdulnew.writeto(outfitsfile)
        print('wrapped fits written as ' + outfitsfile)
        return 1
