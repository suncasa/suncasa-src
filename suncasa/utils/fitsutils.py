import os
import numpy as np
from astropy.io import fits
from sunpy import map as smap
import warnings

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


def header_to_xml(header):
    import xml.etree.ElementTree as ET
    # create the file structure
    tree = ET.ElementTree()
    root = ET.Element('meta')
    elem = ET.Element('fits')
    for k, v in header.items():
        child = ET.Element(k)
        if isinstance(v, bool):
            v = int(v)
        child.text = str(v)
        elem.append(child)
    root.append(elem)
    tree._setroot(root)
    return tree
    # from lxml import etree
    # tree = etree.Element("meta")
    # elem = etree.SubElement(tree,"fits")
    # for k,v in header.items():
    #     child = etree.SubElement(elem, k)
    #     if isinstance(v,bool):
    #         v = int(v)
    #     child.text = str(v)
    # return tree


def write_j2000_image(fname, data, header):
    import glymur
    jp2 = glymur.Jp2k(fname + '.tmp.jp2',
                      ((data - np.min(data)) * 256 / (np.max(data) - np.min(data))).astype(np.uint8), cratios=[20, 10])
    boxes = jp2.box
    header['wavelnth'] = header['crval3']
    header['waveunit'] = header['cunit3']
    xmlobj = header_to_xml(header)
    xmlfile = 'image.xml'
    if os.path.exists(xmlfile):
        os.system('rm -rf {}'.format(xmlfile))
    xmlobj.write(xmlfile)
    xmlbox = glymur.jp2box.XMLBox(filename='image.xml')
    boxes.insert(3, xmlbox)
    jp2_xml = jp2.wrap(fname, boxes=boxes)
    os.system('rm -rf ' + fname + '.tmp.jp2')
    os.system('rm -rf {}'.format(xmlfile))
    return True


def read_compressed_image_fits(fname, verbose=False):
    hdulist = fits.open(fname)
    for i, hdu in enumerate(hdulist):
        try:
            ndim = hdu.data.ndim
            dshape = hdu.data.shape
            header = hdu.header
            slc = [slice(None)] * ndim
            nx = header['NAXIS1']
            ny = header['NAXIS2']
            freqaxis = None
            stokaxis = None
            npol_fits = 1
            nfreq_fits = 1
            for idx in range(ndim):
                v = header['CTYPE{}'.format(idx + 1)]
                if v.startswith('FREQ'):
                    freqaxis = ndim - (idx + 1)
                    nfreq_fits = header['NAXIS{}'.format(idx + 1)]
                if v.startswith('STOKES'):
                    stokaxis = ndim - (idx + 1)
                    npol_fits = header['NAXIS{}'.format(idx + 1)]
            if freqaxis is not None:
                slc[freqaxis] = slice(0, 1)
                rfreqs = (hdu.header['CRVAL{}'.format(ndim - freqaxis)] + hdu.header[
                    'CDELT{}'.format(ndim - freqaxis)] * np.arange(
                    hdu.header['NAXIS{}'.format(ndim - freqaxis)])) / 1e9
            else:
                rfreqs = np.array([hdu.header['RESTFRQ'] / 1e9])
            if stokaxis is not None:
                slc[stokaxis] = slice(0, 1)
            hdu.data[np.isnan(hdu.data)] = 0.0
            rmap = smap.Map(np.squeeze(hdu.data[slc]), hdu.header)
            rdata = hdu.data.copy()
            rheader = hdu.header.copy()
            break
        except Exception as e:
            if verbose:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                print('skipped HDU {}'.format(i))
            rmap, ndim, npol_fits, stokaxis, rfreqs, rdata, rheader = None, None, None, None, None, None, None
    hdulist.close()
    return rmap, ndim, npol_fits, stokaxis, rfreqs, rdata, rheader


def write_compressed_image_fits(fname, data, header, mask=None, fix_invalid=True, filled_value=0.0, **kwargs):
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
        # todo: need to figure out how to keep information of the trimmed axis in the fits
        hdunew = fits.CompImageHDU(data=data, header=header, **kwargs)
        if mask is None:
            hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew])
        else:
            hdumask = fits.CompImageHDU(data=mask.astype(np.uint8), **kwargs)
            hdulnew = fits.HDUList([fits.PrimaryHDU(), hdunew, hdumask])
        hdulnew.writeto(fname, output_verify='fix')
        return 1


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
    df = np.nanmean(np.diff(cfreqs) / np.diff(idx_fits_exist))  ## in case some of the band is missing
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
