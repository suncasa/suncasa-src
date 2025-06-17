__all__ = ['Dspec']

import os
import struct

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from copy import copy

from ..casa_compat import import_casatools, import_casatasks
tasks = import_casatasks('split','hanningsmooth')
split = tasks.get('split')
hanningsmooth = tasks.get('hanningsmooth')

tools = import_casatools(['tbtool', 'mstool', 'qatool'])

tbtool = tools['tbtool']
mstool = tools['mstool']
qatool = tools['qatool']
tb = tbtool()
ms = mstool()
qa = qatool()

from matplotlib.dates import AutoDateFormatter, AutoDateLocator, num2date
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .sources.lwa import rebin1d, rebin2d

stokestype = [
    'Undefined',
    # undefined value = 0
    'I',
    'Q',
    'U',
    # standard stokes parameters
    'V',
    #
    'RR',
    'RL',
    'LR',
    # circular correlation products
    'LL',
    #
    'XX',
    'XY',
    'YX',
    # linear correlation products
    'YY',
    #
    'RX',
    'RY',
    'LX',
    'LY',
    'XR',
    'XL',
    'YR',
    # mixed correlation products
    'YL',
    #
    'PP',
    'PQ',
    'QP',
    # general quasi-orthogonal correlation products
    'QQ',
    #
    'RCircular',
    'LCircular',
    # single dish polarization types
    'Linear',
    # Polarized intensity ((Q^2+U^2+V^2)^(1/2))
    'Ptotal',
    # Linearly Polarized intensity ((Q^2+U^2)^(1/2))
    'Plinear',
    # Polarization Fraction (Ptotal/I)
    'PFtotal',
    # Linear Polarization Fraction (Plinear/I)
    'PFlinear',
    # Linear Polarization Angle (0.5 arctan(U/Q)) (in radians)
    'Pangle']

stokesenum = {}
for k, v in zip(range(len(stokestype)), stokestype):
    stokesenum[k] = v


def calc_bkg_dspec(spec, tim, bkgtim, interp_method='linear'):
    """
    Calculate the background spectrum.

    This function computes the background spectrum from a given spectrogram
    over specified time intervals. The background spectrum can be interpolated
    using various methods.

    :param spec: A 2D array of shape (nfreq, ntime) representing the spectrogram.
    :type spec: ndarray
    :param tim: A time object with the time axis of shape (ntime,).
    :type tim: astropy.time.Time
    :param bkgtim: Background time interval(s). It can be a single time interval (list of two strings)
        or multiple time intervals (list of lists of two strings each).
        Example: ['2024-06-11T19:41:00', '2024-06-11T20:40:58']
                 [['2024-06-11T19:41:00', '2024-06-11T20:40:58'], ['2024-06-11T21:41:00', '2024-06-11T22:40:58']]
    :type bkgtim: list or list of lists
    :param interp_method: The method of interpolation to use. Default is 'linear'.
        Can be any method supported by scipy.interpolate.interp1d.
    :type interp_method: str, optional
    :return: A 2D array of shape (nfreq, ntime) representing the interpolated background spectrum.
    :rtype: ndarray

    :example:
    >>> from astropy.time import Time
    >>> import numpy as np
    >>> spec = np.random.random((100, 200))
    >>> tim = Time(['2024-06-11T19:00:00']*200, scale='utc')
    >>> bkgtim = ['2024-06-11T19:41:00', '2024-06-11T20:40:58']
    >>> calc_bkg_dspec(spec, tim, bkgtim)
    array([[0.52203423, 0.52203423, 0.52203423, ..., 0.52203423, 0.52203423, 0.52203423],
           [0.48498126, 0.48498126, 0.48498126, ..., 0.48498126, 0.48498126, 0.48498126],
           ...,
           [0.49162387, 0.49162387, 0.49162387, ..., 0.49162387, 0.49162387, 0.49162387]])
    """
    from scipy.interpolate import interp1d
    from astropy.time import Time

    # Ensure tim is an array of julian dates
    tim_jd = tim.jd
    bkgtim_jd = Time(bkgtim).jd

    if isinstance(bkgtim[0], str):
        # Single time interval
        tidx_bkg, = np.where((tim_jd >= bkgtim_jd[0]) & (tim_jd <= bkgtim_jd[1]))
        spec_bkg = np.nanmean(spec[:, tidx_bkg], axis=1)[:, np.newaxis]
    else:
        # Multiple time intervals
        all_spec_bkg = []
        for bkgtjd in bkgtim_jd:
            tidx_bkg, = np.where((tim_jd >= bkgtjd[0]) & (tim_jd <= bkgtjd[1]))
            all_spec_bkg.append(np.nanmean(spec[:, tidx_bkg], axis=1))

        # Combine all background spectra
        all_spec_bkg = np.array(all_spec_bkg)  # Shape will be (nfreq, n_intervals)
        bkg_times = np.nanmean(bkgtim_jd,axis=1)

        interp_func = interp1d(bkg_times, all_spec_bkg, kind=interp_method, axis=0, fill_value="extrapolate")
        spec_bkg = interp_func(tim_jd).T  # Shape will be (nfreq, ntime)

    return spec_bkg

class Dspec:
    """
    A class to handle dynamic spectra from radio observations.

    Attributes
    ----------
    data : numpy.ndarray
        The dynamic spectrum data array.
    time_axis : astropy.time.Time
        Time axis of the dynamic spectrum.
    freq_axis : numpy.ndarray
        Frequency axis of the dynamic spectrum in Hz.
    telescope : str
        Name of the telescope used for the observation.
    observatory : str
        Name of the observatory.
    t_label : str
        Label for the time axis.
    f_label : str
        Label for the frequency axis.
    bl : str
        Baseline information.
    uvrange : str
        UV range information.
    pol : str
        Polarization information.
    spec_unit : str, default 'sfu'
        Unit of the dynamic spectrum data ('sfu', 'Jy', or 'K').

    Methods
    -------
    __init__(fname=None, specfile=None, **kwargs):
        Initializes the Dspec object by reading a FITS file or a saved numpy array.

    read(fname, source=None, *args, **kwargs):
        Reads dynamic spectrum data from a file.

    tofits(fitsfile=None, specdata=None, **kwargs):
        Writes the dynamic spectrum data to a FITS file.

    get_dspec(fname=None, **kwargs):
        Extracts dynamic spectrum data from a measurement set.

    wrt_dspec(specfile=None, specdat=None):
        Writes the dynamic spectrum data to a binary file.

    rd_dspec(specdata, **kwargs):
        Reads dynamic spectrum data from a numpy file.

    concat_dspec(specfiles, outfile=None, savespec=False):
        Concatenates multiple dynamic spectrum files along the time axis.

    peek(*args, **kwargs):
        Plots the dynamic spectrum on the current axes.

    plot(pol='I', **kwargs):
        Plots the dynamic spectrum for a given polarization.

    """
    data = None
    time_axis = None
    freq_axis = None
    telescope = None
    observatory = None
    t_label = None
    f_label = None
    bl = None
    uvrange = None
    pol = None
    spec_unit = 'sfu'

    def __init__(self, fname=None, specfile=None, bl='', uvrange='', field='', scan='',
                 datacolumn='data', domedian=False, timeran=None, spw=None, timebin='0s', regridfreq=False,
                 fillnan=None, verbose=False, usetbtool=True, ds_normalised=False):
        """
              Initializes the Dspec object by reading a FITS file or a saved numpy array.

              Parameters
              ----------
              fname : str, optional
                  The file name of the FITS file or the measurement set to be read.
              specfile : str, optional
                  The file name of the saved numpy array (.npz file) that contains the dynamic spectrum data.
              bl : str, optional
                  Baseline selection criteria.
              uvrange : str, optional
                  UV range selection criteria.
              field : str, optional
                  Field selection criteria.
              scan : str, optional
                  Scan selection criteria.
              datacolumn : str, optional
                  Specifies which data column ('data' or 'corrected') to read from the measurement set.
              domedian : bool, optional
                  If True, performs a median over the baseline axis.
              timeran : str or list of str, optional
                  Time range to select data.
              spw : str, optional
                  Spectral window selection criteria.
              timebin : str, optional
                  Time averaging interval.
              regridfreq : bool, optional
                  If True, regrids the frequency axis to have a uniform channel width.
              fillnan : float or int, optional
                  Value to fill in for flagged or NaN values.
              verbose : bool, optional
                  If True, prints detailed information during operation.
              usetbtool : bool, optional
                  If True, uses CASA table tools for data extraction.
              ds_normalised : bool, optional
                  If True, normalizes the dynamic spectrum by the median value.

              """
        if fname:
            if isinstance(fname, str):
                if os.path.exists(fname):
                    if os.path.isdir(fname):
                        self.get_dspec(fname, specfile=specfile, bl=bl, uvrange=uvrange, field=field,
                                       scan=scan, datacolumn=datacolumn, domedian=domedian, timeran=timeran, spw=spw,
                                       timebin=timebin, regridfreq=regridfreq, fillnan=fillnan, verbose=verbose,
                                       usetbtool=usetbtool, ds_normalised=ds_normalised)
                    else:
                        self.read(fname)
                else:
                    raise FileNotFoundError(f"The file {fname} does not exist. Please check your file path.")
            else:
                self.read(fname)

    def read(self, fname, source=None, *args, **kwargs):
        """
        Reads dynamic spectrum data from a file.

        Parameters
        ----------
        fname : str or list of str
            The file name to read the dynamic spectrum data from.
        source : str, optional
            Specifies the data source ('fits', 'suncasa', or 'lwa') to determine the appropriate reader.

        Additional parameters are passed to the specific reader function based on the source.

        """
        ##TODO The existing implementation of the mapping between extensions and instruments requires refinement and sophistication.
        # Explore potential optimization strategies to improve this process.
        from astropy.io import fits

        if type(fname) is str:
            _known_extensions = {
                ('npz'): 'suncasa',
            }
            for extension, readername in _known_extensions.items():
                if fname.lower().endswith(extension):
                    source = readername
            if source is None and (not fname.lower().startswith('eovsa')):
                raise ValueError(f"The filetype provided ({os.path.basename(fname)}) is not supported")

            if source.lower() == 'eovsa' or fname.lower().startswith('eovsa'):
                if fname.endswith('.fts') or fname.endswith('.fits'):
                    from .sources import eovsa
                    s = eovsa.get_dspec(fname, doplot=False)
                    self.data = s['spectrogram']
                    self.time_axis = Time(s['time_axis'], format='mjd')
                    self.freq_axis = s['spectrum_axis'] * 1e9
                    self.spec_unit = 'sfu'
                    self.telescope = 'EOVSA'
                    self.observatory = 'OVRO'


            if source.lower() == 'suncasa':
                spec, tim, freq, bl, pol, spec_unit = self.rd_dspec(fname, spectype='amp', spec_unit='jy')
                self.data = spec
                self.time_axis = Time(tim / 24 / 3600, format='mjd')
                self.freq_axis = freq
                self.bl = bl
                self.pol = pol
                self.spec_unit = spec_unit
                self.telescope = ''
                self.observatory = ''

            if source.lower() == 'lwa':
                if fname.endswith('.fits'):
                    hdu = fits.open(fname) 
                    self.data = hdu[0].data
                    tim = hdu[2].data
                    tmjd = np.array(tim['mjd']) + np.array(tim['time']) / 24. / 3600 / 1000
                    self.time_axis = Time(tmjd, format='mjd')
                    self.freq_axis = hdu[1].data['sfreq'] * 1e9
                    self.pol = [hdu[0].header['POLARIZA']]
                    self.spec_unit = 'sfu'
                    self.telescope = 'LWA'
                    self.observatory = 'OVRO'
                else:
                    from .sources import lwa
                    spec, tim, freq, pol, calfac_x, calfac_y, bkg_flux = lwa.read_data(fname, **kwargs)
                    self.data = spec
                    self.time_axis = Time(tim, format='mjd')
                    self.freq_axis = freq
                    self.pol = pol
                    self.calfac_x = calfac_x # correction factor for X pol, same shape as freq
                    self.calfac_y = calfac_y # correction factor for Y pol, same shape as freq
                    self.bkg_flux = bkg_flux
                    self.spec_unit = 'sfu'
                    self.telescope = 'LWA'
                    self.observatory = 'OVRO'

            if source.lower() =='general' and (fname.endswith('.fits') or fname.endswith('.fts') or fname.endswith('.fit')):
                hdu = fits.open(fname)
                self.data = hdu[0].data
                tim = hdu[2].data
                tmjd = np.array(tim['mjd']) + np.array(tim['time']) / 24. / 3600 / 1000
                self.time_axis = Time(tmjd, format='mjd')
                self.freq_axis = hdu[1].data['sfreq'] * 1e9
                self.pol = [hdu[0].header['POLARIZA']]
                if 'BUNIT' in hdu[0].header:
                    self.spec_unit = hdu[0].header['BUNIT']
                if 'BNAME' in hdu[0].header:
                    self.spec_name = hdu[0].header['BNAME']
                self.telescope = hdu[0].header['TELESCOP']
                self.observatory = ''

        elif source.lower() == 'lwa' and type(fname) is list:
            from .sources import lwa
            spec, tim, freq, pol, calfac_x, calfac_y, bkg_flux = lwa.read_data(fname, **kwargs)
            self.data = spec
            self.time_axis = Time(tim, format='mjd')
            self.freq_axis = freq
            self.pol = pol
            self.calfac_x = calfac_x # correction factor for X pol, same shape as freq
            self.calfac_y = calfac_y # correction factor for Y pol, same shape as freq
            self.bkg_flux = bkg_flux
            self.spec_unit = 'sfu'
            self.telescope = 'LWA'
            self.observatory = 'OVRO'
        elif source.lower()=='ecallisto':
            from .sources import ecallisto
            s = ecallisto.get_dspec(fname, doplot=False)
            self.data = s['spectrogram']
            self.time_axis = Time(s['time_axis'], format='mjd')
            self.freq_axis = s['spectrum_axis'] * 1e6
            self.spec_unit = 'sfu'
            self.telescope = 'e-callisto'
            self.observatory = 'e-callisto'
        else:
            raise ValueError(f"Unsupported data source or type provided: {source}")

    def tofits(self, fitsfile=None, specdata=None, spectype='amp', spec_unit='jy',
               telescope='EOVSA', observatory='Owens Valley Radio Observatory', observer='EOVSA Team'):
        """
        Writes the dynamic spectrum data to a FITS file.

        This method exports the stored dynamic spectrum data into a new FITS file, including relevant metadata and optional input spectrum data.

        Parameters
        ----------
        fitsfile : str, optional
            Path and name of the output FITS file. If not specified, a default name will be generated.
        specdata : dict, optional
            Input dictionary containing dynamic spectrum data and metadata. If not provided, the method uses the data stored in the Dspec object.
        spectype : str, optional
            Specifies the type of the spectrum to be saved ('amp' for amplitude, 'pha' for phase). Default is 'amp'.
        spec_unit : str, optional
            Specifies the unit of the spectrum data ('jy' for Jansky, 'sfu' for Solar Flux Units, 'k' for Kelvin). Default is 'jy'.
        telescope : str, optional
            Name of the telescope with which the data was obtained. Default is 'EOVSA'.
        observatory : str, optional
            Name of the observatory. Default is 'Owens Valley Radio Observatory'.
        observer : str, optional
            Name of the observer or team that made the observation. Default is 'EOVSA Team'.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If the specified `fitsfile` path does not exist.
        ValueError
            If `specdata` is provided but does not contain the required keys.

        Examples
        --------
        >>> dspec = Dspec()
        >>> dspec.tofits(fitsfile='output.fits', specdata=my_specdata, spec_unit='sfu', observer='Sun Observer')

        Note
        ----
        The method can directly utilize the dynamic spectrum data stored within the `Dspec` object if `specdata` is not provided. Ensure the `Dspec` object has been properly initialized with dynamic spectrum data before calling this method without `specdata`.
        """
        from astropy.io import fits

        if hasattr(self, 'data'):
            spec = self.data
            tim = self.time_axis
            freqghz = self.freq_axis / 1e9

        if specdata:
            try:
                spec, tim, freq, bl, pol, spec_unit = self.rd_dspec(specdata, spectype=spectype, spec_unit=spec_unit)
                tim = Time(tim / 24 / 3600, format='mjd')
            except:
                print('Something is wrong with the input specdata.')

        hdu = fits.PrimaryHDU(spec)
        # Set up the extensions: sfreq, ut
        col1 = fits.Column(name='sfreq', format='E', array=freqghz)
        cols1 = fits.ColDefs([col1])
        tbhdu1 = fits.BinTableHDU.from_columns(cols1)
        tbhdu1.name = 'SFREQ'

        # Split up mjd into days and msec
        date_obs = tim[0].isot
        date_end = tim[-1].isot
        ut = tim.mjd
        ut_int = ut.astype(int)
        ut_msec = 1000.0 * 86400.0 * (ut - ut_int)
        ut_ms1 = ut_msec.astype(int)

        # J is the format code for a 32 bit integer, who would have thought
        # http://astropy.readthedocs.org/en/latest/io/fits/usage/table.html
        col3 = fits.Column(name='mjd', format='J', array=ut_int)
        col4 = fits.Column(name='time', format='J', array=ut_ms1)

        cols3 = fits.ColDefs([col3, col4])
        tbhdu3 = fits.BinTableHDU.from_columns(cols3)
        tbhdu3.name = 'UT'

        if hasattr(self, 'calfac_x') and hasattr(self, 'calfac_y') and hasattr(self, 'bkg_flux'):
            # add calfac_x _y to FREQ_CAL
            col5 = fits.Column(name='calfac_x', format='E', array=self.calfac_x)
            col6 = fits.Column(name='calfac_y', format='E', array=self.calfac_y)
            col7 = fits.Column(name='bkg_flux', format='E', array=self.bkg_flux)
            cols5 = fits.ColDefs([col5, col6, col7])
            tbhdu5 = fits.BinTableHDU.from_columns(cols5)
            tbhdu5.name = 'FREQ_CAL'
        else:
            tbhdu5 = fits.BinTableHDU.from_columns([])
            tbhdu5.name = 'FREQ_CAL'
            
        # create an HDUList object to put in header information
        hdulist = fits.HDUList([hdu, tbhdu1, tbhdu3, tbhdu5])

        # primary header
        prihdr = hdulist[0].header
        prihdr.set('FILENAME', os.path.basename(fitsfile))
        prihdr.set('ORIGIN', 'NJIT', 'Location where file was made')
        prihdr.set('DATE', Time.now().isot, 'Date when file was made')
        prihdr.set('OBSERVER', observer, 'Who to appreciate/blame')
        prihdr.set('TELESCOP', telescope, observatory)
        prihdr.set('OBJ_ID', 'SUN', 'Object ID')
        prihdr.set('TYPE', 1, 'Flare Spectrum')
        prihdr.set('DATE_OBS', date_obs, 'Start date/time of observation')
        prihdr.set('DATE_END', date_end, 'End date/time of observation')
        prihdr.set('FREQMIN', min(freqghz), 'Min freq in observation (GHz)')
        prihdr.set('FREQMAX', max(freqghz), 'Max freq in observation (GHz)')
        prihdr.set('XCEN', 0.0, 'Antenna pointing in arcsec from Sun center')
        prihdr.set('YCEN', 0.0, 'Antenna pointing in arcsec from Sun center')
        prihdr.set('POLARIZA', 'I', 'Polarizations present')
        prihdr.set('RESOLUTI', 0.0, 'Resolution value')
        if spec_unit.lower() in ['sfu', 'jy', 'mjy']:
            prihdr.set('BTYPE', 'flx', 'Flux density')
            prihdr.set('BNAME', 'Flux Density')
        elif spec_unit.lower() in ['k']:
            prihdr.set('BTYPE', 'Tb', 'Brightness temperature')
            prihdr.set('BNAME', 'Brightness Temperature')
        prihdr.set('BUNIT', spec_unit, 'Data units')
        # Write the file
        hdulist.writeto(fitsfile, overwrite=True)

    def get_dspec(self, fname=None, specfile=None, bl='', uvrange='', field='', scan='',
                  datacolumn='data',
                  domedian=False, timeran=None, spw=None, timebin='0s', regridfreq=False,
                  hanning=False,
                  applyflag=True, fillnan=None, verbose=False,
                  usetbtool=True, ds_normalised=False):
        if fname.endswith('/'):
            fname = fname[:-1]
        msfile = fname
        if not spw:
            spw = ''
        if not timeran:
            timeran = ''
        if domedian:
            if not uvrange:
                uvrange = '0.2~0.8km'
        if not bl:
            bl = ''
        else:
            uvrange = ''
        # Open the ms and plot dynamic spectrum
        if verbose:
            print('Splitting selected data...')

        if usetbtool:
            if datacolumn.lower() == 'data':
                datacol = 'DATA'
                print('Using DATA column')
            if datacolumn.lower() == 'corrected':
                datacol = 'CORRECTED_DATA'
                print('Using CORRECTED_DATA column')
            if verbose:
                print('using table tool to extract the data')
            try:
                tb.open(fname + '/POLARIZATION')
                corrtype = tb.getcell('CORR_TYPE', 0)
                pol = [stokesenum[p] for p in corrtype]
                tb.close()
            except:
                pol = []

            if hanning:
                hanningsmooth(vis=fname, datacolumn='data', field=field, outputvis=fname + '.tmpms')
                fname = fname + '.tmpms'
            tb.open(fname)
            spwtb = tbtool()
            spwtb.open(fname + '/SPECTRAL_WINDOW')
            ptb = tbtool()
            ptb.open(fname + '/POLARIZATION')
            npol = ptb.getcol('NUM_CORR', 0, 1)[0]
            ptb.close()

            ms.open(fname)
            # ms.selectinit(datadescid=0)
            spwlist = []
            mdata = ms.metadata()
            # determine how many unique spws
            # now it ignores the spw parameter and use all of them
            spws_flds = mdata.spwsforfields()
            spws_ = []
            for k in spws_flds.keys():
                spws_.append(spws_flds[k])
            spws_unq = np.unique(spws_)
            nspw = len(spws_unq)
            nbaselines = mdata.nbaselines()
            nantennas = mdata.nantennas()
            scannumbers = mdata.scannumbers()
            spw_nfrq = []  # List of number of frequencies in each spw
            for i in spws_unq:
                spw_nfrq.append(mdata.nchan(i))
            spw_nfrq = np.array(spw_nfrq)
            nf = np.sum(spw_nfrq)
            smry = mdata.summary()

            scanids = sorted(smry['observationID=0']['arrayID=0'].keys())
            ## todo find a way to determine the target fieldID if multi-fields exist
            smryscan0 = smry['observationID=0']['arrayID=0'][scanids[0]]
            fieldid = ''
            for k in smryscan0.keys():
                if k.startswith('fieldID='):
                    fieldid = k
                    break
            nbl = int(smryscan0[fieldid]['0']['nrows'] / nspw)
            if nbl == nbaselines + nantennas:
                hasautocorr = True
            elif nbl == nbaselines:
                hasautocorr = False
            else:
                raise (ValueError('The baseline number is not correct.'))

            antmask = []
            if uvrange != '' or bl != '':
                ms.open(fname)
                ms.selectinit(datadescid=0, reset=True)
                mdata = ms.metadata()
                antlist = mdata.antennaids()
                smry = mdata.summary()
                mdata.done()
                staql = {'uvdist': uvrange, 'baseline': bl, 'spw': spw, 'field': field, 'scan': scan,
                         'timerange': timeran}
                ### todo the selection only works for uvrange and bl. To make the selection of other items works,
                a = ms.msselect(staql)
                mdata = ms.metadata()
                baselines = mdata.baselines()
                if hasautocorr:
                    for lidx, l in enumerate(antlist):
                        antmask.append(baselines[l][antlist[lidx:]])
                else:
                    for lidx, l in enumerate(antlist):
                        antmask.append(baselines[l][antlist[lidx + 1:]])
                antmask = np.hstack(antmask)
                mdata.done()
                ms.close()
                ms.open(fname)

            scan_ntimes = []  # List of number of times in each scan
            nrows = []
            for s, scanid in enumerate(scanids):
                smryscan = smry['observationID=0']['arrayID=0'][scanid]
                nrow = next(v for k, v in smryscan.items() if 'fieldID=' in k)['nrows']
                nrows.append(nrow)
                scan_ntimes.append(nrow / nspw / nbl)
            scan_ntimes = np.array(scan_ntimes)
            scan_ntimes_integer = scan_ntimes.astype(int)
            if len(np.where(scan_ntimes % scan_ntimes_integer != 0)[0]) != 0:
                # if True:
                scan_ntimes = []  # List of number of times in each scan
                for s, scanid in enumerate(scanids):
                    nrows_scan = []  ## get the nrows for each time. They are not always the SAME!
                    smryscan = smry['observationID=0']['arrayID=0'][scanid]
                    # for k, v in smryscan['fieldID=0'].items():
                    #     if isinstance(v,dict):
                    #         nrows_scan.append(v['nrows'])
                    nrows_scan.append(next(v for k, v in smryscan.items() if 'fieldID=' in k)['nrows'])
                    scan_ntimes.append(nrows[s] / max(set(nrows_scan)))
                    # scan_ntimes.append(
                    #     len(smry['observationID=0']['arrayID=0'][scanids[scannumber]]['fieldID=0'].keys()) - 6)
                scan_ntimes = np.array(scan_ntimes).astype(int)
            else:
                scan_ntimes = scan_ntimes_integer

            nt = np.sum(scan_ntimes)
            times = tb.getcol('TIME')
            if times[nbl] - times[0] != 0:
                # This is frequency/scan sort order
                order = 'f'
            elif times[nbl * nspw - 1] - times[0] != 0:
                # This is time sort order
                order = 't'
            else:
                order = 'f'

            freq = np.zeros(nf, float)
            times = np.zeros(nt, float)
            if verbose:
                print("npol, nf, nt, nbl:", npol, nf, nt, nbl)
            if order == 't':
                specamp = np.zeros((npol, nf, nbl, nt), complex)
                flagf = np.zeros((npol, nf, nbl, nt), int)
                for j in range(nt):
                    fptr = 0
                    # Loop over spw
                    for i, sp in enumerate(spws_unq):
                        # Get channel frequencies for this spw (annoyingly comes out as shape (nf, 1)
                        cfrq = spwtb.getcol('CHAN_FREQ', sp, 1)[:, 0]
                        if j == 0:
                            # Only need this the first time through
                            spwlist += [i] * len(cfrq)
                        if i == 0:
                            times[j] = tb.getcol('TIME', nbl * (i + nspw * j), 1)  # Get the time
                        spec_ = tb.getcol(datacol, nbl * (i + nspw * j), nbl)  # Get complex data for this spw
                        flag = tb.getcol('FLAG', nbl * (i + nspw * j), nbl)  # Get flags for this spw
                        nfrq = len(cfrq)
                        # Apply flags
                        if applyflag:
                            if type(fillnan) in [int, float]:
                                spec_[flag] = float(fillnan)
                            else:
                                spec_[flag] = 0.0
                        # Insert data for this spw into larger array
                        specamp[:, fptr:fptr + nfrq, :, j] = spec_
                        flagf[:, fptr:fptr + nfrq, :, j] = flag
                        freq[fptr:fptr + nfrq] = cfrq
                        fptr += nfrq
            else:
                specf = np.zeros((npol, nf, nt, nbl), complex)  # Array indexes are swapped
                # flagf = np.zeros((npol, nf, nt, nbl), int)  # Array indexes are swapped
                # ant1 = np.zeros((nt, nbl), int)  # Array indexes are swapped
                # ant2 = np.zeros((nt, nbl), int)  # Array indexes are swapped
                iptr = 0
                for j, scanid in enumerate(scanids):
                    # Loop over scans
                    s = scan_ntimes[j]
                    s1 = np.sum(scan_ntimes[:j])  # Start time index
                    s2 = np.sum(scan_ntimes[:j + 1])  # End time index
                    if verbose:
                        print('=======Filling up scan #{0:d} frm:{1:d}-{2:d}======='.format(j, s1, s2))
                    for i, sp in enumerate(spws_unq):
                        # Loop over spectral windows
                        f = spw_nfrq[i]
                        f1 = np.sum(spw_nfrq[:i])  # Start freq index
                        f2 = np.sum(spw_nfrq[:i + 1])  # End freq index
                        if verbose:
                            print('Filling up spw #{0:d} chn:{1:d}--{2:d}'.format(i, f1, f2))
                        spec_ = tb.getcol(datacol, iptr, nbl * s)
                        flag = tb.getcol('FLAG', iptr, nbl * s)
                        if j == 0:
                            cfrq = spwtb.getcol('CHAN_FREQ', sp, 1)[:, 0]
                            freq[f1:f2] = cfrq
                            spwlist += [i] * len(cfrq)
                        times[s1:s2] = tb.getcol('TIME', iptr, nbl * s).reshape(s, nbl)[:, 0]  # Get the times
                        # Apply flags
                        if applyflag:
                            if type(fillnan) in [int, float]:
                                spec_[flag] = float(fillnan)
                            else:
                                spec_[flag] = 0.0
                        # Insert data for this spw into larger array
                        specf[:, f1:f2, s1:s2] = spec_.reshape(npol, f, s, nbl)
                        # flagf[:, f1:f2, s1:s2] = flag.reshape(npol, f, s, nbl)
                        # if i==0:
                        #     ant1[s1:s2] = ant1_.reshape(s, nbl)
                        #     ant2[s1:s2] = ant2_.reshape(s, nbl)
                        # Swap the array indexes back to the desired order
                        # except:
                        #     print('error processing spw {}'.format(i))
                        iptr += nbl * s
                specamp = np.swapaxes(specf, 2, 3)
                # flagf = np.swapaxes(flagf, 2, 3)
            # if applyflag:
            #     specamp = np.ma.masked_array(specamp, flagf)
            tb.close()
            spwtb.close()
            ms.close()
            if len(antmask) > 0:
                specamp = specamp[:, :, np.where(antmask)[0], :]
            (npol, nfreq, nbl, ntim) = specamp.shape
            tim = times
            if hanning:
                os.system('rm -rf {}'.format(fname))
        else:
            # Open the ms and plot dynamic spectrum
            if verbose:
                print('Splitting selected data...')
            vis_spl = './tmpms.splitted'
            if os.path.exists(vis_spl):
                os.system('rm -rf ' + vis_spl)

            # split(vis=msfile, outputvis=vis_spl, timerange=timeran, antenna=bl, field=field, scan=scan, spw=spw,
            #       uvrange=uvrange, timebin=timebin, datacolumn=datacolumn)

            try:
                split(vis=msfile, outputvis=vis_spl, datacolumn=datacolumn, timerange=timeran, spw=spw, antenna=bl,
                      field=field, scan=scan, uvrange=uvrange, timebin=timebin)
            except:
                ms.open(msfile, nomodify=True)
                ms.split(outputms=vis_spl, whichcol=datacolumn, time=timeran, spw=spw, baseline=bl, field=field,
                         scan=scan,
                         uvrange=uvrange, timebin=timebin)
                ms.close()

            # if verbose:
            #   print('Regridding into a single spectral window...')
            # print('Reading data spw by spw')

            try:
                tb.open(vis_spl + '/POLARIZATION')
                corrtype = tb.getcell('CORR_TYPE', 0)
                pol = [stokesenum[p] for p in corrtype]
                tb.close()
            except:
                pol = []

            if regridfreq:
                if verbose:
                    print('Regridding into a single spectral window...')
                ms.open(vis_spl, nomodify=False)
                ms.cvel(outframe='LSRK', mode='frequency', interp='nearest')
                ms.selectinit(datadescid=0, reset=True)
                data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
                specamp = data['amplitude']
                freq = data['axis_info']['freq_axis']['chan_freq']
            else:
                if verbose:
                    print('Concatenating visibility data spw by spw')
                ms.open(vis_spl)
                ms.selectinit(datadescid=0, reset=True)
                spwinfo = ms.getspectralwindowinfo()
                specamp = []
                freq = []
                time = []
                if verbose:
                    print('A total of {0:d} spws to fill'.format(len(spwinfo.keys())))
                for n, descid in enumerate(spwinfo.keys()):
                    ms.selectinit(datadescid=0, reset=True)
                    if verbose:
                        print('filling up spw #{0:d}: {1:s}'.format(n, descid))
                    descid = int(descid)
                    ms.selectinit(datadescid=n)  # , reset=True)
                    data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
                    if verbose:
                        print('shape of this spw', data['amplitude'].shape)
                    specamp_ = data['amplitude']
                    freq_ = data['axis_info']['freq_axis']['chan_freq'].squeeze()
                    if len(freq_.shape) > 1:
                        # chan_freq for each datadecid contains the info for all the spws
                        freq = freq_.transpose().flatten()
                    else:
                        freq.append(freq_)
                    time_ = data['time']
                    if fillnan is not None:
                        flag_ = ms.getdata(['flag', 'time', 'axis_info'], ifraxis=True)['flag']
                        if type(fillnan) in [int, float]:
                            specamp_[flag_] = float(fillnan)
                        else:
                            specamp_[flag_] = 0.0
                    specamp.append(specamp_)
                    time.append(time_)
                specamp = np.concatenate(specamp, axis=1)
                try:
                    # if len(freq.shape) > 1:
                    freq = np.concatenate(freq, axis=0)
                except ValueError:
                    pass
                ms.selectinit(datadescid=0, reset=True)
            ms.close()
            if os.path.exists(vis_spl):
                os.system('rm -rf ' + vis_spl)
            (npol, nfreq, nbl, ntim) = specamp.shape
            freq = freq.reshape(nfreq)

            tim = data['time']

        if verbose:
            print('npol, nfreq, nbl, ntime:', (npol, nfreq, nbl, ntim))
        spec = np.swapaxes(specamp, 2, 1)

        if domedian:
            if verbose:
                print('doing median of all the baselines')
            # mask zero values before median
            # spec_masked = np.ma.masked_where(spec < 1e-9, spec)
            # spec_masked2 = np.ma.masked_invalid(spec)
            # spec_masked = np.ma.masked_array(spec, mask=np.logical_or(spec_masked.mask, spec_masked2.mask))
            # spec_med = np.ma.filled(np.ma.median(spec_masked, axis=1), fill_value=0.)
            spec = np.abs(spec)
            if ds_normalised == False:
                # mask zero values before median
                spec_masked = np.ma.masked_where(spec < 1e-9, spec)
                spec_masked2 = np.ma.masked_invalid(spec)
                spec_masked = np.ma.masked_array(spec, mask=np.logical_or(spec_masked.mask, spec_masked2.mask))
                spec_med = np.ma.filled(np.ma.median(spec_masked, axis=1), fill_value=0.)
                # spec_med = np.nanmedian(spec, axis=1)
                nbl = 1
                ospec = spec_med.reshape((npol, nbl, nfreq, ntim))
            else:
                spec_med_time = np.expand_dims(np.nanmedian(spec, axis=3), axis=3)
                spec_normalised = (spec - spec_med_time) / spec_med_time
                spec_med_bl = np.nanmedian(spec_normalised, axis=1)
                nbl = 1
                ospec = spec_med_bl.reshape((npol, nbl, nfreq, ntim))
                ospec = ospec * 1e4
        else:
            ospec = spec
        # Save the dynamic spectral data
        if not specfile:
            specfile = msfile + '.dspec.npz'
        if os.path.exists(specfile):
            os.system('rm -rf ' + specfile)
        np.savez(specfile, spec=ospec, tim=tim, freq=freq,
                 timeran=timeran, spw=spw, bl=bl, uvrange=uvrange, pol=pol)
        if verbose:
            print('Median dynamic spectrum saved as: ' + specfile)

        self.read(specfile)
        return specfile

    def wrt_dspec(self, specfile=None, specdat=None):
        try:
            specfile
        except NameError:
            print('No input centfile specified for reading. Aborting...')
        if not specdat:
            print('Output file name is not specified, use the default convention')
            specdat = specfile.replace('npz', 'dat')
        specdata = np.load(specfile)
        spec = np.copy(specdata['spec'][:, :, :, :])
        npl, nbl, nf, nt = spec.shape
        print('Dimension of the data cube -- # of pol, # of baseline, # of frequency, # of time:')
        print(npl, nbl, nf, nt)
        nele = npl * nbl * nf * nt
        # need to transpose to nf, nt, npl, nbl
        spec = np.transpose(spec, (2, 3, 0, 1))
        # Turn fmx into a 1d array for writing
        spec = spec.flatten()
        # Pack the data into a byte buffer.  Note the asterisks in the next three lines!
        # If your data includes integers or other data types, use 'i' or other instead of 'f'
        buf = struct.pack(str(nf) + 'f', *specdata['freq'])
        buf += struct.pack(str(nt) + 'd', *specdata['tim'])
        buf += struct.pack(str(nele) + 'f', *spec)
        with open(specdat, 'wb') as f:
            f.write(buf)
        f.close()

    def rd_dspec(self, specdata, spectype='amp', spec_unit='jy'):
        spectype = spectype.lower()
        if type(specdata) is str:
            try:
                specdata = np.load(specdata)
            except:
                specdata = np.load(specdata, encoding='latin1')
        try:
            if np.iscomplexobj(specdata['spec']):
                if spectype == 'amp':
                    if spec_unit.lower() == 'jy':
                        spec = np.abs(specdata['spec']) / 1.e4
                    elif spec_unit.lower() == 'sfu' or spec_unit.lower() == 'k':
                        spec = np.abs(specdata['spec'])
                    else:
                        raise ValueError("Input spec_unit is {}. "
                                         "If spectype = 'amp', spec_unit must be 'jy', 'sfu', or 'k'".format(spec_unit))
                elif spectype == 'pha':
                    spec = np.angle(specdata['spec'])
                else:
                    raise ValueError('spectype must be amp or phase!')
            else:
                if spec_unit.lower() == 'jy':
                    spec = specdata['spec'] / 1.e4
                elif spec_unit.lower() == 'sfu' or spec_unit.lower() == 'k':
                    spec = specdata['spec']
                else:
                    raise ValueError("Input spec_unit is {}. "
                                     "If spectype = 'amp', spec_unit must be 'jy', 'sfu', or 'k'".format(spec_unit))
            if 'bl' in specdata.keys():
                try:
                    bl = specdata['bl'].item()
                except:
                    bl = specdata['bl']
            else:
                bl = ''

            tim = specdata['tim']
            freq = specdata['freq']
            if 'pol' in specdata.keys():
                pol = specdata['pol']
            else:
                pol = []

            return spec, tim, freq, bl, pol, spec_unit
        except:
            raise ValueError('format of specdata not recognized. Check your input')

    def concat_dspec(self, specfiles, outfile=None, savespec=False):
        '''
        concatenate a list of specfiles in time axis
        :param specfiles: a list of specfile to concatenate
        :return: concatenated specdata
        '''
        from tqdm import tqdm
        if isinstance(specfiles, list):
            if len(specfiles) > 1:
                specfiles = sorted(specfiles)
            else:
                print('Abort. Only one specfile is provided.')
                return -1
        else:
            print('Please provide a list of specfiles')
            return -1

        specdata = {}

        specdata_ = np.load(specfiles[0])
        for k, v in specdata_.iteritems():
            specdata[k] = v

        for spfile in tqdm(specfiles[1:]):
            specdata_ = np.load(spfile)
            specdata['spec'] = np.concatenate((specdata['spec'], specdata_['spec']), axis=-1)
            specdata['tim'] = np.hstack((specdata['tim'], specdata_['tim']))

        if savespec:
            if not outfile:
                specfile = 'dspec.npz'
            if os.path.exists(specfile):
                os.system('rm -rf ' + specfile)
            np.savez(specfile, **specdata)
        return specdata



    def peek(self, *args, **kwargs):
        """
        Plot dynamaic spectrum onto current axes.

        Parameters
        ----------
        *args : dict

        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """
        plt.figure()
        ret = self.plot(*args, **kwargs)
        plt.show()
        return ret

    def plot(self, pol='I', vmin=None, vmax=None, norm='log', cmap='viridis', cmap2='viridis', vmin2=None, vmax2=None, figsize=None,
             timerange=None, freqrange=None, bkgtim=None, interp_method='linear',  ignore_gaps=True, freq_unit='GHz', spec_unit=None,
             plot_fast=False, percentile=[1, 99], minmaxpercentile=False, axes=None, **kwargs):
        """
        Plots the dynamic spectrum for a given polarization.

        This method generates a plot of the dynamic spectrum, allowing for customization of the visualization through various parameters such as color maps, normalization, and value ranges.

        :param pol: Polarization for plotting. Default is 'I'.
        :type pol: str, optional
        :param vmin: Minimum intensity value for the color scale. Default is None, which auto-scales.
        :type vmin: float, optional
        :param vmax: Maximum intensity value for the color scale. Default is None, which auto-scales.
        :type vmax: float, optional
        :param norm: Normalization of the color scale. Can be 'linear', 'log', or a custom normalization. Default is 'log'.
        :type norm: str or `matplotlib.colors.Normalize`, optional
        :param cmap: Matplotlib colormap name or instance for the plot. Default is 'viridis'.
        :type cmap: str, optional
        :param cmap2: Matplotlib colormap name or instance for the second polarization plot. Default is 'viridis'.
        :type cmap2: str, optional
        :param vmin2: Minimum intensity value for the color scale of the second polarization. Default is None.
        :type vmin2: float, optional
        :param vmax2: Maximum intensity value for the color scale of the second polarization. Default is None.
        :type vmax2: float, optional
        :param timerange: Time range to plot, formatted as ['start_time', 'end_time']. Times should be in ISO format. Defaults to None, which plots the entire range.
        :type timerange: list of str, optional
        :param freqrange: Frequency range to plot, in units specified by `freq_unit`. Defaults to None, which plots the entire range.
        :type freqrange: list of float, optional
        :param bkgtim: Background time interval(s). It can be a single time interval (list of two strings)
            or multiple time intervals (list of lists of two strings each).
            Example: ['2024-06-11T19:41:00', '2024-06-11T20:40:58']
                     [['2024-06-11T19:41:00', '2024-06-11T20:40:58'], ['2024-06-11T21:41:00', '2024-06-11T22:40:58']]
        :type bkgtim: list or list of lists
        :param interp_method: The method of interpolation to use. Default is 'linear'.
            Can be any method supported by scipy.interpolate.interp1d.
        :type interp_method: str, optional
        :param ignore_gaps: If True, ignores gaps in the frequency axis. Default is True.
        :type ignore_gaps: bool, optional
        :param freq_unit: Unit for the frequency axis ('kHz', 'MHz', 'GHz'). Default is 'GHz'.
        :type freq_unit: str, optional
        :param spec_unit: Unit for the spectrum data ('sfu', 'Jy', or 'K'). If not specified, uses the unit from the Dspec object.
        :type spec_unit: str, optional
        :param plot_fast: If True, uses a faster plotting method which may reduce detail. Default is False.
        :type plot_fast: bool, optional
        :param percentile: Percentile values to use for auto-scaling the color range. Default is [1, 99].
        :type percentile: list of float, optional
        :param minmaxpercentile: If True, uses percentile for vmin and vmax. Default is False.
        :type minmaxpercentile: bool, optional
        :param kwargs: Any additional plot arguments that should be used when plotting.
        :type kwargs: dict

        :return: The matplotlib figure object containing the plot.
        :rtype: `matplotlib.figure.Figure`

        :example:
        >>> dspec = Dspec()
        >>> fig = dspec.plot(pol='I', vmin=0.1, vmax=5, cmap='hot', freqrange=[1, 2], timerange=['2021-01-01T00:00:00', '2021-01-01T01:00:00'])
        >>> plt.show()

        :notes:
        For dual-polarization plots (e.g., 'RRLL'), `cmap2` along with `vmin2` and `vmax2` can be specified to customize the appearance of the second polarization.
        """
        # Set up variables
        import matplotlib.colors as colors
        import matplotlib.pyplot as plt
        from astropy.time import Time

        if pol not in ['RR', 'LL', 'RRLL', 'XX', 'YY', 'XY', 'YX', 'XXYY', 'I', 'V', 'IV', 'IP']:
            print("Please enter 'RR', 'LL', 'RRLL','XX', 'YY', 'XY', 'YX', 'XXYY', 'I', 'V', 'IV', or 'IP' for pol")
            return 0

        spec = self.data

        try:
            cmap = copy(plt.get_cmap(cmap))
        except:
            cmap = copy(plt.get_cmap('viridis'))
        cmap.set_bad(cmap(0.0))

        if norm == 'linear':
            norm = colors.Normalize(vmax=vmax, vmin=vmin)
        if norm == 'log':
            norm = colors.LogNorm(vmax=vmax, vmin=vmin)

        if not isinstance(norm, object):
            print('color normalization is not defined as matplotlib.colors. Use default LogNorm.')
            norm = colors.LogNorm(vmax=vmax, vmin=vmin)

        bl = self.bl
        freq = self.freq_axis
        if spec_unit is None:
            spec_unit = self.spec_unit
        elif spec_unit.lower() != self.spec_unit.lower():
            if spec_unit.lower() == 'sfu' and self.spec_unit.lower() == 'jy':
                spec = np.copy(self.data) / 1e4
            elif spec_unit.lower() == 'jy' and self.spec_unit.lower() == 'sfu':
                spec = np.copy(self.data) * 1e4
            else:
                print('Spectrum unit conversion from {0:s} to {1:s} not supported'.format(self.spec_unit, spec_unit))
                print('Use the original one.')
                spec_unit = self.spec_unit
        # The following is only for the plots
        if spec_unit.lower() == 'jy':
            spec_unit_print = 'Jy'
        if spec_unit.lower() == 'k':
            spec_unit_print = 'K'
        if spec_unit.lower() == 'sfu':
            spec_unit_print = 's.f.u'
        if hasattr(self, 'spec_name'):
            spec_name = self.spec_name
        else:
            spec_name = 'Intensity'

        if spec.ndim == 2:
            nfreq, ntim = len(self.freq_axis), len(self.time_axis)
            npol = 1
            nbl = 1
            polnames = self.pol
        else:
            (npol, nbl, nfreq, ntim) = spec.shape
            polnames = self.pol
            if len(polnames) != npol:
                print('The polarization dimension in the data {0:d} does not match the names {1:d}. Abort.'.format(npol,
                                                                                                                   len(polnames)))

        tim_ = self.time_axis
        tim_plt = tim_.plot_date

        if timerange:
            if isinstance(timerange[0], str):
                timerange = Time(timerange)
            tidx = np.where((tim_ >= timerange[0]) & (tim_ <= timerange[1]))[0]
        else:
            tidx = range(ntim)

        if ignore_gaps:
            df = np.median(freq[1:] - freq[:-1])
            fbreak, = np.where(freq[1:] - freq[:-1] > 2 * df)
            for n, fb in enumerate(fbreak):
                loc = fb + n + 1
                freq = np.concatenate((freq[:loc], np.array([freq[loc - 1] + df]), freq[loc:]))

        if freqrange:
            if freq_unit.lower() == 'ghz':
                fidx = np.where((freq >= freqrange[0] * 1e9) & (freq <= freqrange[1] * 1e9))[0]
            if freq_unit.lower() == 'mhz':
                fidx = np.where((freq >= freqrange[0] * 1e6) & (freq <= freqrange[1] * 1e6))[0]
            if freq_unit.lower() == 'kHz':
                fidx = np.where((freq >= freqrange[0] * 1e3) & (freq <= freqrange[1] * 1e3))[0]
        else:
            fidx = range(nfreq)

        # setup plot parameters
        print('ploting dynamic spectrum...')

        for b in range(nbl):
            if pol not in ['RRLL', 'XXYY', 'IV', 'IP']:
                if spec.ndim == 2:
                    spec_plt = spec
                else:
                    if pol in ['RR', 'XX']:
                        spec_plt = spec[0, b, :, :]
                    elif pol in ['LL', 'YY']:
                        spec_plt = spec[1, b, :, :]
                    elif pol in ['XY']:
                        spec_plt = spec[2, b, :, :]
                    elif pol in ['YX']:
                        spec_plt = spec[3, b, :, :]
                    elif pol == 'I':
                        if ('XX' in polnames) and ('YY' in polnames):
                            spec_plt = spec[0, b, :, :] + spec[1, b, :, :]
                        elif ('RR' in polnames) and ('LL' in polnames):
                            spec_plt = (spec[0, b, :, :] + spec[1, b, :, :]) / 2.
                        elif ('I' in polnames):
                            spec_plt = spec[0, b, :, :]
                    elif pol == 'V':
                        if ('XX' in polnames) and ('YY' in polnames):
                            # TODO: This does not seem to be correct. @Sijie Yu, could you please check?
                            spec_plt = spec[0, b, :, :] - spec[1, b, :, :]
                        elif ('RR' in polnames) and ('LL' in polnames):
                            spec_plt = (spec[0, b, :, :] - spec[1, b, :, :]) / 2.
                        elif ('V' in polnames):
                            spec_plt = spec[1, b, :, :]
                if ignore_gaps:
                    for n, fb in enumerate(fbreak):
                        loc = fb + n + 1
                        spec_plt = np.concatenate((spec_plt[:loc], np.zeros((1, ntim)) + np.nan, spec_plt[loc:]), 0)


                if bkgtim:
                    spec_plt -=calc_bkg_dspec(spec_plt, tim_[tidx], bkgtim, interp_method=interp_method)

                # if plot_fast:


                if figsize is None:
                    figsize = (8, 4)
                if axes:
                    ax = axes
                    fig=None
                else:
                    fig = plt.figure(figsize=figsize, dpi=100)
                    ax = fig.add_subplot(111)

                if freq_unit.lower() == 'ghz':
                    freq_plt = freq / 1e9
                if freq_unit.lower() == 'mhz':
                    freq_plt = freq / 1e6
                if freq_unit.lower() == 'khz':
                    freq_plt = freq / 1e3

                spec_plt = spec_plt[fidx, :][:, tidx]
                tim_plt = tim_plt[tidx]
                freq_plt = freq_plt[fidx]

                # Change the default for Stokes V
                if pol == 'V':
                    cmap = 'gray'
                    if (vmax is None) and (vmin is None):
                        vmax = np.nanmax(np.abs(spec_plt))
                        vmin = -vmax
                    elif (vmax is None) and not (vmin is None):
                        vmax = -vmin
                    elif not (vmax is None) and (vmin is None):
                        vmin = -vmax
                    norm = colors.Normalize(vmax=vmax, vmin=vmin)

                if minmaxpercentile:
                    if percentile[0] > 0 and percentile[1] < 100 and percentile[0] < percentile[1]:
                        norm.vmax = np.nanpercentile(spec_plt, percentile[1])
                        norm.vmin = np.nanpercentile(spec_plt, percentile[0])

                if plot_fast:
                    # rebin the data to speed up plotting
                    ds_shape = spec_plt.shape
                    if ds_shape[1] > 2100:
                        pix_rebin = ds_shape[1] // 2048
                        new_len = (ds_shape[1] // pix_rebin)
                        new_len_total = new_len * pix_rebin

                        spec_plt = rebin2d(spec_plt[:, 0:new_len_total], (ds_shape[0], new_len))
                        # tim_plt = rebin1d(tim_plt[0:new_len_total], 2048)

                    im = ax.imshow(spec_plt, cmap=cmap, norm=norm, aspect='auto', origin='lower',
                                   extent=[tim_plt[0], tim_plt[-1], freq_plt[0], freq_plt[-1]])

                else:
                    im = ax.pcolormesh(tim_plt, freq_plt, spec_plt, cmap=cmap, norm=norm, shading='auto',
                                       rasterized=True)
                    ax.set_xlim(tim_plt[0], tim_plt[-1])
                    ax.set_ylim(freq_plt[0], freq_plt[-1])

                    def format_coord(x, y):
                        col = np.argmin(np.absolute(tim_plt - x))
                        row = np.argmin(np.absolute(freq_plt - y))
                        if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                            timstr = tim_[col].isot
                            flux = spec_plt[row, col]
                            return 'time {0} = {1}, freq {2} = {3:.3f} {4:s}, flux = {5:.2f} {6:s}'.format(col, timstr,
                                                                                                           row,
                                                                                                           y, freq_unit,
                                                                                                           flux,
                                                                                                           spec_unit_print)
                        else:
                            return 'x = {0}, y = {1:.3f}'.format(x, y)

                    ax.format_coord = format_coord

                # make colorbar
                divider = make_axes_locatable(ax)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax, cax=cax_spec)
                clb_spec.set_label(f'{spec_name} [{spec_unit_print}]')


                ax.set_ylabel('Frequency [{0:s}]'.format(freq_unit))
                if bl:
                    ax.set_title('{0:s}-{1:s} dynamic spectrum @ bl {2:s} for pol {3:s}'.
                                 format(self.observatory, self.telescope, bl.split(';')[b], pol))
                else:
                    ax.set_title('{0:s}-{1:s} dynamic spectrum for pol {2:s}'.
                                 format(self.observatory, self.telescope, pol))
                locator = AutoDateLocator(minticks=2)
                ax.xaxis.set_major_locator(locator)
                formatter = AutoDateFormatter(locator)
                formatter.scaled[1 / 24] = '%D %H'
                formatter.scaled[1 / (24 * 60)] = '%H:%M'
                ax.xaxis.set_major_formatter(formatter)
                ax.set_autoscale_on(False)

            else:
                if figsize is None:
                    figsize = (8, 6)
                fig = plt.figure(figsize=figsize, dpi=100)
                if ('RR' in polnames) and ('LL' in polnames):
                    R_plot = np.absolute(spec[0, b, :, :])
                    L_plot = np.absolute(spec[1, b, :, :])
                    I_plot = (R_plot + L_plot) / 2.
                    V_plot = (R_plot - L_plot) / 2.
                if ('I' in polnames) and ('V' in polnames or 'P' in polnames) :
                    I_plot = spec[0, b, :, :]
                    V_plot = spec[1, b, :, :]
                    R_plot = I_plot + V_plot
                    L_plot = I_plot - V_plot
                if pol in ['RRLL', 'XXYY']:
                    spec_plt_1 = R_plot
                    spec_plt_2 = L_plot
                    polstr = [pol[:2], pol[2:]]
                elif pol == 'IV':
                    spec_plt_1 = I_plot
                    spec_plt_2 = V_plot
                    cmap2 = 'RdBu_r'
                    if (vmax2 is None) and (vmin2 is None):
                        vmax2 = np.nanmax(np.abs(spec_plt_2))
                        vmin2 = -vmax2
                    elif (vmax2 is None) and not (vmin2 is None):
                        vmax2 = -vmin2
                    elif not (vmax2 is None) and (vmin2 is None):
                        vmin2 = -vmax2
                    polstr = ['I', 'V']
                elif pol == 'IP':
                    # this is for Stokes I + polarization degree
                    spec_plt_1 = I_plot
                    spec_plt_2 = V_plot / I_plot
                    cmap2 = 'RdBu_r'
                    if (vmax2 is None) and (vmin2 is None):
                        vmax2 = 1.
                        vmin2 = -1.
                    elif (vmax2 is None) and not (vmin2 is None):
                        vmax2 = -vmin2
                    elif not (vmax2 is None) and (vmin2 is None):
                        vmin2 = -vmax2
                    polstr = ['I', 'V/I']

                if ignore_gaps:
                    for n, fb in enumerate(fbreak):
                        loc = fb + n + 1
                        spec_plt_1 = np.concatenate((spec_plt_1[:loc], np.zeros((1, ntim)) + np.nan, spec_plt_1[loc:]),
                                                    0)
                        spec_plt_2 = np.concatenate((spec_plt_2[:loc], np.zeros((1, ntim)) + np.nan, spec_plt_2[loc:]),
                                                    0)


                if bkgtim:
                    spec_plt_1 -= calc_bkg_dspec(spec_plt_1, tim_[tidx], bkgtim, interp_method=interp_method)
                    spec_plt_2 -= calc_bkg_dspec(spec_plt_2, tim_[tidx], bkgtim, interp_method=interp_method)

                ax1 = fig.add_subplot(211)
                if freq_unit.lower() == 'ghz':
                    freq_plt = freq / 1e9
                if freq_unit.lower() == 'mhz':
                    freq_plt = freq / 1e6
                if freq_unit.lower() == 'khz':
                    freq_plt = freq / 1e3

                spec_plt_1 = spec_plt_1[fidx, :][:, tidx]
                spec_plt_2 = spec_plt_2[fidx, :][:, tidx]
                tim_plt = tim_plt[tidx]
                freq_plt = freq_plt[fidx]

                if minmaxpercentile:
                    if percentile[0] > 0 and percentile[1] < 100 and percentile[0] < percentile[1]:
                        norm.vmax = np.nanpercentile(spec_plt_1, percentile[1])
                        norm.vmin = np.nanpercentile(spec_plt_1, percentile[0])

                if plot_fast:
                    # compress in time (idx1)
                    ds_shape = spec_plt_1.shape
                    if ds_shape[1] > 2100:
                        pix_rebin = ds_shape[1] // 2048
                        new_len = (ds_shape[1] // pix_rebin)
                        new_len_total = new_len * pix_rebin

                        spec_plt_1 = rebin2d(spec_plt_1[:, 0:new_len_total], (ds_shape[0], new_len))
                        spec_plt_2 = rebin2d(spec_plt_2[:, 0:new_len_total], (ds_shape[0], new_len))
                        # tim_plt = rebin1d(tim_plt[0:new_len_total], new_len)

                    im = ax1.imshow(spec_plt_1, cmap=cmap, norm=norm, aspect='auto', origin='lower',
                                    extent=[tim_plt[0], tim_plt[-1], freq_plt[0], freq_plt[-1]])
                else:
                    im = ax1.pcolormesh(tim_plt, freq_plt, spec_plt_1, cmap=cmap, norm=norm, shading='auto',
                                        rasterized=True)
                    ax1.set_xlim(tim_plt[0], tim_plt[-1])
                    ax1.set_ylim(freq_plt[0], freq_plt[-1])

                    def format_coord(x, y):
                        col = np.argmin(np.absolute(tim_plt - x))
                        row = np.argmin(np.absolute(freq_plt - y))
                        if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                            timstr = tim_[col].isot
                            flux = spec_plt_1[row, col]
                            return 'time {0} = {1}, freq {2} = {3:.3f} {4:s}, value = {5:.2f} {6:s}'.format(col, timstr,
                                                                                                            row,
                                                                                                            y,
                                                                                                            freq_unit,
                                                                                                            flux,
                                                                                                            spec_unit_print)
                        else:
                            return 'x = {0}, y = {1:.3f}'.format(x, y)

                    ax1.format_coord = format_coord

                divider = make_axes_locatable(ax1)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax1, cax=cax_spec)
                clb_spec.set_label(f'{spec_name} [{spec_unit_print}]')


                ax1.set_ylabel('Frequency [{0:s}]'.format(freq_unit))

                locator = AutoDateLocator(minticks=2)
                ax1.xaxis.set_major_locator(locator)
                # ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
                formatter = AutoDateFormatter(locator)
                formatter.scaled[1 / 24] = '%D %H'
                formatter.scaled[1 / (24 * 60)] = '%H:%M'
                ax1.xaxis.set_major_formatter(formatter)
                if bl:
                    ax1.set_title('{0:s}-{1:s} dynamic spectrum @ bl {2:s} for pol {3:s}'. \
                                  format(self.observatory, self.telescope, bl.split(';')[b], polstr[0]))
                else:
                    ax1.set_title('{0:s}-{1:s} dynamic spectrum for pol {2:s}'. \
                                  format(self.observatory, self.telescope, polstr[0]))

                ax1.set_autoscale_on(False)

                ax2 = fig.add_subplot(212, sharex=ax1, sharey=ax1)
                if cmap2 is None:
                    cmap2 = cmap
                if vmin2 is None:
                    vmin2 = vmin
                if vmax2 is None:
                    vmax2 = vmax
                norm2 = colors.Normalize(vmax=vmax2, vmin=vmin2)



                if plot_fast:
                    im = ax2.imshow(spec_plt_2, cmap=cmap2, norm=norm2, aspect='auto', origin='lower',
                                    extent=[tim_plt[0], tim_plt[-1], freq_plt[0], freq_plt[-1]])
                else:
                    im = ax2.pcolormesh(tim_plt, freq_plt, spec_plt_2, cmap=cmap2, norm=norm2, shading='auto',
                                        rasterized=True)
                    def format_coord(x, y):
                        col = np.argmin(np.absolute(tim_plt - x))
                        row = np.argmin(np.absolute(freq_plt - y))
                        if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                            timstr = tim_[col].isot
                            flux = spec_plt_2[row, col]
                            return 'time {0} = {1}, freq {2} = {3:.3f} {4:s}, value = {5:.2f} {6:s}'.format(col, timstr,
                                                                                                            row,
                                                                                                            y,
                                                                                                            freq_unit,
                                                                                                            flux,
                                                                                                            spec_unit_print)
                        else:
                            return 'x = {0}, y = {1:.3f}'.format(x, y)
                    ax2.format_coord = format_coord

                    ax2.set_xlim(tim_plt[0], tim_plt[-1])
                    ax2.set_ylim(freq_plt[0], freq_plt[-1])

                divider = make_axes_locatable(ax2)
                cax_spec = divider.append_axes('right', size='1.5%', pad=0.05)
                clb_spec = plt.colorbar(im, ax=ax2, cax=cax_spec)
                if pol == 'IP':
                    clb_spec.set_label(polstr[1])
                else:
                    clb_spec.set_label(f'{spec_name} [{spec_unit_print}]')

                locator = AutoDateLocator(minticks=2)
                ax2.xaxis.set_major_locator(locator)
                formatter = AutoDateFormatter(locator)
                formatter.scaled[1 / 24] = '%D %H'
                formatter.scaled[1 / (24 * 60)] = '%H:%M'
                ax2.xaxis.set_major_formatter(formatter)

                ax2.set_ylabel('Frequency [{0:s}]'.format(freq_unit))
                if bl:
                    ax2.set_title('{0:s}-{1:s} dynamic spectrum @ bl {2:s} for pol {3:s}'. \
                                  format(self.observatory, self.telescope, bl.split(';')[b], polstr[1]))
                else:
                    ax2.set_title('{0:s}-{1:s} dynamic spectrum for pol {2:s}'. \
                                  format(self.observatory, self.telescope, polstr[1]))

                ax2.set_autoscale_on(False)

            if fig:
                fig.tight_layout()
        return fig
