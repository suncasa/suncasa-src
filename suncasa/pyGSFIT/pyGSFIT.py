import sys, os
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QPushButton, QLineEdit, QLabel, \
    QTextEdit, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QStatusBar, QProgressBar

from PyQt5.QtCore import QThreadPool
from PyQt5 import QtGui
import pyqtgraph as pg

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas)
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import patches, cm
import gsutils
from suncasa.io import ndfits
from astropy.time import Time
import numpy as np
import sunpy
from sunpy import map as smap
from astropy.io import fits
import astropy.units as u

sys.path.append('./')
from gsutils import ff_emission
import numpy.ma as ma
import warnings

warnings.filterwarnings("ignore")
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(imageAxisOrder='row-major')


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.fname = '<Select or enter a valid EOVSA fits file name>'
        self.aiafname = '<Select or enter a valid AIA fits file name>'
        self.fitsentry = QLineEdit()
        self.title = 'pyGSFIT'
        self.left = 0
        self.top = 0
        self.width = 1200
        self.height = 1000
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self._main = QWidget()
        self.setCentralWidget(self._main)
        self.initUI()
        self.threadpool = QThreadPool()
        self.haseovsamap = False
        self.hasaiamap = False

    def initUI(self):

        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.progressBar = QProgressBar()
        self.progressBar.setGeometry(10, 10, 200, 15)

        layout = QVBoxLayout()
        # Initialize tab screen
        self.tabs = QTabWidget()
        tab_explorer = QWidget()
        tab_fit = QWidget()
        tab_analyzer = QWidget()

        # Add tabs
        self.tabs.addTab(tab_explorer, "Explorer")
        self.tabs.addTab(tab_fit, "Fit")
        self.tabs.addTab(tab_analyzer, "Analyzer")

        # Each tab's user interface is complex, so this splits them into separate functions.
        self.initUItab_explorer()
        self.initUItab_fit()
        self.initUItab_analyzer()

        #self.tabs.currentChanged.connect(self.tabChanged)

        # Add tabs to widget
        layout.addWidget(self.tabs)
        self._main.setLayout(layout)

        self.show()

    # Explorer Tab User Interface
    #
    def initUItab_explorer(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()

        # Create Data Select tab
        upperbox = QVBoxLayout()  # Upper box has two hboxes: leftbox and quicklook plot box
        topbox = QHBoxLayout()  # Leftbox is the left of the upper box for file selection and fits header display

        # EOVSA FITS Filename entry
        eo_selection_box = QHBoxLayout()
        topbox.addLayout(eo_selection_box)
        # Create Browse button
        eo_selection_button = QPushButton("Load EOVSA")
        eo_selection_button.clicked.connect(self.eofile_select)
        eo_selection_box.addWidget(eo_selection_button)

        # Create LineEdit widget for FITS filename
        self.fitsentry.resize(8 * len(self.fname), 20)
        eo_selection_box.addWidget(self.fitsentry)

        # AIA FITS Filename entry
        aia_selection_box = QHBoxLayout()
        topbox.addLayout(aia_selection_box)
        # Create Browse button
        aia_selection_button = QPushButton("Load AIA")
        aia_selection_button.clicked.connect(self.aiafile_select)
        aia_selection_box.addWidget(aia_selection_button)
        # Create LineEdit widget for AIA FITS file
        self.aiafitsentry = QLineEdit()
        self.aiafitsentry.resize(8 * len(self.fname), 20)
        aia_selection_box.addWidget(self.aiafitsentry)

        # Add left box of the upper box
        upperbox.addLayout(topbox)

        # Create label and TextEdit widget for FITS header information
        # fitsinfobox = QVBoxLayout()
        # topbox.addLayout(fitsinfobox)
        # self.infoEdit = QTextEdit()
        # self.infoEdit.setReadOnly(True)
        # self.infoEdit.setMinimumHeight(100)
        # self.infoEdit.setMaximumHeight(400)
        # self.infoEdit.setMinimumWidth(300)
        # self.infoEdit.setMaximumWidth(500)
        # fitsinfobox.addWidget(QLabel("FITS Information"))
        # fitsinfobox.addWidget(self.infoEdit)

        # Add quicklook plotting area
        qlookarea = QVBoxLayout()
        self.qlookcanvas = FigureCanvas(Figure(figsize=(8, 4)))
        self.qlooktoolbar = NavigationToolbar(self.qlookcanvas, self)
        qlookarea.addWidget(self.qlooktoolbar)
        qlookarea.addWidget(self.qlookcanvas)
        self.qlook_axs = self.qlookcanvas.figure.subplots(nrows=1, ncols=4)

        upperbox.addLayout(qlookarea)
        mainlayout.addLayout(upperbox)

        # lowerbox
        lowerbox = QGridLayout()  # lower box has two hboxes: left for multi-panel display and right for spectrum

        # Add plotting area for multi-panel EOVSA images
        # self.meocanvas = FigureCanvas(Figure(figsize=(6, 6)))
        self.meocanvas = pg.ImageView(name='EOVSA Explorer')
        # self.meotoolbar = NavigationToolbar(self.meocanvas, self)
        # meoplotarea.addWidget(self.meotoolbar)
        lowerbox.addWidget(self.meocanvas, 0, 0)
        lowerbox.setColumnStretch(0, 1.5)
        # self.meo_axs = self.meocanvas.figure.subplots(nrows=2, ncols=2)
        # self.plot_axs.clear()
        # im = self.dspec_ax.pcolormesh(pd, fghz,np.clip(spec, minval, maxval))
        # self.dspec_ax.xaxis_date()
        # self.dspec_ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        # self.dspec_ax.set_ylabel('Frequency [GHz]')
        # self.dspec_ax.set_xlabel('Time [UT]')
        # lowerbox.addLayout(meoplotarea)

        # right box for spectral plots
        specplotarea = QVBoxLayout()
        # self.speccanvas = FigureCanvas(Figure(figsize=(4, 6)))
        self.speccanvas = pg.PlotWidget()
        # self.spectoolbar = NavigationToolbar(self.speccanvas, self)
        # specplotarea.addWidget(self.spectoolbar)
        specplotarea.addWidget(self.speccanvas)
        # self.spec_axs = self.speccanvas.figure.subplots(nrows=1, ncols=1)
        lowerbox.addLayout(specplotarea, 0, 1)
        lowerbox.setColumnStretch(1, 1)

        mainlayout.addLayout(lowerbox)
        self.tabs.widget(0).setLayout(mainlayout)

    # Fit Tab User Interface
    #
    def initUItab_fit(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the tab for analyzing fit results"))
        self.tabs.widget(1).setLayout(mainlayout)

    # Analyzer Tab User Interface
    #
    def initUItab_analyzer(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the tab for analyzing fit results"))
        self.tabs.widget(1).setLayout(mainlayout)

    def eofile_select_return(self):
        ''' Called when the FITS filename LineEdit widget gets a carriage-return.
            Trys to read FITS header and return header info (no data read at this time)
        '''
        self.fname = self.fitsentry.text()
        self.fitsdata = None
        try:
            rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rfdelts = ndfits.read(self.fname)
            # self.infoEdit.setPlainText(repr(rheader))
        except:
            self.statusBar.showMessage('Filename is not a valid FITS file', 2000)
            self.fname = '<Select or enter a valid fits filename>'
            self.fitsentry.setText(self.fname)
            # self.infoEdit.setPlainText('')

        self.plot_qlookmap()
        # self.initspecplot()
        self.plot_pg_eovsamap()

    def aiafile_select_return(self):
        ''' Called when the FITS filename LineEdit widget gets a carriage-return.
            Trys to read FITS header and return header info (no data read at this time)
        '''
        self.aiafname = self.aiafitsentry.text()
        self.fitsdata = None
        try:
            hdu = fits.open(self.aiafname)
            # self.infoEdit.setPlainText(repr(hdu[1].header))
        except:
            self.statusBar.showMessage('Filename is not a valid FITS file', 2000)
            self.aiafname = '<Select or enter a valid fits filename>'
            # self.aiafitsentry.setText(self.fname)
            # self.infoEdit.setPlainText('')

        self.plot_qlookmap()

    def eofile_select(self):
        """ Handle Browse button for EOVSA FITS file"""
        # self.fname = QFileDialog.getExistingDirectory(self, 'Select FITS File', './', QFileDialog.ShowDirsOnly)
        self.fname, _file_filter = QFileDialog.getOpenFileName(self, 'Select EOVSA FITS File to open', './',
                                                               "FITS Images (*.fits *.fit *.ft)")
        self.fitsentry.setText(self.fname)
        self.eofile_select_return()

    def aiafile_select(self):
        """ Handle Browse button for AIA FITS file """
        self.aiafname, _file_filter = QFileDialog.getOpenFileName(self, 'Select AIA FITS File to open', './',
                                                                  "FITS Images (*.fits *.fit *.ft)")
        self.aiafitsentry.setText(self.aiafname)
        self.aiafile_select_return()

    def initspecplot(self):
        """ Initial Spectral Plot if no data has been loaded (not used for pyqtgraph method)"""
        errobjs = []
        ax = self.spec_axs
        # for ax in self.spec_axs.reshape(-1):
        errobjs.append(ax.errorbar([], [], yerr=[], linestyle='', marker='o', mfc='none', mec='k', alpha=1.0))

        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim([1, 20])
        ax.set_ylim([0.1, 1000])
        ax.set_xticks([1, 5, 10, 20])
        ax.set_xticklabels([1, 5, 10, 20])
        ax.set_xticks([1, 5, 10, 20])
        ax.set_yticks([])
        ax.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
        ax.set_ylabel('T$_b$ [MK]')
        ax.set_xlabel('Frequency [GHz]')

        x = np.linspace(1, 20, 10)
        for ll in [-1, 0, 1, 2, 3, 4]:
            y = 10. ** (-2 * np.log10(x) + ll)
            ax.plot(x, y, 'k--', alpha=0.1)
            # y2 = 10. ** (-4 * np.log10(x) + ll)
            # y3 = 10. ** (-8 * np.log10(x) + ll)
            # ax_eospec.plot(x, y, 'k--', x, y2, 'k:', x, y3, 'k-.', alpha=0.1)

        self.speccanvas.figure.subplots_adjust(left=0.15, right=0.95,
                                               bottom=0.15, top=0.95)
        self.speccanvas.draw()
        self.errobjs = errobjs

    def plot_qlookmap(self):
        """Quicklook plot in the upper box using matplotlib.pyplot and sunpy.map"""
        from suncasa.utils import plot_mapX as pmX
        # Plot a quicklook map
        # self.qlook_ax.clear()
        ax0 = self.qlook_axs[0]
        bds = [5, 15, 25]
        if os.path.exists(self.fname):
            rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rfdelts = ndfits.read(self.fname)
            self.rmap = rmap
            self.rdata = rdata
            self.rfreqs = rfreqs
            self.rheader = rheader
            self.bdinfo = bdinfo = ndfits.get_bdinfo(rfreqs, rfdelts)
            self.cfreqs = cfreqs = bdinfo['cfreqs']
            self.cfreqs_all = cfreqs_all = bdinfo['cfreqs_all']
            self.freq_dist = lambda fq: (fq - cfreqs_all[0]) / (cfreqs_all[-1] - cfreqs_all[0])
            self.opencontour = True
            self.clevels = np.array([0.3, 0.7])
            self.calpha = 1.
            self.haseovsamap = True
            nspw = len(self.rfreqs)
            eodate = Time(self.rmap.date.mjd + self.rmap.exposure_time.value / 2. / 24 / 3600, format='mjd')
            rsun_obs = sunpy.coordinates.sun.angular_radius(eodate).value
            solar_limb = patches.Circle((0, 0), radius=rsun_obs, fill=False, color='k', lw=1, linestyle=':')
            ax0.add_patch(solar_limb)
            ny, nx = self.rmap.data.shape
            x0, x1 = (np.array([1, self.rmap.meta['NAXIS1']]) - self.rmap.meta['CRPIX1']) * self.rmap.meta['CDELT1'] + \
                     self.rmap.meta['CRVAL1']
            y0, y1 = (np.array([1, self.rmap.meta['NAXIS2']]) - self.rmap.meta['CRPIX2']) * self.rmap.meta['CDELT2'] + \
                     self.rmap.meta['CRVAL2']
            rect = patches.Rectangle((x0, y0), x1 - x0, y1 - y0, color='k', alpha=0.7, lw=1, fill=False)
            ax0.add_patch(rect)
            mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)
            icmap = plt.get_cmap('RdYlBu')
            # now plot 3 panels of EOVSA images
            # plot the images
            eocmap = plt.get_cmap('viridis')
            for n, bd in enumerate(bds):
                ax = self.qlook_axs[n + 1]
                cfreq = self.cfreqs[bd]
                eomap_ = smap.Map(self.rdata[bd], self.rheader)
                eomap = pmX.Sunmap(eomap_)
                eomap.imshow(axes=ax, cmap=eocmap)
                eomap.draw_grid(axes=ax)
                eomap.draw_limb(axes=ax)
                ax.set_xlim()
                ax.set_ylim()
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title('')
                ax.text(0.01, 0.98, '{0:.1f} GHz'.format(cfreq), ha='left', va='top',
                        fontsize=10, color='w', transform=ax.transAxes)
                ax.set_aspect('equal')
        else:
            self.statusBar.showMessage('EOVSA FITS file does not exist', 2000)
            self.fname = '<Select or enter a valid fits filename>'
            self.fitsentry.setText(self.fname)
            # self.infoEdit.setPlainText('')
            self.haseovsamap = False

        if os.path.exists(self.aiafname):
            try:
                aiacmap = plt.get_cmap('gray_r')
                aiamap = smap.Map(self.aiafname)
                aiamap.plot(axes=ax0, cmap=aiacmap)
                self.hasaiamap = True
            except:
                self.statusBar.showMessage('Something is wrong with the provided AIA FITS file', 2000)
        else:
            self.statusBar.showMessage('AIA FITS file does not exist', 2000)
            self.aiafname = '<Select or enter a valid fits filename>'
            self.aiafitsentry.setText(self.aiafname)
            # self.infoEdit.setPlainText('')
            self.hasaiamap = False
        cts = []
        if self.hasaiamap:
            aiamap.plot(axes=ax0, cmap=aiacmap)
        if self.haseovsamap:
            for s, sp in enumerate(self.rfreqs):
                data = self.rdata[s, ...]
                clvls = self.clevels * np.nanmax(data)
                rcmap = [icmap(self.freq_dist(self.cfreqs[s]))] * len(clvls)
                if self.opencontour:
                    cts.append(ax0.contour(mapx, mapy, data, levels=clvls,
                                           colors=rcmap,
                                           alpha=self.calpha))
                else:
                    cts.append(ax0.contourf(mapx, mapy, data, levels=clvls,
                                            colors=rcmap,
                                            alpha=self.calpha))

        ax0.set_xlim([-1200, 1200])
        ax0.set_ylim([-1200, 1200])
        ax0.set_xlabel('Solar-X [arcsec]')
        ax0.set_ylabel('Solar-y [arcsec]')
        # ax0.set_title('')
        ax0.set_aspect('equal')
        self.qlookcanvas.figure.subplots_adjust(left=0.05, right=0.95,
                                                bottom=0.06, top=0.95,
                                                hspace=0.1, wspace=0.1)
        self.qlookcanvas.draw()

    def plot_pg_eovsamap(self):
        """This is to plot the eovsamap with the pyqtgraph's ImageView Widget"""
        if os.path.exists(self.fname):
            rmap, rdata, rheader, ndim, npol_fits, stokaxis, rfreqs, rfdelts = ndfits.read(self.fname)
            self.rmap = rmap
            self.rdata = rdata
            self.rfreqs = rfreqs
            self.rheader = rheader
            self.bdinfo = bdinfo = ndfits.get_bdinfo(rfreqs, rfdelts)
            self.cfreqs = cfreqs = bdinfo['cfreqs']
            self.cfreqs_all = cfreqs_all = bdinfo['cfreqs_all']
            self.freq_dist = lambda fq: (fq - cfreqs_all[0]) / (cfreqs_all[-1] - cfreqs_all[0])
            self.opencontour = True
            self.clevels = np.array([0.3, 0.7])
            self.calpha = 1.
            self.haseovsamap = True
            nspw = len(self.rfreqs)
            eodate = Time(self.rmap.date.mjd + self.rmap.exposure_time.value / 2. / 24 / 3600, format='mjd')
            ny, nx = self.rmap.data.shape
            x0, x1 = (np.array([1, self.rmap.meta['NAXIS1']]) - self.rmap.meta['CRPIX1']) * self.rmap.meta['CDELT1'] + \
                     self.rmap.meta['CRVAL1']
            y0, y1 = (np.array([1, self.rmap.meta['NAXIS2']]) - self.rmap.meta['CRPIX2']) * self.rmap.meta['CDELT2'] + \
                     self.rmap.meta['CRVAL2']
            dx = self.rmap.meta['CDELT1']
            dy = self.rmap.meta['CDELT2']
            mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)

            # plot the images
            # need to reverse the y axis to emulate matplotlib.pyplot.imshow's origin='lower' option
            self.meocanvas.setImage(self.rdata[:, ::-1, :], xvals=self.cfreqs)
            ## Set a custom color map
            colors = [
                (0, 0, 0),
                (45, 5, 61),
                (84, 42, 55),
                (150, 87, 60),
                (208, 171, 141),
                (255, 255, 255)
            ]
            cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=colors)
            self.meocanvas.setColorMap(cmap)
            nf, nx, ny = self.rdata.shape
            self.roi = pg.RectROI([nx / 2, ny / 2], [20, 20], pen=(0, 9))
            # self.roi.sigRegionChanged.connect(self.update_spec)
            self.roi.sigRegionChanged.connect(self.update_pgspec)
            self.meocanvas.addItem(self.roi)

    def update_pgspec(self):
        """Use Pyqtgraph's PlotWidget for the spectral plot"""
        self.freqghz_bound = [1.0, 18.0]
        self.subim = self.roi.getArrayRegion(self.rdata[:, ::-1, :], self.meocanvas.getImageItem(), axes=(2, 1))
        # print(self.subim.shape)
        tb = np.nanmax(self.subim, axis=(1, 2))
        tb_ma = ma.masked_less_equal(tb, 0)
        freqghz = self.cfreqs
        freqghz_ma = ma.masked_outside(freqghz, self.freqghz_bound[0], self.freqghz_bound[1])
        mask_fit = np.logical_or(freqghz_ma.mask, tb_ma.mask)
        freqghz_ma = ma.masked_array(freqghz, mask_fit)
        tb_ma = ma.masked_array(tb, mask_fit)
        logtb_ma = np.log10(tb_ma)

        tb_err = 0.1 * tb
        logtb_err = tb_err / tb / np.log(10.)
        tb_err_ma = ma.masked_array(tb_err, tb_ma.mask)
        logtb_err_ma = ma.masked_array(logtb_err, tb_ma.mask)
        errobjs = []
        # plt = pg.plot()
        # ax.clear()
        # x = np.linspace(1, 20, 10)
        # for ll in [-1, 0, 1, 2, 3, 4]:
        #    y = 10. ** (-2 * np.log10(x) + ll)
        #    plt.plot(x, y, 'k--', alpha=0.1)

        self.speccanvas.clear()
        dataplot = self.speccanvas.plot(x=freqghz_ma, y=logtb_ma, symbol='o')
        errplot = pg.ErrorBarItem(x=freqghz_ma, y=logtb_ma, top=logtb_err_ma,
                                  bottom=logtb_err_ma, beam=0.5)
        # dataplot.setLogMode(True, True)
        self.speccanvas.addItem(dataplot)
        self.speccanvas.addItem(errplot)

    def update_spec(self):
        """Use Matplotlib.pyplot for the spectral plot"""
        self.freqghz_bound = [1.0, 18.0]
        self.subim = self.roi.getArrayRegion(self.rdata[:, ::-1, :], self.meocanvas.getImageItem(), axes=(2, 1))
        # print(self.subim.shape)
        tb = np.nanmax(self.subim, axis=(1, 2))
        tb_ma = ma.masked_less_equal(tb, 0)
        freqghz = self.cfreqs
        freqghz_ma = ma.masked_outside(freqghz, self.freqghz_bound[0], self.freqghz_bound[1])
        mask_fit = np.logical_or(freqghz_ma.mask, tb_ma.mask)
        freqghz_ma = ma.masked_array(freqghz, mask_fit)
        tb_ma = ma.masked_array(tb, mask_fit)

        tb_err = tb * 0.0
        tb_err[:] = 1.e6
        tb_err_ma = ma.masked_array(tb_err, tb_ma.mask)
        errobjs = []
        ax = self.spec_axs
        ax.clear()
        x = np.linspace(1, 20, 10)
        for ll in [-1, 0, 1, 2, 3, 4]:
            y = 10. ** (-2 * np.log10(x) + ll)
            ax.plot(x, y, 'k--', alpha=0.1)

        self.speccanvas.figure.subplots_adjust(left=0.15, right=0.95,
                                               bottom=0.15, top=0.95)
        self.spec_axs.errorbar(freqghz_ma, tb_ma / 1.e6, yerr=tb_err_ma / 1.e6, marker='.', ms=1,
                               linestyle='', c='k')
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim([1, 20])
        ax.set_ylim([0.1, 1000])
        ax.set_xticks([1, 5, 10, 20])
        ax.set_xticklabels([1, 5, 10, 20])
        ax.set_xticks([1, 5, 10, 20])
        ax.set_yticks([])
        ax.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
        ax.set_ylabel('T$_b$ [MK]')
        ax.set_xlabel('Frequency [GHz]')

        self.speccanvas.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
