import sys, os
# from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QPushButton, QSlider, QLineEdit, QLabel, \
#    QTextEdit, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QStatusBar, QProgressBar

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import pyqtgraph as pg

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas)
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import patches, cm
from suncasa.io import ndfits
from astropy.time import Time
import numpy as np
import sunpy
from sunpy import map as smap
from astropy.io import fits
import astropy.units as u
import lmfit
from gstools import GSCostFunctions as gscf

# from suncasa.pyGSFIT import gsutils
# from gsutils import ff_emission
import numpy.ma as ma
import warnings

sys.path.append('./')
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
        self.width = 1600
        self.height = 1100
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self._main = QWidget()
        self.setCentralWidget(self._main)
        self.fit_method = 'nelder'
        self.fit_params = lmfit.Parameters()
        self.fit_params.add_many(('Bx100G', 2., True, 0.1, 100., None, None),
                                 ('log_nnth', 7., True, 3., 11, None, None),
                                 ('delta', 4., True, 1., 30., None, None),
                                 ('Emin_keV', 10., True, 1., 100., None, None),
                                 ('Emax_MeV', 10., False, 0.05, 100., None, None),
                                 ('theta', 45., False, 0.01, 89.9, None, None),
                                 ('log_nth', 10, True, 4., 13., None, None),
                                 ('T_MK', 10., False, 0.1, 100, None, None),
                                 ('depth_Mm', 10., False, 1., 100., None, None))
        self.fit_function = gscf.SinglePowerLawMinimizerOneSrc
        self.threadpool = QThreadPool()
        self.has_eovsamap = False
        self.has_aiamap = False
        self.has_rms = False
        self.has_rois = False
        self.fbar = None
        self.tb_upperbound = 5e9  # lower bound of brightness temperature to consider
        self.tb_lowerbound = 1e3  # lower bound of brightness temperature to consider
        self.freq_bound = [1.0, 18.0]
        self.rois = []
        self.nroi = 0
        self.roi_select_idx = 0
        self.initUI()



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
        # self.tabs.addTab(tab_fit, "Fit")
        self.tabs.addTab(tab_analyzer, "Analyzer")

        # Each tab's user interface is complex, so this splits them into separate functions.
        self.initUItab_explorer()
        # self.initUItab_fit()
        self.initUItab_analyzer()

        # self.tabs.currentChanged.connect(self.tabChanged)

        # Add tabs to widget
        layout.addWidget(self.tabs)
        self._main.setLayout(layout)

        self.show()

    # Explorer Tab User Interface
    #
    def initUItab_explorer(self):
        # Create main layout (a Vertical Layout)
        main_layout = QHBoxLayout()

        # Creat Data Display and Fit Tab
        data_layout = QVBoxLayout()
        fit_layout = QVBoxLayout()

        ###### The Following is for datalayout ######
        # Upper box of the data layout has two hboxes: top box and quicklook plot box
        data_layout_upperbox = QVBoxLayout()
        # Top of the upper box of the data layout is for file selection and fits header display
        file_selection_box = QHBoxLayout()

        # EOVSA FITS Filename entry
        eo_selection_box = QHBoxLayout()
        file_selection_box.addLayout(eo_selection_box)
        # Create Browse button
        eo_selection_button = QPushButton("Load EOVSA")
        eo_selection_button.clicked.connect(self.eofile_select)
        eo_selection_box.addWidget(eo_selection_button)

        # Create LineEdit widget for FITS filename
        self.fitsentry.resize(8 * len(self.fname), 20)
        eo_selection_box.addWidget(self.fitsentry)

        # AIA FITS Filename entry
        aia_selection_box = QHBoxLayout()
        file_selection_box.addLayout(aia_selection_box)
        # Create Browse button
        aia_selection_button = QPushButton("Load AIA")
        aia_selection_button.clicked.connect(self.aiafile_select)
        aia_selection_box.addWidget(aia_selection_button)
        # Create LineEdit widget for AIA FITS file
        self.aiafitsentry = QLineEdit()
        self.aiafitsentry.resize(8 * len(self.fname), 20)
        aia_selection_box.addWidget(self.aiafitsentry)

        # Add top box of the upper box in data layout
        data_layout_upperbox.addLayout(file_selection_box)

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

        # Bottom of the upper box in data layout: quick look plotting area
        qlookarea = QVBoxLayout()
        self.qlookcanvas = FigureCanvas(Figure(figsize=(8, 4)))
        self.qlooktoolbar = NavigationToolbar(self.qlookcanvas, self)
        qlookarea.addWidget(self.qlooktoolbar)
        qlookarea.addWidget(self.qlookcanvas)
        # self.qlook_axs = self.qlookcanvas.figure.subplots(nrows=1, ncols=4)
        gs = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[1, 3], left=0.08, right=0.95, wspace=0.15)
        gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0], wspace=0)
        gs2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.05)
        self.qlook_axs = []
        self.qlook_axs.append(self.qlookcanvas.figure.add_subplot(gs1[0]))
        self.qlook_axs.append(self.qlookcanvas.figure.add_subplot(gs2[0]))
        self.qlook_axs.append(self.qlookcanvas.figure.add_subplot(gs2[1], sharex=self.qlook_axs[-1],
                                                                  sharey=self.qlook_axs[-1]))
        self.qlook_axs.append(self.qlookcanvas.figure.add_subplot(gs2[2], sharex=self.qlook_axs[-1],
                                                                  sharey=self.qlook_axs[-1]))

        data_layout_upperbox.addLayout(qlookarea)
        data_layout.addLayout(data_layout_upperbox)
        # data_layout.setRowStretch(0, 1.)

        # lowerbox
        data_layout_lowerbox = QGridLayout()  # lower box has two hboxes: left for multi-panel display and right for spectrum
        pgimgbox = QVBoxLayout()  # left box of lower box for pg ImageView and associated buttons
        pgbuttonbox = QHBoxLayout()  # bottom of pgimgbox for a bunch or buttons

        # RMS Region Selection Button
        self.rms_selection_button = QPushButton("Select RMS Rgn")
        self.rms_selection_button.setStyleSheet("background-color : lightgrey")
        self.rms_selection_button.setCheckable(True)
        self.rms_selection_button.clicked.connect(self.rms_rgn_select)
        pgbuttonbox.addWidget(self.rms_selection_button)

        # ADD ROI Region for Obtaining Spectra
        self.add_roi_button = QPushButton("Add ROI")
        self.add_roi_button.setStyleSheet("background-color : lightgrey")
        pgbuttonbox.addWidget(self.add_roi_button)
        self.add_roi_button.clicked.connect(self.roi_select)

        ## todo: add this button later
        self.add_roi_grid_button = QPushButton("Add ROI Grid")
        self.add_roi_grid_button.setStyleSheet("background-color : lightgrey")

        pgbuttonbox.addWidget(self.add_roi_grid_button)

        # Add plotting area for multi-panel EOVSA images
        self.pg_img_canvas = pg.ImageView(name='EOVSA Explorer')

        pgimgbox.addWidget(self.pg_img_canvas)
        pgimgbox.addLayout(pgbuttonbox)
        data_layout_lowerbox.addLayout(pgimgbox, 0, 0)
        data_layout_lowerbox.setColumnStretch(0, 1.5)
        self.pg_img_canvas.sigTimeChanged.connect(self.update_fbar)

        # right box for spectral plots
        specplotarea = QVBoxLayout()
        # self.speccanvas = FigureCanvas(Figure(figsize=(4, 6)))
        self.speccanvas = pg.PlotWidget()
        # self.spectoolbar = NavigationToolbar(self.speccanvas, self)
        # specplotarea.addWidget(self.spectoolbar)
        specplotarea.addWidget(self.speccanvas)
        # self.spec_axs = self.speccanvas.figure.subplots(nrows=1, ncols=1)
        data_layout_lowerbox.addLayout(specplotarea, 0, 1)
        data_layout_lowerbox.setColumnStretch(1, 1)

        data_layout.addLayout(data_layout_lowerbox)
        main_layout.addLayout(data_layout)

        ####### The following is for fit layout on the right of the main window ######
        # Group 1: ROI Selection Group
        roi_select_group = QGroupBox("ROI Selection")
        roi_group_box = QVBoxLayout()
        roi_select_group.setLayout(roi_group_box)
        roi_select_box = QHBoxLayout()

        # Slider for ROI selection
        self.roi_select_slider = QSlider(Qt.Horizontal)
        self.roi_select_slider.setMinimum(0)
        self.roi_select_slider.setMaximum(1)
        self.roi_select_slider.setSingleStep(1)
        self.roi_select_slider.setTickPosition(QSlider.TicksBelow)
        self.roi_select_slider.setTickInterval(1)
        self.roi_select_slider.valueChanged.connect(self.roi_slider_valuechange)
        roi_select_box.addWidget(QLabel("ROI ID #"))
        roi_select_box.addWidget(self.roi_select_slider)
        roi_group_box.addLayout(roi_select_box)

        # Text Box for ROI information
        self.roi_info = QTextEdit()
        self.roi_info.setReadOnly(True)
        self.roi_info.setMinimumHeight(50)
        self.roi_info.setMaximumHeight(100)
        # self.roi_info.setMinimumWidth(200)
        # self.roi_info.setMaximumWidth(200)
        # roi_group_box.addWidget(QLabel("ROI Information"))
        roi_group_box.addWidget(self.roi_info)

        # Do fit button for selected ROI
        do_spec_fit_button = QPushButton("Fit Selected ROI")
        do_spec_fit_button.clicked.connect(self.do_spec_fit)
        roi_group_box.addWidget(do_spec_fit_button)

        fit_layout.addWidget(roi_select_group)


        ####### Group 2: Fit Settings Group
        fit_setting_group = QGroupBox("Fit Settings")
        fit_setting_box = QVBoxLayout()
        fit_setting_group.setLayout(fit_setting_box)
        freq_bound_box = QHBoxLayout()

        # Group 2: Lower Frequency Bounds for fit
        self.freq_lowbound_slider = QSlider(Qt.Horizontal)
        self.freq_lowbound_slider.setMinimum(1)
        self.freq_lowbound_slider.setMaximum(18)
        self.freq_lowbound_slider.setValue(1)
        self.freq_lowbound_slider.setSingleStep(0.2)
        self.freq_lowbound_slider.setTickPosition(QSlider.TicksBelow)
        self.freq_lowbound_slider.setTickInterval(1)
        self.freq_lowbound_slider.valueChanged.connect(self.freq_lowbound_valuechange)
        freq_bound_box.addWidget(QLabel("Min Freq"))
        freq_bound_box.addWidget(self.freq_lowbound_slider)

        # Group 2: Upper Frequency Bounds for fit
        self.freq_hibound_slider = QSlider(Qt.Horizontal)
        self.freq_hibound_slider.setMinimum(1)
        self.freq_hibound_slider.setMaximum(18)
        self.freq_hibound_slider.setValue(18)
        self.freq_hibound_slider.setSingleStep(0.2)
        self.freq_hibound_slider.setTickPosition(QSlider.TicksBelow)
        self.freq_hibound_slider.setTickInterval(1)
        self.freq_hibound_slider.valueChanged.connect(self.freq_hibound_valuechange)
        freq_bound_box.addWidget(QLabel("Max Freq"))
        freq_bound_box.addWidget(self.freq_hibound_slider)

        # Group 2: Fit Method
        fit_method_box = QHBoxLayout()
        self.fit_method_selector_widget = QComboBox()
        self.fit_method_selector_widget.addItems(["nelder", "leastsq", "differential_evolution", "basinhopping"])
        self.fit_method_selector_widget.currentIndexChanged.connect(self.fit_method_selector)
        fit_method_box.addWidget(QLabel("Fit Method"))
        fit_method_box.addWidget(self.fit_method_selector_widget)

        # Group 2: Electron Function
        ele_function_box = QHBoxLayout()
        self.ele_function_selector_widget = QComboBox()
        self.ele_function_selector_widget.addItems(
            ["powerlaw", "double_Powerlaw", "thermal f-f + gyrores", "thermal f-f"])
        self.ele_function_selector_widget.currentIndexChanged.connect(self.ele_function_selector)
        ele_function_box.addWidget(QLabel("Electron Dist. Function"))
        ele_function_box.addWidget(self.ele_function_selector_widget)

        fit_setting_box.addLayout(freq_bound_box)
        fit_setting_box.addLayout(fit_method_box)
        fit_setting_box.addLayout(ele_function_box)
        fit_layout.addWidget(fit_setting_group)

        ##### Group 3: Fit Parameters Group
        fit_param_group = QGroupBox("Fit Parameters")
        fit_param_box = QGridLayout()
        fit_param_group.setLayout(fit_param_box)
        nparam = len(self.fit_params)

        self.param_init_value_widgets = []
        self.param_vary_widgets = []
        self.param_min_widgets = []
        self.param_max_widgets = []
        self.param_fit_value_widgets = []

        fit_param_box.addWidget(QLabel('Name'), 0, 0)
        fit_param_box.addWidget(QLabel('Initial Guess'), 0, 1)
        fit_param_box.addWidget(QLabel('Vary'), 0, 2)
        fit_param_box.addWidget(QLabel('Minimum'), 0, 3)
        fit_param_box.addWidget(QLabel('Maximum'), 0, 4)
        fit_param_box.addWidget(QLabel('Fit Results'), 0, 5)
        for n, key in enumerate(self.fit_params):
            #param_layout = QHBoxLayout()
            #param_layout.addWidget(QLabel(key))
            param_init_value_widget = QDoubleSpinBox()
            param_init_value_widget.setDecimals(1)
            param_init_value_widget.setValue(self.fit_params[key].init_value)
            param_init_value_widget.valueChanged.connect(self.update_params)

            param_vary_widget = QCheckBox()
            param_vary_widget.setChecked(self.fit_params[key].vary)
            param_vary_widget.toggled.connect(self.update_params)

            param_min_widget = QDoubleSpinBox()
            param_min_widget.setDecimals(1)
            param_min_widget.setValue(self.fit_params[key].min)
            param_min_widget.valueChanged.connect(self.update_params)

            param_max_widget = QDoubleSpinBox()
            param_max_widget.setDecimals(1)
            param_max_widget.setValue(self.fit_params[key].max)
            param_max_widget.valueChanged.connect(self.update_params)

            param_fit_value_widget = QDoubleSpinBox()
            param_fit_value_widget.setDecimals(1)
            param_fit_value_widget.setValue(self.fit_params[key].value)
            param_fit_value_widget.valueChanged.connect(self.update_params)

            fit_param_box.addWidget(QLabel(key), n+1, 0)
            fit_param_box.addWidget(param_init_value_widget, n+1, 1)
            fit_param_box.addWidget(param_vary_widget, n+1, 2)
            fit_param_box.addWidget(param_min_widget, n+1, 3)
            fit_param_box.addWidget(param_max_widget, n+1, 4)
            fit_param_box.addWidget(param_fit_value_widget, n+1, 5)

            self.param_init_value_widgets.append(param_init_value_widget)
            self.param_vary_widgets.append(param_vary_widget)
            self.param_min_widgets.append(param_min_widget)
            self.param_max_widgets.append(param_max_widget)
            self.param_fit_value_widgets.append(param_fit_value_widget)

        fit_layout.addWidget(fit_param_group)
        main_layout.addLayout(fit_layout)
        # data_layout.setRowStretch(1, 1.2)
        self.tabs.widget(0).setLayout(main_layout)

    # Fit Tab User Interface
    #
    # def initUItab_fit(self):
    #    # Create main layout (a Vertical Layout)
    #    mainlayout = QVBoxLayout()
    #    mainlayout.addWidget(QLabel("This is the tab for analyzing fit results"))
    #    self.tabs.widget(1).setLayout(mainlayout)

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
        self.init_pgspecplot()

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

    def init_pgspecplot(self):
        """ Initial Spectral Plot if no data has been loaded """
        xticksv = list(range(0, 4))
        self.xticks = []
        xticksv_minor = []
        for v in xticksv:
            vp = 10 ** v
            xticksv_minor += list(np.arange(2 * vp, 10 * vp, 1 * vp))
            self.xticks.append((v, '{:.0f}'.format(vp)))
        self.xticks_minor = []
        for v in xticksv_minor:
            self.xticks_minor.append((np.log10(v), ''))
        yticksv = list(range(1, 15))
        yticksv_minor = []
        for v in yticksv:
            vp = 10 ** v
            yticksv_minor += list(np.arange(2 * vp, 10 * vp, 1 * vp))
        self.yticks_minor = []
        for v in yticksv_minor:
            self.yticks_minor.append((np.log10(v), ''))
        self.yticks = []
        for v in yticksv:
            if v >= 6:
                self.yticks.append([v, r'{:.0f}'.format(10 ** v / 1e6)])
            if v == 1:
                self.yticks.append([v, r'{:.5f}'.format(10 ** v / 1e6)])
            elif v == 2:
                self.yticks.append([v, r'{:.4f}'.format(10 ** v / 1e6)])
            elif v == 3:
                self.yticks.append([v, r'{:.3f}'.format(10 ** v / 1e6)])
            elif v == 4:
                self.yticks.append([v, r'{:.2f}'.format(10 ** v / 1e6)])
            elif v == 5:
                self.yticks.append([v, r'{:.1f}'.format(10 ** v / 1e6)])
            else:
                self.yticks.append([v, r'{:.0f}'.format(10 ** v / 1e6)])

        self.update_pgspec()

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
            self.has_eovsamap = True
            nspw = len(self.rfreqs)
            bds = np.linspace(0, nspw, 5)[1:4].astype(np.int)
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
                ax.set_xlabel('Solar X [arcsec]')
                ax.set_title('')
                if n == 0:
                    ax.set_ylabel('Solar Y [arcsec]')
                else:
                    ax.set_ylabel('')
                    ax.set_yticklabels([])
                ax.text(0.01, 0.98, '{0:.1f} GHz'.format(cfreq), ha='left', va='top',
                        fontsize=10, color='w', transform=ax.transAxes)
                ax.set_aspect('equal')
        else:
            self.statusBar.showMessage('EOVSA FITS file does not exist', 2000)
            self.fname = '<Select or enter a valid fits filename>'
            self.fitsentry.setText(self.fname)
            # self.infoEdit.setPlainText('')
            self.has_eovsamap = False

        if os.path.exists(self.aiafname):
            try:
                aiacmap = plt.get_cmap('gray_r')
                aiamap = smap.Map(self.aiafname)
                aiamap.plot(axes=ax0, cmap=aiacmap)
                self.has_aiamap = True
            except:
                self.statusBar.showMessage('Something is wrong with the provided AIA FITS file', 2000)
        else:
            self.statusBar.showMessage('AIA FITS file does not exist', 2000)
            self.aiafname = '<Select or enter a valid fits filename>'
            self.aiafitsentry.setText(self.aiafname)
            # self.infoEdit.setPlainText('')
            self.has_aiamap = False
        cts = []
        if self.has_aiamap:
            aiamap.plot(axes=ax0, cmap=aiacmap)
        if self.has_eovsamap:
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
        ax0.set_xlabel('Solar X [arcsec]')
        ax0.set_ylabel('Solar Y [arcsec]')
        # ax0.set_title('')
        ax0.set_aspect('equal')
        # self.qlookcanvas.figure.subplots_adjust(left=0.05, right=0.95,
        #                                        bottom=0.06, top=0.95,
        #                                        hspace=0.1, wspace=0.1)
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
            self.has_eovsamap = True
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
            self.pg_img_canvas.setImage(self.rdata[:, ::-1, :], xvals=self.cfreqs)
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
            self.pg_img_canvas.setColorMap(cmap)
            # nf, nx, ny = self.rdata.shape
            # self.roi = pg.RectROI([nx / 2 - 10, ny / 2 - 10], [20, 20], pen=(0, 9))
            # self.roi.sigRegionChanged.connect(self.update_spec)
            # self.roi.sigRegionChanged.connect(self.update_pgspec)
            # self.meocanvas.addItem(self.roi)

    def get_roi_specs(self):
        for n, roi in enumerate(self.rois):
            subim = roi.getArrayRegion(self.rdata[:, ::-1, :], self.pg_img_canvas.getImageItem(), axes=(2, 1))
            # print(self.subim.shape)
            roi.freqghz = self.cfreqs
            roi.tb_max = np.nanmax(subim, axis=(1, 2))

        self.update_fitmask()

    def update_fitmask(self):
        for n, roi in enumerate(self.rois):
            tb_max_tofit0 = ma.masked_less_equal(roi.tb_max, self.tb_lowerbound)
            freqghz_tofit0 = ma.masked_outside(self.cfreqs, self.freq_bound[0], self.freq_bound[1])
            roi.mask_tofit = mask_tofit = np.logical_or(freqghz_tofit0.mask, tb_max_tofit0.mask)
            roi.freqghz_tofit = ma.masked_array(self.cfreqs, mask_tofit)
            roi.tb_max_tofit = ma.masked_array(roi.tb_max, mask_tofit)
            if self.has_rms:
                roi.tb_rms_tofit = ma.masked_array(self.rms, roi.mask_tofit)

    def plot_pgspec(self):
        self.update_fitmask()
        if self.has_rms:
            rmsplot = self.speccanvas.plot(x=np.log10(self.cfreqs), y=np.log10(self.rms), pen='k',
                                           symbol='d', symbolPen='k', symbolBrush=None)
            self.speccanvas.addItem(rmsplot)
        self.spec_dataplots = []
        self.spec_dataplots_tofit = []
        for n, roi in enumerate(self.rois):
            if n == self.roi_select_idx:
                symbolfill = (n, 9)
            else:
                symbolfill = None
            spec_dataplot = self.speccanvas.plot(x=np.log10(roi.freqghz), y=np.log10(roi.tb_max), pen=None,
                                                 symbol='o', symbolPen=(n, 9), symbolBrush=None)
            spec_dataplot_tofit = self.speccanvas.plot(x=np.log10(roi.freqghz_tofit), y=np.log10(roi.tb_max_tofit),
                                                       pen=None,
                                                       symbol='o', symbolPen=(n, 9), symbolBrush=symbolfill)
            self.speccanvas.addItem(spec_dataplot)
            self.speccanvas.addItem(spec_dataplot_tofit)
            self.spec_dataplots.append(spec_dataplot)
            self.spec_dataplots_tofit.append(spec_dataplot_tofit)

            # Add errorbar if rms is defined
            if self.has_rms:
                tb_rms = self.rms
                tb_bounds_min = np.maximum(roi.tb_max - tb_rms, np.ones_like(roi.tb_max) * self.tb_lowerbound)
                tb_bounds_max = roi.tb_max + tb_rms
                tb_bounds_min_ma = np.maximum(roi.tb_max_tofit - roi.tb_rms_tofit,
                                              np.ones_like(roi.tb_max_tofit) * self.tb_lowerbound)
                tb_bounds_max_ma = roi.tb_max_tofit + roi.tb_rms_tofit
                errplot = pg.ErrorBarItem(x=np.log10(roi.freqghz), y=np.log10(roi.tb_max),
                                          top=np.log10(tb_bounds_max / roi.tb_max),
                                          bottom=np.log10(roi.tb_max / tb_bounds_min), beam=0.025, pen=(n, 9))
                self.speccanvas.addItem(errplot)

        self.fbar = self.speccanvas.plot(x=np.log10([self.pg_img_canvas.timeLine.getXPos()] * 2), y=[1, 15], pen='k')
        self.speccanvas.addItem(self.fbar)
        self.speccanvas.setLimits(yMin=np.log10(self.tb_lowerbound), yMax=np.log10(self.tb_upperbound))
        xax = self.speccanvas.getAxis('bottom')
        yax = self.speccanvas.getAxis('left')
        xax.setLabel("Frequency [GHz]")
        yax.setLabel("Brightness Temperature [MK]")
        xax.setTicks([self.xticks, self.xticks_minor])
        yax.setTicks([self.yticks, self.yticks_minor])

    def pgspec_add_boundbox(self):
        # add frequency bound
        self.freq_bound_rgn = pg.LinearRegionItem(brush=(20, 50, 200, 20))
        self.freq_bound_rgn.setZValue(10)
        self.freq_bound_rgn.setRegion([np.log10(self.freq_bound[0]), np.log10(self.freq_bound[1])])
        self.speccanvas.addItem(self.freq_bound_rgn)
        self.freq_bound_rgn.sigRegionChanged.connect(self.update_freq_bound_rgn)

    def update_pgspec(self):
        """Use Pyqtgraph's PlotWidget for the spectral plot"""
        self.speccanvas.clear()
        self.plot_pgspec()
        self.pgspec_add_boundbox()

    def update_fbar(self):
        if self.fbar is not None:
            try:
                self.speccanvas.removeItem(self.fbar)
                self.fbar = self.speccanvas.plot(x=np.log10([self.pg_img_canvas.timeLine.getXPos()] * 2), y=[1, 15],
                                                 pen='k')
                self.speccanvas.addItem(self.fbar)
            except:
                pass

    def rms_rgn_select(self):
        """Select a region to calculate rms on spectra"""
        if self.rms_selection_button.isChecked():
            self.rms_selection_button.setStyleSheet("background-color : lightblue")
            self.rms_roi = pg.RectROI([0, 0], [40, 40], pen='w')
            self.pg_img_canvas.addItem(self.rms_roi)
            self.rms_rgn_return()
            self.rms_roi.sigRegionChanged.connect(self.rms_rgn_return)
        else:
            # if rms already exists, remove the ROI box
            self.rms_selection_button.setStyleSheet("background-color : lightgrey")
            self.pg_img_canvas.removeItem(self.rms_roi)

    def rms_rgn_return(self):
        """Select a region to calculate rms on spectra"""
        self.rms = np.std(self.rms_roi.getArrayRegion(self.rdata[:, ::-1, :],
                                                      self.pg_img_canvas.getImageItem(), axes=(2, 1)), axis=(1, 2))
        self.has_rms = True
        self.update_pgspec()

    def roi_select(self):
        """Add a ROI region to the selection"""
        nf, nx, ny = self.rdata.shape
        newroi = pg.RectROI([nx / 2 - 10, ny / 2 - 10], [10, 10], pen=(len(self.rois), 9))
        self.rois.append(newroi)
        self.pg_img_canvas.addItem(self.rois[-1])
        self.roi_select_return()
        newroi.sigRegionChanged.connect(self.roi_select_return)
        self.nroi = len(self.rois)
        self.roi_slider_rangechange()

    def roi_select_return(self):
        self.get_roi_specs()
        self.has_rois = True
        self.update_pgspec()

    def roi_slider_rangechange(self):
        self.roi_select_slider.setMaximum(self.nroi - 1)

    def roi_slider_valuechange(self):
        self.roi_select_idx = self.roi_select_slider.value()
        self.roi_info.setPlainText('Selected ROI {}'.format(self.roi_select_idx))
        self.update_pgspec()

    def freq_lowbound_valuechange(self):
        self.freq_bound[0] = self.freq_lowbound_slider.value()
        self.update_pgspec()

    def freq_hibound_valuechange(self):
        self.freq_bound[1] = self.freq_hibound_slider.value()
        self.update_pgspec()

    def update_freq_bound_rgn(self):
        self.freq_bound_rgn.setZValue(10)
        min_logfreq, max_logfreq = self.freq_bound_rgn.getRegion()
        self.freq_bound = [10. ** min_logfreq, 10. ** max_logfreq]
        self.update_fitmask()
        for spec_dataplot, spec_dataplot_tofit in zip(self.spec_dataplots, self.spec_dataplots_tofit):
            self.speccanvas.removeItem(spec_dataplot)
            self.speccanvas.removeItem(spec_dataplot_tofit)
        self.plot_pgspec()

    def fit_method_selector(self):
        print("Selected Fit Method is: {}".format(self.fit_method_selector_widget.currentText()))
        self.fit_method = self.fit_method_selector_widget.currentText()

    def ele_function_selector(self):
        print("Selected Electron Distribution Function is: {}".format(self.ele_function_selector_widget.currentText()))
        self.ele_dist = self.ele_function_selector_widget.currentText()
        self.init_params()


    def init_params(self):
        if self.ele_dist == 'powerlaw':
            self.fit_params = lmfit.Parameters()
            self.fit_params.add_many(('Bx100G', 2., True, 0.1, 100., None, None),
                                     ('log_nnth', 7., True, 3., 11, None, None),
                                     ('delta', 4., True, 1., 30., None, None),
                                     ('Emin_keV', 10., True, 1., 100., None, None),
                                     ('Emax_MeV', 10., False, 0.05, 100., None, None),
                                     ('theta', 45., True, 0.01, 89.9, None, None),
                                     ('log_nth', 10, True, 4., 13., None, None),
                                     ('T_MK', 10., False, 0.1, 100, None, None),
                                     ('depth_Mm', 10., False, 1., 100., None, None))
            self.fit_function = gscf.SinglePowerLawMinimizerOneSrc
        #print(self.fit_params)

    def update_params(self):
        print('==========Parameters Updated To the Following=======')
        for n, key in enumerate(self.fit_params):
            self.fit_params[key].init_value = self.param_init_value_widgets[n].value()
            self.fit_params[key].vary = self.param_vary_widgets[n].isChecked()
            self.fit_params[key].min = self.param_min_widgets[n].value()
            self.fit_params[key].max = self.param_max_widgets[n].value()
            self.fit_params[key].value = self.param_init_value_widgets[n].value()
            print(key, ':', self.fit_params[key], self.fit_params[key].vary)

    def do_spec_fit(self):
        if hasattr(self, 'spec_fitplot'):
            self.speccanvas.removeItem(self.spec_fitplot)
        freqghz_tofit = self.rois[self.roi_select_idx].freqghz_tofit.compressed()
        tb_max_tofit = self.rois[self.roi_select_idx].tb_max_tofit.compressed()
        tb_rms_tofit = self.rois[self.roi_select_idx].tb_rms_tofit.compressed()

        mini = lmfit.Minimizer(gscf.SinglePowerLawMinimizerOneSrc, self.fit_params,
                               fcn_args=(freqghz_tofit,),
                               fcn_kws={'tb': tb_max_tofit, 'tb_err': tb_rms_tofit},
                               nan_policy='omit')
        method = self.fit_method
        mi = mini.minimize(method=method)
        print(method + ' minimization results')
        print(lmfit.fit_report(mi.params))
        self.fit_params_res = mi.params
        print('==========Update Fit Parameter Results=======')
        for n, key in enumerate(self.fit_params_res):
            self.param_fit_value_widgets[n].setValue(self.fit_params_res[key].value)

        freqghz_toplot = np.logspace(0, np.log10(20.), 100)
        tb_max_fit_res = self.fit_function(mi.params, freqghz_toplot)
        self.spec_fitplot = self.speccanvas.plot(x=np.log10(freqghz_toplot), y=np.log10(tb_max_fit_res),
                                                 pen=dict(color=pg.mkColor(self.roi_select_idx), width=4),
                                             symbol=None, symbolBrush=None)
        self.speccanvas.addItem(self.spec_fitplot)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
