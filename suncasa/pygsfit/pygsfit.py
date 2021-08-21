import sys, os, copy
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
import numpy.ma as ma
import sunpy
from sunpy import map as smap
import astropy
from astropy.io import fits
import astropy.units as u
from astropy import wcs
import lmfit
filedir = os.path.dirname(os.path.realpath(__file__))
print(filedir)
sys.path.append(filedir)
import gstools
import roi_utils
import warnings

warnings.filterwarnings("ignore")
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(imageAxisOrder='row-major')

SMALL_SIZE = 8
MEDIUM_SIZE = 9
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fit_param_text = {'Bx100G': 'B [\u00d7100 G]',
                  'log_nnth': 'log(n<sub>nth</sub>) [cm<sup>-3</sup>]',
                  'delta': '\u03b4',
                  'Emin_keV': 'E<sub>min</sub> [keV]',
                  'Emax_MeV': 'E<sub>max</sub> [MeV]',
                  'theta': '\u03b8 [degree]',
                  'log_nth': 'log(n<sub>th</sub>) [cm<sup>-3</sup>]',
                  'T_MK': 'T [MK]',
                  'depth_asec': 'depth [arcsec]',
                  'area_asec2': 'area [arcsec<sup>2</sup>]'}


class LineSegmentROIX(pg.ROI):
    r"""
    ROI subclass with two freely-moving handles defining a line.

    ============== =============================================================
    **Arguments**
    positions      (list of two length-2 sequences) The endpoints of the line
                   segment. Note that, unlike the handle positions specified in
                   other ROIs, these positions must be expressed in the normal
                   coordinate system of the ROI, rather than (0 to 1) relative
                   to the size of the ROI.
    \**args        All extra keyword arguments are passed to ROI()
    ============== =============================================================
    """

    def __init__(self, positions=(None, None), pos=None, handles=(None, None), **args):
        if pos is None:
            pos = [0, 0]

        pg.ROI.__init__(self, pos, [1, 1], **args)
        if len(positions) > 2:
            raise Exception("LineSegmentROI must be defined by exactly 2 positions. For more points, use PolyLineROI.")

        for i, p in enumerate(positions):
            self.addFreeHandle(p, item=handles[i])

    @property
    def endpoints(self):
        # must not be cached because self.handles may change.
        return [h['item'] for h in self.handles]

    def listPoints(self):
        return [p['item'].pos() for p in self.handles]

    def getState(self):
        state = pg.ROI.getState(self)
        state['points'] = [pg.Point(h.pos()) for h in self.getHandles()]
        return state

    def saveState(self):
        state = pg.ROI.saveState(self)
        state['points'] = [tuple(h.pos()) for h in self.getHandles()]
        return state

    def setState(self, state):
        pg.ROI.setState(self, state)
        p1 = [state['points'][0][0] + state['pos'][0], state['points'][0][1] + state['pos'][1]]
        p2 = [state['points'][1][0] + state['pos'][0], state['points'][1][1] + state['pos'][1]]
        self.movePoint(self.getHandles()[0], p1, finish=False)
        self.movePoint(self.getHandles()[1], p2)

    def paint(self, p, *args):
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(self.currentPen)
        h1 = self.endpoints[0].pos()
        h2 = self.endpoints[1].pos()
        p.drawLine(h1, h2)

    def boundingRect(self):
        return self.shape().boundingRect()

    def shape(self):
        p = QPainterPath()

        h1 = self.endpoints[0].pos()
        h2 = self.endpoints[1].pos()
        dh = h2 - h1
        if dh.length() == 0:
            return p
        pxv = self.pixelVectors(dh)[1]
        if pxv is None:
            return p

        pxv *= 4

        p.moveTo(h1 + pxv)
        p.lineTo(h2 + pxv)
        p.lineTo(h2 - pxv)
        p.lineTo(h1 - pxv)
        p.lineTo(h1 + pxv)

        return p

    def getArrayRegion(self, data, img, axes=(0, 1), order=1, returnMappedCoords=False, **kwds):
        """
        Use the position of this ROI relative to an imageItem to pull a slice
        from an array.

        Since this pulls 1D data from a 2D coordinate system, the return value
        will have ndim = data.ndim-1

        See :meth:`~pytqgraph.ROI.getArrayRegion` for a description of the
        arguments.
        """
        print([h.pos() for h in self.endpoints])
        imgPts = [self.mapToItem(img, h.pos()) for h in self.endpoints]
        rgns = []
        coords = []

        d = pg.Point(imgPts[1] - imgPts[0])
        o = pg.Point(imgPts[0])
        rgn = pg.functions.affineSlice(data, shape=(int(d.length()),), vectors=[pg.Point(d.norm())], origin=o,
                                       axes=axes,
                                       order=order, returnCoords=returnMappedCoords, **kwds)
        print(rgn)

        return rgn


class PolyLineROIX(pg.PolyLineROI):
    # def __init__(self, positions, closed=False, pos=None, **args):
    #     if pos is None:
    #         pos = [0, 0]
    #
    #     self.closed = closed
    #     self.segments = []
    #     pg.ROI.__init__(self, pos, size=[1,1], **args)
    #
    #     self.setPoints(positions)

    @property
    def endpoints(self):
        # must not be cached because self.handles may change.
        return [h['item'] for h in self.handles]

    def listPoints(self):
        return [p['item'].pos() for p in self.handles]

    def getArrayRegion(self, data, img, axes=(0, 1), order=1, returnMappedCoords=False, **kwds):
        """
        Use the position of this ROI relative to an imageItem to pull a slice
        from an array.

        Since this pulls 1D data from a 2D coordinate system, the return value
        will have ndim = data.ndim-1

        See :meth:`~pytqgraph.ROI.getArrayRegion` for a description of the
        arguments.
        """
        # print([h.pos() for h in self.endpoints])
        imgPts = [self.mapToItem(img, h.pos()) for h in self.endpoints]
        rgns = []
        coords = []

        rgns = []
        for idx, imgPt in enumerate(imgPts[:-1]):
            d = pg.Point(imgPts[idx + 1] - imgPts[idx])
            o = pg.Point(imgPts[idx])
            rgn = pg.functions.affineSlice(data, shape=(int(d.length()),), vectors=[pg.Point(d.norm())], origin=o,
                                           axes=axes,
                                           order=order, returnCoords=returnMappedCoords, **kwds)
            rgns.append(rgn)

        return np.hstack(rgns)


class App(QMainWindow):

    def __init__(self):
        super().__init__()

        self.eoimg_fname = '<Select or enter a valid EOVSA image fits file name>'
        self.eodspec_fname = '<Select or enter a valid EOVSA spectrogram fits file name>'
        self.aiafname = '<Select or enter a valid AIA fits file name>'
        # self.eoimg_fitsentry = QLineEdit()
        # self.eodspec_fitsentry = QLineEdit()
        self.title = 'pygsfit'
        self.left = 0
        self.top = 0
        self.width = 1600
        self.height = 1000
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self._main = QWidget()
        self.setCentralWidget(self._main)
        self.fit_method = 'nelder'
        self.fit_params = lmfit.Parameters()
        self.fit_params.add_many(('Bx100G', 2., True, 0.1, 100., None, None),
                                 ('log_nnth', 5., True, 3., 11, None, None),
                                 ('delta', 4., True, 1., 30., None, None),
                                 ('Emin_keV', 10., False, 1., 100., None, None),
                                 ('Emax_MeV', 10., False, 0.05, 100., None, None),
                                 ('theta', 45., True, 0.01, 89.9, None, None),
                                 ('log_nth', 10, True, 4., 13., None, None),
                                 ('T_MK', 1., False, 0.1, 100, None, None),
                                 ('depth_asec', 5., False, 1., 100., None, None))
        self.fit_params_nvarys = 5
        self.fit_kws = {'maxiter': 2000, 'xatol': 0.01, 'fatol': 0.01}
        self.fit_function = gstools.GSCostFunctions.SinglePowerLawMinimizerOneSrc
        self.threadpool = QThreadPool()
        self.has_eovsamap = False
        self.has_dspec = False
        self.has_stokes = False
        self.has_aiamap = False
        self.has_bkg = False
        self.has_rois = False
        self.fbar = None
        self.data_freq_bound = [1.0, 18.0]  # Initial Frequency Bound of the instrument
        self.tb_spec_bound = [1e3, 5e9]  # Bound of brightness temperature; the lower bound is set to the fit limit
        self.flx_spec_bound = [1e-4, 1e5]  # Bound of flux density; the lower bound is set to the fit limit
        self.fit_freq_bound = [1.0, 18.0]
        self.roi_freq_bound = [1.0, 18.0]
        self.spec_frac_err = 0.1  # fractional error of the intensity (assumed to be due to flux calibration error)
        self.spec_rmsplots = []
        self.spec_dataplots = []
        self.spec_dataplots_tofit = []
        self.roi_grid_size = 2
        self.rois = [[]]
        self.roi_group_idx = 0
        self.nroi_current_group = 0
        self.current_roi_idx = 0
        self.distSpecCanvasSet = {}
        self.pol_select_idx = 0
        self.spec_in_tb = True
        self.is_calibrated_tp = True
        self.qlookimg_axs = None
        self.qlookdspec_ax = None
        # some controls for qlookplot
        self.opencontour = True
        self.clevels = np.array([0.3, 0.7])
        self.calpha = 1.
        self.pgcmap = self._create_pgcmap(cmap='viridis', ncolorstop=6)
        # initialize the window
        # self.initUI()
        self.initUItab_explorer()
        # ## quick input for debug --------------
        # self.eoimg_file_select()
        # ## quick input for debug --------------

    def _create_pgcmap(self, cmap='viridis', ncolorstop=6):
        """This is to create the cmap for pyqtgraph's ImageView Widgets"""
        from matplotlib import cm
        mpl_cmap = cm.get_cmap(cmap, ncolorstop)
        mpl_colors = mpl_cmap(np.linspace(0, 1, ncolorstop)) * 255
        colors = []
        for s in range(ncolorstop):
            colors.append((int(mpl_colors[s, 0]), int(mpl_colors[s, 1]), int(mpl_colors[s, 2])))

        cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, ncolorstop), color=colors)
        return cmap

    def _createToolBars(self):
        # Using a title
        # fileToolBar = self.addToolBar("File")
        iconsize = QSize(50, 50)
        selectToolBar = self.addToolBar("ADD")
        self.addSelectAction = QAction(QIcon("{}/resources/add-button.svg".format(filedir)), "&Add ROI", self)
        self.addSelectAction.setText('Add ROI')
        self.addSelectAction.triggered.connect(self.add_new_roi)
        selectToolBar.addAction(self.addSelectAction)

        ## define ROI toolbar
        roiToolBar = self.addToolBar("ROI")

        self.rectButton = QToolButton(self, text="RectROI", checkable=True)
        self.rectButton.setIcon(QIcon("{}/resources/roi-rect.svg".format(filedir)))
        self.rectButton.setChecked(True)
        self.rectButton.setIconSize(iconsize)
        self.rectButton.setToolTip('Rectangle ROI tool')
        self.rectButton.setStyleSheet("QToolButton::hover"
                                      "{"
                                      "background-color : #555555;"
                                      "}")

        self.elpsButton = QToolButton(self, text="EllipseROI", checkable=True)
        self.elpsButton.setIcon(QIcon("{}/resources/roi-ellipse.svg".format(filedir)))
        self.elpsButton.setIconSize(iconsize)
        self.elpsButton.setToolTip('Ellipse ROI tool')
        self.elpsButton.setStyleSheet("QToolButton::hover"
                                      "{"
                                      "background-color : #555555;"
                                      "}")

        self.polygonButton = QToolButton(self, text="PolygonROI", checkable=True)
        self.polygonButton.setIcon(QIcon("{}/resources/roi-polygon.svg".format(filedir)))
        self.polygonButton.setIconSize(iconsize)
        self.polygonButton.setToolTip('Polygon ROI tool')
        self.polygonButton.setStyleSheet("QToolButton::hover"
                                         "{"
                                         "background-color : #555555;"
                                         "}")
        ## define Slice toolbar
        sliceToolBar = self.addToolBar("sliceROI")

        self.lineButton = QToolButton(self, text="LineSegmentROI", checkable=True)
        self.lineButton.setIcon(QIcon("{}/resources/slice-line.svg".format(filedir)))
        self.lineButton.setIconSize(iconsize)
        self.lineButton.setToolTip('Line slice ROI tool')
        self.lineButton.setStyleSheet("QToolButton::hover"
                                      "{"
                                      "background-color : #555555;"
                                      "}")

        self.polyLineButton = QToolButton(self, text="PolyLineROI", checkable=True)
        self.polyLineButton.setIcon(QIcon("{}/resources/slice-polyline.svg".format(filedir)))
        self.polyLineButton.setIconSize(iconsize)
        self.polyLineButton.setToolTip('PolyLine slice ROI tool')
        self.polyLineButton.setStyleSheet("QToolButton::hover"
                                          "{"
                                          "background-color : #555555;"
                                          "}")

        # self.multiRectButton = QToolButton(self, text="MultiRectROI", checkable=True)
        # self.multiRectButton.setIcon(QIcon("{}/resources/slice-line.svg".format(filedir)))
        # self.multiRectButton.setIconSize(iconsize)
        # self.multiRectButton.setToolTip('PolyLine slice tool')
        # self.multiRectButton.setStyleSheet("QToolButton::hover"
        #                                    "{"
        #                                    "background-color : #555555;"
        #                                    "}")

        self.toolBarButtonGroup = QButtonGroup(self, exclusive=True)
        for button in [
            self.rectButton,
            self.elpsButton,
            self.polygonButton
        ]:
            roiToolBar.addWidget(button)
            self.toolBarButtonGroup.addButton(button)
        # for button in [self.lineButton, self.polyLineButton, self.multiRectButton]:
        for button in [self.lineButton, self.polyLineButton]:
            sliceToolBar.addWidget(button)
            self.toolBarButtonGroup.addButton(button)

        self.add2slice = QCheckBox('Add ROI to slice?')
        self.add2slice.setChecked(False)
        self.add2slice.setToolTip(
            'If True, the selected ROI will be added to the last slice ROI. If there no slice ROI exists, use the sliceROI tools on the right to create one.')
        roiToolBar.addWidget(self.add2slice)
        # self.add2slice.toggled.connect(self.is_calibrated_tp_state)

        # self.addSelectAction.set(True)
        # self.rectAction = QAction(QIcon("{}/resources/roi-rect.svg".format(filedir)), "&Rectangle", self)
        # self.rectAction.setCheckable(True)
        # self.circAction = QAction(QIcon("{}/resources/roi-circ.svg".format(filedir)), "&Circle", self)
        # self.circAction.setCheckable(True)
        # self.polylineAction = QAction(QIcon("{}/resources/roi-polyline.svg".format(filedir)), "&PolyLine", self)
        # self.polylineAction.setCheckable(True)
        # roiToolBar.addAction(self.rectAction)
        # roiToolBar.addAction(self.circAction)
        # roiToolBar.addAction(self.polylineAction)
        # # # Using a QToolBar object
        # # editToolBar = QToolBar("Edit", self)
        # # self.addToolBar(editToolBar)
        # # # Using a QToolBar object and a toolbar area
        # # helpToolBar = QToolBar("Image", self)
        # # self.addToolBar(Qt.TopToolBarArea, helpToolBar)

    def _createMenuBar(self):
        menubar = self.menuBar()
        ## PYQT5 does not support native menubar on MacOS
        menubar.setNativeMenuBar(False)
        # fileMenu = QMenu("&File", self)
        # menubar.addMenu(fileMenu)
        actionFile = menubar.addMenu("File")

        action_loadEOVSAimage = QAction("Load EOVSA Image", self)
        action_loadEOVSAimage.triggered.connect(self.eoimg_file_select)
        actionFile.addAction(action_loadEOVSAimage)
        action_loadAIA = QAction("Load AIA", self)
        action_loadAIA.triggered.connect(self.aiafile_select)
        actionFile.addAction(action_loadAIA)
        action_loadEOVSAspectrogram = QAction("Load EOVSA Spectrogram", self)
        action_loadEOVSAspectrogram.triggered.connect(self.eodspec_file_select)
        actionFile.addAction(action_loadEOVSAspectrogram)
        # actionFile.addAction("Open AIA")
        # actionFile.addAction("Open EOVSA Spectrogram")
        actionFile.addSeparator()
        actionFile.addAction("Quit")

        # Creating menus using a title
        editMenu = menubar.addMenu("&Edit")
        # helpMenu = menuBar.addMenu(" &Help")
        # self.menu_layout.addWidget(menuBar)

    # def initUI(self):
    #     self.statusBar = QStatusBar()
    #     self.setStatusBar(self.statusBar)
    #     self.progressBar = QProgressBar()
    #     self.progressBar.setGeometry(10, 10, 200, 15)
    #
    #     layout = QVBoxLayout()
    #     # Initialize tab screen
    #     self.tabs = QTabWidget()
    #     tab_explorer = QWidget()
    #     tab_fit = QWidget()
    #     tab_analyzer = QWidget()
    #
    #     # Add tabs
    #     self.tabs.addTab(tab_explorer, "Explorer")
    #     # self.tabs.addTab(tab_fit, "Fit")
    #     self.tabs.addTab(tab_analyzer, "Analyzer")
    #
    #     # Each tab's user interface is complex, so this splits them into separate functions.
    #     self.initUItab_explorer()
    #     # self.initUItab_fit()
    #     self.initUItab_analyzer()
    #
    #     # self.tabs.currentChanged.connect(self.tabChanged)
    #
    #     # Add tabs to widget
    #     layout.addWidget(self.tabs)
    #     self._main.setLayout(layout)
    #
    #     self.show()

    # Explorer Tab User Interface
    #
    def initUItab_explorer(self):
        # Create main layout (a Vertical Layout)
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.progressBar = QProgressBar()
        self.progressBar.setGeometry(10, 10, 200, 15)
        self._createMenuBar()
        self._createToolBars()
        main_layout = QHBoxLayout()

        # Creat Data Display and Fit Tab
        data_layout = QVBoxLayout()
        fit_layout = QVBoxLayout()

        ###### The Following is for datalayout ######
        # Upper box of the data layout has two hboxes: top box and quicklook plot box
        data_layout_upperbox = QVBoxLayout()
        # # Top of the upper box of the data layout is for file selection and fits header display
        # file_selection_box = QGridLayout()
        #
        # # EOVSA Image FITS Filename entry
        # eoimg_selection_box = QHBoxLayout()
        # file_selection_box.addLayout(eoimg_selection_box, 0, 0)
        # # Create Browse button
        # eoimg_selection_button = QPushButton("Load EOVSA Image")
        # eoimg_selection_button.clicked.connect(self.eoimg_file_select)
        # eoimg_selection_box.addWidget(eoimg_selection_button)
        # # Create LineEdit widget for FITS filename
        # self.eoimg_fitsentry.resize(8 * len(self.eoimg_fname), 15)
        # eoimg_selection_box.addWidget(self.eoimg_fitsentry)
        #
        # # EOVSA Spectrogram FITS Filename entry
        # eodspec_selection_box = QHBoxLayout()
        # file_selection_box.addLayout(eodspec_selection_box, 0, 1)
        # # Create Browse button
        # eodspec_selection_button = QPushButton("Load EOVSA Spectrogram")
        # eodspec_selection_button.clicked.connect(self.eodspec_file_select)
        # eodspec_selection_box.addWidget(eodspec_selection_button)
        # # Create LineEdit widget for FITS filename
        # self.eodspec_fitsentry.resize(8 * len(self.eoimg_fname), 15)
        # eodspec_selection_box.addWidget(self.eodspec_fitsentry)
        #
        # self.is_calibrated_tp_button = QRadioButton('Is Calibrated TP?')
        # self.is_calibrated_tp_button.setChecked(self.is_calibrated_tp)
        # self.is_calibrated_tp_button.toggled.connect(self.is_calibrated_tp_state)
        # eodspec_selection_box.addWidget(self.is_calibrated_tp_button)
        #
        # # AIA FITS Filename entry
        # aia_selection_box = QHBoxLayout()
        # file_selection_box.addLayout(aia_selection_box, 1, 0)
        # # Create Browse button
        # aia_selection_button = QPushButton("Load AIA")
        # aia_selection_button.clicked.connect(self.aiafile_select)
        # aia_selection_box.addWidget(aia_selection_button)
        # # Create LineEdit widget for AIA FITS file
        # self.aiaimg_fitsentry = QLineEdit()
        # self.aiaimg_fitsentry.resize(8 * len(self.eoimg_fname), 15)
        # aia_selection_box.addWidget(self.aiaimg_fitsentry)
        #
        # # Add top box of the upper box in data layout
        # data_layout_upperbox.addLayout(file_selection_box)

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
        data_layout.addLayout(data_layout_upperbox)

        # # toolbox
        # ## ADD ROI type selection
        # self.roi_type_button_box = QHBoxLayout()
        # self.roi_type_selector_widget = QComboBox()
        # self.roi_type_selector_widget.addItems(["Rect", "Line", "differential_evolution"])
        # # self.roi_type_selector_widget.currentIndexChanged.connect(self.fit_method_selector)
        # self.roi_type_button_box.addWidget(QLabel("placeholder for toolbar"))
        # self.roi_type_button_box.addWidget(self.roi_type_selector_widget)
        #
        # # middle box of the
        # data_layout_middlebox = QVBoxLayout()
        # data_layout_middlebox.addLayout(self.roi_type_button_box)
        # data_layout.addLayout(data_layout_middlebox)

        # lowerbox
        data_layout_lowerbox = QGridLayout()  # lower box has two hboxes: left for multi-panel display and right for spectrum
        pgimgbox = QVBoxLayout()  # left box of lower box for pg ImageView and associated buttons

        # Status box
        pg_img_status_box = QVBoxLayout()
        self.pg_img_mouse_pos_widget = QLabel("")
        pg_img_status_box.addWidget(self.pg_img_mouse_pos_widget)
        self.pg_img_roi_info_widget = QLabel("")
        # pg_img_status_box.addWidget(self.pg_img_roi_info_widget)
        self.pg_img_bkg_roi_info_widget = QLabel("")
        pg_img_status_box.addWidget(self.pg_img_bkg_roi_info_widget)

        # Add plotting area for multi-panel EOVSA images
        self.pg_img_plot = pg.PlotItem(labels={'bottom': ('Solar X [arcsec]', ''), 'left': ('Solar Y [arcsec]', '')})
        self.pg_img_canvas = pg.ImageView(name='EOVSA Explorer', view=self.pg_img_plot)

        pgimgbox.addWidget(self.pg_img_canvas)
        # pgimgbox.addLayout(pgbuttonbox)
        pgimgbox.addLayout(pg_img_status_box)
        data_layout_lowerbox.addLayout(pgimgbox, 0, 0)
        data_layout_lowerbox.setColumnStretch(0, 2)
        self.pg_img_canvas.sigTimeChanged.connect(self.update_fbar)

        # right box for spectral plots
        self.specplotarea = QVBoxLayout()
        # self.speccanvas = FigureCanvas(Figure(figsize=(4, 6)))
        self.speccanvas = pg.PlotWidget()
        # self.spectoolbar = NavigationToolbar(self.speccanvas, self)
        # self.specplotarea.addWidget(self.spectoolbar)
        self.specplotarea.addWidget(self.speccanvas)
        # self.spec_axs = self.speccanvas.figure.subplots(nrows=1, ncols=1)
        data_layout_lowerbox.addLayout(self.specplotarea, 0, 1)
        data_layout_lowerbox.setColumnStretch(1, 1)

        # add a toggle between Tb and flux density
        specplotmode_box = QVBoxLayout()
        tb_flx_button_group = QButtonGroup(self)
        self.plot_tb_button = QRadioButton("Brightness Temperature", self)
        self.plot_tb_button.toggled.connect(self.tb_flx_btnstate)
        self.plot_flx_button = QRadioButton("Flux Density", self)
        self.plot_flx_button.toggled.connect(self.tb_flx_btnstate)
        tb_flx_button_group.addButton(self.plot_tb_button)
        tb_flx_button_group.addButton(self.plot_flx_button)
        self.plot_tb_button.setChecked(True)
        # self.plot_flx_button.setChecked(False)
        specplotmode_box.addWidget(self.plot_tb_button)
        specplotmode_box.addWidget(self.plot_flx_button)
        self.specplotarea.addLayout(specplotmode_box)

        data_layout.addLayout(data_layout_lowerbox)

        data_layout_lowerbox2 = QVBoxLayout()
        # Create a button to toggle the qlookimg box.
        qlookimglabel = QLabel("solar Image")
        self.qlookimgbutton = QToolButton()
        self.qlookimgbutton.setArrowType(Qt.RightArrow)
        self.qlookimgbutton.setToolButtonStyle(Qt.ToolButtonIconOnly)
        self.qlookimgbutton.setFixedSize(20, 20)
        self.qlookimgbutton.setCheckable(True)
        self.qlookimgbutton.toggled.connect(self.showqlookimg)

        # Create a button to toggle the qlookspec box.
        qlookdspeclabel = QLabel("Spectrogram")
        self.qlookdspecbutton = QToolButton()
        self.qlookdspecbutton.setArrowType(Qt.RightArrow)
        self.qlookdspecbutton.setToolButtonStyle(Qt.ToolButtonIconOnly)
        self.qlookdspecbutton.setFixedSize(20, 20)
        self.qlookdspecbutton.setCheckable(True)
        self.qlookdspecbutton.setChecked(False)
        self.qlookdspecbutton.toggled.connect(self.showqlookdspec)

        self.is_calibrated_tp_button = QCheckBox('Is Calibrated TP?')
        self.is_calibrated_tp_button.setChecked(self.is_calibrated_tp)
        self.is_calibrated_tp_button.toggled.connect(self.is_calibrated_tp_state)

        qlookbuttonbox = QHBoxLayout()
        qlookbuttonbox_l = QHBoxLayout()
        qlookbuttonbox_l.addWidget(self.qlookimgbutton)
        qlookbuttonbox_l.addWidget(qlookimglabel)
        qlookbuttonbox_r = QHBoxLayout()
        qlookbuttonbox_r.addWidget(self.qlookdspecbutton)
        qlookbuttonbox_r.addWidget(qlookdspeclabel)
        qlookbuttonbox_r.addWidget(self.is_calibrated_tp_button)
        qlookbuttonbox.addLayout(qlookbuttonbox_l)
        qlookbuttonbox.addLayout(qlookbuttonbox_r)
        qlookbuttonbox.setStretch(0, 1)
        qlookbuttonbox.setStretch(1, 2)
        data_layout_lowerbox2.addLayout(qlookbuttonbox)

        # Bottom of the upper box in data layout: quick look plotting area
        qlookarea = QHBoxLayout()

        self.qlookimgbox = QVBoxLayout()
        self.qlookdspecbox = QVBoxLayout()
        self.qlookdummybox = QVBoxLayout()

        dummy_spacer = QLabel('')
        self.qlookdummybox.addWidget(dummy_spacer)
        qlookarea.addLayout(self.qlookimgbox)
        qlookarea.addLayout(self.qlookdspecbox)
        qlookarea.addLayout(self.qlookdummybox)

        data_layout_lowerbox2.addLayout(qlookarea)
        data_layout.addLayout(data_layout_lowerbox2)
        # data_layout.setRowStretch(0, 1.)

        # data_layout.addLayout(data_layout_bottombox)
        main_layout.addLayout(data_layout)

        ####### The following is for fit layout on the right of the main window ######
        # Group 1: ROI Definition Group
        roi_definition_group = QGroupBox("ROI Definition")
        roi_definition_group_box = QVBoxLayout()
        roi_definition_group.setLayout(roi_definition_group_box)

        # ADD ROI Region for Obtaining Spectra
        roi_button_box = QHBoxLayout()
        self.add_roi_button = QPushButton("Add New ROI")
        self.add_roi_button.setStyleSheet("background-color : lightgrey")
        roi_button_box.addWidget(self.add_roi_button)
        self.add_roi_button.clicked.connect(self.add_new_roi)

        self.add_to_roigroup_button = QToolButton(self)
        self.add_to_roigroup_button.setText("0")
        self.add_to_roigroup_button.setPopupMode(QToolButton.MenuButtonPopup)
        self.add_to_roigroup_button.setMenu(QMenu(self.add_to_roigroup_button))
        self.add_to_roigroup_widget = QListWidget()
        self.add_to_roigroup_widget.addItems(['0', '1', '2'])  ##todo: update ROI group # on the fly
        self.add_to_roigroup_widget.itemClicked.connect(self.add_to_roigroup_selection)
        action = QWidgetAction(self.add_to_roigroup_button)
        action.setDefaultWidget(self.add_to_roigroup_widget)
        self.add_to_roigroup_button.menu().addAction(action)
        roi_button_box.addWidget(QLabel("to Group"))
        roi_button_box.addWidget(self.add_to_roigroup_button)

        self.roi_freq_lowbound_selector = QDoubleSpinBox()
        self.roi_freq_lowbound_selector.setDecimals(1)
        self.roi_freq_lowbound_selector.setRange(self.data_freq_bound[0], self.data_freq_bound[1])
        self.roi_freq_lowbound_selector.setSingleStep(0.1)
        self.roi_freq_lowbound_selector.setValue(self.data_freq_bound[0])
        self.roi_freq_lowbound_selector.valueChanged.connect(self.roi_freq_lowbound_valuechange)
        roi_button_box.addWidget(QLabel("Min Freq (GHz)"))
        roi_button_box.addWidget(self.roi_freq_lowbound_selector)

        self.roi_freq_hibound_selector = QDoubleSpinBox()
        self.roi_freq_hibound_selector.setDecimals(1)
        self.roi_freq_hibound_selector.setRange(self.data_freq_bound[0], self.data_freq_bound[1])
        self.roi_freq_hibound_selector.setSingleStep(0.1)
        self.roi_freq_hibound_selector.setValue(self.data_freq_bound[1])
        self.roi_freq_hibound_selector.valueChanged.connect(self.roi_freq_hibound_valuechange)
        roi_button_box.addWidget(QLabel("Max Freq (GHz)"))
        roi_button_box.addWidget(self.roi_freq_hibound_selector)

        # ADD presets selection box
        self.add_manual_rois_button = QToolButton()
        self.add_manual_rois_button.setText('Define ROIs')
        self.add_manual_rois_button.clicked.connect(self.add_manually_defined_rois)
        self.add_manual_rois_button.setPopupMode(QToolButton.MenuButtonPopup)
        self.group_roi_op_menu = QMenu()
        self.group_roi_op_menu.addAction('Save Group', self.group_roi_op_selector)
        self.group_roi_op_menu.addAction('Load Group', self.group_roi_op_selector)
        self.add_manual_rois_button.setMenu(self.group_roi_op_menu)
        roi_button_box.addWidget(self.add_manual_rois_button)
        roi_definition_group_box.addLayout(roi_button_box)

        roi_grid_box = QHBoxLayout()
        ## todo: add this button to draw a box and define a group of ROIs with a given size
        self.add_roi_grid_button = QPushButton("Add ROI Grid")
        self.add_roi_grid_button.setStyleSheet("background-color : lightgrey")
        roi_grid_box.addWidget(self.add_roi_grid_button)

        self.roi_grid_size_selector = QSpinBox()
        self.roi_grid_size_selector.setRange(1, 1000)
        self.roi_grid_size_selector.setSingleStep(1)
        self.roi_grid_size_selector.setValue(2)
        self.roi_grid_size_selector.valueChanged.connect(self.roi_grid_size_valuechange)
        roi_grid_box.addWidget(QLabel("Grid Size (pix)"))
        roi_grid_box.addWidget(self.roi_grid_size_selector)
        roi_definition_group_box.addLayout(roi_grid_box)

        fit_layout.addWidget(roi_definition_group)
        # Group 2: ROI Selection Group
        roi_selection_group = QGroupBox("ROI Selection")
        roi_selection_group_box = QVBoxLayout()
        roi_selection_group.setLayout(roi_selection_group_box)

        roi_selection_button_box = QHBoxLayout()

        # ROI GROUP Selection Button
        self.roigroup_selection_button = QToolButton(self)
        self.roigroup_selection_button.setText('0')
        self.roigroup_selection_button.setPopupMode(QToolButton.MenuButtonPopup)
        self.roigroup_selection_button.setMenu(QMenu(self.roigroup_selection_button))
        self.roigroup_selection_widget = QListWidget()
        self.roigroup_selection_widget.addItems(['0', '1', '2'])  ##todo: update ROI group list on the fly
        self.roigroup_selection_widget.itemClicked.connect(self.roigroup_selection)
        action = QWidgetAction(self.roigroup_selection_button)
        action.setDefaultWidget(self.roigroup_selection_widget)
        self.roigroup_selection_button.menu().addAction(action)
        roi_selection_button_box.addWidget(QLabel("Select ROI Group"))
        roi_selection_button_box.addWidget(self.roigroup_selection_button)

        # ROI Selection Button
        self.roi_selection_button = QToolButton(self)
        self.roi_selection_button.setText(str(self.current_roi_idx))
        self.roi_selection_button.setPopupMode(QToolButton.MenuButtonPopup)
        self.roi_selection_button.setMenu(QMenu(self.roi_selection_button))
        self.roi_selection_widget = QListWidget()
        self.roi_selection_widget.addItems([str(i) for i in range(self.nroi_current_group)])
        self.roi_selection_widget.itemClicked.connect(self.roi_selection_action)
        action = QWidgetAction(self.roi_selection_button)
        action.setDefaultWidget(self.roi_selection_widget)
        self.roi_selection_button.menu().addAction(action)
        roi_selection_button_box.addWidget(QLabel("Select ROI"))
        roi_selection_button_box.addWidget(self.roi_selection_button)
        # Slider for ROI selection
        # self.roi_select_slider = QSlider(Qt.Horizontal)
        # self.roi_select_slider.setMinimum(0)
        # self.roi_select_slider.setMaximum(1)
        # self.roi_select_slider.setSingleStep(1)
        # self.roi_select_slider.setTickPosition(QSlider.TicksBelow)
        # self.roi_select_slider.setTickInterval(1)
        # self.roi_select_slider.valueChanged.connect(self.roi_slider_valuechange)
        # roi_select_box.addWidget(QLabel("ROI ID #"))
        # roi_select_box.addWidget(self.roi_select_slider)
        roi_selection_group_box.addLayout(roi_selection_button_box)

        # Text Box for ROI information
        self.roi_info = QTextEdit()
        self.roi_info.setReadOnly(True)
        self.roi_info.setMinimumHeight(50)
        self.roi_info.setMaximumHeight(100)
        # self.roi_info.setMinimumWidth(200)
        # self.roi_info.setMaximumWidth(200)
        # roi_group_box.addWidget(QLabel("ROI Information"))
        roi_selection_group_box.addWidget(self.roi_info)

        # Layout for computing and applying total power calibration factor
        tpcal_correction_box = QHBoxLayout()
        # Button for combining flux from a ROI group
        combine_flux_button = QPushButton("Combine Flux of Current ROI Group")
        combine_flux_button.clicked.connect(self.combine_roi_group_flux)
        tpcal_correction_box.addWidget(combine_flux_button)

        # Button for calculating total power calibration factor
        calc_tpcal_factor_button = QPushButton("Compute TP Cal factor")
        calc_tpcal_factor_button.clicked.connect(self.calc_tpcal_factor)
        tpcal_correction_box.addWidget(calc_tpcal_factor_button)

        # Button for applying total power calibration factor
        self.apply_tpcal_factor_button = QRadioButton("Apply?")
        self.apply_tpcal_factor_button.clicked.connect(self.apply_tpcal_factor)
        tpcal_correction_box.addWidget(self.apply_tpcal_factor_button)

        roi_selection_group_box.addLayout(tpcal_correction_box)

        ## Add buttons for doing spectral fit
        fit_button_box = QHBoxLayout()
        # Background selection button
        # self.bkg_selection_button = QPushButton("Select Background")
        # self.bkg_selection_button.setStyleSheet("background-color : lightgrey")
        # self.bkg_selection_button.setCheckable(True)
        # self.bkg_selection_button.clicked.connect(self.bkg_rgn_select)
        # fit_button_box.addWidget(self.bkg_selection_button)

        # Do Fit Button
        do_spec_fit_button = QPushButton("Fit Selected ROI")
        do_spec_fit_button.clicked.connect(self.do_spec_fit)
        fit_button_box.addWidget(do_spec_fit_button)

        # Add fit button box to roi group
        roi_selection_group_box.addLayout(fit_button_box)

        fit_layout.addWidget(roi_selection_group)

        ####### Group 2: Fit Settings Group
        fit_setting_group = QGroupBox("Fit Settings")
        fit_setting_box = QVBoxLayout()
        fit_setting_group.setLayout(fit_setting_box)
        spec_adjust_box = QHBoxLayout()

        # Group 2: Lower Frequency Bounds for fit
        self.freq_lowbound_selector = QDoubleSpinBox()
        self.freq_lowbound_selector.setDecimals(1)
        self.freq_lowbound_selector.setRange(self.data_freq_bound[0], self.data_freq_bound[1])
        self.freq_lowbound_selector.setSingleStep(0.1)
        self.freq_lowbound_selector.setValue(self.data_freq_bound[0])
        self.freq_lowbound_selector.valueChanged.connect(self.freq_lowbound_valuechange)
        spec_adjust_box.addWidget(QLabel("Min Freq (GHz)"))
        spec_adjust_box.addWidget(self.freq_lowbound_selector)

        # Group 2: Upper Frequency Bounds for fit
        self.freq_hibound_selector = QDoubleSpinBox()
        self.freq_hibound_selector.setDecimals(1)
        self.freq_hibound_selector.setRange(self.data_freq_bound[0], self.data_freq_bound[1])
        self.freq_hibound_selector.setSingleStep(0.1)
        self.freq_hibound_selector.setValue(self.data_freq_bound[1])
        self.freq_hibound_selector.valueChanged.connect(self.freq_hibound_valuechange)
        spec_adjust_box.addWidget(QLabel("Max Freq (GHz)"))
        spec_adjust_box.addWidget(self.freq_hibound_selector)

        # Group 2: Fractional intensity error for fit
        self.spec_frac_err_selector = QDoubleSpinBox()
        self.spec_frac_err_selector.setDecimals(2)
        self.spec_frac_err_selector.setRange(0.0, 1.0)
        self.spec_frac_err_selector.setSingleStep(0.02)
        self.spec_frac_err_selector.setValue(self.spec_frac_err)
        self.spec_frac_err_selector.valueChanged.connect(self.spec_frac_err_valuechange)
        spec_adjust_box.addWidget(QLabel("Frac. Intensity Error"))
        spec_adjust_box.addWidget(self.spec_frac_err_selector)

        fit_setting_box.addLayout(spec_adjust_box)

        # Group 2: Fit Method
        self.fit_method_box = QHBoxLayout()
        self.fit_method_selector_widget = QComboBox()
        self.fit_method_selector_widget.addItems(["nelder", "basinhopping", "differential_evolution"])
        self.fit_method_selector_widget.currentIndexChanged.connect(self.fit_method_selector)
        self.fit_method_box.addWidget(QLabel("Fit Method"))
        self.fit_method_box.addWidget(self.fit_method_selector_widget)

        # Group 2: Fit Method Dependent Keywords
        self.fit_kws_box = QHBoxLayout()
        self.update_fit_kws_widgets()

        # Group 2: Electron Function
        ele_function_box = QHBoxLayout()
        self.ele_function_selector_widget = QComboBox()
        self.ele_function_selector_widget.addItems(
            ["powerlaw", "double_Powerlaw", "thermal f-f + gyrores", "thermal f-f"])
        self.ele_function_selector_widget.currentIndexChanged.connect(self.ele_function_selector)
        ele_function_box.addWidget(QLabel("Electron Dist. Function"))
        ele_function_box.addWidget(self.ele_function_selector_widget)

        fit_setting_box.addLayout(self.fit_method_box)
        fit_setting_box.addLayout(self.fit_kws_box)
        fit_setting_box.addLayout(ele_function_box)
        fit_layout.addWidget(fit_setting_group)

        ##### Group 3: Fit Parameters Group
        fit_param_group = QGroupBox("Fit Parameters")
        self.fit_param_box = QGridLayout()
        fit_param_group.setLayout(self.fit_param_box)
        self.update_fit_param_widgets()

        fit_layout.addWidget(fit_param_group)
        main_layout.addLayout(fit_layout)

        self._main.setLayout(main_layout)
        self.show()

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

    def eoimg_file_select(self):
        """ Handle Browse button for EOVSA FITS file"""
        # self.fname = QFileDialog.getExistingDirectory(self, 'Select FITS File', './', QFileDialog.ShowDirsOnly)
        ## quick input for debug -------------
        self.eoimg_fname, _file_filter = QFileDialog.getOpenFileName(self, 'Select EOVSA Spectral Image FITS File',
                                                                     './', 'FITS Images (*.fits *.fit *.ft *.fts)')
        # self.eoimg_fname = 'EOVSA_20210507T190135.000000.outim.image.allbd.fits'
        ## quick input for debug -------------
        # self.eoimg_fitsentry.setText(self.eoimg_fname)
        self.eoimg_file_select_return()

    def eoimg_file_select_return(self):
        ''' Called when the FITS filename LineEdit widget gets a carriage-return.
            Trys to read FITS header and return header info (no data read at this time)
        '''

        # print('self.qlookimg_axs is None:',self.qlookimg_axs is None)
        if self.qlookimg_axs is not None:
            for ax in self.qlookimg_axs:
                ax.cla()
        #     self.showqlookimg(showplt=True)
        # self.eoimg_fname = self.eoimg_fitsentry.text()
        self.eoimg_fitsdata = None
        self.has_eovsamap = False

        self.pg_img_canvas.clear()
        if self.has_bkg:
            self.pg_img_canvas.removeItem(self.bkg_roi)
            self.has_bkg = False
        if self.has_rois:
            for roi_group in self.rois:
                for roi in roi_group:
                    self.pg_img_canvas.removeItem(roi)
            self.rois = [[]]
            self.has_rois = False

        try:
            meta, data = ndfits.read(self.eoimg_fname)
            if meta['naxis'] < 3:
                print('Input fits file must have at least 3 dimensions. Abort..')
            elif meta['naxis'] == 3:
                print('Input fits file does not have stokes axis. Assume Stokes I.')
                data = np.expand_dims(data, axis=0)
                self.pol_axis = 0
                self.pol_names = ['I']
                self.has_stokes = False
            elif meta['naxis'] == 4:
                print('Input fits file has stokes axis, located at index {0:d} of the data cube'.
                      format(meta['pol_axis']))
                self.has_stokes = True
                self.pol_axis = meta['pol_axis']
                self.pol_names = meta['pol_names']
            self.meta = meta
            self.data = data
            self.cfreqs = meta['ref_cfreqs'] / 1e9  # convert to GHz
            self.freqdelts = meta['ref_freqdelts'] / 1e9  # convert to GHz
            self.freq_dist = lambda fq: (fq - self.cfreqs[0]) / (self.cfreqs[-1] - self.cfreqs[0])
            self.x0, self.x1 = (np.array([1, meta['nx']]) - meta['header']['CRPIX1']) * meta['header']['CDELT1'] + \
                               meta['header']['CRVAL1']
            self.y0, self.y1 = (np.array([1, meta['ny']]) - meta['header']['CRPIX2']) * meta['header']['CDELT2'] + \
                               meta['header']['CRVAL2']
            self.dx, self.dy = self.meta['header']['CDELT1'], self.meta['header']['CDELT2']
            self.xcen, self.ycen = [(self.x0 + self.x1) / 2.0, (self.y0 + self.y1) / 2.0]
            self.xsiz, self.ysiz = [self.meta['nx'] * self.dx,
                                    self.meta['ny'] * self.dy]
            self.mapx, self.mapy = np.linspace(self.x0, self.x1, meta['nx']), np.linspace(self.y0, self.y1, meta['ny'])
            self.tp_cal_factor = np.ones_like(self.cfreqs)
            self.has_eovsamap = True
            with fits.open(self.eoimg_fname, mode='readonly') as wcs_hdul:
                self.eo_wcs = wcs.WCS(wcs_hdul[0].header)
            # self.infoEdit.setPlainText(repr(rheader))
        except:
            self.statusBar.showMessage('Filename is not a valid FITS file', 2000)
            self.eoimg_fname = '<Select or enter a valid fits filename>'
            # self.eoimg_fitsentry.setText(self.eoimg_fname)
            # self.infoEdit.setPlainText('')

        if not self.qlookimgbutton.isChecked():
            self.qlookimgbutton.setChecked(True)
        else:
            self.showqlookimg()
        # Clean up all existing plots

        # self.plot_qlookmap()
        self.init_pgspecplot()
        self.plot_pg_eovsamap()

    def eodspec_file_select(self):
        """ Handle Browse button for EOVSA FITS file"""
        # self.fname = QFileDialog.getExistingDirectory(self, 'Select FITS File', './', QFileDialog.ShowDirsOnly)
        self.eodspec_fname, _file_filter = QFileDialog.getOpenFileName(self, 'Select EOVSA Dynamic Spectrum FITS File',
                                                                       './', 'FITS Images (*.fits *.fit *.ft *.fts)')
        # self.eodspec_fitsentry.setText(self.eodspec_fname)
        self.eodspec_file_select_return()

    def eodspec_file_select_return(self):
        try:
            hdulist = fits.open(self.eodspec_fname)
            dspec = hdulist[0].data
            header = hdulist[0].header
            observatory = header['telescop']
            pol = header['polariza']
            fghz = np.array(astropy.table.Table(hdulist[1].data)['sfreq'])
            tim_ = astropy.table.Table(hdulist[2].data)
            tmjd = np.array(tim_['mjd']) + np.array(tim_['time']) / 24. / 3600 / 1000
            tim = Time(tmjd, format='mjd')
            self.dspec = {'dspec': dspec, 'time_axis': tim, 'freq_axis': fghz, 'observatory': observatory, 'pol': pol}
            self.has_dspec = True
            self.eodspec_fname = '<Select or enter a valid fits filename>'
            # self.eodspec_fitsentry.setText(self.eodspec_fname)
        except:
            self.statusBar.showMessage('{} is not a valid dynamic spectrum FITS file'.format(self.eodspec_fname))

        if not self.qlookdspecbutton.isChecked():
            self.qlookdspecbutton.setChecked(True)
        else:
            self.showqlookdspec()

    def is_calibrated_tp_state(self):
        if self.is_calibrated_tp_button.isChecked() == True:
            self.statusBar.showMessage('Loaded spectrogram is calibrated total power dynamic spectrum.')
            self.is_calibrated_tp = True
        else:
            self.statusBar.showMessage('Loaded spectrogram is *not* calibrated total power dynamic spectrum.')
            self.is_calibrated_tp = False

    def aiafile_select(self):
        """ Handle Browse button for AIA FITS file """
        self.aiafname, _file_filter = QFileDialog.getOpenFileName(self, 'Select AIA FITS File to open', './',
                                                                  "FITS Images (*.fits *.fit *.ft)")
        # self.aiaimg_fitsentry.setText(self.aiafname)
        self.aiafile_select_return()

    def aiafile_select_return(self):
        ''' Called when the FITS filename LineEdit widget gets a carriage-return.
            Trys to read FITS header and return header info (no data read at this time)
        '''
        # self.aiafname = self.aiaimg_fitsentry.text()
        self.eoimg_fitsdata = None
        try:
            hdu = fits.open(self.aiafname)
            # self.infoEdit.setPlainText(repr(hdu[1].header))
        except:
            self.statusBar.showMessage('Filename is not a valid FITS file', 2000)
            self.aiafname = '<Select or enter a valid fits filename>'
            # self.aiafitsentry.setText(self.fname)
            # self.infoEdit.setPlainText('')
        if not self.qlookimgbutton.isChecked():
            self.qlookimgbutton.setChecked(True)
        else:
            self.showqlookimg()
        # self.plot_qlookmap()

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

        if self.spec_in_tb:
            yticksv = list(range(1, 15))
        else:
            yticksv = list(range(-3, 3))
        yticksv_minor = []
        for v in yticksv:
            vp = 10 ** v
            yticksv_minor += list(np.arange(2 * vp, 10 * vp, 1 * vp))
        self.yticks_minor = []
        for v in yticksv_minor:
            self.yticks_minor.append((np.log10(v), ''))
        self.yticks = []
        if self.spec_in_tb:
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
        else:
            for v in yticksv:
                if v >= 0:
                    self.yticks.append([v, r'{:.0f}'.format(10 ** v)])
                if v == -3:
                    self.yticks.append([v, r'{:.3f}'.format(10 ** v)])
                if v == -2:
                    self.yticks.append([v, r'{:.2f}'.format(10 ** v)])
                if v == -1:
                    self.yticks.append([v, r'{:.1f}'.format(10 ** v)])

        self.update_pgspec()

    def plot_qlookmap(self):
        """Quicklook plot in the upper box using matplotlib.pyplot and sunpy.map"""
        # Plot a quicklook map
        # self.qlook_ax.clear()
        ax0 = self.qlookimg_axs[0]

        if self.has_eovsamap:
            nspw = self.meta['nfreq']
            self.eoimg_date = eoimg_date = Time(self.meta['refmap'].date.mjd +
                                                self.meta['refmap'].exposure_time.value / 2. / 24 / 3600, format='mjd')
            eotimestr = eoimg_date.isot[:-4]
            rsun_obs = sunpy.coordinates.sun.angular_radius(eoimg_date).value
            solar_limb = patches.Circle((0, 0), radius=rsun_obs, fill=False, color='k', lw=1, linestyle=':')
            ax0.add_patch(solar_limb)
            rect = patches.Rectangle((self.x0, self.y0), self.x1 - self.x0, self.y1 - self.y0,
                                     color='k', alpha=0.7, lw=1, fill=False)
            ax0.add_patch(rect)
            icmap = plt.get_cmap('RdYlBu')

            self.qlookimg_canvas.figure.suptitle('EOVSA at {}'.format(eotimestr))
        else:
            self.statusBar.showMessage('EOVSA FITS file does not exist', 2000)
            self.eoimg_fname = '<Select or enter a valid fits filename>'
            # self.eoimg_fitsentry.setText(self.eoimg_fname)
            # self.infoEdit.setPlainText('')
            self.has_eovsamap = False

        if os.path.exists(self.aiafname):
            try:
                aiamap = smap.Map(self.aiafname)
                self.has_aiamap = True
            except:
                self.statusBar.showMessage('Something is wrong with the provided AIA FITS file', 2000)
        else:
            self.statusBar.showMessage('AIA FITS file does not exist', 2000)
            self.aiafname = '<Select or enter a valid fits filename>'
            # self.aiaimg_fitsentry.setText(self.aiafname)
            # self.infoEdit.setPlainText('')
            self.has_aiamap = False
        cts = []
        if self.has_aiamap:
            aiacmap = plt.get_cmap('gray_r')
            aiamap.plot(axes=ax0, cmap=aiacmap)
            ax0.set_title('')
            aia_tit_str = 'AIA {0:.0f} at {1:s}'.format(aiamap.wavelength.value, aiamap.date.isot[:19])
            ax0.text(0.02, 0.98, aia_tit_str, ha='left', va='top', transform=ax0.transAxes, fontsize=10)
        if self.has_eovsamap:
            for s, sp in enumerate(self.cfreqs):
                data = self.data[self.pol_select_idx, s, ...]
                clvls = self.clevels * np.nanmax(data)
                rcmap = [icmap(self.freq_dist(self.cfreqs[s]))] * len(clvls)
                if self.opencontour:
                    cts.append(ax0.contour(self.mapx, self.mapy, data, levels=clvls,
                                           colors=rcmap,
                                           alpha=self.calpha))
                else:
                    cts.append(ax0.contourf(self.mapx, self.mapy, data, levels=clvls,
                                            colors=rcmap,
                                            alpha=self.calpha))

        ax0.set_xlim([-1200, 1200])
        ax0.set_ylim([-1200, 1200])
        ax0.set_xlabel('Solar X [arcsec]')
        ax0.set_ylabel('Solar Y [arcsec]')
        # ax0.set_title('')
        ax0.set_aspect('equal')
        # self.qlookimg_canvas.figure.subplots_adjust(left=0.10, right=0.95,
        #                                        bottom=0.10, top=0.95,
        #                                        hspace=0, wspace=0.35)
        self.qlookimg_canvas.draw()

    def plot_dspec(self, cmap='viridis', vmin=None, vmax=None):
        from matplotlib import dates as mdates
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        if not vmin:
            vmin = np.nanmin(vmin)
        if not vmax:
            vmax = np.nanmax(vmax)
        ax = self.qlookdspec_ax
        tim = self.dspec['time_axis']
        tim_plt = self.dspec['time_axis'].plot_date
        fghz = self.dspec['freq_axis']
        spec = self.dspec['dspec']
        observatory = self.dspec['observatory']
        pol = self.dspec['pol']
        ##todo: now I am only using the first polarization and first baseline
        if spec.ndim < 2:
            print('Dynamic spectrum needs at least 2 dimensions. We have {0:d} here.'.format(spec.ndim))
            return
        elif spec.ndim == 2:
            nfreq, ntim = len(fghz), len(tim_plt)
            npol = 1
            nbl = 1
            spec_plt = spec
        else:
            (npol, nbl, nfreq, ntim) = spec.shape
            print('Dynamic spectrum has more than 2 dimensions {0:d}. '
                  'I am only using the first polarization and first baseline'.format(spec.ndim))
            spec_plt = spec[0, 0]

        im_spec = ax.pcolormesh(tim_plt, fghz, spec_plt, cmap=cmap,
                                vmin=vmin, vmax=vmax, rasterized=True)
        ax.set_title("{0:s} Stokes {1:s} Spectrogram on {2:s}".format(observatory, pol, tim[0].iso[:10]), fontsize=10)
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.set_xlim(tim_plt[0], tim_plt[-1])
        ax.set_ylim(fghz[0], fghz[-1])
        ax.set_xlabel('Time [UT]', fontsize=9)
        ax.set_ylabel('Frequency [GHz]')
        for xlabel in ax.get_xmajorticklabels():
            xlabel.set_rotation(30)
            xlabel.set_horizontalalignment("right")

        # add vertical bar to indicate the time of the EOVSA image
        if hasattr(self, 'eoimg_date'):
            ax.plot([self.eoimg_date.plot_date] * 2, [1, 20], color='w', lw=1)

        divider = make_axes_locatable(ax)
        cax_spec = divider.append_axes('right', size='3.0%', pad=0.05)
        cax_spec.tick_params(direction='out')
        clb_spec = plt.colorbar(im_spec, ax=ax, cax=cax_spec)
        clb_spec.set_label('Flux [sfu]')
        self.qlookdspec_canvas.figure.subplots_adjust(left=0.1, right=0.85,
                                                      bottom=0.20, top=0.92,
                                                      hspace=0, wspace=0)
        self.qlookdspec_canvas.draw()

    def plot_pg_eovsamap(self):
        """This is to plot the eovsamap with the pyqtgraph's ImageView Widget"""
        if self.has_eovsamap:
            # plot the images
            # need to invert the y axis to put the origin to the lower left (tried invertY but labels are screwed up)
            # self.pg_img_plot.setLimits(xMin=self.x0, xMax=self.x1, yMin=self.y0, yMax=self.y1)
            # self.pgdata = self.data[self.pol_select_idx, :, ::-1, :].reshape(
            #     (self.meta['nfreq'], self.meta['ny'], self.meta['nx']))
            self.pgdata = self.data[self.pol_select_idx, :, :, :].reshape(
                (self.meta['nfreq'], self.meta['ny'], self.meta['nx']))
            # self.pg_img_canvas.setImage(self.pgdata, xvals=self.cfreqs)
            self.pg_img_canvas.setImage(self.pgdata, xvals=self.cfreqs, pos=[self.x0, self.y0],
                                        scale=[self.meta['header']['CDELT1'], self.meta['header']['CDELT2']])
            # self.pg_img_canvas.setImage(self.data[self.pol_select_idx], xvals=self.cfreqs)
            self.pg_img_canvas.getView().invertY(False)

            self.pg_img_canvas.setColorMap(self.pgcmap)
            self.pg_freq_current = self.pg_img_canvas.timeLine.getXPos()
            self.pg_freq_idx = np.argmin(np.abs(self.cfreqs - self.pg_freq_current))
            self.pg_img_canvas.getImageItem().hoverEvent = self.pg_map_hover_event

            # define the initial background ROI region
            if not self.has_bkg:
                self.bkg_roi = pg.RectROI([self.x0 + self.dx / 2.0, self.y1 + self.dx / 2.0 - self.ysiz / 5],
                                          [self.xsiz / 5, self.ysiz / 5], pen='w')
                self.bkg_roi.addScaleHandle([1.0, 0.0], [0.0, 1.0])
                self.bkg_roi.addRotateHandle([1.0, 0.5], [0.5, 0.5])
                self.pg_img_canvas.addItem(self.bkg_roi)
                self.has_bkg = True
                self.bkg_roi_label = pg.TextItem("Background", anchor=(0, 0), color='w')
                self.bkg_roi_label.setParentItem(self.bkg_roi)
                self.bkg_rgn_update()
                self.bkg_roi.sigRegionChanged.connect(self.bkg_rgn_update)

    def pg_map_hover_event(self, event):
        """Show the position, pixel, and value under the mouse cursor.
        """
        vbox = self.pg_img_canvas.getView()
        if event.isExit():
            # self.pg_img_hover_label.setText("")
            self.pg_img_mouse_pos_widget.setText("")
            return

        pos = event.pos()
        i, j = pos.y(), pos.x()
        i = int(np.clip(i, 0, self.meta['ny'] - 1))
        j = int(np.clip(j, 0, self.meta['nx'] - 1))
        self.pg_freq_current = self.pg_img_canvas.timeLine.getXPos()
        self.pg_freq_idx = np.argmin(np.abs(self.cfreqs - self.pg_freq_current))
        val = self.pgdata[self.pg_freq_idx, i, j]
        # ppos = self.pg_img_canvas.getImageItem().mapToParent(pos)
        solx, soly = self.mapx[j], self.mapy[::-1][i]
        text_to_display = '[Cursor] x: {0:6.1f}", y: {1:6.1f}", freq: {3:4.1f} GHz, ' \
                          'T<sub>B</sub>={4:6.2f} MK'.format(solx, soly, self.pg_freq_idx,
                                                             self.pg_freq_current, val / 1e6)
        self.pg_img_mouse_pos_widget.setText(text_to_display)
        # print(text_to_display)
        # self.pg_img_hover_label.setText(text_to_display, color='w')

    def update_fitmask(self):
        # update fit mask and the values for currently selected roi group
        if self.has_rois:
            for roi in self.rois[self.roi_group_idx]:
                # define error in spectrum
                if self.spec_in_tb:
                    spec = roi.tb_max
                    spec_bound = self.tb_spec_bound
                    spec_rms = self.bkg_roi.tb_rms
                else:
                    spec = roi.total_flux
                    spec_bound = self.flx_spec_bound
                    spec_rms = gstools.sfu2tb(roi.freqghz * 1e9 * u.Hz, self.bkg_roi.tb_rms * u.K,
                                              area=roi.total_area * u.arcsec ** 2, reverse=True).value
                # add fractional err in quadrature
                spec_err = np.sqrt(spec_rms ** 2. + (self.spec_frac_err * spec) ** 2.)

                spec_tofit0 = ma.masked_less_equal(spec, spec_bound[0])
                freqghz_tofit0 = ma.masked_outside(self.cfreqs, self.fit_freq_bound[0], self.fit_freq_bound[1])
                roi.mask_tofit = mask_tofit = np.logical_or(np.logical_or(freqghz_tofit0.mask, spec_tofit0.mask),
                                                            roi.freq_mask)
                roi.freqghz_tofit = ma.masked_array(self.cfreqs, mask_tofit)
                roi.spec_tofit = ma.masked_array(spec, mask_tofit)
                roi.tb_max_tofit = ma.masked_array(roi.tb_max, mask_tofit)
                roi.total_flux_tofit = ma.masked_array(roi.total_flux, mask_tofit)

                if self.has_bkg:
                    roi.spec_err_tofit = ma.masked_array(spec_err, roi.mask_tofit)
                    roi.tb_rms_tofit = ma.masked_array(self.bkg_roi.tb_rms, roi.mask_tofit)

    def init_pgdistspec_widget(self):
        """Use Pyqtgraph's PlotWidget for the distance-spectral plot"""
        plot = pg.PlotItem(labels={'bottom': ('Frequency [GHz]', ''), 'left': ('Distance [arcsec]', '')})
        self.distSpecCanvasSet[self.new_roi.roi_id] = pg.ImageView(view=plot)
        self.specplotarea.insertWidget(0, self.distSpecCanvasSet[self.new_roi.roi_id])
        for i in range(self.specplotarea.count()):
            self.specplotarea.setStretch(i, 1)

    def update_pgspec(self):
        """Use Pyqtgraph's PlotWidget for the spectral plot"""
        self.speccanvas.clear()
        if self.spec_in_tb:
            spec_bound = self.tb_spec_bound
        else:
            spec_bound = self.flx_spec_bound
        self.speccanvas.setLimits(yMin=np.log10(spec_bound[0]), yMax=np.log10(spec_bound[1]))
        self.plot_pgspec()
        self.pgspec_add_boundbox()

    def plot_pgspec(self):
        self.update_fitmask()
        if self.spec_in_tb:
            spec_bound = self.tb_spec_bound
        else:
            spec_bound = self.flx_spec_bound

        self.spec_dataplots = []
        self.spec_dataplots_tofit = []
        self.spec_rmsplots = []
        current_roi_group = self.rois[self.roi_group_idx]
        for n, roi in enumerate(current_roi_group):
            if roi.type == "sliceROI":
                if self.spec_in_tb:
                    spectrogram = roi.tb_im
                else:
                    spectrogram = roi.tb_im_flux
                im_pix_siz = np.nanmean([self.meta['header']['CDELT1'], self.meta['header']['CDELT2']])
                self.distSpecCanvasSet[roi.roi_id].clear()
                self.distSpecCanvasSet[roi.roi_id].setImage(spectrogram.T, pos=[self.cfreqs[0], 0],
                                                            scale=[np.nanmean(np.diff(self.cfreqs)), im_pix_siz])
                self.distSpecCanvasSet[roi.roi_id].setColorMap(self.pgcmap)
                view = self.distSpecCanvasSet[roi.roi_id].getView()
                view.invertY(False)
                view.setAspectLocked(False)
                nf_sub, nd_sub = spectrogram.shape
                view.setLimits(xMin=self.cfreqs[0], xMax=self.cfreqs[-1], yMin=0, yMax=nd_sub * im_pix_siz)
            else:
                if n == self.current_roi_idx or (n == len(current_roi_group) + self.current_roi_idx):
                    symbolfill = (n, 9)
                else:
                    symbolfill = None
                if self.spec_in_tb:
                    spec = roi.tb_max
                    spec_bound = self.tb_spec_bound
                else:
                    spec = roi.total_flux
                    spec_bound = self.flx_spec_bound
                if ma.is_masked(roi.freqghz):
                    log_freqghz = np.log10(roi.freqghz.compressed())
                    log_spec = np.log10(spec.compressed())
                else:
                    log_freqghz = np.log10(roi.freqghz)
                    log_spec = np.log10(spec)
                spec_dataplot = self.speccanvas.plot(x=log_freqghz, y=log_spec, pen=None,
                                                     symbol='o', symbolPen=(n, 9), symbolBrush=None)
                spec_dataplot_tofit = self.speccanvas.plot(x=np.log10(roi.freqghz_tofit.compressed()),
                                                           y=np.log10(roi.spec_tofit.compressed()),
                                                           pen=None,
                                                           symbol='o', symbolPen=(n, 9), symbolBrush=symbolfill)

                self.speccanvas.addItem(spec_dataplot)
                self.speccanvas.addItem(spec_dataplot_tofit)
                self.spec_dataplots.append(spec_dataplot)
                self.spec_dataplots_tofit.append(spec_dataplot_tofit)
                # print('ROI', n)
                # print(roi.freqghz_tofit.compressed())
                # print(roi.spec_tofit.compressed())

                # Add errorbar if rms is defined
                if self.has_bkg:
                    # define error in spectrum
                    if self.spec_in_tb:
                        spec_rms = self.bkg_roi.tb_rms
                    else:
                        spec_rms = gstools.sfu2tb(roi.freqghz * 1e9 * u.Hz, self.bkg_roi.tb_rms * u.K,
                                                  area=roi.total_area * u.arcsec ** 2, reverse=True).value
                    # add fractional err in quadrature
                    spec_err = np.sqrt(spec_rms ** 2. + (self.spec_frac_err * spec) ** 2.)
                    spec_rmsplot = self.speccanvas.plot(x=np.log10(self.cfreqs), y=np.log10(spec_rms), pen='k',
                                                        symbol='d', symbolPen='k', symbolBrush=None)
                    self.speccanvas.addItem(spec_rmsplot)
                    self.spec_rmsplots.append(spec_rmsplot)
                else:
                    spec_err = (self.spec_frac_err * spec) ** 2.
                err_bounds_min = np.maximum(spec - spec_err, np.ones_like(spec) * spec_bound[0])
                err_bounds_max = spec + spec_err
                errplot = pg.ErrorBarItem(x=np.log10(roi.freqghz), y=np.log10(spec),
                                          top=np.log10(err_bounds_max) - np.log10(spec),
                                          bottom=np.log10(spec) - np.log10(err_bounds_min), beam=0.025, pen=(n, 9))

                self.speccanvas.addItem(errplot)

            self.fbar = self.speccanvas.plot(x=np.log10([self.pg_img_canvas.timeLine.getXPos()] * 2),
                                             y=[np.log10(spec_bound[0]), np.log10(spec_bound[1])], pen='k')
            self.speccanvas.addItem(self.fbar)
            self.speccanvas.setLimits(yMin=np.log10(spec_bound[0]), yMax=np.log10(spec_bound[1]))
            xax = self.speccanvas.getAxis('bottom')
            yax = self.speccanvas.getAxis('left')
            xax.setLabel("Frequency [GHz]")
            if self.spec_in_tb:
                yax.setLabel("Brightness Temperature [MK]")
            else:
                yax.setLabel("Flux Density [sfu]")
            xax.setTicks([self.xticks, self.xticks_minor])
            yax.setTicks([self.yticks, self.yticks_minor])

    def pgspec_add_boundbox(self):
        # add frequency bound
        self.rgn_freq_bound = pg.LinearRegionItem(brush=(20, 50, 200, 20))
        self.rgn_freq_bound.setZValue(10)
        self.rgn_freq_bound.setRegion([np.log10(self.fit_freq_bound[0]), np.log10(self.fit_freq_bound[1])])
        self.speccanvas.addItem(self.rgn_freq_bound)
        self.rgn_freq_bound.sigRegionChangeFinished.connect(self.update_freq_bound_rgn)

    def update_fbar(self):
        if self.fbar is not None:
            try:
                self.speccanvas.removeItem(self.fbar)
                self.fbar = self.speccanvas.plot(x=np.log10([self.pg_img_canvas.timeLine.getXPos()] * 2), y=[1, 15],
                                                 pen='k')
                self.speccanvas.addItem(self.fbar)
            except:
                pass

    #    def bkg_rgn_select(self):
    #        """Select a region to calculate rms on spectra"""
    #        if self.bkg_selection_button.isChecked():
    #            self.bkg_selection_button.setStyleSheet("background-color : lightblue")
    #            #self.bkg_roi = pg.RectROI([0, 0], [40, 40], pen='w')
    #            #self.pg_img_canvas.addItem(self.bkg_roi)
    #            #self.bkg_rgn_update()
    #            self.bkg_roi.sigRegionChanged.connect(self.bkg_rgn_update)
    #        else:
    #            # if unChecked remove the ROI box from the plot
    #            self.bkg_selection_button.setStyleSheet("background-color : lightgrey")
    #            self.pg_img_canvas.removeItem(self.bkg_roi)

    def bkg_rgn_update(self):
        """Select a region to calculate rms on spectra"""
        bkg_subim = self.bkg_roi.getArrayRegion(self.pgdata, self.pg_img_canvas.getImageItem(), axes=(1, 2))
        nf_bkg, ny_bkg, nx_bkg = bkg_subim.shape
        if self.bkg_roi.pos()[1] < self.y0 + self.ysiz * 0.1:
            self.bkg_roi_label.setAnchor((0, 1))
        else:
            self.bkg_roi_label.setAnchor((0, 0))
        # if self.bkg_roi.pos()[1] < self.y0 + self.ysiz * 0.1:
        self.bkg_roi.freqghz = self.cfreqs
        self.bkg_roi.tb_mean = np.nanmean(bkg_subim, axis=(1, 2))
        self.bkg_roi.tb_rms = np.std(bkg_subim, axis=(1, 2))
        self.bkg_roi.tb_max = np.nanmax(bkg_subim, axis=(1, 2))
        self.bkg_roi.total_pix = ny_bkg * nx_bkg
        # Total area of the ROI in arcsec^2
        self.bkg_roi.total_area = self.bkg_roi.total_pix * self.meta['header']['CDELT1'] * self.meta['header']['CDELT2']
        # Total flux of the ROI in sfu
        self.bkg_roi.total_flux = gstools.sfu2tb(self.bkg_roi.freqghz * 1e9 * u.Hz, self.bkg_roi.tb_mean * u.K,
                                                 area=self.bkg_roi.total_area * u.arcsec ** 2, reverse=True).value
        self.pg_img_bkg_roi_info_widget.setText('[Background] freq: {0:.2f} GHz, '
                                                'T<sub>B</sub><sup>rms</sup>: {1:.2f} MK, '
                                                'T<sub>B</sub><sup>mean</sup>: {2:.2f} MK, Flux: {3:.2f} sfu'.
                                                format(self.pg_freq_current,
                                                       self.bkg_roi.tb_rms[self.pg_freq_idx] / 1e6,
                                                       self.bkg_roi.tb_mean[self.pg_freq_idx] / 1e6,
                                                       self.bkg_roi.total_flux[self.pg_freq_idx]))
        # self.bkg_roi_label = pg.TextItem("Background", anchor=(0, 0), color='w')
        # self.bkg_roi_label.setParentItem(self.bkg_roi)
        self.update_pgspec()

    def findlastsliceROI(self):
        roi_type_list = [roi.type for roi in self.rois[self.roi_group_idx]][::-1]
        if "sliceROI" in roi_type_list:
            return self.rois[self.roi_group_idx][::-1][roi_type_list.index("sliceROI")]
        else:
            return None

    def add_new_roi(self):
        """Add a ROI region to the selection"""
        self.roi_type = self.toolBarButtonGroup.checkedButton().text()
        ischildROI = False
        if self.add2slice.isChecked():
            if self.roi_type in ["LineSegmentROI", "PolyLineROI"]:
                roi_xcen = self.xcen
                roi_ycen = self.ycen
                roi_size = np.array([self.xsiz, self.ysiz]) / 20
            elif self.roi_type in ["PolygonROI"]:
                self.statusBar.showMessage(
                    '"add ROI to Slice" is not yet supported for PolygonROI.')
                return
            else:
                lastsliceROI = self.findlastsliceROI()
                if lastsliceROI is None:
                    self.statusBar.showMessage(
                        'There is no sliceROI existed in the current ROI group, use the sliceROI tools to add one first.')
                    return
                else:
                    ischildROI = True
                    # print(lastsliceROI)
                    childROI = lastsliceROI.childROI
                    nchildROI = len(childROI.keys())

                    pts = [h.pos() for h in lastsliceROI.endpoints]
                    # print(pts)
                    # print(childROI.keys())
                    dist = [0]
                    pts_orig = lastsliceROI.pos()
                    pts_x = [pts[0].x()]
                    pts_y = [pts[0].y()]
                    for idx, pt in enumerate(pts[:-1]):
                        d = pg.Point(pts[idx + 1] - pts[idx])
                        dist.append(d.length())
                        pts_x.append(pts[idx + 1].x())
                        pts_y.append(pts[idx + 1].y())
                    dist = np.cumsum(dist)
                    dnew = np.linspace(dist[0], dist[-1], nchildROI + 1)
                    pts_x = np.array(pts_x) + pts_orig.x()
                    pts_y = np.array(pts_y) + pts_orig.y()
                    pts_xnew = np.interp(dnew, dist, pts_x)
                    pts_ynew = np.interp(dnew, dist, pts_y)
                    roi_xcen = pts_xnew[-1]
                    roi_ycen = pts_ynew[-1]
                    if nchildROI > 0:
                        roi_size = childROI[0].size()
                    else:
                        roi_size = np.array([self.xsiz, self.ysiz]) / 20

                    for idx, (k, v) in enumerate(childROI.items()):
                        v.sigRegionChanged.disconnect(self.calc_roi_spec)
                        v.setPos((pts_xnew[idx] - roi_size[0] / 2.0, pts_ynew[idx] - roi_size[0] / 2.0))
                        v.setSize(roi_size)
                    for idx, (k, v) in enumerate(childROI.items()):
                        v.sigRegionChanged.connect(self.calc_roi_spec)
                    # if nchildROI == 0:
                    #     pass
                    # else:
                    #     pass

                    # if there is no slice roi. return. infobar output: add a slice first.
                    # else:
                    #   if there is childROI
                    #       read the last childROI's size and distribute this childROI
                    #   else
                    #       use default size and place this childROI in the center of the slice.
        else:
            roi_xcen = self.xcen
            roi_ycen = self.ycen
            roi_size = np.array([self.xsiz, self.ysiz]) / 20

        # if ischildROI:
        #     colorid = len(lastsliceROI.childROI.keys())
        # else:
        colorid = len(self.rois[self.roi_group_idx])
        if self.roi_type == "RectROI":
            self.new_roi = pg.RectROI([roi_xcen - roi_size[0] / 2.0,
                                       roi_ycen - roi_size[1] / 2.0], roi_size,
                                      pen=(colorid, 9),
                                      removable=True, centered=True)
            # self.new_roi.setSize(self.new_roi.size(), center=(0.5, 0.5))
            # self.new_roi.setAngle(self.new_roi.angle(), center=(0.5, 0.5))
            self.new_roi.addRotateHandle([1.0, 0.5], [0.5, 0.5])
        elif self.roi_type == "EllipseROI":
            self.new_roi = pg.EllipseROI([roi_xcen - roi_size[0] / 2.0,
                                          roi_ycen - roi_size[1] / 2.0], roi_size,
                                         pen=(colorid, 9), removable=True)
            # self.new_roi.setSize(self.new_roi.size(), center=(0.5, 0.5))
            # self.new_roi.setAngle(self.new_roi.angle(), center=(0.5, 0.5))
        elif self.roi_type == "PolygonROI":
            self.new_roi = pg.PolyLineROI([[roi_xcen - roi_size[0] * 1,
                                            roi_ycen - roi_size[1] * 0.5],
                                           [roi_xcen + roi_size[0] * 0,
                                            roi_ycen - roi_size[1] * 0.5],
                                           [roi_xcen + roi_size[0] * 1,
                                            roi_ycen + roi_size[1] * 0.5],
                                           [roi_xcen - roi_size[0] * 0,
                                            roi_ycen + roi_size[1] * 0.5]],
                                          pen=(colorid, 9), removable=True, closed=True)
            self.new_roi.addRotateHandle([1.0, 0.5], [0.5, 0.5])
        elif self.roi_type == "LineSegmentROI":
            self.new_roi = pg.LineSegmentROI([[roi_xcen - roi_size[0] / 2.0,
                                               roi_ycen - roi_size[1] / 2.0],
                                              [roi_xcen + roi_size[0] / 2.0,
                                               roi_ycen + roi_size[1] / 2.0]],
                                             pen=(colorid, 9),
                                             removable=True)
        elif self.roi_type == "PolyLineROI":
            self.new_roi = PolyLineROIX([[roi_xcen - roi_size[0] / 2.0,
                                          roi_ycen - roi_size[1] / 2.0],
                                         [roi_xcen + roi_size[0] / 2.0,
                                          roi_ycen + roi_size[1] / 2.0]],
                                        pen=(colorid, 9),
                                        removable=True, closed=False)
        else:
            self.new_roi = pg.RectROI([roi_xcen - roi_size[0] / 2.0,
                                       roi_ycen - roi_size[1] / 2.0], roi_size,
                                      pen=(colorid, 9),
                                      removable=True)
            self.new_roi.addRotateHandle([1.0, 0.5], [0.5, 0.5])

        self.new_roi.ischildROI = ischildROI
        self.pg_img_canvas.addItem(self.new_roi)
        self.new_roi.freq_mask = np.ones_like(self.cfreqs) * False
        self.new_roi.sigRegionChanged.connect(self.calc_roi_spec)
        self.new_roi.sigRemoveRequested.connect(self.remove_ROI)

        # choose which group to add
        self.add_to_roigroup_selection()
        if self.roi_type in ["LineSegmentROI", "PolyLineROI"]:
            self.new_roi.type = "sliceROI"
            self.new_roi.childROI = {}
        else:
            self.new_roi.type = "ROI"
            if ischildROI:
                self.new_roi.colorid = colorid
                self.new_roi.parentROI = lastsliceROI
                self.new_roi.childROIkey = nchildROI
                self.new_roi.parentROI.childROI[self.new_roi.childROIkey] = self.new_roi

        self.rois[self.roi_group_idx].append(self.new_roi)
        self.nroi_current_group = len(self.rois[self.roi_group_idx])
        self.roi_selection_widget.clear()
        self.roi_selection_widget.addItems([str(i) for i in range(self.nroi_current_group)])
        self.current_roi_idx = self.nroi_current_group - 1
        self.new_roi.roi_id = self.current_roi_idx
        self.has_rois = True
        # self.roi_freq_lowbound_selector.setValue(self.data_freq_bound[0])
        # self.roi_freq_hibound_selector.setValue(self.data_freq_bound[1])
        if self.new_roi.type == "sliceROI":
            self.init_pgdistspec_widget()
        self.calc_roi_spec(None)

        # self.roi_slider_rangechange()

    def remove_ROI(self, evt):
        if evt.ischildROI:
            evt.parentROI.childROI.pop(evt.childROIkey)
            chilROI = {}
            for idx, (k, v) in enumerate(evt.parentROI.childROI.items()):
                chilROI[idx] = v
                v.childROIkey = idx
            del evt.parentROI.childROI
            evt.parentROI.childROI = chilROI
        if evt.type in ["sliceROI"]:
            self.distSpecCanvasSet[evt.roi_id].deleteLater()
            for idx, (k, v) in enumerate(evt.childROI.items()):
                self.rois[self.roi_group_idx].remove(v)
                self.pg_img_canvas.removeItem(v)
        self.rois[self.roi_group_idx].remove(evt)

        self.nroi_current_group = len(self.rois[self.roi_group_idx])
        self.roi_selection_widget.clear()
        self.roi_selection_widget.addItems([str(i) for i in range(self.nroi_current_group)])
        self.current_roi_idx = self.nroi_current_group - 1
        for idx, roi in enumerate(self.rois[self.roi_group_idx]):
            if roi.type != "sliceROI":
                roi.setPen((idx, 9))
        self.pg_img_canvas.removeItem(evt)

        if np.count_nonzero(self.rois) < 1:
            self.has_rois = False
            self.speccanvas.clear()
        else:
            self.calc_roi_spec(None)
        del evt

    # self.view.scene().removeItem(evt)

    # ROI = None

    def add_to_roigroup_selection(self):
        items = self.add_to_roigroup_widget.selectedItems()
        if len(items) > 0:
            self.add_to_roigroup_button.setText(items[0].text())
            self.roi_group_idx = int(items[0].text())
        else:
            self.roi_group_idx = 0

    def roigroup_selection(self):
        items = self.roigroup_selection_widget.selectedItems()
        if len(items) > 0:
            self.roigroup_selection_button.setText(items[0].text())
            self.roi_group_idx = int(items[0].text())
        else:
            self.roi_group_idx = 0

    def roi_selection_action(self):
        items = self.roi_selection_widget.selectedItems()
        if len(items) > 0:
            self.roi_selection_button.setText(items[0].text())
            self.current_roi_idx = int(items[0].text())
        else:
            self.current_roi_idx = len(self.rois[self.roi_group_idx]) - 1
        self.update_pgspec()
        # print(self.current_roi_idx)

    def group_roi_op_selector(self):
        cur_action = self.sender()
        cur_op = cur_action.text()
        if cur_op == 'Save Group':
            roi_utils.save_roi_group(self)

    def exec_customized_rois_window(self):
        try:
            self.customized_rois_Form = QDialog()
            ui = roi_utils.roi_dialog(img_size=[self.meta['nx'], self.meta['ny']], cfreqs = self.cfreqs)
            ui.setupUi(self.customized_rois_Form)
            self.customized_rois_Form.show()
            cur_result = self.customized_rois_Form.exec()
            crtf_list = ui.getResult(self.customized_rois_Form)
            print('cur_result: ',crtf_list)
            return (crtf_list,cur_result)
        except:
            msg_box = QMessageBox(QMessageBox.Warning, 'No EOVSA Image Loaded', 'Load EOVSA Image first!')
            msg_box.exec_()

    def add_manually_defined_rois(self):
        dialog_output = self.exec_customized_rois_window()
        if not dialog_output[1] or len(dialog_output[0]) == 0:
            print('No ROI is added!')
        else:
            roi_utils.add_md_rois(self, inp_str_list =dialog_output[0])

    def calc_roi_spec(self, evt):
        # print('=================Update ROI SPEC===============')
        # roi = self.rois[self.roi_group_idx][self.current_roi_idx]
        # roi = self.new_roi
        # print('calc_roi_spec', evt)

        ## if the signal emitter is a child ROI of a sliceROI
        if evt is not None:
            rois2update = [evt]
            if evt.ischildROI:
                ROIsize = evt.size()
                ROIangle = evt.angle()
                ROIpos = evt.pos()
                # print(evt.parentROI.childROI)
                for idx, (k, v) in enumerate(evt.parentROI.childROI.items()):
                    update = 0
                    v.sigRegionChanged.disconnect(self.calc_roi_spec)
                    if ROIsize != v.size():
                        update = 1
                        v.setSize(ROIsize, center=(0.5, 0.5))
                    if ROIangle != v.angle():
                        update = 1
                        v.setAngle(ROIangle, center=(0.5, 0.5))
                    if ROIpos != v.pos():
                        update = 0
                    # print(update,v)
                    if update:
                        rois2update.append(v)
                for idx, (k, v) in enumerate(evt.parentROI.childROI.items()):
                    v.sigRegionChanged.connect(self.calc_roi_spec)
        else:
            rois2update = self.rois[self.roi_group_idx]

        for roi in rois2update:
            ## todo: investigate why using axes = (2, 1) returns entirely different (wrong!) results

            subim = roi.getArrayRegion(self.pgdata,
                                       self.pg_img_canvas.getImageItem(), axes=(1, 2))

            roi.freqghz = self.cfreqs
            roi.freq_bound = copy.copy(self.roi_freq_bound)
            tb2flx_1d = gstools.sfu2tb(roi.freqghz * 1e9 * u.Hz, np.ones_like(roi.freqghz) * u.K,
                                       area=self.meta['header']['CDELT1'] * self.meta['header'][
                                           'CDELT2'] * u.arcsec ** 2, reverse=True).value

            if roi.type == "sliceROI":
                nf_sub, nd_sub = subim.shape
                roi.total_pix = nd_sub

                roi.tb_max = np.nanmax(subim, axis=1)
                roi.tb_mean = np.nanmean(subim, axis=1)

                roi.tb_im = subim

                tb2flx_2d = np.tile(tb2flx_1d, nd_sub).reshape(nd_sub, nf_sub).transpose()
                roi.tb_im_flux = roi.tb_im * tb2flx_2d
            else:
                nf_sub, ny_sub, nx_sub = subim.shape
                roi.total_pix = ny_sub * nx_sub

                roi.tb_max = np.nanmax(subim, axis=(1, 2))
                roi.tb_mean = np.nanmean(subim, axis=(1, 2))

            # Total area of the ROI in arcsec^2
            roi.total_area = roi.total_pix * self.meta['header']['CDELT1'] * self.meta['header']['CDELT2']
            # Total flux of the ROI in sfu
            roi.total_flux = gstools.sfu2tb(roi.freqghz * 1e9 * u.Hz, roi.tb_mean * u.K,
                                            area=roi.total_area * u.arcsec ** 2, reverse=True).value
            # print(tb2flx_1d * roi.tb_mean * roi.total_pix /roi.total_flux)
            if np.sum(roi.freq_mask) > 0:
                roi.freqghz = ma.masked_array(roi.freqghz, roi.freq_mask)
                roi.tb_max = ma.masked_array(roi.tb_max, roi.freq_mask)
                roi.tb_mean = ma.masked_array(roi.tb_mean, roi.freq_mask)
                roi.total_flux = ma.masked_array(roi.total_flux, roi.freq_mask)

        if 'roi' in vars():
            self.roi_info.setText('[Current ROI] x: {0:6.1f}", y: {1:6.1f}", freq: {2:4.1f} GHz, '
                                  'T<sub>B</sub><sup>max</sup>: {3:5.1f} MK, '
                                  'T<sub>B</sub><sup>mean</sup>: {4:5.1f} MK, Flux: {5:5.1f} sfu'.
                                  format(roi.pos()[0] + roi.size()[0] / 2, roi.pos()[1] + roi.size()[1] / 2,
                                         self.pg_freq_current,
                                         roi.tb_max[self.pg_freq_idx] / 1e6,
                                         roi.tb_mean[self.pg_freq_idx] / 1e6,
                                         roi.total_flux[self.pg_freq_idx]))
        # self.update_fitmask()
        self.update_pgspec()

    def combine_roi_group_flux(self):
        """
        Function to combine flux density of a selected ROI group
        @return:
            self.img_tot_fghz: a masked array of frequencies of the selected ROI group
            self.img_tot_flux: a masked array of total flux density of the selected ROI group
            self.img_tot_flux_plot: a handler for pyqtgraph plot of the total flux
        """
        self.spec_in_tb = False
        self.update_pgspec()
        rois_ = self.rois[self.roi_group_idx]
        nroi = len(rois_)
        if nroi > 1:
            for n in range(nroi):
                roi = rois_[n]
                if n == 0:
                    img_tot_flux_arr = rois_[n].total_flux
                    freq_masks = rois_[n].freq_mask
                else:
                    img_tot_flux_arr = ma.vstack((img_tot_flux_arr, rois_[n].total_flux))
                    freq_masks = ma.vstack((freq_masks, rois_[n].freq_mask))
            self.img_tot_flux = ma.max(img_tot_flux_arr, axis=0)
            img_tot_freq_mask = np.logical_and.reduce(freq_masks)
            self.img_tot_fghz = ma.masked_array(self.cfreqs, img_tot_freq_mask)
        else:
            self.img_tot_flux = rois_[0].total_flux
            self.img_tot_fghz = rois_[0].freqghz
        if self.spec_in_tb == False:
            if hasattr(self, 'img_tot_flux_plot'):
                self.speccanvas.removeItem(self.img_tot_flux_plot)
            else:
                self.img_tot_flux_plot = self.speccanvas.plot(x=np.log10(self.img_tot_fghz),
                                                              y=np.log10(self.img_tot_flux),
                                                              pen=dict(color='k', width=4),
                                                              symbol=None, symbolBrush=None)
                self.speccanvas.addItem(self.img_tot_flux_plot)

    def calc_tpcal_factor(self):
        # find out and plot total power spectrum at the given time
        self.spec_in_tb = False
        if hasattr(self, 'eoimg_date') and hasattr(self, 'dspec') and self.is_calibrated_tp:
            t_idx = np.argmin(np.abs(self.dspec['time_axis'] - self.eoimg_date))
            self.tp_spec = self.dspec['dspec'][:, t_idx]
            self.tp_fghz = self.dspec['freq_axis']
            if hasattr(self, 'tp_flux_plot'):
                self.speccanvas.removeItem(self.tp_flux_plot)
            else:
                self.tp_flux_plot = self.speccanvas.plot(x=np.log10(self.tp_fghz), y=np.log10(self.tp_spec),
                                                         pen=dict(color='k', width=8),
                                                         symbol=None, symbolBrush=None)
                self.speccanvas.addItem(self.tp_flux_plot)
            if hasattr(self, 'img_tot_flux') and hasattr(self, 'img_tot_fghz'):
                self.tp_cal_factor = np.ones_like(self.img_tot_flux)
                for n, fghz in enumerate(self.img_tot_fghz):
                    fidx_tp = np.argmin(np.abs(self.tp_fghz - self.img_tot_fghz[n]))
                    self.tp_cal_factor[n] = self.img_tot_flux[n] / self.tp_spec[fidx_tp]
                print('Total Power Calibration Factor updated')
                print(self.tp_cal_factor)
        else:
            print('Either image time or calibrated total power dynamic spectrum does not exist.')

    def apply_tpcal_factor(self):
        if self.apply_tpcal_factor_button.isChecked() == True:
            self.statusBar.showMessage('Apply total power correction factor to data.')
            self.data[self.pol_select_idx] /= self.tp_cal_factor[:, None, None]
            self.tpcal_factor_applied = True
            self.calc_roi_spec(None)
        else:
            self.statusBar.showMessage('Unapply total power correction factor to data.')
            self.data[self.pol_select_idx] *= self.tp_cal_factor[:, None, None]
            self.tpcal_factor_applied = False
            self.calc_roi_spec(None)

    # def roi_slider_rangechange(self):
    #    self.roi_select_slider.setMaximum(self.nroi - 1)

    # def roi_slider_valuechange(self):
    #    self.current_roi_idx = self.roi_select_slider.value()
    #    self.roi_info.setPlainText('Selected ROI {}'.format(self.current_roi_idx))
    #    self.update_pgspec()

    def roi_freq_lowbound_valuechange(self):
        self.roi_freq_bound[0] = self.roi_freq_lowbound_selector.value()
        self.statusBar.showMessage('Selected Lower Frequency Bound for ROI is {0:.1f} GHz'.
                                   format(self.roi_freq_bound[0]))
        self.update_freq_mask(self.new_roi)
        self.update_pgspec()

    def roi_freq_hibound_valuechange(self):
        self.roi_freq_bound[1] = self.roi_freq_hibound_selector.value()
        self.statusBar.showMessage('Selected Higher Frequency Bound for ROI is {0:.1f} GHz'.
                                   format(self.roi_freq_bound[1]))
        self.update_freq_mask(self.new_roi)
        self.update_pgspec()

    def update_freq_mask(self, roi):
        # unmask everything first
        if ma.is_masked(roi.freqghz):
            roi.freqghz.mask = ma.nomask
            roi.tb_max.mask = ma.nomask
            roi.tb_mean.mask = ma.nomask
            roi.total_flux.mask = ma.nomask

        # update the masks
        roi.freq_bound = self.roi_freq_bound
        roi.freqghz = ma.masked_outside(self.cfreqs, self.roi_freq_bound[0], self.roi_freq_bound[1])
        roi.freq_mask = roi.freqghz.mask
        roi.tb_max = ma.masked_array(roi.tb_max, roi.freq_mask)
        roi.tb_mean = ma.masked_array(roi.tb_mean, roi.freq_mask)
        roi.total_flux = ma.masked_array(roi.total_flux, roi.freq_mask)

    def roi_grid_size_valuechange(self):
        self.roi_grid_size = self.roi_grid_size_selector.value()

    def freq_lowbound_valuechange(self):
        self.fit_freq_bound[0] = self.freq_lowbound_selector.value()
        self.update_pgspec()

    def freq_hibound_valuechange(self):
        self.fit_freq_bound[1] = self.freq_hibound_selector.value()
        self.update_pgspec()

    def spec_frac_err_valuechange(self):
        self.spec_frac_err = self.spec_frac_err_selector.value()
        self.update_pgspec()

    def update_freq_bound_rgn(self):
        self.rgn_freq_bound.setZValue(10)
        min_logfreq, max_logfreq = self.rgn_freq_bound.getRegion()
        self.fit_freq_bound = [10. ** min_logfreq, 10. ** max_logfreq]
        self.update_fitmask()
        self.freq_lowbound_selector.setValue(self.fit_freq_bound[0])
        self.freq_hibound_selector.setValue(self.fit_freq_bound[1])
        for spec_dataplot, spec_dataplot_tofit, spec_rmsplot in \
                zip(self.spec_dataplots, self.spec_dataplots_tofit, self.spec_rmsplots):
            self.speccanvas.removeItem(spec_dataplot)
            self.speccanvas.removeItem(spec_dataplot_tofit)
            self.speccanvas.removeItem(spec_rmsplot)
        self.plot_pgspec()

    def fit_method_selector(self):
        print("Selected Fit Method is: {}".format(self.fit_method_selector_widget.currentText()))
        self.fit_method = self.fit_method_selector_widget.currentText()
        self.init_fit_kws()
        self.update_fit_kws_widgets()

    def ele_function_selector(self):
        print("Selected Electron Distribution Function is: {}".format(self.ele_function_selector_widget.currentText()))
        self.ele_dist = self.ele_function_selector_widget.currentText()
        self.init_params()

    def init_params(self):
        if self.ele_dist == 'powerlaw':
            self.fit_params = lmfit.Parameters()
            self.fit_params.add_many(('Bx100G', 2., True, 0.1, 100., None, None),
                                     ('log_nnth', 5., True, 3., 11, None, None),
                                     ('delta', 4., True, 1., 30., None, None),
                                     ('Emin_keV', 10., False, 1., 100., None, None),
                                     ('Emax_MeV', 10., False, 0.05, 100., None, None),
                                     ('theta', 45., True, 0.01, 89.9, None, None),
                                     ('log_nth', 10, True, 4., 13., None, None),
                                     ('T_MK', 1., False, 0.1, 100, None, None),
                                     ('depth_asec', 5., False, 1., 100., None, None))
            self.fit_params_nvarys = 0
            for key, par in self.fit_params.items():
                if par.vary:
                    self.fit_params_nvarys += 1

            self.fit_function = gstools.GSCostFunctions.SinglePowerLawMinimizerOneSrc
            self.update_fit_param_widgets()

        if self.ele_dist == 'thermal f-f':
            ## todo: thermal free-free cost function to be added
            self.fit_params = lmfit.Parameters()
            self.fit_params.add_many(('theta', 45., True, 0.01, 89.9, None, None),
                                     ('log_nth', 10, True, 4., 13., None, None),
                                     ('T_MK', 1., False, 0.1, 100, None, None),
                                     ('depth_asec', 5., False, 1., 100., None, None),
                                     ('area_asec2', 25., False, 1., 10000., None, None))
            self.fit_params_nvarys = 0
            for key, par in self.fit_params.items():
                if par.vary:
                    self.fit_params_nvarys += 1

            self.update_fit_param_widgets()


        if self.ele_dist == 'thermal f-f + gyrores':
            self.fit_params = lmfit.Parameters()
            self.fit_params.add_many(('Bx100G', 2., True, 0.1, 100., None, None),
                                     ('log_nnth', 5., False, 3., 11, None, None),
                                     ('delta', 4., False, 1., 30., None, None),
                                     ('Emin_keV', 10., False, 1., 100., None, None),
                                     ('Emax_MeV', 10., False, 0.05, 100., None, None),
                                     ('theta', 45., True, 0.01, 89.9, None, None),
                                     ('log_nth', 10, True, 4., 13., None, None),
                                     ('T_MK', 1., True, 0.1, 100, None, None),
                                     ('depth_asec', 5., False, 1., 100., None, None))
            self.fit_params_nvarys = 0
            for key, par in self.fit_params.items():
                if par.vary:
                    self.fit_params_nvarys += 1

            self.fit_function = gstools.GSCostFunctions.Ff_Gyroresonance_MinimizerOneSrc
            self.update_fit_param_widgets()

    def init_fit_kws(self):
        # first refresh the widgets
        if self.fit_method == 'nelder':
            self.fit_kws = {'maxiter': 2000, 'xatol': 0.01, 'fatol': 0.01}
        if self.fit_method == 'basinhopping':
            self.fit_kws = {'niter': 50, 'T': 90., 'stepsize': 0.8,
                            'interval': 25}

    def update_fit_kws_widgets(self):
        # first delete every widget for the fit keywords
        if self.fit_kws_box.count() > 0:
            for n in reversed(range(self.fit_kws_box.count())):
                self.fit_kws_box.itemAt(n).widget().deleteLater()

        self.fit_kws_key_widgets = []
        self.fit_kws_value_widgets = []
        for n, key in enumerate(self.fit_kws):
            fit_kws_key_widget = QLabel(key)
            self.fit_kws_box.addWidget(fit_kws_key_widget)
            if type(self.fit_kws[key]) == int:
                fit_kws_value_widget = QSpinBox()
                fit_kws_value_widget.setRange(0, 10000)
                fit_kws_value_widget.setValue(self.fit_kws[key])
                fit_kws_value_widget.valueChanged.connect(self.update_fit_kws)
            if type(self.fit_kws[key]) == float:
                fit_kws_value_widget = QDoubleSpinBox()
                fit_kws_value_widget.setRange(0, 10000)
                fit_kws_value_widget.setValue(self.fit_kws[key])
                fit_kws_value_widget.setDecimals(2)
                fit_kws_value_widget.valueChanged.connect(self.update_fit_kws)
            if type(self.fit_kws[key]) == str:
                fit_kws_value_widget = QTextEdit()
                fit_kws_value_widget.setText(self.fit_kws[key])
                fit_kws_value_widget.valueChanged.connect(self.update_fit_kws)
            self.fit_kws_key_widgets.append(fit_kws_key_widget)
            self.fit_kws_value_widgets.append(fit_kws_value_widget)
            self.fit_kws_box.addWidget(fit_kws_value_widget)

    def update_fit_param_widgets(self):
        # first delete every widget for the fit parameters
        if self.fit_param_box.count() > 0:
            for n in reversed(range(self.fit_param_box.count())):
                self.fit_param_box.itemAt(n).widget().deleteLater()

        self.param_init_value_widgets = []
        self.param_vary_widgets = []
        self.param_min_widgets = []
        self.param_max_widgets = []
        self.param_fit_value_widgets = []
        self.fit_param_box.addWidget(QLabel('Name'), 0, 0)
        self.fit_param_box.addWidget(QLabel('Initial Guess'), 0, 1)
        self.fit_param_box.addWidget(QLabel('Vary'), 0, 2)
        self.fit_param_box.addWidget(QLabel('Minimum'), 0, 3)
        self.fit_param_box.addWidget(QLabel('Maximum'), 0, 4)
        self.fit_param_box.addWidget(QLabel('Fit Results'), 0, 5)
        for n, key in enumerate(self.fit_params):
            # param_layout = QHBoxLayout()
            # param_layout.addWidget(QLabel(key))
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

            self.fit_param_box.addWidget(QLabel(fit_param_text[key]), n + 1, 0)
            self.fit_param_box.addWidget(param_init_value_widget, n + 1, 1)
            self.fit_param_box.addWidget(param_vary_widget, n + 1, 2)
            self.fit_param_box.addWidget(param_min_widget, n + 1, 3)
            self.fit_param_box.addWidget(param_max_widget, n + 1, 4)
            self.fit_param_box.addWidget(param_fit_value_widget, n + 1, 5)

            self.param_init_value_widgets.append(param_init_value_widget)
            self.param_vary_widgets.append(param_vary_widget)
            self.param_min_widgets.append(param_min_widget)
            self.param_max_widgets.append(param_max_widget)
            self.param_fit_value_widgets.append(param_fit_value_widget)

    def update_fit_kws(self):
        # print('==========Parameters Updated To the Following=======')
        for n, key in enumerate(self.fit_kws):
            if (isinstance(self.fit_kws_value_widgets[n], QSpinBox)
                    or isinstance(self.fit_kws_value_widgets[n], QDoubleSpinBox)):
                self.fit_kws[key] = self.fit_kws_value_widgets[n].value()
            if isinstance(self.fit_kws_value_widgets[n], QLineEdit):
                self.fit_kws[key] = self.fit_kws_value_widgets[n].toPlainText()
        # print(self.fit_kws)

    def update_params(self):
        # print('==========Parameters Updated To the Following=======')
        self.fit_params_nvarys = 0
        for n, key in enumerate(self.fit_params):
            self.fit_params[key].init_value = self.param_init_value_widgets[n].value()
            self.fit_params[key].vary = self.param_vary_widgets[n].isChecked()
            self.fit_params[key].min = self.param_min_widgets[n].value()
            self.fit_params[key].max = self.param_max_widgets[n].value()
            self.fit_params[key].value = self.param_init_value_widgets[n].value()
            if self.fit_params[key].vary:
                self.fit_params_nvarys += 1

    def tb_flx_btnstate(self):
        if self.plot_tb_button.isChecked() == True:
            self.statusBar.showMessage('Plot Brightness Temperature')
            self.spec_in_tb = True
            if 'area_asec2' in self.fit_params.keys():
                del self.fit_params['area_asec2']
                self.update_fit_param_widgets()
        else:
            self.statusBar.showMessage('Plot Flux Density')
            self.spec_in_tb = False
            if not 'area_asec2' in self.fit_params.keys():
                self.fit_params.add_many(('area_asec2', 25., False, 1., 10000., None, None))
                self.update_fit_param_widgets()
        self.init_pgspecplot()
        self.update_pgspec()

    def showqlookimg(self):
        for i in range(self.qlookimgbox.count()):
            self.qlookimgbox.itemAt(i).widget().deleteLater()
        if self.qlookimgbutton.isChecked():
            self.qlookimg_canvas = FigureCanvas(Figure(figsize=(6, 4)))
            self.qlookimg_toolbar = NavigationToolbar(self.qlookimg_canvas, self)
            self.qlookimg_axs = [self.qlookimg_canvas.figure.subplots(1, 1)]
            self.qlookimgbox.addWidget(self.qlookimg_canvas)
            self.qlookimgbox.addWidget(self.qlookimg_toolbar)
            self.qlookimgbutton.setArrowType(Qt.DownArrow)
            self.plot_qlookmap()
            if self.qlookdspecbutton.isChecked():
                self.qlookimgbox.parent().setStretch(0, 1)
                self.qlookimgbox.parent().setStretch(1, 2)
                self.qlookimgbox.parent().setStretch(2, 0)
            else:
                self.qlookimgbox.parent().setStretch(0, 0)
                self.qlookimgbox.parent().setStretch(0, 0)
                self.qlookimgbox.parent().setStretch(1, 0)
        else:
            self.qlookimgbutton.setArrowType(Qt.RightArrow)
            self.qlookimgbox.parent().setStretch(0, 0)
            self.qlookimgbox.parent().setStretch(0, 0)
            self.qlookimgbox.parent().setStretch(1, 0)

    def showqlookdspec(self):
        for i in range(self.qlookdspecbox.count()):
            self.qlookdspecbox.itemAt(i).widget().deleteLater()
        if self.qlookdspecbutton.isChecked():
            self.qlookdspec_canvas = FigureCanvas(Figure(figsize=(6, 4)))
            self.qlookdspec_toolbar = NavigationToolbar(self.qlookdspec_canvas, self)
            self.qlookdspec_ax = self.qlookdspec_canvas.figure.subplots(1, 1)
            self.qlookdspecbox.addWidget(self.qlookdspec_canvas)
            self.qlookdspecbox.addWidget(self.qlookdspec_toolbar)
            self.qlookdspecbutton.setArrowType(Qt.DownArrow)
            if self.qlookimgbutton.isChecked():
                self.qlookdspecbox.parent().setStretch(0, 1)
                self.qlookdspecbox.parent().setStretch(1, 2)
                self.qlookdspecbox.parent().setStretch(2, 0)
            else:
                self.qlookdspecbox.parent().setStretch(0, 0)
                self.qlookdspecbox.parent().setStretch(1, 0)
                self.qlookdspecbox.parent().setStretch(2, 0)
            if self.has_dspec:
                self.plot_dspec()
        else:
            self.qlookdspecbutton.setArrowType(Qt.RightArrow)
            self.qlookdspecbox.parent().setStretch(0, 0)
            self.qlookdspecbox.parent().setStretch(1, 0)
            self.qlookdspecbox.parent().setStretch(2, 0)

    def do_spec_fit(self):
        roi = self.rois[self.roi_group_idx][self.current_roi_idx]
        freqghz_tofit = roi.freqghz_tofit.compressed()
        spec_tofit = roi.spec_tofit.compressed()
        spec_err_tofit = roi.spec_err_tofit.compressed()
        max_nfev = 1000
        ## Set up fit keywords
        fit_kws = self.fit_kws
        if self.fit_method == 'basinhopping':
            fit_kws['minimizer_kwargs'] = {'method': 'Nelder-Mead'}
            max_nfev *= self.fit_kws['niter'] * (self.fit_params_nvarys + 1)
        if self.fit_method == 'nelder':
            fit_kws = {'options': self.fit_kws}
        print(fit_kws)

        if hasattr(self, 'spec_fitplot'):
            self.speccanvas.removeItem(self.spec_fitplot)

        mini = lmfit.Minimizer(gstools.GSCostFunctions.SinglePowerLawMinimizerOneSrc, self.fit_params,
                               fcn_args=(freqghz_tofit,),
                               fcn_kws={'spec': spec_tofit, 'spec_err': spec_err_tofit, 'spec_in_tb': self.spec_in_tb},
                               max_nfev=max_nfev, nan_policy='omit')
        method = self.fit_method
        mi = mini.minimize(method=method, **fit_kws)
        print(method + ' minimization results')
        print(lmfit.fit_report(mi, show_correl=True))
        self.fit_params_res = mi.params
        print('==========Fit Parameters Updated=======')
        for n, key in enumerate(self.fit_params_res):
            self.param_fit_value_widgets[n].setValue(self.fit_params_res[key].value)

        freqghz_toplot = np.logspace(0, np.log10(20.), 100)
        spec_fit_res = self.fit_function(mi.params, freqghz_toplot, spec_in_tb=self.spec_in_tb)
        self.spec_fitplot = self.speccanvas.plot(x=np.log10(freqghz_toplot), y=np.log10(spec_fit_res),
                                                 pen=dict(color=pg.mkColor(self.current_roi_idx), width=4),
                                                 symbol=None, symbolBrush=None)
        self.speccanvas.addItem(self.spec_fitplot)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
