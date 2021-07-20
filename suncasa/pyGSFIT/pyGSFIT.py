import sys, os, traceback, glob
from time import sleep
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QPushButton, QCheckBox, QLineEdit, QLabel, \
    QTextEdit, QListWidgetItem, QWidget, QAction, QTabWidget, QTableView, QVBoxLayout, QHBoxLayout, QGridLayout, \
    QStatusBar, QProgressBar, QFrame, QFormLayout, QGroupBox, QButtonGroup, QRadioButton, QComboBox, QSizePolicy

from PyQt5.QtGui import QIcon, QFont, QColor
from PyQt5.QtCore import pyqtSlot, Qt, QAbstractTableModel, pyqtSignal, QRunnable, QObject, QThreadPool
from casatools import ms as mstool

from astropy.time import Time
from casatools import table
from astropy.io import fits
from copy import deepcopy
from functools import partial
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from matplotlib.dates import DateFormatter
import numpy as np

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.fname = '<Select or enter a valid fits file name>'
        self.fitsentry = QLineEdit()
        self.title = 'pyGSFIT'
        self.left = 0
        self.top = 0
        self.width = 1200
        self.height = 900
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self._main = QWidget()
        self.setCentralWidget(self._main)
        self.initUI()
        self.threadpool = QThreadPool()


    def initUI(self):

        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.progressBar = QProgressBar()
        self.progressBar.setGeometry(10,10,200,15)

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

        self.tabs.currentChanged.connect(self.tabChanged)

        # Add tabs to widget
        layout.addWidget(self.tabs)
        self._main.setLayout(layout)

        self.show()

    def tabChanged(self, i):
        if i == 2:
            self.update_params()

    # Explorer Tab User Interface
    #
    def initUItab_explorer(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()

        # Create Data Select tab
        upperbox = QHBoxLayout()  # Upper box has two hboxes: leftbox and quicklook plot box
        leftbox = QVBoxLayout()   # Has a gridlayout: filenamebox; and groupBox: selectBox

        # Filename entry
        filenamebox = QHBoxLayout()
        leftbox.addLayout(filenamebox)
        # Create Browse button
        myButton1 = QPushButton("Browse")
        myButton1.clicked.connect(self.on_click)
        filenamebox.addWidget(myButton1)
        # Create LineEdit widget for ms filename
        self.fitsentry.resize(8*len(self.fname), 20)
        filenamebox.addWidget(QLabel("FITS Filename"))
        filenamebox.addWidget(self.fitsentry)

        upperbox.addLayout(leftbox)

        # Create label and TextEdit widget for FITS information
        fitsinfobox = QVBoxLayout()
        leftbox.addLayout(fitsinfobox)
        self.infoEdit = QTextEdit()
        self.infoEdit.setReadOnly(True)
        self.infoEdit.setMinimumHeight(100)
        self.infoEdit.setMaximumHeight(200)
        self.infoEdit.setMinimumWidth(300)
        self.infoEdit.setMaximumWidth(500)
        #fitsinfobox.addWidget(QLabel("FITS Information"))
        fitsinfobox.addWidget(self.infoEdit)


        # Add quicklook plotting area
        qlookarea = QVBoxLayout()
        self.qlookcanvas = FigureCanvas(Figure(figsize=(6, 6)))
        qlookarea.addWidget(self.qlookcanvas)
        self.plot_axs = self.qlookcanvas.figure.subplots(nrows=1, ncols=1)
        upperbox.addLayout(qlookarea)

        mainlayout.addLayout(upperbox)

        # Add more plotting area
        plotarea = QHBoxLayout()
        self.plotcanvas = FigureCanvas(Figure(figsize=(8, 4)))
        plotarea.addWidget(self.plotcanvas)
        self.plot_axs = self.plotcanvas.figure.subplots(nrows=1, ncols=2)
        #self.plot_axs.clear()
        #im = self.dspec_ax.pcolormesh(pd, fghz,np.clip(spec, minval, maxval))
        #self.dspec_ax.xaxis_date()
        #self.dspec_ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        #self.dspec_ax.set_ylabel('Frequency [GHz]')
        #self.dspec_ax.set_xlabel('Time [UT]')

        mainlayout.addLayout(plotarea)
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

    def on_return(self):
        ''' Called when the FITS filename LineEdit widget gets a carriage-return.
            Trys to read FITS header and return header info (no data read at this time)
        '''
        self.fname = self.fitsentry.text()
        self.fitsdata = None
        try:
            print(self.fname)
        except:
            self.statusBar.showMessage('Filename is not a valid FITS file', 2000)
            self.fname = '<Select or enter a valid fits filename>'
            self.fitsentry.setText(self.fname)
            self.infoEdit.setPlainText('')

    def on_click(self):
        ''' Handle Browse button '''
        #self.fname = QFileDialog.getExistingDirectory(self, 'Select FITS File', './', QFileDialog.ShowDirsOnly)
        self.fname, _file_filter = QFileDialog.getOpenFileName(self, 'Select FITS File to open', './', "FITS Images (*.fits *.fit *.ft)")
        self.fitsentry.setText(self.fname)
        self.on_return()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
