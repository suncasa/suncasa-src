import sys, os, traceback, glob
from time import sleep
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QPushButton, QCheckBox, QLineEdit, QLabel, QTextEdit, QListWidgetItem, QWidget, QAction, QTabWidget, QTableView, QVBoxLayout, QHBoxLayout, QGridLayout, QStatusBar, QProgressBar, QFrame, QFormLayout, QGroupBox, QButtonGroup, QRadioButton, QComboBox, QSizePolicy

from PyQt5.QtGui import QIcon, QFont, QColor
from PyQt5.QtCore import pyqtSlot, Qt, QAbstractTableModel, pyqtSignal, QRunnable, QObject, QThreadPool
from pr_summary import pr_summary
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


class WorkerSignals(QObject):
    '''
    Defines the finished signal available from a running worker thread.
    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)


class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        #self.kwargs['progress_callback'] = self.signals.progress

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            #print(self.fn)
            #print(*self.args)
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        finally:
            self.signals.finished.emit()  # Done

class TableModel(QAbstractTableModel):
    ''' Defines the action of the tclean table
    '''
    def __init__(self, data):
        super(TableModel, self).__init__()
        self._data = data
        
    def data(self, index, role):
        if role == Qt.DisplayRole:
            if index.column() > 1: return ""
            return self._data[index.row()][index.column()]
        if role == Qt.BackgroundRole:
            if self._data[index.row()][0][0] != ' ':
                return QColor('#ffccaa')
        if role == Qt.EditRole:
            return self._data[index.row()][index.column()]
        if role == Qt.ToolTipRole:
            if index.column() == 0:
                if len(self._data[index.row()]) == 4:
                    return self._data[index.row()][3]
            
    def setData(self, index, value, role):
        if role == Qt.EditRole:
            #value = self.checkdata(index, value)
            self._data[index.row()][index.column()] = value
            self.dataChanged.emit(index, index)
            return True
            
    def rowCount(self, index):
        return len(self._data)
        
    def columnCount(self, index):
        return len(self._data[0])

    def flags(self, index):
        if index.column() !=3:
            flag_item = Qt.ItemIsSelectable|Qt.ItemIsEnabled
        if index.column() == 1 and len(self._data[index.row()][2]) == 0: 
            flag_item |= Qt.ItemIsEditable
        return flag_item 

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'EOVSA Imager'
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
        tab1 = QWidget()
        tab2 = QWidget()
        tab3 = QWidget()
        tab4 = QWidget()
        tab5 = QWidget()
        tab6 = QWidget()

        # Add tabs
        self.tabs.addTab(tab1,"Data Select")
        self.tabs.addTab(tab2,"Viewer")
        self.tabs.addTab(tab3,"Imager")
        self.tabs.addTab(tab4,"Selfcal")
        self.tabs.addTab(tab5,"Production")
        self.tabs.addTab(tab6,"Export")

        # Each tab's user interface is complex, so this splits them into separate functions.
        self.initUItab1()
        self.initUItab2()
        self.initUItab3()
        self.initUItab4()
        self.initUItab5()
        self.initUItab6()

        self.tabs.currentChanged.connect(self.tabChanged)
        
        # Add tabs to widget
        layout.addWidget(self.tabs)
        self._main.setLayout(layout)

        self.show()

    def tabChanged(self, i):
        if i == 2:
            self.update_params()

#
# Data Select Tab User Interface        
#
    def initUItab1(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()

        # Create Data Select tab
        upperbox = QHBoxLayout()  # Has two hboxes: leftbox and msinfobox
        leftbox = QVBoxLayout()   # Has a gridlayout: filenamebox; and groupBox: selectBox
        # Filename entry
        filenamebox = QGridLayout()
        leftbox.addLayout(filenamebox)
        # Create LineEdit widget for ms filename
        self.msentry = QLineEdit()
        self.fname = '<Select or enter a valid ms filename>'
        self.msentry.resize(8*len(self.fname),20)
        self.ms = None       # No ms yet
        self.msdata = None  # No amplitude data read from ms yet
        self.msentry.setText(self.fname)
        self.msentry.returnPressed.connect(self.on_return)
        filenamebox.addWidget(QLabel("MS Filename"),0,0)
        filenamebox.addWidget(self.msentry,1,0,1,4)
        # Create Browse button
        myButton1 = QPushButton("Browse")
        myButton1.clicked.connect(self.on_click)
        filenamebox.addWidget(myButton1,1,5)
        upperbox.addLayout(leftbox)

        # Create label and TextEdit widget for ms information
        msinfobox = QVBoxLayout()
        upperbox.addLayout(msinfobox)
        self.infoEdit = QTextEdit()
        #f = QFont("Courier",9)
        #self.infoEdit.setCurrentFont(f)
        self.infoEdit.setReadOnly(True)
        self.infoEdit.setMinimumHeight(300)
        self.infoEdit.setMinimumWidth(550)
        msinfobox.addWidget(QLabel("MS Information"))
        msinfobox.addWidget(self.infoEdit)
        mainlayout.addLayout(upperbox)

        # Data Selection
        selectBox = QGroupBox("Data Selection Criteria")
        selectarea = QFormLayout()
        self.trangeEdit = QLineEdit()
        selectarea.addRow("Timerange", self.trangeEdit)
        self.spwEdit = QLineEdit()
        selectarea.addRow("Sp. Window", self.spwEdit)
        self.baselineEdit = QLineEdit()
        selectarea.addRow("Baseline", self.baselineEdit)
        self.stokesEdit = QLineEdit()
        selectarea.addRow("Stokes", self.stokesEdit)
        self.uvrangeEdit = QLineEdit()
        selectarea.addRow("UV range", self.uvrangeEdit)
        selectBox.setLayout(selectarea)
        leftbox.addWidget(selectBox)
        drawselect = QPushButton('Draw Selection on Plot')
        drawselect.clicked.connect(self.drawsel)
        leftbox.addWidget(drawselect)
        leftbox.addStretch(1)
        
        #hbox.addStretch(1)
        
        playarea = QHBoxLayout()
        playarea.addStretch(1)
        xferarea = QVBoxLayout()
        xferDownRight = QPushButton('Use Selection for Plot Limits')
        xferDownRight.setIcon(QIcon('icons/down_right.png'))
        xferDownRight.clicked.connect(self.xferDR)
        xferarea.addWidget(xferDownRight)
        xferUpLeft = QPushButton('Use Plot Limits for Selection')
        xferUpLeft.setIcon(QIcon('icons/up_left.png'))
        xferUpLeft.clicked.connect(self.xferUL)
        xferarea.addWidget(xferUpLeft)
        xferarea.addStretch(1)
        playarea.addLayout(xferarea)
        # Add a figure to the canvas.
        plotarea = QVBoxLayout()
        self.speccanvas = FigureCanvas(Figure(figsize=(8, 3)))
        plotarea.addWidget(self.speccanvas)
        xferarea.addWidget(NavigationToolbar(self.speccanvas, self))
        maxmin = QHBoxLayout()
        minlabel = QLabel('Min[sfu]')
        maxlabel = QLabel('Max[sfu]')
        self.minentry = QLineEdit()
        self.maxentry = QLineEdit()
        maxmin.addStretch(1)
        maxmin.addWidget(minlabel)
        maxmin.addWidget(self.minentry)
        maxmin.addWidget(maxlabel)
        maxmin.addWidget(self.maxentry)
        maxmin.addStretch(1)
        self.ignore_gaps = QCheckBox('Ignore Gaps')
        maxmin.addWidget(self.ignore_gaps)
        maxmin.addStretch(1)
        self.minentry.returnPressed.connect(self.plot_data)
        self.maxentry.returnPressed.connect(self.plot_data)
        plotarea.addLayout(maxmin)
        playarea.addLayout(plotarea)
#        self.addToolBar(Qt.BottomToolBarArea, NavigationToolbar(self.speccanvas, self))
        self.dspec_ax = self.speccanvas.figure.subplots()
        baseline_layout = QGridLayout()
        baseline_layout.setSpacing(1)
        ant1 = []
        ant2 = []
        nant = 13
        bl2ord = get_bl2ord(nant)     # Lookup table for 13-ant baseline matrix
        nbl = int(nant*(nant+1)/2)
        self.blcheck = [0]*nbl        # This will hold the check box widgets for the baselines
        self.blchecked = [False]*nbl  # This is the state of the buttons
        self.prev = None
        for i in range(nant):
            ant1.append(QPushButton('{:2d}'.format(i+1)))
            ant1[-1].setCheckable(True)
            ant1[-1].setMaximumSize(12,12)   # Ant button size in pixels
            ant1[-1].setStyleSheet("font : 8px;")
            ant1[-1].clicked.connect(partial(self.ant_click,ant1[-1]))
            baseline_layout.addWidget(ant1[-1],0,i)    # Grid location of ant buttons along top
            ant2.append(QPushButton('{:2d}'.format(i+1)))
            ant2[-1].setCheckable(True)
            ant2[-1].setMaximumSize(12,12)   # Ant button size in pixels
            ant2[-1].setStyleSheet("font : 8px;")
            ant2[-1].clicked.connect(partial(self.ant_click,ant2[-1]))
            baseline_layout.addWidget(ant2[-1],i+1,nant)  # Grid location of ant buttons along side
            for j in range(i,nant):
                if i == j:
                    self.blcheck[bl2ord[i,j]] = QCheckBox("")
                    button = self.blcheck[bl2ord[i,j]]  # Just saves typing
                    button.setStyleSheet("background-color : #8888ff;")
                else:
                    self.blcheck[bl2ord[i,j]] = QCheckBox("")
                    button = self.blcheck[bl2ord[i,j]]  # Just saves typing
                    button.setStyleSheet("background-color : #ff8888;")
                button.setMaximumSize(12,12)
                button.clicked.connect(self.baselineClicked)
                baseline_layout.addWidget(button,i+1,j)   # Grid location of baseline buttons
        self.ant1buttons = ant1
        self.ant2buttons = ant2
        unsel = QPushButton('Unselect All')
        unsel.setStyleSheet("font : 10px;")
        unsel.clicked.connect(partial(self.sel_clicked, False))
        sel = QPushButton('Select All')
        sel.setStyleSheet("font : 10px;")
        sel.clicked.connect(partial(self.sel_clicked, True))
        baseline_layout.addWidget(unsel,9,0,2,7)
        baseline_layout.addWidget(sel,11,0,2,7)
        space = QVBoxLayout()
        # Add Polarization Button Group
        self.polGroup = QButtonGroup()
        self.polGroup.buttonClicked.connect(self.get_data)
        xxButton = QRadioButton("XX")
        xxButton.setChecked(True)
        yyButton = QRadioButton("YY")
        xyButton = QRadioButton("XY")
        yxButton = QRadioButton("YX")
        rrButton = QRadioButton("RR")
        llButton = QRadioButton("LL")
        self.polGroup.addButton(xxButton)
        self.polGroup.addButton(yyButton)
        self.polGroup.addButton(xyButton)
        self.polGroup.addButton(yxButton)
        self.polGroup.addButton(rrButton)
        self.polGroup.addButton(llButton)
        pspace = QHBoxLayout()
        pspace.addWidget(xxButton)
        pspace.addWidget(yyButton)
        pspace.addWidget(xyButton)
        pspace.addWidget(yxButton)
        pspace.addWidget(rrButton)
        pspace.addWidget(llButton)
        space.addStretch(1)
        space.addLayout(pspace)
        space.addWidget(QLabel('Baseline Selection Map'))
        space.addLayout(baseline_layout)
        space.addStretch(1)
        gobutton = QPushButton('Go')
        gobutton.clicked.connect(self.get_data)
        gobutton.setStyleSheet("background-color : #ffff88;")
        space.addWidget(gobutton)
        playarea.addLayout(space)
        mainlayout.addLayout(playarea)
        self.tabs.widget(0).setLayout(mainlayout)

#
# Viewer Tab User Interface        
#
    def initUItab2(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the viewer tab"))
        self.tabs.widget(1).setLayout(mainlayout)

#
# Imager Tab User Interface        
#
    def initUItab3(self):
        # Create main layout (a Horizontal Layout)
        mainlayout = QHBoxLayout()
        # Create a left and right side (both Vertical Layouts)
        leftlayout = QVBoxLayout()
        rightlayout = QVBoxLayout()
        
        # Create a table for interacting with parameters
        self.table = QTableView()
        self.table.doubleClicked.connect(self.tblRowClicked)
        # For playing with, I will make a static table of tclean parameters
        tbl = [
                  ["Data Selection","",[]],
                  ["  vis","''",[]],
                  #["  field","''",[]],
                  ["  spw","''",[],"Spectral Window:\n default: ''=all; examples:\n spw='0~2,4'; spectral windows 0,1,2,4 (all channels)\n spw='0:5~61'; spw 0, channels 5 to 61\n spw='<2';   spectral windows less than 2 (i.e. 0,1)" ],
                  ["  timerange","''",[],"Range of time to select from data\n default: '' (all); examples:\n timerange = 'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'\n Note: if YYYY/MM/DD is missing date defaults to first day in data set\n timerange='09:14:0~09:54:0' picks 40 min on first day\n timerange='25:00:00~27:30:00' picks 1 hr to 3 hr 30 min on NEXT day\n timerange='09:44:00' pick data within one integration of time\n timerange='> 10:24:00' data after this time"],
                  ["  uvrange","''",[],"Select data within uvrange (default unit is meters) [default: '' (all)]\n examples:\n uvrange='0~1000klambda'; uvrange from 0-1000 kilo-lambda\n uvrange='> 4klambda';uvranges greater than 4 kilo lambda"],
                  ["  antenna","''",[],"Select data on 0-based antenna/baseline index [default: '' (all)]\n examples:\n antenna='0~5,7~12&0~5,7~12'; all baselines not including antenna index 6\n antenna='5&6;7&8'; baselines 5-6 and 7-8\n antenna='5'; all baselines with antenna index 5\n antenna='5,6,9'; all baselines with antenna index numbers 5,6,9"],
                  ["  datacolumn","'data'",["'data'","'corrected'","'model'"]],
                  ["Image Definition","",[]],
                  ["  imagename","''",[]],
                  ["  imsize","[128,128]",[]],
                  ["  cellsize","'2arcsec'",[]],
                  ["  phaseshift","[0, 0]",[],"X, Y offset of center of map from Sun Center"],
                  ["  stokes","'XX'",["'XX'","'YY'","'I'","'V'","'IV'","'RR'","'LL'"]],
                  ["  startmodel", "''",[]],
                  ["  specmode","'mfs'",["'mfs'","'cubedata'"]],
                  ["Deconvolution Options","",[]],
                  ["  deconvolver","'multiscale'",["'hogbom'","'clark'","'clarkstokes'","'multiscale'","'mem'"]],
                  ["  scales", "[1,5,10]",[]],
                  ["  restoringbeam","''",[]],
                  ["  pbcor","False",['True','False']],
                  ["Weighting","",[]],
                  ["  weighting","'briggs'",["'natural'","'uniform'","'briggs'"]],
                  ["  robust","0.5",[]],
                  ["  uvtaper","''",[]],
                  ["Other Options","",[]],
                  ["  niter","0",[]],
                  ["  gain","0.1",[]],
                  ["  threshold","0",[]],
                  ["  interactive","False",['True','False']],
                  ["  mask","''",[]]]
        
        self.table.verticalHeader().setDefaultSectionSize(14)  # Sets height of cells to 14px
        self.table.setModel(TableModel(tbl))  # The TableModel class is defined at the top of this file
        
        # For parameters that have to be only a limited set of fixed values, create a combobox dropdown
        # for editing them.
        for idx in range(len(tbl)):
            if len(tbl[idx]) > 2:
                if len(tbl[idx][2]) > 1:
                    i = self.table.model().index(idx,2)
                    c = QComboBox()
                    for item in tbl[idx][2]:
                        c.addItem(item)
                    c.currentTextChanged.connect(self.handleCombo)
                    c.index = self.table.model().index(idx,1)
                    self.table.setIndexWidget(i,c)

        self.table.resizeColumnsToContents()
        self.table.model().dataChanged.connect(self.update_view)

        # Determine the header rows in the table (indicated by NOT starting with a blank space).
        self.headerrows = []
        for i,tblrow in enumerate(tbl):
            if tbl[i][0][0] != ' ':
                self.headerrows.append(i)
                self.table.setSpan(i,0,1,2)
        self.headerrows.append(len(tbl))  #Add length of table, for finding length of last section
        
        titlelayout = QHBoxLayout()
        titlelayout.addWidget(QLabel("TCLEAN Parameters"))
        titlelayout.addSpacing(100)
        updateButton = QPushButton("Update Parameters")
        titlelayout.addWidget(updateButton)
        updateButton.clicked.connect(self.update_params)
        titlelayout.addSpacing(100)
        scriptButton = QPushButton("Save to CASA Script")
        titlelayout.addWidget(scriptButton)
        scriptButton.clicked.connect(self.save2CASAscript)
        titlelayout.addStretch(1)
        tablelayout = QHBoxLayout()
        self.table.setMinimumSize(600,300)
        tablelayout.addWidget(self.table)
        self.nofits = QCheckBox('Skip conversion to FITS')
        self.nocleanup = QCheckBox('Keep CASA image files')
        tablelayout.addStretch(1)
        leftlayout.addLayout(titlelayout)
        leftlayout.addLayout(tablelayout)
        leftlayout.addWidget(self.nofits)
        leftlayout.addWidget(self.nocleanup)
        self.scriptEdit = QTextEdit()
        #f = QFont("Courier",9)
        #self.infoEdit.setCurrentFont(f)
        self.scriptEdit.setReadOnly(True)
        self.scriptEdit.setMinimumHeight(300)
        self.scriptEdit.setMinimumWidth(550)
        leftlayout.addWidget(QLabel("Generated Script"))
        leftlayout.addWidget(self.scriptEdit)
        execButton = QPushButton('Execute Script')
        execButton.clicked.connect(self.execscript)
        eblayout = QHBoxLayout()
        eblayout.addWidget(execButton)
        eblayout.addStretch(1)
        leftlayout.addLayout(eblayout)
        leftlayout.addStretch(1)
        mainlayout.addLayout(leftlayout)
        
        self.imgcanvas = FigureCanvas(Figure(figsize=(7, 6)))
        rightlayout.addWidget(self.imgcanvas)
        rightlayout.addWidget(NavigationToolbar(self.imgcanvas, self))
        self.img_ax = self.imgcanvas.figure.subplots()

        mainlayout.addLayout(rightlayout)
        
        self.tabs.widget(2).setLayout(mainlayout)
        
    def save2CASAscript(self):
        script = ['from casatasks import tclean']
        tbl = self.table.model()._data
        for k,row in enumerate(tbl):
            if row[0][0] == ' ':
                # This is a parameter row.  Most such rows can be defined directly, but some require translation
                param = row[0].lstrip()
                if param == 'phaseshift':
                    # Determine appropriate phase center and set accordingly
                    try:
                        pshift = [float(i) for i in row[1].replace('[','').replace(']','').split(',')]
                    except:
                        self.statusBar.showMessage('Error translating phaseshift parameter--using zero',2000)
                        pshift = [0.0,0.0]
                    script.append("phasecenter='"+pshift2pcenter(pshift,self.ms)+"'")
                elif param == 'imagename':
                    if row[1] == "''":
                        msstem = os.path.basename(os.path.splitext(self.fname)[0])
                        script.append("imagename='images/"+msstem+"'")
                    else:
                        script.append(param+"="+row[1])
                else:
                    if row[1] == "''":
                        script.append(param+"=''")
                    else:
                        script.append(param+"="+row[1])

        script.append("tclean(vis=vis, selectdata=True, field='', spw=spw, timerange=timerange, uvrange=uvrange, antenna=antenna, scan='', observation='', intent='', datacolumn=datacolumn, imagename=imagename, imsize=imsize, cell=cellsize, phasecenter=phasecenter, stokes=stokes, projection='SIN', startmodel='', specmode=specmode, reffreq='', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='standard', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=True, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver=deconvolver, scales=scales, nterms=2, smallscalebias=0.0, restoration=True, restoringbeam=restoringbeam, pbcor=pbcor, outlierfile='', weighting=weighting, robust=robust, noise='1.0Jy', npixels=0, uvtaper=uvtaper, niter=niter, gain=gain, threshold=threshold, nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=interactive, usemask='user', mask=mask, pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)")
        if not self.nofits.isChecked():
            script.append("import helioim2fits as hf")
            script.append("fitsfile = hf.imreg(vis=vis,imagefile=imagename+'.image',fitsdir='fits',timerange=timerange)")
        if not self.nocleanup.isChecked():
            script.append("import glob, shutil")
            script.append("for file in glob.glob(imagename+'.*'): shutil.rmtree(file)")
            if self.nofits.isChecked():
                self.statusBar.showMessage('Warning! You have selected no FITs output and do not keep CASA images, so you will get no output!',2000)
        self.script = script
        self.scriptEdit.setPlainText('\n'.join(script))

    def execscript(self):
        #print('Before call:')
        #print(exec)
        #print(['\n'.join(self.script)])
        worker = Worker(exec, '\n'.join(self.script))
        worker.signals.finished.connect(self.thread_complete)
        self.threadpool.start(worker)
        print('The thread is started--just waiting for the signal.')
        
    def thread_complete(self):
        print('The script execution thread is done!')
        print('The fits file should have been completed, and can be read now...')
        sleep(1) # Make sure file is closed...
        list_of_files = glob.glob('fits/*')
        latest_file = max(list_of_files, key=os.path.getctime)
        img, h = fits.getdata(latest_file, header=True)
        img.shape = (h['NAXIS2'],h['NAXIS1'])
        self.img_ax.cla()
        xval = np.linspace((h['CRVAL1']-h['NAXIS1'])*h['CDELT1']/2., (h['CRVAL1']+h['NAXIS1'])*h['CDELT1']/2.,h['NAXIS1'])
        yval = np.linspace((h['CRVAL2']-h['NAXIS2'])*h['CDELT2']/2., (h['CRVAL2']+h['NAXIS2'])*h['CDELT2']/2.,h['NAXIS2'])
        im = self.img_ax.pcolormesh(xval, yval, img)
        
        self.img_ax.set_aspect('equal')
        self.img_ax.set_ylabel('Y [arcsec]')
        self.img_ax.set_xlabel('X [arcsec]')
        self.img_ax.set_title(latest_file)
        self.imgcanvas.draw()
        
    def update_params(self):
        ''' Take entries from Data Select tab and update the tclean parameter table
        '''
        tbl = self.table.model()._data
        params = ['vis', 'spw', 'timerange', 'uvrange', 'antenna', 'stokes']
        curvalues = [self.fname, self.spwEdit.text(), self.trangeEdit.text(), self.uvrangeEdit.text(), 
                                    self.baselineEdit.text(), self.stokesEdit.text()]
        for k,row in enumerate(tbl):
            r = row[0][2:]  # parameter name (after removing two blank spaces)
            for idx, p in enumerate(params):
                if r == p:
                    i = self.table.model().index(k,1)   # Table index of parameter value 
                    # Write the new value into the table (all of these are strings, so include quotes)
                    self.table.model().setData(i, "'"+curvalues[idx]+"'", Qt.EditRole)
                    if r == 'stokes':
                        # Stokes has a combobox, so set it according to the updated parameter
                        i2 = self.table.model().index(k,2)  # Table index of combobox widget
                        c = self.table.indexWidget(i2)
                        loc = c.findText("'"+curvalues[idx]+"'")
                        if loc == -1:
                            # If curvalue is not found in the combobox, indicate stokes as unset.
                            self.table.model().setData(i, '<unset>', Qt.EditRole)
                        else: 
                            c.setCurrentIndex(loc)

    def update_view(self):
        self.table.resizeColumnsToContents()
    
    def handleCombo(self):
        sender = self.sender()
        text = sender.currentText()
        self.table.model().setData(sender.index, text, Qt.EditRole)
        
    def tblRowClicked(self, index):
        ''' If a header row is double-clicked, show or hide the rows corresponding to that section
        '''
        hrows = self.headerrows   # Row numbers corresponding to headers   
        for n,row in enumerate(hrows):
            if index.row() == row:
                # The row of the clicked cell is a header row, so identify the lines between it
                # and the next header and toggle hidden
                k1 = hrows[n] + 1
                k2 = hrows[n+1]
                if self.table.isRowHidden(k1):
                    # Rows are aleady hidden, so show them
                    for i in range(k1,k2):
                        self.table.showRow(i)
                else:
                    # Rows are visible, so hide them
                    for i in range(k1,k2):
                        self.table.hideRow(i)
                break
#     
# Selfcal Tab User Interface        
#
    def initUItab4(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the selfcal tab"))
        self.tabs.widget(3).setLayout(mainlayout)

#
# Production Tab User Interface        
#
    def initUItab5(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the production tab"))
        self.tabs.widget(4).setLayout(mainlayout)

#
# Export Tab User Interface        
#
    def initUItab6(self):
        # Create main layout (a Vertical Layout)
        mainlayout = QVBoxLayout()
        mainlayout.addWidget(QLabel("This is the export tab"))
        self.tabs.widget(5).setLayout(mainlayout)

    def baselineClicked(self):
        button = self.sender()
        if button.isChecked():
            if not self.prev is None:
                if self.prev.isChecked():
                    self.prev.setChecked(False)
            self.prev = button
            self.get_data()
        
    def on_return(self):
        ''' Called when the ms filename LineEdit widget gets a carriage-return.
            Trys to connect to the ms and return some info (no data read at this time)
        '''
        self.fname = self.msentry.text()
        self.msdata = None
        try:
            if self.ms:
                self.ms.close()
            else:
                self.ms = mstool()
            self.ms.open(self.fname)
            self.msentry.setText(self.fname)
            lines = pr_summary(self.ms)
            self.infoEdit.setPlainText('\n'.join(lines))
        except:
            self.statusBar.showMessage('Filename is not a valid ms',2000)
            self.fname = '<Select or enter a valid ms filename>'
            self.msentry.setText(self.fname)
            self.infoEdit.setPlainText('')
        try:
            t1str, t2str = lines[2].replace('-','/').split(': ')[1:3]
            trange = t1str[0:10]+'/'+t1str[11:21]+'~'+t2str[0:10]+'/'+t2str[11:21]
        except:
            trange = '<unset>'
        try:
            spwstr = lines[-1].strip().split(' ')[0]
            spw = '0~'+spwstr
        except:
            spw = '<unset>'
        self.trangeEdit.setText(trange)
        self.spwEdit.setText(spw)
        self.baselineEdit.setText('0~12&0~12')
        self.stokesEdit.setText('XX')
        self.uvrangeEdit.setText('0~2km')
        if spwstr == '30':
            self.ignore_gaps.setChecked(True)
        else:
            self.ignore_gaps.setChecked(False)

    def on_click(self):
        ''' Handle Browse button '''
        self.fname = QFileDialog.getExistingDirectory(self, 'Select MS','./',QFileDialog.ShowDirsOnly)
        self.msentry.setText(self.fname)
        self.on_return()

    def ant_click(self, button):
        ''' An antenna button was clicked, so check the checkboxes corresponding
            to the currently checked antennas.  button is a dictionary with keys
            'button' (actual button handle) and 'iant' (ordinal location in list
            of buttons).  
            Each antenna has two buttons, so we have to ensure that both are toggled.  
            iant is an integer from 0-25.  If even, toggle buttons iant and iant+1.  
            If odd, toggle buttons iant and iant-1.
        '''
        k = 0
        ant = -1
        for b in self.ant1buttons:
            if b == button:
                self.ant2buttons[k].setChecked(button.isChecked())
                ant = k
                break
            else:
                k+=1
        k = 0
        for b in self.ant2buttons:
            if b == button:
                self.ant1buttons[k].setChecked(button.isChecked())
                ant = k
                break
            else:
                k+=1
        self.statusBar.showMessage('Antenna '+str(ant+1)+' hit',1000)
        nant = 13
        bl2ord = get_bl2ord(nant)
        # Reflect antenna state on baselines
        for i in range(nant):
            for j in range(i,nant):
                if i != j:
                    if i == ant:
                        self.blcheck[bl2ord[i,j]].setChecked(button.isChecked())
                    if j == ant:
                        self.blcheck[bl2ord[i,j]].setChecked(button.isChecked())
        # Have to go back and set any baselines whose antennas are checked.
        for i in range(nant):
            if self.ant1buttons[i].isChecked():
                for k in range(i):
                    self.blcheck[bl2ord[k,i]].setChecked(True)
                for j in range(i+1,nant):
                    self.blcheck[bl2ord[i,j]].setChecked(True)
                    
    def sel_clicked(self, state):
        nant = 13
        bl2ord = get_bl2ord(nant)
        for i in range(nant):
            self.ant1buttons[i].setChecked(state)
            self.ant2buttons[i].setChecked(state)
            for j in range(i,nant):
                self.blcheck[bl2ord[i,j]].setChecked(state)

    def get_data(self):
        # Find which baselines are checked:
        bl = []
        blid = []
        isauto = []
        ord13 = get_bl2ord(13)
        ord16 = get_bl2ord(16)
        for i in range(13):
            for j in range(i,13):
                if self.blcheck[ord13[i,j]].isChecked():
                    bl.append(ord16[i,j])
                    isauto.append(i == j)
                    blid.append([i,j])
        if len(bl) == 0:
            self.statusBar.showMessage('No baselines have been selected.',2000)
            return
        if len(bl) != 1:
            self.statusBar.showMessage('Only one baseline can be shown.',2000)
            return

        # If no data have been read, read the data
        if self.msdata is None:
            self.read_msdata(fillnan=np.nan)
        npol, nf, nbl, nt = self.msdata['data'].shape
        self.plotbl = bl[0]

        # Determine which polarization state is selected by radiobuttons
        polnum = {'XX':0, 'YY':1, 'XY':2, 'YX':3, 'RR':4, 'LL':5}
        ipol = polnum[self.polGroup.checkedButton().text()]
        
        if ipol < 4:
            # The selected polarization is just one of the native polarizations (XX, YY, XY, or YX)
            spec = np.abs(self.msdata['data'][ipol,:,self.plotbl])
        elif ipol == 4:
            # The selected polarization is RR (not implemented because we only have amplitudes...)
            xx = self.msdata['data'][0,:,self.plotbl]
            yy = self.msdata['data'][1,:,self.plotbl]
            xy = self.msdata['data'][2,:,self.plotbl]
            yx = self.msdata['data'][3,:,self.plotbl]
            rr = (xx + yy + 1j*(xy - yx))/2
            spec = np.abs(rr)  
        elif ipol == 5:
            # The selected polarization is LL (not implemented because we only have amplitudes...)
            xx = self.msdata['data'][0,:,self.plotbl]
            yy = self.msdata['data'][1,:,self.plotbl]
            xy = self.msdata['data'][2,:,self.plotbl]
            yx = self.msdata['data'][3,:,self.plotbl]
            ll = (xx + yy - 1j*(xy - yx))/2
            spec = np.abs(ll)  

        i, j = blid[0] 
        if isauto[0]:
            self.ptitle = 'Auto-Correlation for Ant {}'.format(i+1)
        else:
            self.ptitle = 'Cross-Correlation for Baseline {}-{}'.format(i+1,j+1)
        spec1d = spec.flatten()
        spec1dg = np.sort(spec1d[~np.isnan(spec1d)])
        if len(spec1dg) == 0:
            self.statusBar.showMessage('This baseline is all NaN',2000)
            self.dspec_ax.set_title('This baseline is all NaN')
            self.speccanvas.draw()
        else:
            n95 = int(len(spec1dg)*0.95)
            ampmax = spec1dg[n95]
            ampmin = spec1dg[0]
            self.minentry.setText('{:10.2f}'.format(ampmin/10000.))
            self.maxentry.setText('{:10.2f}'.format(ampmax/10000.))
            self.plot_data()

    def plot_data(self):
        ''' Plot the data selected by get_data().
        '''
        # Determine which polarization state is selected by radiobuttons
        polnum = {'XX':0, 'YY':1, 'XY':2, 'YX':3, 'RR':4, 'LL':5}
        ipol = polnum[self.polGroup.checkedButton().text()]
        
        if ipol < 4:
            # The selected polarizatin is just one of the native polarizations (XX, YY, XY, or YX)
            spec = np.abs(self.msdata['data'][ipol,:,self.plotbl])
        elif ipol == 4:
            # The selected polarization is RR
            xx = self.msdata['data'][0,:,self.plotbl]
            yy = self.msdata['data'][1,:,self.plotbl]
            xy = self.msdata['data'][2,:,self.plotbl]
            yx = self.msdata['data'][3,:,self.plotbl]
            rr = (xx + yy + 1j*(xy - yx))/2
            spec = np.abs(rr)  
        elif ipol == 5:
            # The selected polarization is LL
            xx = self.msdata['data'][0,:,self.plotbl]
            yy = self.msdata['data'][1,:,self.plotbl]
            xy = self.msdata['data'][2,:,self.plotbl]
            yx = self.msdata['data'][3,:,self.plotbl]
            ll = (xx + yy - 1j*(xy - yx))/2
            spec = np.abs(ll)  

        try:
            minval = float(self.minentry.text())*10000.
        except:
            self.statusBar.showMessage('Could not interpret minimum sfu value',2000)
            self.minentry.setText('{:10.2f}'.format(np.nanmin(spec)/10000.))
            minval = float(self.minentry.text())*10000.
            
        try:
            maxval = float(self.maxentry.text())*10000.
        except:
            self.statusBar.showMessage('Could not interpret maximum sfu value',2000)
            self.maxentry.setText('{:10.2f}'.format(np.nanmax(spec)/10000.))
            maxval = float(self.maxentry.text())*10000.
            
        for i in range(len(self.dspec_ax.images)):
            # Clear any existing images without changing the axis limits
            self.dspec_ax.images[i].remove()
            
        fghz = self.msdata['freq']/1e9
        times = self.msdata['time']
        nf, nt = spec.shape
        if not self.ignore_gaps.isChecked():
            # Insert nan rows and columns in gaps
            df = np.median(fghz[1:] - fghz[:-1])
            fbreak, = np.where(fghz[1:] - fghz[:-1] > 2*df)
            for n,fb in enumerate(fbreak):
                loc = fb + n + 1
                fghz = np.concatenate((fghz[:loc],np.array([fghz[loc-1]+df]),fghz[loc:]))
                spec = np.concatenate((spec[:loc],np.zeros((1,nt))+np.nan,spec[loc:]),0)
            nf, nt = spec.shape
            dt = np.median(times[1:] - times[:-1])
            tbreak, = np.where(times[1:] - times[:-1] > 2*dt)
            for n,tb in enumerate(tbreak):
                loc = tb + n + 1
                times = np.concatenate((times[:loc],np.array([times[loc-1]+dt]),times[loc:]))
                spec = np.concatenate((spec[:,:loc],np.zeros((nf,1))+np.nan,spec[:,loc:]),1)
            nf, nt = spec.shape

        pd = Time(times/86400.,format='mjd').plot_date
        x0, x1, y0, y1 = self.dspec_ax.axis()   # Save current plot extent
        self.dspec_ax.clear()
        im = self.dspec_ax.pcolormesh(pd,fghz,np.clip(spec, minval, maxval))
        self.dspec_ax.xaxis_date()
        self.dspec_ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        self.dspec_ax.set_ylabel('Frequency [GHz]')
        self.dspec_ax.set_xlabel('Time [UT]')
        if x0 == 0. and x1 == 1.:
            # Plot extent is 0,1 => cleared plot, so do not update plot extent
            pass
        else:
            # Restore saved plot extent
            self.dspec_ax.set_xlim(x0,x1)
            self.dspec_ax.set_ylim(y0,y1)
#        self.dspec_ax.imshow(np.clip(spec, minval, maxval), interpolation='nearest', aspect='auto', origin='lower')
        self.dspec_ax.set_title(self.ptitle)
        self.speccanvas.draw()
        
    def drawsel(self):
        ''' Uses timerange, spw, and stokes selection to draw a box on the plot.
        '''
        #if len(self.dspec_ax.get_images()) == 0:
        #    # No plot exists, so exit
        #    self.statusBar.showMessage("No plot exists yet. Please select a baseline to plot first.", 2000)
        #    return
        pd1, pd2 = trange2pd(self.trangeEdit.text())
        if pd1 is None:
            self.statusBar.showMessage(pd2, 2000)
            return
        f1, f2 = spw2frange(self.spwEdit.text(), self.msdata)
        if f1 is None:
            self.statusBar.showMessage(f2, 2000)
            return
        # Remove any existing lines before plotting a new one
        for line in self.dspec_ax.get_lines():
            line.remove()
        self.dspec_ax.plot([pd1, pd2, pd2, pd1, pd1],[f1, f1, f2, f2, f1],color='white')
        self.speccanvas.draw()
        
    def xferUL(self):
        ''' Uses the plot limits to define the timerange and spw selections.
        '''
        trange = pd2trange(self.dspec_ax.get_xlim())
        self.trangeEdit.setText(trange)
        f1, f2 = self.dspec_ax.get_ylim()
        idx1, = np.where(self.msdata['freq'] >= f1*1.e9)
        idx2, = np.where(self.msdata['freq'] <= f2*1.e9)
        spw = str(self.msdata['spwlist'][idx1[0]])+'~'+str(self.msdata['spwlist'][idx2[-1]])
        self.spwEdit.setText(spw)
    
    def xferDR(self):
        ''' Uses timerange, spw, and stokes selection to define the plot limits.
        '''
        pd = trange2pd(self.trangeEdit.text())
        if pd[0] is None:
            self.statusBar.showMessage(pd[1], 2000)
            return
        frange = spw2frange(self.spwEdit.text(), self.msdata)
        if frange[0] is None:
            self.statusBar.showMessage(frange[1], 2000)
            return
        self.dspec_ax.set_xlim(pd)
        self.dspec_ax.set_ylim(frange)            
        self.speccanvas.draw()
            
    def read_msdata(self, fillnan=None):
        ''' Read amplitudes for ALL data in the measurement set, returning
            a dictionary of the data, list of times, and list of frequencies.
            Fill flagged data with fillnan value if not None.
            
            Several changes to make this a LOT faster.  But nbl is hard-coded...
        '''
        ms = self.ms
        #if type(ms) != casatools.ms.ms:
        #    self.statusBar.showMessage('Not connected to a valid ms!',2000)
        #    return
            
        self.statusBar.addWidget(self.progressBar)
        self.progressBar.setFormat('Reading file...%p% done.')
        self.progressBar.setValue(0)
        self.progressBar.show()
        
#        ms.selectinit(datadescid=0, reset=True)
        spwinfo = ms.getspectralwindowinfo()
        nspw = len(spwinfo.keys())
#        spec = []
#        freq = []
        spwlist = []
#        for descid in range(len(spwinfo.keys())):
#            ms.selectinit(datadescid=0, reset=True)
#            ms.selectinit(datadescid=descid)
#            data = ms.getdata(['data', 'time', 'axis_info'], ifraxis=True)
#            spec_ = data['data']
#            freq_ = data['axis_info']['freq_axis']['chan_freq']
#            spwlist += [descid]*len(freq_)
#            if fillnan is not None:
#                flag_ = ms.getdata(['flag', 'time', 'axis_info'], ifraxis=True)['flag']
#                if type(fillnan) in [int, float]:
#                    spec_[flag_] = float(fillnan)
#                else:
#                    spec_[flag_] = 0.0
#            spec.append(spec_)
#            freq.append(freq_)
#            self.progressBar.setValue(100*(descid+1)/nspw)

        tb = table()
        tb.open(self.fname)
        spwtb = table()
        spwtb.open(self.fname+'/SPECTRAL_WINDOW')
        ptb = table()
        ptb.open(self.fname+'/POLARIZATION')
        
        mdata = ms.metadata()
        nspw = mdata.nspw()
        nbl = mdata.nbaselines() + mdata.nantennas()
        nscans = mdata.nscans()
        spw_nfrq = []     # List of number of frequencies in each spw
        for i in range(nspw):
            spw_nfrq.append(mdata.nchan(i))
        spw_nfrq = np.array(spw_nfrq)
        nf = np.sum(spw_nfrq)
        smry = mdata.summary()
        scan_ntimes = []  # List of number of times in each scan
        for iscan in range(nscans):
            scan_ntimes.append(len(smry['observationID=0']['arrayID=0']['scan='+str(iscan)]['fieldID=0'].keys()) - 6)
        scan_ntimes = np.array(scan_ntimes)
        nt = np.sum(scan_ntimes)
        times = tb.getcol('TIME')
        if times[nbl] - times[0] != 0:
            # This is frequency/scan sort order
            order = 'f'
        elif times[nbl*nspw] - times[0] !=0:
            # This is time sort order
            order = 't'
        npol = ptb.getcol('NUM_CORR',0,1)[0]
        ptb.close()
        freq = np.zeros(nf, float)
        times = np.zeros(nt, float)
        if order == 't':
            spec = np.zeros((npol, nf, nbl, nt), np.complex)
            for j in range(nt):
                fptr = 0
                # Loop over spw
                for i in range(nspw):
                    cfrq = spwtb.getcol('CHAN_FREQ',i,1)[:,0]     # Get channel frequencies for this spw (annoyingly comes out as shape (nf, 1)
                    if j == 0:
                        # Only need this the first time through
                        spwlist += [i]*len(cfrq)
                    if i == 0:
                        times[j] = tb.getcol('TIME',nbl*(i+nspw*j),1)  # Get the time
                    spec_ = tb.getcol('DATA',nbl*(i+nspw*j),nbl)  # Get complex data for this spw
                    flag = tb.getcol('FLAG',nbl*(i+nspw*j),nbl)   # Get flags for this spw
                    nfrq = len(cfrq)
                    # Apply flags
                    if type(fillnan) in [int, float]:
                        spec_[flag] = float(fillnan)
                    else:
                        spec_[flag] = 0.0
                    # Insert data for this spw into larger array
                    spec[:, fptr:fptr+nfrq, :, j] = spec_
                    freq[fptr:fptr+nfrq] = cfrq
                    fptr += nfrq
                self.progressBar.setValue(100*j/nt)
        else:
            specf = np.zeros((npol, nf, nt, nbl), np.complex)  # Array indexes are swapped
            iptr = 0
            for j in range(nscans):
                # Loop over scans
                for i in range(nspw):
                    #Loop over spectral windows
                    s = scan_ntimes[j]
                    f = spw_nfrq[i]
                    s1 = np.sum(scan_ntimes[:j])     # Start time index
                    s2 = np.sum(scan_ntimes[:j+1])   # End time index
                    f1 = np.sum(spw_nfrq[:i])         # Start freq index
                    f2 = np.sum(spw_nfrq[:i+1])       # End freq index
                    spec_ = tb.getcol('DATA',iptr,nbl*s)
                    flag  = tb.getcol('FLAG',iptr,nbl*s)
                    if j == 0:
                        cfrq = spwtb.getcol('CHAN_FREQ',i,1)[:,0]
                        freq[f1:f2] = cfrq
                        spwlist += [i]*len(cfrq)
                    times[s1:s2] = tb.getcol('TIME', iptr, nbl*s).reshape(s, nbl)[:,0]  # Get the times
                    iptr += nbl*s
                    # Apply flags
                    if type(fillnan) in [int, float]:
                        spec_[flag] = float(fillnan)
                    else:
                        spec_[flag] = 0.0
                    # Insert data for this spw into larger array
                    specf[:, f1:f2, s1:s2] = spec_.reshape(npol,f,s,nbl)
                    # Swap the array indexes back to the desired order
                    spec = np.swapaxes(specf,2,3)
                self.progressBar.setValue(100*j/nscans)
        tb.close()
        spwtb.close()

        self.statusBar.removeWidget(self.progressBar)
        self.dspec_ax.clear()  # Clear plot axis in preparation for a new dynamic spectrum
        self.msdata = {'data':spec, 'freq':freq, 'time':times, 'name':ms.name, 'spwlist':np.array(spwlist)}
                        
def get_bl2ord(nant=16):
    bl2ord = np.zeros((nant,nant),dtype='int')
    k = 0
    for i in range(nant):
        for j in range(i,nant):
            bl2ord[i,j] = k
            bl2ord[j,i] = k
            k += 1
    return bl2ord
    
def trange2pd(trange):
    ''' Convert a standard time range string to a pair of plot_date values.
        
        An error returns [None, msg]
    '''
    try:
        tstr1, tstr2 = trange.replace('/','-').split('~')
        return Time([tstr1[:10]+' '+tstr1[11:], tstr2[:10]+' '+tstr2[11:]]).plot_date
    except:
        return [None, 'Error interpreting timerange string']
        
def pd2trange(pd):
    ''' Convert a pair of plot_date values to a standard time trange string.
    
        An error returns None
    '''
    try:
        t1str = Time(pd[0],format='plot_date').iso.replace('-','/').replace(' ','/')
        t2str = Time(pd[1],format='plot_date').iso.replace('-','/').replace(' ','/')
        return t1str+'~'+t2str
    except:
        return None

def spw2frange(spw,msdata):
    ''' Convert a standard spectral window range to a pair of frequency values.
    '''
    try:
        spwrange = spw.split('~')
        if len(spwrange) == 1:
            spwrange *= 2
        spwlist = msdata['spwlist']
        fghz = msdata['freq']/1.e9
        idx1, = np.where(spwlist == int(spwrange[0]))
        idx2, = np.where(spwlist == int(spwrange[1]))
        return [fghz[idx1[0]], fghz[idx2[-1]]]
    except:
        return [None, 'Error interpreting spw string']

def pshift2pcenter(pshift,ms):
    ''' Convert phaseshift in heliocentric X,Y to RA, Dec using phase center in ms file

        pshift is a 2-element list of floats, e.g. [0.0, 0.0]

        output is the formatted RA Dec of the phase center, for use as phasecenter input in TCLEAN.
    '''
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from sunpy.coordinates import sun
    asec2rad = np.pi/(3600.*180.)
    
    mdata = ms.metadata()
    pcenter = mdata.phasecenter()
    # Phase center in radians
    ra0, dec0 = (pcenter['m0']['value'], pcenter['m1']['value'])
    trange = mdata.timerangeforobs(0)
    tstart = Time(trange['begin']['m0']['value'],format='mjd')  # Start time of file
    p0rad = sun.P(tstart).rad
    # Convert shifts in heliocentric to geocentric (i.e. rotate for solar P-angle)
    xgeo = (pshift[0]*np.cos(p0rad) - pshift[1]*np.sin(p0rad))*asec2rad
    ygeo = (pshift[0]*np.sin(p0rad) + pshift[1]*np.cos(p0rad))*asec2rad
    ra = ra0 - xgeo/np.cos(dec0)
    dec = dec0 + ygeo
    c = SkyCoord(ra=ra*u.rad, dec=dec*u.rad)
    rax = c.ra.hms
    decx = c.dec.dms
    phasecenter='J2000 {:02d}:{:02d}:{:07.4f} {:3d}.{:02d}.{:06.3f}'.format(int(rax.h),int(rax.m),rax.s,
                                                                          int(decx.d),int(abs(decx.m)),abs(decx.s))
    return phasecenter
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
