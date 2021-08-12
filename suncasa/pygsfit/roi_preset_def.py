import pyqtgraph as pg
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os
import pickle
from regions import Regions
from PyQt5.QtGui import QDesktopServices


def add_preset_roi_selection(self, preset=None):
    if preset == 'Img_Inte_ROIs':
        # def get_num_of_rois():
        selections = [str(x) for x in np.arange(2, 11).tolist()]
        seleted_item, okPressed = QInputDialog.getItem(None, "Number of ROIs", "Number of ROIs", selections, 0, False)
        num_of_rois = int(seleted_item)
        if okPressed:
            print(num_of_rois, ' ROIs are created')
            add_new_group(self)
            freq_boundry = np.linspace(0, len(self.cfreqs) - 1, num_of_rois + 1, dtype=int)
            size_list = [int(max(100.0 * self.cfreqs[0] / self.cfreqs[freq_ind], 5.)) for freq_ind in freq_boundry[:-1]]
            for nri in range(num_of_rois):
                add_customized_roi(self, start_point_shift=np.ones((2), dtype=int) * int(size_list[nri] / 2),
                                   roi_size=np.ones((2), dtype=int) * size_list[nri],
                                   freq_index_range=[freq_boundry[nri], freq_boundry[nri + 1]])
        else:
            preset = 'Presets'
    elif preset == 'Save Group':
        selections = []
        for group_index, cur_group in enumerate(self.rois):
            selections.append('Group {}, {} ROI(s)'.format(group_index, len(cur_group)))
        seleted_item, okPressed = QInputDialog.getItem(None, "Which Group?", "List of Group", selections, 0, False)
        if okPressed:
            print(seleted_item, ' is selected')
            seleted_group_index = int(seleted_item.split(' ')[1][:-1])
            fileName, ok2 = QFileDialog.getSaveFileName(None,
                                                        "Save the selected group ()",
                                                        os.getcwd(),
                                                        "All Files (*);;pickle save (*.p)")
            if ok2:
                with open(fileName, 'wb') as handle:
                    pickle.dump([cur_roi.saveState() for cur_roi in self.rois[seleted_group_index]], handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
                    print('The selected group has been saved to ', fileName)
            else:
                preset = 'Presets'
        else:
            preset = 'Presets'
    if preset == 'Presets':
        print('No preset is selected')
        # self.roi_selection_presets_widget.setCurrentIndex(0)


def add_customized_roi(self, crtf_str):
    """Add a ROI region to the selection"""
    pg_roi_obj, freq_range = crtf_to_pgroi(crtf_str=crtf_str, eo_wcs=self.eo_wcs, pen_arg=(len(self.rois[self.roi_group_idx]), 9))
    self.new_roi = pg_roi_obj
    self.pg_img_canvas.addItem(self.new_roi)
    self.new_roi.freq_mask = np.ones_like(self.cfreqs) * False
    self.new_roi.sigRegionChanged.connect(self.calc_roi_spec)
    # choose which group to add
    self.add_to_roigroup_selection()
    self.rois[self.roi_group_idx].append(self.new_roi)
    self.nroi_current_group = len(self.rois[self.roi_group_idx])
    self.roi_selection_widget.clear()
    self.roi_selection_widget.addItems([str(i) for i in range(self.nroi_current_group)])
    self.current_roi_idx = self.nroi_current_group - 1
    self.has_rois = True
    self.calc_roi_spec()
    self.roi_freq_lowbound_selector.setValue(freq_range[0])
    self.roi_freq_hibound_selector.setValue(freq_range[1])


def add_new_group(self):
    if self.has_rois:
        self.rois.append([])
        self.roi_group_idx += 1
        if len(self.rois) > self.add_to_roigroup_widget.count():
            self.add_to_roigroup_widget.addItem(str(self.roi_group_idx))
            self.roigroup_selection_widget.addItem(str(self.roi_group_idx))
    else:
        self.roi_group_idx = 0
    self.add_to_roigroup_widget.setCurrentRow(self.roi_group_idx)
    self.add_to_roigroup_button.setText(self.add_to_roigroup_widget.currentItem().text())
    self.roigroup_selection_widget.setCurrentRow(self.roi_group_idx)
    self.roigroup_selection_button.setText(self.roigroup_selection_widget.currentItem().text())


class roi_dialog(object):
    def __init__(self, img_size, world_converter=None):
        self.world = False
        self.img_size = img_size
        self.display_list = []
        self.input_example = '[[{0}, {1}],[{2}, {3}]],[1.0, 18.0]'.format(int(img_size[0] * 0.3),
                                                                          int(img_size[1] * 0.3),
                                                                          int(img_size[0] * 0.7),
                                                                          int(img_size[1] * 0.7))

    def setupUi(self, Dialog):
        Dialog.setObjectName("Manually Define ROI(s)")
        Dialog.resize(400, 302)
        self.ok_cancel = QDialogButtonBox(Dialog)
        self.ok_cancel.setGeometry(QRect(300, 40, 81, 71))
        self.ok_cancel.setOrientation(Qt.Vertical)
        self.ok_cancel.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        self.ok_cancel.setObjectName("ok_cancel")
        self.rois_list_view = QListView(Dialog)
        self.rois_list_view.setGeometry(QRect(10, 10, 281, 141))
        self.rois_list_view.setObjectName("rois_list_view")
        self.slm = QStringListModel()
        self.slm.setStringList(self.display_list)
        # listView.clicked.connect(self.checkItem)
        self.rois_list_view.setModel(self.slm)
        self.lineEdit = QLineEdit(Dialog)
        self.lineEdit.setGeometry(QRect(10, 220, 361, 21))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit.setText(self.input_example)
        #self.label = QLabel(Dialog)
        #self.label.setGeometry(QRect(10, 160, 381, 21))
        #self.label.setObjectName("label")
        self.pushButton_url = QPushButton(Dialog)
        self.pushButton_url.setGeometry(QRect(10, 160, 113, 32))
        self.pushButton_url.setObjectName("pushButton_url")
        self.pushButton_url.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://casaguides.nrao.edu/index.php/CASA_Region_Format#Global_definitions")))
        #self.checkBox = QCheckBox(Dialog)
        #self.checkBox.setGeometry(QRect(10, 250, 87, 20))
        #self.checkBox.setObjectName("checkBox")
        self.pushButton = QPushButton(Dialog)
        self.pushButton.setGeometry(QRect(270, 250, 113, 32))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.save_to_list)
        self.comboBox = QComboBox(Dialog)
        self.comboBox.setGeometry(QRect(130, 250, 104, 26))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.label_2 = QLabel(Dialog)
        self.label_2.setGeometry(QRect(10, 190, 221, 16))
        self.label_2.setObjectName("label_2")

        self.retranslateUi(Dialog)
        self.ok_cancel.accepted.connect(Dialog.accept)
        self.ok_cancel.rejected.connect(Dialog.reject)
        QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Manually Define ROI(s)"))
        #self.label.setText(
        #    _translate("Dialog", "FOV in pixel/arcsec, freq_boundry in GHz"))
        #self.checkBox.setText(_translate("Dialog", "World"))
        self.pushButton_url.setText(_translate("Dialog", "CASA Region Format"))
        self.pushButton.setText(_translate("Dialog", "Add to List"))
        self.comboBox.setItemText(0, _translate("Dialog", "box"))
        self.comboBox.setItemText(1, _translate("Dialog", "centerbox"))
        self.comboBox.setItemText(2, _translate("Dialog", "rotbox"))
        self.comboBox.setItemText(3, _translate("Dialog", "ellipse"))
        self.comboBox.setItemText(4, _translate("Dialog", "circle"))
        self.comboBox.setItemText(5, _translate("Dialog", "poly"))
        self.comboBox.setItemText(6, _translate("Dialog", "annulus"))
        self.label_2.setText(_translate("Dialog", "e.g. [[x1,y1],[x2,y2]],[freq_lo,freq_hi]"))

    def save_to_list(self):
        self.input_string = self.lineEdit.text()
        self.cur_info_string = 'FOV, Freq_range, ' + self.input_string
        self.display_list.append(self.cur_info_string)
        self.slm.setStringList(self.display_list)
        self.rois_list_view.setModel(self.slm)


'''
class pygsfit_region(Regions):
    def __init__(self, cur_map_):
        super().__init__()
        self.crtf_str=''
    def crtf_2_pyqtroi(self, crtf_str, eo_wcs):
        try:
            crtf_region = self.parse(crtf_str, format='crtf')[0]
        except:
            print('check your input string format: https://casaguides.nrao.edu/index.php/CASA_Region_Format#Global_definitions')
        if 'Sky' in str(type(crtf_region)):
            crtf_region = crtf_region.to_pixel(eo_wcs)
        if 
'''


def crtf_to_pgroi(crtf_str, eo_wcs, pen_arg):
    try:
        crtf_region = Regions.parse(crtf_str, format='crtf')[0]
    except:
        print(
            'check your input string format: https://casaguides.nrao.edu/index.php/CASA_Region_Format#Global_definitions')
    if 'Sky' in str(type(crtf_region)):
        crtf_region = crtf_region.to_pixel(eo_wcs)

    def get_corner_rotate(center_xy, width, height, rot_angle):
        distance = np.sqrt(width ** 2 + height ** 2)
        cur_theta = np.pi / 2. - np.arctan(width / height) - rot_angle
        cor_x = center_xy[0] - distance * np.sin(cur_theta)
        cor_y = center_xy[1] - distance * np.cos(cur_theta)
        return [cor_x, cor_y]

    if 'CirclePixel' in str(type(crtf_region)):
        pg_roi = pg.CircleROI(pos=np.asarray(crtf_region.center.xy) - np.ones(2) * crtf_region.radius, radius=crtf_region.radius, pen = pen_arg)
    elif 'RectPixel' in str(type(crtf_region)):
        pg_roi = pg.RectROI(
            pos=get_corner_rotate(center_xy=crtf_region.center.xy, width=crtf_region.width, height=crtf_region.height,
                                  rot_angle=crtf_region.angle.to_value(unit='radian')),
            size=[crtf_region.width, crtf_region.height], angle=crtf_region.angle.to_value(unit='degree'), pen = pen_arg)
    elif 'EllipsePixel' in str(type(crtf_region)):
        pg_roi = pg.EllipseROI(pos=get_corner_rotate(center_xy=crtf_region.center.xy, width=crtf_region.width,
                                                height=crtf_region.height,
                                                rot_angle=crtf_region.angle.to_value(unit='radian')),
                                                size=[crtf_region.width, crtf_region.height],
                                                angle=crtf_region.angle.to_value(unit='degree'), pen = pen_arg)
    #todo: add line ROI
    else:
        raise NameError('Sorry, only rectangle, circle, and ellipse are supported for this moment.')
    freq_range = [1.0,18.0]
    if 'range' in crtf_region.meta:
        if str(crtf_region.meta['range'][0].unit) == 'GHz':
            freq_range = [crtf_region.meta['range'][0].value, crtf_region.meta['range'][1].value]
            print('Set freq range to {}'.format(freq_range))
    return (pg_roi, freq_range)


