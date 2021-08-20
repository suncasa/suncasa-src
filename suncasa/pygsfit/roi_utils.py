import pyqtgraph as pg
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os
import pickle
from regions import Regions, CRTFRegionParserError
from PyQt5.QtGui import QDesktopServices

class roi_dialog(object):
    def __init__(self, img_size, cfreqs):
        self.world = False
        self.img_size = img_size
        self.cfreqs = cfreqs
        self.display_list = []
        # self.input_example = '[[{0}, {1}],[{2}, {3}]],[1.0, 18.0]'.format(int(img_size[0] * 0.3),
        #                                                                  int(img_size[1] * 0.3),
        #                                                                  int(img_size[0] * 0.7),
        #                                                                  int(img_size[1] * 0.7))
        PixSzeLt = [int(self.img_size[0] * 0.1), int(self.img_size[1] * 0.1), int(self.img_size[0] * 0.3),
                    int(self.img_size[1] * 0.3), int(self.img_size[0] * 0.5), int(self.img_size[1] * 0.5)]
        self.input_examples = []
        self.input_examples.append(
            'box[[{0}pix, {1}pix], [{2}pix, {3}pix]], range=[1.0GHz, 18.0GHz]'.format(PixSzeLt[2], PixSzeLt[3],
                                                                                      PixSzeLt[0], PixSzeLt[1]))
        self.input_examples.append(
            'centerbox[[{0}pix, {1}pix], [{2}pix, {3}pix]], range=[1.0GHz, 18.0GHz]'.format(PixSzeLt[4], PixSzeLt[5],
                                                                                            PixSzeLt[0], PixSzeLt[1]))
        self.input_examples.append(
            'rotbox[[{0}pix, {1}pix], [{2}pix, {3}pix], 15deg], range=[1.0GHz, 18.0GHz]'.format(PixSzeLt[2],
                                                                                                PixSzeLt[3],
                                                                                                PixSzeLt[0],
                                                                                                PixSzeLt[1]))
        self.input_examples.append(
            'ellipse[[{0}pix, {1}pix], [{2}pix, {3}pix], 15deg], range=[1.0GHz, 18.0GHz]'.format(PixSzeLt[4],
                                                                                                 PixSzeLt[5],
                                                                                                 PixSzeLt[0],
                                                                                                 PixSzeLt[0] * 0.75))
        self.input_examples.append(
            'circle[[{0}pix, {1}pix], {2}pix], range=[1.0GHz, 18.0GHz]'.format(PixSzeLt[4], PixSzeLt[5], PixSzeLt[0]))
        self.input_example = self.input_examples[0]

    def setupUi(self, Dialog):
        Dialog.setObjectName("Manually Defined ROI(s)")
        Dialog.resize(600, 302)
        self.ok_cancel = QDialogButtonBox(Dialog)
        self.ok_cancel.setGeometry(QRect(500, 40, 81, 71))
        self.ok_cancel.setOrientation(Qt.Vertical)
        self.ok_cancel.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        self.ok_cancel.setObjectName("ok_cancel")
        self.rois_list_view = QListView(Dialog)
        self.rois_list_view.setGeometry(QRect(10, 10, 481, 200))
        self.rois_list_view.setObjectName("rois_list_view")
        self.slm = QStringListModel()
        self.slm.setStringList(self.display_list)
        self.selected_index = [cur_item.row() for cur_item in self.rois_list_view.selectedIndexes()]

        # listView.clicked.connect(self.checkItem)
        self.rois_list_view.setModel(self.slm)
        self.lineEdit = QLineEdit(Dialog)
        self.lineEdit.setGeometry(QRect(10, 220, 561, 21))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit.setText(self.input_example)
        # self.label = QLabel(Dialog)
        # self.label.setGeometry(QRect(10, 160, 381, 21))
        # self.label.setObjectName("label")
        self.pushButton_url = QPushButton(Dialog)
        self.pushButton_url.setGeometry(QRect(500, 180, 81, 32))
        self.pushButton_url.setObjectName("pushButton_url")
        self.pushButton_url.clicked.connect(lambda: QDesktopServices.openUrl(
            QUrl("https://casaguides.nrao.edu/index.php/CASA_Region_Format#Global_definitions")))
        self.pushButton_delete = QPushButton(Dialog)
        self.pushButton_delete.setGeometry(QRect(500, 140, 81, 32))
        self.pushButton_delete.setObjectName("pushButton_delete")
        self.pushButton_delete.clicked.connect(self.delete_from_list)
        # self.checkBox = QCheckBox(Dialog)
        # self.checkBox.setGeometry(QRect(10, 250, 87, 20))
        # self.checkBox.setObjectName("checkBox")
        self.pushButton_add_to_list = QPushButton(Dialog)
        self.pushButton_add_to_list.setGeometry(QRect(470, 250, 113, 32))
        self.pushButton_add_to_list.setObjectName("pushButton")
        self.pushButton_add_to_list.clicked.connect(self.add_to_list)
        self.pushButton_img_flx_preset = QPushButton(Dialog)
        self.pushButton_img_flx_preset.setGeometry(QRect(10, 250, 120, 32))
        self.pushButton_img_flx_preset.setObjectName("img_flx_roi")
        self.pushButton_img_flx_preset.clicked.connect(self.rois_for_cal_flux)
        self.pushButton_load_file = QPushButton(Dialog)
        self.pushButton_load_file.setGeometry(QRect(200, 250, 140, 32))
        self.pushButton_load_file.setObjectName("load_roi_file")
        self.pushButton_load_file.clicked.connect(self.roi_file_select)
        self.shape_comboBox = QComboBox(Dialog)
        self.shape_comboBox.setGeometry(QRect(350, 250, 104, 26))
        self.shape_comboBox.setObjectName("comboBox")
        self.shape_comboBox.addItems([''] * 7)
        # self.label_2 = QLabel(Dialog)
        # self.label_2.setGeometry(QRect(10, 190, 221, 16))
        # self.label_2.setObjectName("label_2")

        self.retranslateUi(Dialog)
        self.ok_cancel.accepted.connect(Dialog.accept)
        self.ok_cancel.rejected.connect(Dialog.reject)
        QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Manually Define ROI(s)"))
        # self.label.setText(
        #    _translate("Dialog", "FOV in pixel/arcsec, freq_boundry in GHz"))
        # self.checkBox.setText(_translate("Dialog", "World"))
        self.pushButton_url.setText(_translate("Dialog", "Doc"))
        self.pushButton_img_flx_preset.setText(_translate("Dialog", "EOVSA Presets"))
        self.pushButton_load_file.setText(_translate("Dialog", "Load ROI File"))
        self.pushButton_add_to_list.setText(_translate("Dialog", "Add to List"))
        self.pushButton_delete.setText(_translate("Dialog", "Delete"))
        self.shape_comboBox.setItemText(0, _translate("Dialog", "box"))
        self.shape_comboBox.setItemText(1, _translate("Dialog", "centerbox"))
        self.shape_comboBox.setItemText(2, _translate("Dialog", "rotbox"))
        self.shape_comboBox.setItemText(3, _translate("Dialog", "ellipse"))
        self.shape_comboBox.setItemText(4, _translate("Dialog", "circle"))
        self.shape_comboBox.setItemText(5, _translate("Dialog", "poly"))
        self.shape_comboBox.setItemText(6, _translate("Dialog", "annulus"))
        self.shape_comboBox.currentIndexChanged.connect(self.change_example)
        # self.ok_cancel.accepted.connect(Dialog.accept)
        # self.ok_cancel.rejected.connect(Dialog.reject)
        # self.label_2.setText(_translate("Dialog", "e.g. [[x1,y1],[x2,y2]],[freq_lo,freq_hi]"))

    def add_to_list(self):
        self.input_string = self.lineEdit.text()
        self.cur_info_string = self.input_string
        try:
            tmp_crtf_region = Regions.parse(self.cur_info_string, format='crtf')[0]
            self.display_list.append(self.cur_info_string)
            self.slm.setStringList(self.display_list)
            self.rois_list_view.setModel(self.slm)
        except:
            msg_box = QMessageBox(QMessageBox.Warning, 'Invalid Input!', 'The input is not a vliad CRTF string!')
            msg_box.exec_()

    def delete_from_list(self):
        self.selected_sting = [cur_item.row() for cur_item in self.rois_list_view.selectedIndexes()]
        del self.display_list[self.selected_sting[0]]
        self.slm.setStringList(self.display_list)
        self.rois_list_view.setModel(self.slm)

    def change_example(self):
        self.input_example = self.input_examples[self.shape_comboBox.currentIndex()]
        self.lineEdit.setText(self.input_example)

    # @staticmethod
    def getResult(self, Dialog):
        # dialog = roi_dialog(parent)
        # cur_result = dialog.exec_()
        return self.display_list

    def rois_for_cal_flux(self):
        selections = [str(x) for x in np.arange(2, 6).tolist()]
        seleted_item, okPressed = QInputDialog.getItem(None, "Number of ROIs", "Number of ROIs", selections, 0, False)
        num_of_rois = int(seleted_item)
        print(self.cfreqs)
        if okPressed:
            print(num_of_rois, ' ROIs are created')
            freq_boundry = np.linspace(0, len(self.cfreqs) - 1, num_of_rois + 1, dtype=int)
            size_list = [int(max(100.0 * self.cfreqs[0] / self.cfreqs[freq_ind], 5.)) for freq_ind in
                         freq_boundry[:-1]]
            crtf_str_list = []
            for nri in range(num_of_rois):
                crtf_str_list.append(
                    'centerbox[[{}pix, {}pix], [{}pix, {}pix]], range=[{}GHz, {}GHz]'.format(int(self.img_size[0] / 2),
                                                                                             int(self.img_size[1] / 2),
                                                                                             size_list[nri],
                                                                                             size_list[nri],
                                                                                             self.cfreqs[freq_boundry[nri]],
                                                                                             self.cfreqs[freq_boundry[nri + 1]]))
            for cur_str in crtf_str_list:
                self.display_list.append(cur_str)
            self.slm.setStringList(self.display_list)
            self.rois_list_view.setModel(self.slm)
        else:
            return
    def roi_file_select(self, Dialog):
        cur_file_name, _file_filter = QFileDialog.getOpenFileName(None, 'Select Region save',
                                                                     './', 'Region save (*.crtf *.reg *.ds9 *.fits)')
        #self.fname = 'EOVSA_20210507T190205.000000.outim.image.allbd.fits'
        cur_format = cur_file_name.split('.')[-1]
        if cur_format == 'reg':
            cur_format = 'ds9'
        try:
            cur_region = Regions.read(cur_file_name, format=cur_file_name.split('.')[-1])
        except:
            msg_box = QMessageBox(QMessageBox.Warning, 'Invalid Input!', 'The input can not be converted!')
            msg_box.exec_()
        self.display_list.append(cur_region.serialize(format='crtf'))
        self.slm.setStringList(self.display_list)
        self.rois_list_view.setModel(self.slm)
        return



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

    print(str(type(crtf_region)))
    if 'CirclePixel' in str(type(crtf_region)):
        pg_roi = pg.CircleROI(pos=np.asarray(crtf_region.center.xy) - np.ones(2) * crtf_region.radius,
                              radius=crtf_region.radius, pen=pen_arg)
    elif 'RectanglePixel' in str(type(crtf_region)):
        pg_roi = pg.RectROI(
            pos=get_corner_rotate(center_xy=crtf_region.center.xy, width=crtf_region.width, height=crtf_region.height,
                                  rot_angle=crtf_region.angle.to_value(unit='radian')),
            size=[crtf_region.width, crtf_region.height], angle=crtf_region.angle.to_value(unit='degree'), pen=pen_arg)
    elif 'EllipsePixel' in str(type(crtf_region)):
        pg_roi = pg.EllipseROI(pos=get_corner_rotate(center_xy=crtf_region.center.xy, width=crtf_region.width,
                                                     height=crtf_region.height,
                                                     rot_angle=crtf_region.angle.to_value(unit='radian')),
                               size=[crtf_region.width, crtf_region.height],
                               angle=crtf_region.angle.to_value(unit='degree'), pen=pen_arg)
    # todo: add line ROI
    else:
        raise NameError('Sorry, only rectangle, circle, and ellipse are supported for this moment.')
    freq_range = [1.0, 18.0]
    if 'range' in crtf_region.meta:
        if str(crtf_region.meta['range'][0].unit) == 'GHz':
            freq_range = [crtf_region.meta['range'][0].value, crtf_region.meta['range'][1].value]
            print('Set freq range to {}'.format(freq_range))
    return (pg_roi, freq_range)

def rois_for_cal_flux(self1):
    selections = [str(x) for x in np.arange(2, 6).tolist()]
    seleted_item, okPressed = QInputDialog.getItem(None, "Number of ROIs", "Number of ROIs", selections, 0, False)
    num_of_rois = int(seleted_item)
    if okPressed:
        print(num_of_rois, ' ROIs are created')
        add_new_group(self1)
        freq_boundry = np.linspace(0, len(self1.cfreqs) - 1, num_of_rois + 1, dtype=int)
        size_list = [int(max(100.0 * self1.cfreqs[0] / self1.cfreqs[freq_ind], 5.)) for freq_ind in freq_boundry[:-1]]
        crtf_str_list = []
        for nri in range(num_of_rois):
            crtf_str_list.append(
                'centerbox[[{}pix, {}pix], [{}pix, {}pix]], range=[{}GHz, {}GHz]'.format(int(self1.meta['nx'] / 2),
                                                                                         int(self1.meta['ny'] / 2),
                                                                                         size_list[nri], size_list[nri],
                                                                                         freq_boundry[nri],
                                                                                         freq_boundry[nri + 1]))
        add_md_rois(self1,crtf_str_list)
    else:
        print('No ROI is defined.')


def save_roi_group(self1):
    data_saved = False
    selections = []
    for group_index, cur_group in enumerate(self1.rois):
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
                pickle.dump([cur_roi.saveState() for cur_roi in self1.rois[seleted_group_index]], handle,
                            protocol=pickle.HIGHEST_PROTOCOL)
                print('The selected group has been saved to ', fileName)
                data_saved = True
    if not data_saved:
        print('No ROI is saved.')


def add_md_rois(self1, inp_str_list):
    add_new_group(self1)
    for si, cur_str in enumerate(inp_str_list):
        add_customized_roi(self1=self1, crtf_str=cur_str)


def add_customized_roi(self1, crtf_str):
    """Add a ROI region to the selection"""
    pg_roi_obj, freq_range = crtf_to_pgroi(crtf_str=crtf_str, eo_wcs=self1.eo_wcs,
                                           pen_arg=(len(self1.rois[self1.roi_group_idx]), 9))
    self1.new_roi = pg_roi_obj
    self1.pg_img_canvas.addItem(self1.new_roi)
    self1.new_roi.freq_mask = np.ones_like(self1.cfreqs) * False
    self1.new_roi.sigRegionChanged.connect(self1.calc_roi_spec)
    # choose which group to add
    self1.add_to_roigroup_selection()
    self1.rois[self1.roi_group_idx].append(self1.new_roi)
    self1.nroi_current_group = len(self1.rois[self1.roi_group_idx])
    self1.roi_selection_widget.clear()
    self1.roi_selection_widget.addItems([str(i) for i in range(self1.nroi_current_group)])
    self1.current_roi_idx = self1.nroi_current_group - 1
    self1.has_rois = True
    self1.calc_roi_spec()
    self1.roi_freq_lowbound_selector.setValue(freq_range[0])
    self1.roi_freq_hibound_selector.setValue(freq_range[1])


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