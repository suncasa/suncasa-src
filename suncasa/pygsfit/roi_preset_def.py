import pyqtgraph as pg
import numpy as np
from PyQt5.QtWidgets import QInputDialog, QFileDialog
import os
import pickle


def add_preset_roi_selection(self, preset=None):
    if preset == 'Img_Inte_ROIs':
        #def get_num_of_rois():
        selections = [str(x) for x in np.arange(2, 11).tolist()]
        seleted_item, okPressed = QInputDialog.getItem(None, "Number of ROIs", "Number of ROIs", selections, 0, False)
        num_of_rois = int(seleted_item)
        if okPressed:
            print(num_of_rois , ' ROIs are created')
            add_new_group(self)
            freq_boundry = np.linspace(0, len(self.cfreqs)-1, num_of_rois + 1, dtype=int)
            size_list = [int(max(100.0 * self.cfreqs[0] / self.cfreqs[freq_ind], 5.)) for freq_ind in freq_boundry[:-1]]
            for nri in range(num_of_rois):
                add_customized_roi(self, start_point_shift=np.ones((2),dtype=int) * int(size_list[nri] / 2),
                                   roi_size=np.ones((2),dtype=int) * size_list[nri],
                                   freq_index_range=[freq_boundry[nri], freq_boundry[nri + 1]])
        else:
            preset = 'Presets'
    elif preset == 'Save Group':
        selections = []
        for group_index, cur_group in enumerate(self.rois):
            selections.append('Group {}, {} ROI(s)'.format(group_index,len(cur_group)))
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
                    pickle.dump([cur_roi.saveState() for cur_roi in self.rois[seleted_group_index]], handle, protocol=pickle.HIGHEST_PROTOCOL)
                    print('The selected group has been saved to ', fileName)
            else:
                preset = 'Presets'
        else:
            preset = 'Presets'
    if preset == 'Presets':
        print('No preset is selected')
        self.roi_selection_presets_widget.setCurrentIndex(0)


def add_customized_roi(self, start_point_shift, roi_size, freq_index_range):
    """Add a ROI region to the selection"""
    self.new_roi = pg.RectROI([self.meta['nx'] / 2 - start_point_shift[0], self.meta['nx'] / 2 - start_point_shift[1]],
                              roi_size, pen=(len(self.rois[self.roi_group_idx]), 9))
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
    self.roi_freq_lowbound_selector.setValue(self.cfreqs[freq_index_range[0]])
    self.roi_freq_hibound_selector.setValue(self.cfreqs[freq_index_range[1]])


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
