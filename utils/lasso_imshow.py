#!/usr/bin/env python

"""
Function for interactively selecting part of an array displayed as an image with matplotlib.
"""

import matplotlib.pyplot as plt
from matplotlib import is_interactive
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector, RectangleSelector
import numpy as np


def path_bbox(p):
    """
    Return rectangular bounding box of given path.
    Parameters
    ----------
    p : array_like
        Array of vertices with shape Nx2.
    Returns
    -------
    bbox : array_like
        Array of bounding box vertices with shape 4x2.
    """

    assert p.ndim == 2
    assert p.shape[1] == 2

    ix_min = p[:, 0].argmin()
    ix_max = p[:, 0].argmax()
    iy_min = p[:, 1].argmin()
    iy_max = p[:, 1].argmax()

    return np.array([[p[ix_min, 0], p[iy_min, 1]], [p[ix_min, 0], p[iy_max, 1]], [p[ix_max, 0], p[iy_max, 1]], [p[ix_max, 0], p[iy_min, 1]]])


def imshow_select(data, selected=None, selector='lasso', bbox=False, sav_img=False, contour_color='w', plotonly=False):
    """
    Display array as image with region selector.

    Parameters
    ----------
    data : array_like
        Array to display.
    selector : str
        Region selector. For `lasso`, use `LassoSelector`; for `rectangle`,
        use `RectangleSelector`.
    bbox : bool
        If True, only return array within rectangular bounding box of selected region.
        Otherwise, return array with same dimensions as `data` such that selected region
        contains the corresponding values from `data` and the remainder contains 0.
    Returns
    -------
    region : array_like
        Data for selected region.
    mask : array_like
        Boolean mask with same shape of `data` for selecting the returned region from `data`.
    """

    interactive = is_interactive()
    if not interactive:
        plt.ion()
    fig, axs = plt.subplots(ncols=2, nrows=1, sharex=True, sharey=True, subplot_kw={'adjustable': 'box-forced'})
    axs = axs.ravel()
    axs[0].imshow(data, origin='lower')
    axs[1].imshow(np.zeros_like(data), origin='lower')
    plt.subplots_adjust()

    x, y = np.meshgrid(np.arange(data.shape[1], dtype=int), np.arange(data.shape[0], dtype=int))
    pix = np.vstack((x.flatten(), y.flatten())).T

    # Store data in dict value to permit overwriting by nested
    # functions in Python 2.7:
    collects = {}
    if selected is None:
        selected = {}
    else:
        if 'mask' in selected.keys():
            collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
            if 'data' not in selected.keys():
                selected['data'] = np.zeros_like(data)
                selected['data'][selected['mask']] = data[selected['mask']]
            axs[1].imshow(selected['data'], origin='lower')
            fig.canvas.draw_idle()
        else:
            selected = {}

    def _onselect_lasso(verts):
        verts = np.array(verts)
        p = Path(verts)
        ind = p.contains_points(pix, radius=0.5)
        if 'data' in selected.keys() and (not bbox):
            selected['data'].flat[ind] = data.flat[ind]
            selected['mask'].flat[ind] = True
            if 'contour_ax_orig' in collects.keys():
                for coll in collects['contour_ax_orig'].collections:
                    axs[0].collections.remove(coll)
                collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
        else:
            selected['data'] = np.zeros_like(data)
            selected['mask'] = np.tile(False, data.shape)
            selected['data'].flat[ind] = data.flat[ind]
            selected['mask'].flat[ind] = True
            collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
        axs[1].imshow(selected['data'], origin='lower')
        fig.canvas.draw_idle()
        if bbox:
            b = path_bbox(verts)
            selected['data'] = selected['data'][int(min(b[:, 1])):int(max(b[:, 1])), int(min(b[:, 0])):int(max(b[:, 0]))]

    def _deselect_lasso(verts):
        verts = np.array(verts)
        p = Path(verts)
        ind = p.contains_points(pix, radius=0.5)
        if 'data' in selected.keys():
            selected['data'].flat[ind] = 0.
            selected['mask'].flat[ind] = False
            if 'contour_ax_orig' in collects.keys():
                for coll in collects['contour_ax_orig'].collections:
                    axs[0].collections.remove(coll)
                collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
            axs[1].imshow(selected['data'], origin='lower')
            fig.canvas.draw_idle()

    def _onselect_rectangle(start, end):
        verts = np.array([[start.xdata, start.ydata], [start.xdata, end.ydata], [end.xdata, end.ydata], [end.xdata, start.ydata]], int)
        p = Path(verts)
        ind = p.contains_points(pix, radius=0.5)
        if 'data' in selected.keys() and (not bbox):
            selected['data'].flat[ind] = data.flat[ind]
            selected['mask'].flat[ind] = True
            if 'contour_ax_orig' in collects.keys():
                for coll in collects['contour_ax_orig'].collections:
                    axs[0].collections.remove(coll)
                collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
        else:
            selected['data'] = np.zeros_like(data)
            selected['mask'] = np.tile(False, data.shape)
            selected['data'].flat[ind] = data.flat[ind]
            selected['mask'].flat[ind] = True
            collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
        axs[1].imshow(selected['data'], origin='lower')
        fig.canvas.draw_idle()
        if bbox:
            b = path_bbox(verts)
            selected['data'] = selected['data'][min(b[:, 1]):max(b[:, 1]), min(b[:, 0]):max(b[:, 0])]

    def _deselect_rectangle(start, end):
        verts = np.array([[start.xdata, start.ydata], [start.xdata, end.ydata], [end.xdata, end.ydata], [end.xdata, start.ydata]], int)
        p = Path(verts)
        ind = p.contains_points(pix, radius=0.5)
        if 'data' in selected.keys():
            selected['data'].flat[ind] = 0.
            selected['mask'].flat[ind] = False
            if 'contour_ax_orig' in collects.keys():
                for coll in collects['contour_ax_orig'].collections:
                    axs[0].collections.remove(coll)
                collects['contour_ax_orig'] = axs[0].contour(selected['mask'].astype(np.int), levels=[0.5], linewidths=0.75, colors=contour_color)
            axs[1].imshow(selected['data'], origin='lower')
            fig.canvas.draw_idle()

    name_to_selector = {'lasso': LassoSelector, 'rectangle': RectangleSelector}
    selector = name_to_selector[selector]
    onselect_dict = {LassoSelector: _onselect_lasso, RectangleSelector: _onselect_rectangle}
    deselect_dict = {LassoSelector: _deselect_lasso, RectangleSelector: _deselect_rectangle}
    kwargs_dict = {LassoSelector: {}, RectangleSelector: {'interactive': True}}

    lasso = selector(axs[0], onselect_dict[selector], **kwargs_dict[selector])
    lasso2 = selector(axs[1], deselect_dict[selector], **kwargs_dict[selector])
    if not plotonly:
        raw_input('Press Enter when done')
    lasso.disconnect_events()
    lasso2.disconnect_events()
    if not interactive:
        plt.ioff()
    if sav_img:
        if type(sav_img) is str:
            fig.savefig(sav_img)
        else:
            fig.savefig('selected.png')
    return selected['data'], selected['mask']

# if __name__ == '__main__':
#     from skimage.data import coins
#
#     data = coins()
#     selected, mask = imshow_select(data, 'lasso', True)
#     plt.imsave('selected.png', selected)
#     plt.imsave('mask.png', mask)
