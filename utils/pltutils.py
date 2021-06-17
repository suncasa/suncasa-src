def multicolor_text(axes, x, y, textin, cmap=None,wratio=0.5, bbox={}, **kw):
    import matplotlib.pyplot as plt
    from matplotlib import transforms
    """
        Take a list of strings ``ls`` and colors ``lc`` and place them next to each
        other, with text ls[i] being shown in color lc[i].

        This example shows how to do both vertical and horizontal text, and will
        pass all keyword arguments to plt.text, so you can set the font size,
        family, etc.
        """
    fig = plt.gcf()
    t = axes.transAxes

    # horizontal version
    if cmap is None:
        cmap = plt.cm.RdYlBu

    textin = list(textin)
    ntextin = len(textin)
    for idx, s in enumerate(textin):
        # s_ = ['_']*ntextin
        # s_[idx] = s
        c = cmap(float(idx) / (ntextin - 1))
        # tx = axes.text(x, y, ''.join(s_), color=c, transform=t, **kw)
        tx = axes.text(x, y, s, color=c, transform=t, **kw)
        if bbox:
            tx.set_bbox(bbox)
        tx.draw(fig.canvas.get_renderer())
        ex = tx.get_window_extent()
        t = transforms.offset_copy(tx._transform,fig=fig, x=ex.width*wratio, units='points')


def align_marker(marker, halign='center', valign='middle', ):
    from matplotlib import markers
    from matplotlib.path import Path
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """

    if isinstance(halign, str):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, str):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)