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