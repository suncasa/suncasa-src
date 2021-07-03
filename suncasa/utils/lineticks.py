import numpy as np
from matplotlib import transforms


def get_perp_vec(u1, u2, direction=1):
    """Return the unit vector perpendicular to the vector u2-u1."""

    x1, y1 = u1
    x2, y2 = u2
    vx, vy = x2 - x1, y2 - y1
    v = np.linalg.norm((vx, vy))
    wx, wy = -vy / v * direction, vx / v * direction
    return wx, wy


def get_av_vec(u1, u2):
    """Return the average unit vector between u1 and u2."""

    u1x, u1y = u1
    u2x, u2y = u2
    dx, dy = u1x + u2x, u1y + u2y
    dlen = np.linalg.norm((dx, dy))
    return dx / dlen, dy / dlen


class LineTicks:
    def __init__(self, line, idx, tick_length, tick_shift=3.0, direction=1, label=None, label_color='k',
                 fontsize='xx-small', **kwargs):
        self.line = line
        self.idx = idx
        self.tick_length = tick_length
        self.direction = direction
        self.label = label
        self.ticks = []
        self.tick_labels = []
        self.label_color = label_color
        self.fontsize = fontsize
        self.tick_shift = tick_shift
        self.tick_styles = kwargs
        # If no colour is specified for the ticks, set it to the line colour
        if not set(('c', 'color')).intersection(kwargs.keys()):
            self.tick_styles['color'] = self.line.get_color()

        self.ax = line.axes
        self.ax.callbacks.connect('xlim_changed', self.on_change_lims)
        self.ax.callbacks.connect('ylim_changed', self.on_change_lims)
        cid = self.ax.figure.canvas.mpl_connect('resize_event', self.on_resize)

        # get scale tick_length in display coordinates to data coordinates
        xy1, xy2 = self.ax.transData.inverted().transform([[0, 0], [self.tick_length, 0]])
        self.tick_length = xy2[0] - xy1[0]

        self.add_ticks(self.ax)

    def add_ticks(self, ax):
        ax.set_autoscale_on(False)  # Otherwise, infinite loop
        # Transform to  display coordinates
        # z =ax.transData.transform(np.array(self.line.get_data()).T)
        z = np.array(self.line.get_data()).T
        x, y = zip(*z)

        # Remove existing ticks
        for tick in self.ticks:
            ax.lines.remove(tick)
        # Remove references to the ticks so they can be garbage-collected
        self.ticks = []

        # Remove any existing tick labels
        for ticklabel in self.tick_labels:
            ax.texts.remove(ticklabel)
        self.tick_labels = []

        for j, i in enumerate(self.idx):
            if i == 0:
                # The first tick is perpendicular to the line between the
                # first two points
                tx, ty = get_perp_vec((x[0], y[0]), (x[1], y[1]),
                                      self.direction)
            elif i == len(x) - 1:
                # The last tick is perpendicular to the line between the
                # last two points
                tx, ty = get_perp_vec((x[-2], y[-2]), (x[-1], y[-1]),
                                      self.direction)
            else:
                # General tick marks bisect the incoming and outgoing line
                # segments
                u1 = get_perp_vec((x[i - 1], y[i - 1]), (x[i], y[i]),
                                  self.direction)
                u2 = get_perp_vec((x[i], y[i]), (x[i + 1], y[i + 1]),
                                  self.direction)
                tx, ty = get_av_vec(u1, u2)
            tx, ty = self.tick_length * tx, self.tick_length * ty
            this_tick, = ax.plot((x[i], x[i] + tx), (y[i], y[i] + ty), **self.tick_styles)
            self.ticks.append(this_tick)

            if self.label:
                this_ticklabel = ax.text(x[i] + tx * self.tick_shift, y[i] + ty * self.tick_shift, self.label[j],
                                         ha='center', va='center', clip_on=True, color=self.label_color,
                                         fontsize=self.fontsize)
                self.tick_labels.append(this_ticklabel)

    def on_change_lims(self, ax):
        self.add_ticks(ax)

    def on_resize(self, event):
        self.add_ticks(self.ax)
