# -*- coding: utf-8 -*-
"""

Conventions:

"constrain_x" means to constrain the variable with either
another kiwisolver variable, or a float.  i.e. `constrain_width(0.2)`
will set a constraint that the width has to be 0.2 and this constraint is
permanent - i.e. it will not be removed if it becomes obsolete.

"edit_x" means to set x to a value (just a float), and that this value can
change.  So `edit_width(0.2)` will set width to be 0.2, but `edit_width(0.3)`
will allow it to change to 0.3 later.  Note that these values are still just
"suggestions" in `kiwisolver` parlance, and could be over-ridden by
other constrains.

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools
import kiwisolver as kiwi
import logging
import numpy as np
import warnings

import matplotlib

_log = logging.getLogger(__name__)


# renderers can be complicated
def get_renderer(fig):
    if fig._cachedRenderer:
        renderer = fig._cachedRenderer
    else:
        canvas = fig.canvas
        if canvas and hasattr(canvas, "get_renderer"):
            renderer = canvas.get_renderer()
        else:
            # not sure if this can happen
            # seems to with PDF...
            _log.info("constrained_layout : falling back to Agg renderer")
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            canvas = FigureCanvasAgg(fig)
            renderer = canvas.get_renderer()

    return renderer


class LayoutBox(object):
    """
    Basic rectangle representation using kiwi solver variables
    """

    def __init__(self, parent=None, name='', tightwidth=False,
                 tightheight=False, artist=None,
                 lower_left=(0, 0), upper_right=(1, 1), pos=False,
                 subplot=False, h_pad=None, w_pad=None):
        Variable = kiwi.Variable
        self.parent = parent
        self.name = name
        sn = self.name + '_'
        if parent is None:
            self.solver = kiwi.Solver()
            self.constrained_layout_called = 0
        else:
            self.solver = parent.solver
            self.constrained_layout_called = None
            # parent wants to know about this child!
            parent.add_child(self)
        # keep track of artist associated w/ this layout.  Can be none
        self.artist = artist
        # keep track if this box is supposed to be a pos that is constrained
        # by the parent.
        self.pos = pos
        # keep track of whether we need to match this subplot up with others.
        self.subplot = subplot

        # we need the str below for Py 2 which complains the string is unicode
        self.top = Variable(str(sn + 'top'))
        self.bottom = Variable(str(sn + 'bottom'))
        self.left = Variable(str(sn + 'left'))
        self.right = Variable(str(sn + 'right'))

        self.width = Variable(str(sn + 'width'))
        self.height = Variable(str(sn + 'height'))
        self.h_center = Variable(str(sn + 'h_center'))
        self.v_center = Variable(str(sn + 'v_center'))

        self.min_width = Variable(str(sn + 'min_width'))
        self.min_height = Variable(str(sn + 'min_height'))
        self.pref_width = Variable(str(sn + 'pref_width'))
        self.pref_height = Variable(str(sn + 'pref_height'))
        # margis are only used for axes-position layout boxes.  maybe should
        # be a separate subclass:
        self.left_margin = Variable(str(sn + 'left_margin'))
        self.right_margin = Variable(str(sn + 'right_margin'))
        self.bottom_margin = Variable(str(sn + 'bottom_margin'))
        self.top_margin = Variable(str(sn + 'top_margin'))
        # mins
        self.left_margin_min = Variable(str(sn + 'left_margin_min'))
        self.right_margin_min = Variable(str(sn + 'right_margin_min'))
        self.bottom_margin_min = Variable(str(sn + 'bottom_margin_min'))
        self.top_margin_min = Variable(str(sn + 'top_margin_min'))

        right, top = upper_right
        left, bottom = lower_left
        self.tightheight = tightheight
        self.tightwidth = tightwidth
        self.add_constraints()
        self.children = []
        self.subplotspec = None
        if self.pos:
            self.constrain_margins()
        self.h_pad = h_pad
        self.w_pad = w_pad

    def constrain_margins(self):
        """
        Only do this for pos.  This sets a variable distance
        margin between the position of the axes and the outer edge of
        the axes.

        Margins are variable because they change with the fogure size.

        Margin minimums are set to make room for axes decorations.  However,
        the margins can be larger if we are mathicng the position size to
        otehr axes.
        """
        sol = self.solver

        # left
        if not sol.hasEditVariable(self.left_margin_min):
            sol.addEditVariable(self.left_margin_min, 'strong')
            sol.suggestValue(self.left_margin_min, 0.0001)
        c = (self.left_margin == self.left - self.parent.left)
        self.solver.addConstraint(c | 'required')
        c = (self.left_margin >= self.left_margin_min)
        self.solver.addConstraint(c | 'strong')

        # right
        if not sol.hasEditVariable(self.right_margin_min):
            sol.addEditVariable(self.right_margin_min, 'strong')
            sol.suggestValue(self.right_margin_min, 0.0001)
        c = (self.right_margin == self.parent.right - self.right)
        self.solver.addConstraint(c | 'required')
        c = (self.right_margin >= self.right_margin_min)
        self.solver.addConstraint(c | 'required')
        # bottom
        if not sol.hasEditVariable(self.bottom_margin_min):
            sol.addEditVariable(self.bottom_margin_min, 'strong')
            sol.suggestValue(self.bottom_margin_min, 0.0001)
        c = (self.bottom_margin == self.bottom - self.parent.bottom)
        self.solver.addConstraint(c | 'required')
        c = (self.bottom_margin >= self.bottom_margin_min)
        self.solver.addConstraint(c | 'required')
        # top
        if not sol.hasEditVariable(self.top_margin_min):
            sol.addEditVariable(self.top_margin_min, 'strong')
            sol.suggestValue(self.top_margin_min, 0.0001)
        c = (self.top_margin == self.parent.top - self.top)
        self.solver.addConstraint(c | 'required')
        c = (self.top_margin >= self.top_margin_min)
        self.solver.addConstraint(c | 'required')

    def add_child(self, child):
        self.children += [child]

    def remove_child(self, child):
        try:
            self.children.remove(child)
        except ValueError:
            _log.info("Tried to remove child that doesn't belong to parent")

    def add_constraints(self):
        sol = self.solver
        # never let width and height go negative.
        for i in [self.min_width, self.min_height]:
            sol.addEditVariable(i, 1e9)
            sol.suggestValue(i, 0.0)
        # define relation ships between things thing width and right and left
        self.hard_constraints()
        # self.soft_constraints()
        if self.parent:
            self.parent_constrain()
        # sol.updateVariables()

    def parent_constrain(self):
        parent = self.parent
        hc = [self.left >= parent.left,
              self.bottom >= parent.bottom,
              self.top <= parent.top,
              self.right <= parent.right]
        for c in hc:
            self.solver.addConstraint(c | 'required')

    def hard_constraints(self):
        hc = [self.width == self.right - self.left,
              self.height == self.top - self.bottom,
              self.h_center == (self.left + self.right) * 0.5,
              self.v_center == (self.top + self.bottom) * 0.5,
              self.width >= self.min_width,
              self.height >= self.min_height]
        for c in hc:
            self.solver.addConstraint(c | 'required')

    def soft_constraints(self):
        sol = self.solver
        if self.tightwidth:
            suggest = 0.
        else:
            suggest = 20.
        c = (self.pref_width == suggest)
        for i in c:
            sol.addConstraint(i | 'required')
        if self.tightheight:
            suggest = 0.
        else:
            suggest = 20.
        c = (self.pref_height == suggest)
        for i in c:
            sol.addConstraint(i | 'required')

        c = [(self.width >= suggest),
             (self.height >= suggest)]
        for i in c:
            sol.addConstraint(i | 150000)

    def set_parent(self, parent):
        ''' replace the parent of this with the new parent
        '''
        self.parent = parent
        self.parent_constrain()

    def constrain_geometry(self, left, bottom, right, top, strength='strong'):
        hc = [self.left == left,
              self.right == right,
              self.bottom == bottom,
              self.top == top]
        for c in hc:
            self.solver.addConstraint((c | strength))
        # self.solver.updateVariables()

    def constrain_same(self, other, strength='strong'):
        """
        Make the layoutbox have same position as other layoutbox
        """
        hc = [self.left == other.left,
              self.right == other.right,
              self.bottom == other.bottom,
              self.top == other.top]
        for c in hc:
            self.solver.addConstraint((c | strength))

    def constrain_left_margin(self, margin, strength='strong'):
        c = (self.left == self.parent.left + margin)
        self.solver.addConstraint(c | strength)

    def edit_left_margin_min(self, margin):
        self.solver.suggestValue(self.left_margin_min, margin)

    def constrain_right_margin(self, margin, strength='strong'):
        c = (self.right == self.parent.right - margin)
        self.solver.addConstraint(c | strength)

    def edit_right_margin_min(self, margin):
        self.solver.suggestValue(self.right_margin_min, margin)

    def constrain_bottom_margin(self, margin, strength='strong'):
        c = (self.bottom == self.parent.bottom + margin)
        self.solver.addConstraint(c | strength)

    def edit_bottom_margin_min(self, margin):
        self.solver.suggestValue(self.bottom_margin_min, margin)

    def constrain_top_margin(self, margin, strength='strong'):
        c = (self.top == self.parent.top - margin)
        self.solver.addConstraint(c | strength)

    def edit_top_margin_min(self, margin):
        self.solver.suggestValue(self.top_margin_min, margin)

    def get_rect(self):
        return (self.left.value(), self.bottom.value(),
                self.width.value(), self.height.value())

    def update_variables(self):
        '''
        Update *all* the variables that are part of the solver this LayoutBox
        is created with
        '''
        self.solver.updateVariables()

    def edit_height(self, height, strength='strong'):
        '''
        Set the height of the layout box.

        This is done as an editable variable so that the value can change
        due to resizing.
        '''
        sol = self.solver
        for i in [self.height]:
            if not sol.hasEditVariable(i):
                sol.addEditVariable(i, strength)
        sol.suggestValue(self.height, height)

    def constrain_height(self, height, strength='strong'):
        '''
        Constrain the height of the layout box.  height is
        either a float or a layoutbox.height.
        '''
        c = (self.height == height)
        self.solver.addConstraint(c | strength)

    def constrain_height_min(self, height, strength='strong'):
        c = (self.height >= height)
        self.solver.addConstraint(c | strength)

    def edit_width(self, width, strength='strong'):
        sol = self.solver
        for i in [self.width]:
            if not sol.hasEditVariable(i):
                sol.addEditVariable(i, strength)
        sol.suggestValue(self.width, width)

    def constrain_width(self, width, strength='strong'):
        '''
        Constrain the width of the layout box.  `width` is
        either a float or a layoutbox.width.
        '''
        c = (self.width == width)
        self.solver.addConstraint(c | strength)

    def constrain_width_min(self, width, strength='strong'):
        c = (self.width >= width)
        self.solver.addConstraint(c | strength)

    def constrain_left(self, left,  strength='strong'):
        c = (self.left == left)
        self.solver.addConstraint(c | strength)

    def constrain_bottom(self, bottom, strength='strong'):
        c = (self.bottom == bottom)
        self.solver.addConstraint(c | strength)

    def constrain_right(self, right, strength='strong'):
        c = (self.right == right)
        self.solver.addConstraint(c | strength)

    def constrain_top(self, top, strength='strong'):
        c = (self.top == top)
        self.solver.addConstraint(c | strength)

    def _is_subplotspec_layoutbox(self):
        '''
        Helper to check if this layoutbox is the layoutbox of a
        subplotspec
        '''
        name = (self.name).split('.')[-1][:-3]
        if name == 'ss':
            return True
        return False

    def _is_gridspec_layoutbox(self):
        '''
        Helper to check if this layoutbox is the layoutbox of a
        gridspec
        '''
        name = (self.name).split('.')[-1][:-3]
        if name == 'gridspec':
            return True
        return False

    def find_child_subplots(self):
        '''
        Find children of this layout box that are subplots.  We want to line
        poss up, and this is an easy way to find them all.
        '''
        if self.subplot:
            subplots = [self]
        else:
            subplots = []
        for child in self.children:
            subplots += child.find_child_subplots()
        return subplots

    def layout_from_subplotspec(self, subspec,
                                name='', artist=None, pos=False):
        '''  Make a layout box from a subplotspec. The layout box is
        constrained to be a fraction of the width/height of the parent,
        and be a fraction of the parent width/height from the left/bottom
        of the parent.  Therefore the parent can move around and the
        layout for the subplot spec should move with it.

        The parent is *usually* the gridspec that made the subplotspec.??
        '''
        lb = LayoutBox(parent=self, name=name, artist=artist, pos=pos)
        gs = subspec.get_gridspec()
        nrows, ncols = gs.get_geometry()
        parent = self.parent

        # OK, now, we want to set the position of this subplotspec
        # based on its subplotspec parameters.  The new gridspec will inherit.

        # from gridspec.  prob should be new method in gridspec
        left = 0.0
        right = 1.0
        bottom = 0.0
        top = 1.0
        totWidth = right-left
        totHeight = top-bottom
        hspace = 0.
        wspace = 0.

        # calculate accumulated heights of columns
        cellH = totHeight / (nrows + hspace * (nrows - 1))
        sepH = hspace*cellH

        if gs._row_height_ratios is not None:
            netHeight = cellH * nrows
            tr = float(sum(gs._row_height_ratios))
            cellHeights = [netHeight*r/tr for r in gs._row_height_ratios]
        else:
            cellHeights = [cellH] * nrows

        sepHeights = [0] + ([sepH] * (nrows - 1))
        cellHs = np.add.accumulate(np.ravel(
                list(zip(sepHeights, cellHeights))))

        # calculate accumulated widths of rows
        cellW = totWidth/(ncols + wspace * (ncols - 1))
        sepW = wspace*cellW

        if gs._col_width_ratios is not None:
            netWidth = cellW * ncols
            tr = float(sum(gs._col_width_ratios))
            cellWidths = [netWidth * r / tr for r in gs._col_width_ratios]
        else:
            cellWidths = [cellW] * ncols

        sepWidths = [0] + ([sepW] * (ncols - 1))
        cellWs = np.add.accumulate(np.ravel(list(zip(sepWidths, cellWidths))))

        figTops = [top - cellHs[2 * rowNum] for rowNum in range(nrows)]
        figBottoms = [top - cellHs[2 * rowNum + 1] for rowNum in range(nrows)]
        figLefts = [left + cellWs[2 * colNum] for colNum in range(ncols)]
        figRights = [left + cellWs[2 * colNum + 1] for colNum in range(ncols)]

        rowNum, colNum = divmod(subspec.num1, ncols)
        figBottom = figBottoms[rowNum]
        figTop = figTops[rowNum]
        figLeft = figLefts[colNum]
        figRight = figRights[colNum]

        if subspec.num2 is not None:

            rowNum2, colNum2 = divmod(subspec.num2, ncols)
            figBottom2 = figBottoms[rowNum2]
            figTop2 = figTops[rowNum2]
            figLeft2 = figLefts[colNum2]
            figRight2 = figRights[colNum2]

            figBottom = min(figBottom, figBottom2)
            figLeft = min(figLeft, figLeft2)
            figTop = max(figTop, figTop2)
            figRight = max(figRight, figRight2)
        # These are numbers relative to 0,0,1,1.  Need to constrain
        # relative to parent.

        width = figRight - figLeft
        height = figTop - figBottom
        parent = self.parent
        cs = [self.left == parent.left + parent.width * figLeft,
              self.bottom == parent.bottom + parent.height * figBottom,
              self.width == parent.width * width,
              self.height == parent.height * height]
        for c in cs:
            self.solver.addConstraint((c | 'required'))

        return lb

    def __repr__(self):
        args = (self.name, self.left.value(), self.bottom.value(),
                self.right.value(), self.top.value())
        return ('LayoutBox: %25s, (left: %1.3f) (bot: %1.3f) '
               '(right: %1.3f)  (top: %1.3f) ') % args


# Utility functions that act on layoutboxes...
def hstack(boxes, padding=0, strength='strong'):
    '''
    Stack LayoutBox instances from left to right.
    `padding` is in figure-relative units.
    '''

    for i in range(1, len(boxes)):
        c = (boxes[i-1].right + padding <= boxes[i].left)
        boxes[i].solver.addConstraint(c | strength)


def hpack(boxes, padding=0, strength='strong'):
    '''
    Stack LayoutBox instances from left to right.
    '''

    for i in range(1, len(boxes)):
        c = (boxes[i-1].right + padding == boxes[i].left)
        boxes[i].solver.addConstraint(c | strength)


def vstack(boxes, padding=0, strength='strong'):
    '''
    Stack LayoutBox instances from top to bottom
    '''

    for i in range(1, len(boxes)):
        c = (boxes[i-1].bottom - padding >= boxes[i].top)
        boxes[i].solver.addConstraint(c | strength)


def vpack(boxes, padding=0, strength='strong'):
    '''
    Stack LayoutBox instances from top to bottom
    '''

    for i in range(1, len(boxes)):
        c = (boxes[i-1].bottom - padding >= boxes[i].top)
        boxes[i].solver.addConstraint(c | strength)


def match_heights(boxes, height_ratios=None, strength='medium'):
    '''
    Stack LayoutBox instances from top to bottom
    '''

    if height_ratios is None:
        height_ratios = np.ones(len(boxes))
    for i in range(1, len(boxes)):
        c = (boxes[i-1].height ==
             boxes[i].height*height_ratios[i-1]/height_ratios[i])
        boxes[i].solver.addConstraint(c | strength)


def match_widths(boxes, width_ratios=None, strength='medium'):
    '''
    Stack LayoutBox instances from top to bottom
    '''

    if width_ratios is None:
        width_ratios = np.ones(len(boxes))
    for i in range(1, len(boxes)):
        c = (boxes[i-1].width ==
             boxes[i].width*width_ratios[i-1]/width_ratios[i])
        boxes[i].solver.addConstraint(c | strength)


def vstackeq(boxes, padding=0, height_ratios=None):
    vstack(boxes, padding=padding)
    match_heights(boxes, height_ratios=height_ratios)


def hstackeq(boxes, padding=0, width_ratios=None):
    hstack(boxes, padding=padding)
    match_widths(boxes, width_ratios=width_ratios)


def align(boxes, attr, strength='strong'):
    cons = []
    for box in boxes[1:]:
        cons = (getattr(boxes[0], attr) == getattr(box, attr))
        boxes[0].solver.addConstraint(cons | strength)


def match_top_margins(boxes, levels=1):
    box0 = boxes[0]
    top0 = box0
    for n in range(levels):
        top0 = top0.parent
    for box in boxes[1:]:
        topb = box
        for n in range(levels):
            topb = topb.parent
        c = (box0.top-top0.top == box.top-topb.top)
        box0.solver.addConstraint(c | 'strong')


def match_bottom_margins(boxes, levels=1):
    box0 = boxes[0]
    top0 = box0
    for n in range(levels):
        top0 = top0.parent
    for box in boxes[1:]:
        topb = box
        for n in range(levels):
            topb = topb.parent
        c = (box0.bottom-top0.bottom == box.bottom-topb.bottom)
        box0.solver.addConstraint(c | 'strong')


def match_left_margins(boxes, levels=1):
    box0 = boxes[0]
    top0 = box0
    for n in range(levels):
        top0 = top0.parent
    for box in boxes[1:]:
        topb = box
        for n in range(levels):
            topb = topb.parent
        c = (box0.left-top0.left == box.left-topb.left)
        box0.solver.addConstraint(c | 'strong')


def match_right_margins(boxes, levels=1):
    box0 = boxes[0]
    top0 = box0
    for n in range(levels):
        top0 = top0.parent
    for box in boxes[1:]:
        topb = box
        for n in range(levels):
            topb = topb.parent
        c = (box0.right-top0.right == box.right-topb.right)
        box0.solver.addConstraint(c | 'strong')


def match_width_margins(boxes, levels=1):
    match_left_margins(boxes, levels=levels)
    match_right_margins(boxes, levels=levels)


def match_height_margins(boxes, levels=1):
    match_top_margins(boxes, levels=levels)
    match_bottom_margins(boxes, levels=levels)


def match_margins(boxes, levels=1):
    match_width_margins(boxes, levels=levels)
    match_height_margins(boxes, levels=levels)


_layoutboxobjnum = itertools.count()


def seq_id():
    '''
    Generate a short sequential id for layoutbox objects...
    '''

    global _layoutboxobjnum

    return ('%03d' % (next(_layoutboxobjnum)))


def print_children(lb):
    '''
    Print the children of the layoutbox
    '''
    print(lb)
    for child in lb.children:
        print_children(child)


def nonetree(lb):
    '''
    Make all elements in this tree none...  This signals not to do any more
    layout.
    '''
    if lb is not None:
        if lb.parent is None:
            # Clear the solver.  Hopefully this garbage collects.
            lb.solver.reset()
            nonechildren(lb)
        else:
            nonetree(lb.parent)


def nonechildren(lb):
    for child in lb.children:
        nonechildren(child)
    lb.artist._layoutbox = None
    lb = None


def print_tree(lb):
    '''
    Print the tree of layoutboxes
    '''

    if lb.parent is None:
        print('LayoutBox Tree\n')
        print('==============\n')
        print_children(lb)
        print('\n')
    else:
        print_tree(lb.parent)


def plot_children(fig, box, level=0, printit=True):
    '''
    Simple plotting to show where boxes are
    '''
    import matplotlib
    import matplotlib.pyplot as plt

    if isinstance(fig, matplotlib.figure.Figure):
        ax = fig.add_axes([0., 0., 1., 1.])
        ax.set_facecolor([1., 1., 1., 0.7])
        ax.set_alpha(0.3)
        fig.draw(fig.canvas.get_renderer())
    else:
        ax = fig

    import matplotlib.patches as patches
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    if printit:
        print("Level:", level)
    for child in box.children:
        rect = child.get_rect()
        if printit:
            print(child)
        ax.add_patch(
            patches.Rectangle(
                (child.left.value(), child.bottom.value()),   # (x,y)
                child.width.value(),          # width
                child.height.value(),          # height
                fc='none',
                alpha=0.8,
                ec=colors[level]
            )
        )
        if level > 0:
            name = child.name.split('.')[-1]
            if level % 2 == 0:
                ax.text(child.left.value(), child.bottom.value(), name,
                        size=12-level, color=colors[level])
            else:
                ax.text(child.right.value(), child.top.value(), name,
                        ha='right', va='top', size=12-level,
                        color=colors[level])

        plot_children(ax, child, level=level+1, printit=printit)
