"""
Place a table below the x-axis at location loc.

The table consists of a grid of cells.

The grid need not be rectangular and can have holes.

Cells are added by specifying their row and column.

For the purposes of positioning the cell at (0, 0) is
assumed to be at the top left and the cell at (max_row, max_col)
is assumed to be at bottom right.

You can add additional cells outside this range to have convenient
ways of positioning more interesting grids.

Author    : John Gill <jng@europe.renre.com>
Copyright : 2004 John Gill and John Hunter
License   : matplotlib license

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

import warnings

from . import artist
from .artist import Artist, allow_rasterization
from .patches import Rectangle
from matplotlib import docstring
from .text import Text
from .transforms import Bbox
from matplotlib.path import Path


class Cell(Rectangle):
    """
    A cell is a Rectangle with some associated text.

    """
    PAD = 0.1  # padding between text and rectangle

    def __init__(self, xy, width, height,
                 edgecolor='k', facecolor='w',
                 fill=True,
                 text='',
                 loc=None,
                 fontproperties=None
                 ):

        # Call base
        Rectangle.__init__(self, xy, width=width, height=height,
                           edgecolor=edgecolor, facecolor=facecolor)
        self.set_clip_on(False)

        # Create text object
        if loc is None:
            loc = 'right'
        self._loc = loc
        self._text = Text(x=xy[0], y=xy[1], text=text,
                          fontproperties=fontproperties)
        self._text.set_clip_on(False)

    def set_transform(self, trans):
        Rectangle.set_transform(self, trans)
        # the text does not get the transform!
        self.stale = True

    def set_figure(self, fig):
        Rectangle.set_figure(self, fig)
        self._text.set_figure(fig)

    def get_text(self):
        'Return the cell Text intance'
        return self._text

    def set_fontsize(self, size):
        self._text.set_fontsize(size)
        self.stale = True

    def get_fontsize(self):
        'Return the cell fontsize'
        return self._text.get_fontsize()

    def auto_set_font_size(self, renderer):
        """ Shrink font size until text fits. """
        fontsize = self.get_fontsize()
        required = self.get_required_width(renderer)
        while fontsize > 1 and required > self.get_width():
            fontsize -= 1
            self.set_fontsize(fontsize)
            required = self.get_required_width(renderer)

        return fontsize

    @allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        # draw the rectangle
        Rectangle.draw(self, renderer)

        # position the text
        self._set_text_position(renderer)
        self._text.draw(renderer)
        self.stale = False

    def _set_text_position(self, renderer):
        """ Set text up so it draws in the right place.

        Currently support 'left', 'center' and 'right'
        """
        bbox = self.get_window_extent(renderer)
        l, b, w, h = bbox.bounds

        # draw in center vertically
        self._text.set_verticalalignment('center')
        y = b + (h / 2.0)

        # now position horizontally
        if self._loc == 'center':
            self._text.set_horizontalalignment('center')
            x = l + (w / 2.0)
        elif self._loc == 'left':
            self._text.set_horizontalalignment('left')
            x = l + (w * self.PAD)
        else:
            self._text.set_horizontalalignment('right')
            x = l + (w * (1.0 - self.PAD))

        self._text.set_position((x, y))

    def get_text_bounds(self, renderer):
        """ Get text bounds in axes co-ordinates. """
        bbox = self._text.get_window_extent(renderer)
        bboxa = bbox.inverse_transformed(self.get_data_transform())
        return bboxa.bounds

    def get_required_width(self, renderer):
        """ Get width required for this cell. """
        l, b, w, h = self.get_text_bounds(renderer)
        return w * (1.0 + (2.0 * self.PAD))

    def set_text_props(self, **kwargs):
        'update the text properties with kwargs'
        self._text.update(kwargs)
        self.stale = True


class CustomCell(Cell):
    """
    A subclass of Cell where the sides may be visibly toggled.

    """

    _edges = 'BRTL'
    _edge_aliases = {'open':         '',
                     'closed':       _edges,  # default
                     'horizontal':   'BT',
                     'vertical':     'RL'
                     }

    def __init__(self, *args, **kwargs):
        visible_edges = kwargs.pop('visible_edges')
        Cell.__init__(self, *args, **kwargs)
        self.visible_edges = visible_edges

    @property
    def visible_edges(self):
        return self._visible_edges

    @visible_edges.setter
    def visible_edges(self, value):
        if value is None:
            self._visible_edges = self._edges
        elif value in self._edge_aliases:
            self._visible_edges = self._edge_aliases[value]
        else:
            for edge in value:
                if edge not in self._edges:
                    raise ValueError('Invalid edge param {}, must only be one '
                                     'of {} or string of {}'.format(
                                         value,
                                         ", ".join(self._edge_aliases),
                                         ", ".join(self._edges)))
            self._visible_edges = value
        self.stale = True

    def get_path(self):
        'Return a path where the edges specified by _visible_edges are drawn'

        codes = [Path.MOVETO]

        for edge in self._edges:
            if edge in self._visible_edges:
                codes.append(Path.LINETO)
            else:
                codes.append(Path.MOVETO)

        if Path.MOVETO not in codes[1:]:  # All sides are visible
            codes[-1] = Path.CLOSEPOLY

        return Path(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            codes,
            readonly=True
            )


class Table(Artist):
    """
    Create a table of cells.

    Table can have (optional) row and column headers.

    Each entry in the table can be either text or patches.

    Column widths and row heights for the table can be specified.

    Return value is a sequence of text, line and patch instances that make
    up the table
    """
    codes = {'best': 0,
             'upper right':  1,  # default
             'upper left':   2,
             'lower left':   3,
             'lower right':  4,
             'center left':  5,
             'center right': 6,
             'lower center': 7,
             'upper center': 8,
             'center':       9,
             'top right':    10,
             'top left':     11,
             'bottom left':  12,
             'bottom right': 13,
             'right':        14,
             'left':         15,
             'top':          16,
             'bottom':       17,
             }

    FONTSIZE = 10
    AXESPAD = 0.02    # the border between the axes and table edge

    def __init__(self, ax, loc=None, bbox=None, **kwargs):

        Artist.__init__(self)

        if isinstance(loc, six.string_types) and loc not in self.codes:
            warnings.warn('Unrecognized location %s. Falling back on '
                          'bottom; valid locations are\n%s\t' %
                          (loc, '\n\t'.join(self.codes)))
            loc = 'bottom'
        if isinstance(loc, six.string_types):
            loc = self.codes.get(loc, 1)
        self.set_figure(ax.figure)
        self._axes = ax
        self._loc = loc
        self._bbox = bbox

        # use axes coords
        self.set_transform(ax.transAxes)

        self._texts = []
        self._cells = {}
        self._edges = None
        self._autoRows = []
        self._autoColumns = []
        self._autoFontsize = True
        self.update(kwargs)

        self.set_clip_on(False)

    def add_cell(self, row, col, *args, **kwargs):
        """
        Add a cell to the table.

        Parameters
        ----------
        row : int
            Row index
        col : int
            Column index

        Returns
        -------
        `CustomCell`: Automatically created cell

        """
        xy = (0, 0)
        cell = CustomCell(xy, visible_edges=self.edges, *args, **kwargs)
        self[row, col] = cell
        return cell

    def __setitem__(self, position, cell):
        """
        Set a customcell in a given position
        """
        if not isinstance(cell, CustomCell):
            raise TypeError('Table only accepts CustomCell')
        try:
            row, col = position[0], position[1]
        except Exception:
            raise KeyError('Only tuples length 2 are accepted as coordinates')
        cell.set_figure(self.figure)
        cell.set_transform(self.get_transform())
        cell.set_clip_on(False)
        self._cells[row, col] = cell
        self.stale = True

    def __getitem__(self, position):
        """
        Retreive a custom cell from a given position
        """
        try:
            row, col = position[0], position[1]
        except Exception:
            raise KeyError('Only tuples length 2 are accepted as coordinates')
        return self._cells[row, col]

    @property
    def edges(self):
        return self._edges

    @edges.setter
    def edges(self, value):
        self._edges = value
        self.stale = True

    def _approx_text_height(self):
        return (self.FONTSIZE / 72.0 * self.figure.dpi /
                self._axes.bbox.height * 1.2)

    @allow_rasterization
    def draw(self, renderer):
        # Need a renderer to do hit tests on mouseevent; assume the last one
        # will do
        if renderer is None:
            renderer = self.figure._cachedRenderer
        if renderer is None:
            raise RuntimeError('No renderer defined')

        if not self.get_visible():
            return
        renderer.open_group('table')
        self._update_positions(renderer)

        for key in sorted(self._cells):
            self._cells[key].draw(renderer)

        renderer.close_group('table')
        self.stale = False

    def _get_grid_bbox(self, renderer):
        """Get a bbox, in axes co-ordinates for the cells.

        Only include those in the range (0,0) to (maxRow, maxCol)"""
        boxes = [cell.get_window_extent(renderer)
                 for (row, col), cell in six.iteritems(self._cells)
                 if row >= 0 and col >= 0]
        bbox = Bbox.union(boxes)
        return bbox.inverse_transformed(self.get_transform())

    def contains(self, mouseevent):
        """Test whether the mouse event occurred in the table.

        Returns T/F, {}
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)

        # TODO: Return index of the cell containing the cursor so that the user
        # doesn't have to bind to each one individually.
        renderer = self.figure._cachedRenderer
        if renderer is not None:
            boxes = [cell.get_window_extent(renderer)
                     for (row, col), cell in six.iteritems(self._cells)
                     if row >= 0 and col >= 0]
            bbox = Bbox.union(boxes)
            return bbox.contains(mouseevent.x, mouseevent.y), {}
        else:
            return False, {}

    def get_children(self):
        'Return the Artists contained by the table'
        return list(six.itervalues(self._cells))
    get_child_artists = get_children  # backward compatibility

    def get_window_extent(self, renderer):
        'Return the bounding box of the table in window coords'
        boxes = [cell.get_window_extent(renderer)
                 for cell in six.itervalues(self._cells)]
        return Bbox.union(boxes)

    def _do_cell_alignment(self):
        """ Calculate row heights and column widths.

        Position cells accordingly.
        """
        # Calculate row/column widths
        widths = {}
        heights = {}
        for (row, col), cell in six.iteritems(self._cells):
            height = heights.setdefault(row, 0.0)
            heights[row] = max(height, cell.get_height())
            width = widths.setdefault(col, 0.0)
            widths[col] = max(width, cell.get_width())

        # work out left position for each column
        xpos = 0
        lefts = {}
        for col in sorted(widths):
            lefts[col] = xpos
            xpos += widths[col]

        ypos = 0
        bottoms = {}
        for row in sorted(heights, reverse=True):
            bottoms[row] = ypos
            ypos += heights[row]

        # set cell positions
        for (row, col), cell in six.iteritems(self._cells):
            cell.set_x(lefts[col])
            cell.set_y(bottoms[row])

    def auto_set_column_width(self, col):
        """ Given column indexs in either List, Tuple or int. Will be able to
        automatically set the columns into optimal sizes.

        Here is the example of the input, which triger automatic adjustment on
        columns to optimal size by given index numbers.
        -1: the row labling
        0: the 1st column
        1: the 2nd column

        Args:
            col(List): list of indexs
            >>>table.auto_set_column_width([-1,0,1])

            col(Tuple): tuple of indexs
            >>>table.auto_set_column_width((-1,0,1))

            col(int): index integer
            >>>table.auto_set_column_width(-1)
            >>>table.auto_set_column_width(0)
            >>>table.auto_set_column_width(1)
        """
        # check for col possibility on iteration
        try:
            iter(col)
        except (TypeError, AttributeError):
            self._autoColumns.append(col)
        else:
            for cell in col:
                self._autoColumns.append(cell)

        self.stale = True

    def _auto_set_column_width(self, col, renderer):
        """ Automagically set width for column.
        """
        cells = [key for key in self._cells if key[1] == col]

        # find max width
        width = 0
        for cell in cells:
            c = self._cells[cell]
            width = max(c.get_required_width(renderer), width)

        # Now set the widths
        for cell in cells:
            self._cells[cell].set_width(width)

    def auto_set_font_size(self, value=True):
        """ Automatically set font size. """
        self._autoFontsize = value
        self.stale = True

    def _auto_set_font_size(self, renderer):

        if len(self._cells) == 0:
            return
        fontsize = list(six.itervalues(self._cells))[0].get_fontsize()
        cells = []
        for key, cell in six.iteritems(self._cells):
            # ignore auto-sized columns
            if key[1] in self._autoColumns:
                continue
            size = cell.auto_set_font_size(renderer)
            fontsize = min(fontsize, size)
            cells.append(cell)

        # now set all fontsizes equal
        for cell in six.itervalues(self._cells):
            cell.set_fontsize(fontsize)

    def scale(self, xscale, yscale):
        """ Scale column widths by xscale and row heights by yscale. """
        for c in six.itervalues(self._cells):
            c.set_width(c.get_width() * xscale)
            c.set_height(c.get_height() * yscale)

    def set_fontsize(self, size):
        """
        Set the fontsize of the cell text

        ACCEPTS: a float in points
        """

        for cell in six.itervalues(self._cells):
            cell.set_fontsize(size)
        self.stale = True

    def _offset(self, ox, oy):
        'Move all the artists by ox,oy (axes coords)'

        for c in six.itervalues(self._cells):
            x, y = c.get_x(), c.get_y()
            c.set_x(x + ox)
            c.set_y(y + oy)

    def _update_positions(self, renderer):
        # called from renderer to allow more precise estimates of
        # widths and heights with get_window_extent

        # Do any auto width setting
        for col in self._autoColumns:
            self._auto_set_column_width(col, renderer)

        if self._autoFontsize:
            self._auto_set_font_size(renderer)

        # Align all the cells
        self._do_cell_alignment()

        bbox = self._get_grid_bbox(renderer)
        l, b, w, h = bbox.bounds

        if self._bbox is not None:
            # Position according to bbox
            rl, rb, rw, rh = self._bbox
            self.scale(rw / w, rh / h)
            ox = rl - l
            oy = rb - b
            self._do_cell_alignment()
        else:
            # Position using loc
            (BEST, UR, UL, LL, LR, CL, CR, LC, UC, C,
             TR, TL, BL, BR, R, L, T, B) = xrange(len(self.codes))
            # defaults for center
            ox = (0.5 - w / 2) - l
            oy = (0.5 - h / 2) - b
            if self._loc in (UL, LL, CL):   # left
                ox = self.AXESPAD - l
            if self._loc in (BEST, UR, LR, R, CR):  # right
                ox = 1 - (l + w + self.AXESPAD)
            if self._loc in (BEST, UR, UL, UC):     # upper
                oy = 1 - (b + h + self.AXESPAD)
            if self._loc in (LL, LR, LC):           # lower
                oy = self.AXESPAD - b
            if self._loc in (LC, UC, C):            # center x
                ox = (0.5 - w / 2) - l
            if self._loc in (CL, CR, C):            # center y
                oy = (0.5 - h / 2) - b

            if self._loc in (TL, BL, L):            # out left
                ox = - (l + w)
            if self._loc in (TR, BR, R):            # out right
                ox = 1.0 - l
            if self._loc in (TR, TL, T):            # out top
                oy = 1.0 - b
            if self._loc in (BL, BR, B):           # out bottom
                oy = - (b + h)

        self._offset(ox, oy)

    def get_celld(self):
        'return a dict of cells in the table'
        return self._cells


def table(ax,
          cellText=None, cellColours=None,
          cellLoc='right', colWidths=None,
          rowLabels=None, rowColours=None, rowLoc='left',
          colLabels=None, colColours=None, colLoc='center',
          loc='bottom', bbox=None, edges='closed',
          **kwargs):
    """
    TABLE(cellText=None, cellColours=None,
          cellLoc='right', colWidths=None,
          rowLabels=None, rowColours=None, rowLoc='left',
          colLabels=None, colColours=None, colLoc='center',
          loc='bottom', bbox=None, edges='closed')

    Factory function to generate a Table instance.

    Thanks to John Gill for providing the class and table.
    """

    if cellColours is None and cellText is None:
        raise ValueError('At least one argument from "cellColours" or '
                         '"cellText" must be provided to create a table.')

    # Check we have some cellText
    if cellText is None:
        # assume just colours are needed
        rows = len(cellColours)
        cols = len(cellColours[0])
        cellText = [[''] * cols] * rows

    rows = len(cellText)
    cols = len(cellText[0])
    for row in cellText:
        if len(row) != cols:
            raise ValueError("Each row in 'cellText' must have {} columns"
                             .format(cols))

    if cellColours is not None:
        if len(cellColours) != rows:
            raise ValueError("'cellColours' must have {} rows".format(rows))
        for row in cellColours:
            if len(row) != cols:
                raise ValueError("Each row in 'cellColours' must have {} "
                                 "columns".format(cols))
    else:
        cellColours = ['w' * cols] * rows

    # Set colwidths if not given
    if colWidths is None:
        colWidths = [1.0 / cols] * cols

    # Fill in missing information for column
    # and row labels
    rowLabelWidth = 0
    if rowLabels is None:
        if rowColours is not None:
            rowLabels = [''] * rows
            rowLabelWidth = colWidths[0]
    elif rowColours is None:
        rowColours = 'w' * rows

    if rowLabels is not None:
        if len(rowLabels) != rows:
            raise ValueError("'rowLabels' must be of length {0}".format(rows))

    # If we have column labels, need to shift
    # the text and colour arrays down 1 row
    offset = 1
    if colLabels is None:
        if colColours is not None:
            colLabels = [''] * cols
        else:
            offset = 0
    elif colColours is None:
        colColours = 'w' * cols

    # Set up cell colours if not given
    if cellColours is None:
        cellColours = ['w' * cols] * rows

    # Now create the table
    table = Table(ax, loc, bbox, **kwargs)
    table.edges = edges
    height = table._approx_text_height()

    # Add the cells
    for row in xrange(rows):
        for col in xrange(cols):
            table.add_cell(row + offset, col,
                           width=colWidths[col], height=height,
                           text=cellText[row][col],
                           facecolor=cellColours[row][col],
                           loc=cellLoc)
    # Do column labels
    if colLabels is not None:
        for col in xrange(cols):
            table.add_cell(0, col,
                           width=colWidths[col], height=height,
                           text=colLabels[col], facecolor=colColours[col],
                           loc=colLoc)

    # Do row labels
    if rowLabels is not None:
        for row in xrange(rows):
            table.add_cell(row + offset, -1,
                           width=rowLabelWidth or 1e-15, height=height,
                           text=rowLabels[row], facecolor=rowColours[row],
                           loc=rowLoc)
        if rowLabelWidth == 0:
            table.auto_set_column_width(-1)

    ax.add_table(table)
    return table


docstring.interpd.update(Table=artist.kwdoc(Table))
