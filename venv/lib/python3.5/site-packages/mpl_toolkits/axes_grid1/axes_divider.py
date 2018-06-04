"""
The axes_divider module provides helper classes to adjust the positions of
multiple axes at drawing time.

 Divider: this is the class that is used to calculate the axes
    position. It divides the given rectangular area into several sub
    rectangles. You initialize the divider by setting the horizontal
    and vertical lists of sizes that the division will be based on. You
    then use the new_locator method, whose return value is a callable
    object that can be used to set the axes_locator of the axes.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import map

import matplotlib.transforms as mtransforms

from matplotlib.axes import SubplotBase

from . import axes_size as Size


class Divider(object):
    """
    This class calculates the axes position. It
    divides the given rectangular area into several
    sub-rectangles. You initialize the divider by setting the
    horizontal and vertical lists of sizes
    (:mod:`mpl_toolkits.axes_grid.axes_size`) that the division will
    be based on. You then use the new_locator method to create a
    callable object that can be used as the axes_locator of the
    axes.
    """

    def __init__(self, fig, pos, horizontal, vertical,
                 aspect=None, anchor="C"):
        """
        Parameters
        ----------
        fig : Figure
        pos : tuple of 4 floats
            position of the rectangle that will be divided
        horizontal : list of :mod:`~mpl_toolkits.axes_grid.axes_size`
            sizes for horizontal division
        vertical : list of :mod:`~mpl_toolkits.axes_grid.axes_size`
            sizes for vertical division
        aspect : bool
            if True, the overall rectangular area is reduced
            so that the relative part of the horizontal and
            vertical scales have the same scale.
        anchor : {'C', 'SW', 'S', 'SE', 'E', 'NE', 'N', 'NW', 'W'}
            placement of the reduced rectangle when *aspect* is True
        """

        self._fig = fig
        self._pos = pos
        self._horizontal = horizontal
        self._vertical = vertical
        self._anchor = anchor
        self._aspect = aspect
        self._xrefindex = 0
        self._yrefindex = 0
        self._locator = None

    def get_horizontal_sizes(self, renderer):
        return [s.get_size(renderer) for s in self.get_horizontal()]

    def get_vertical_sizes(self, renderer):
        return [s.get_size(renderer) for s in self.get_vertical()]

    def get_vsize_hsize(self):

        from .axes_size import AddList

        vsize = AddList(self.get_vertical())
        hsize = AddList(self.get_horizontal())

        return vsize, hsize

    @staticmethod
    def _calc_k(l, total_size):

        rs_sum, as_sum = 0., 0.

        for _rs, _as in l:
            rs_sum += _rs
            as_sum += _as

        if rs_sum != 0.:
            k = (total_size - as_sum) / rs_sum
            return k
        else:
            return 0.

    @staticmethod
    def _calc_offsets(l, k):

        offsets = [0.]

        #for s in l:
        for _rs, _as in l:
            #_rs, _as = s.get_size(renderer)
            offsets.append(offsets[-1] + _rs*k + _as)

        return offsets

    def set_position(self, pos):
        """
        set the position of the rectangle.

        Parameters
        ----------
        pos : tuple of 4 floats
            position of the rectangle that will be divided
        """
        self._pos = pos

    def get_position(self):
        "return the position of the rectangle."
        return self._pos

    def set_anchor(self, anchor):
        """
        Parameters
        ----------
        anchor : {'C', 'SW', 'S', 'SE', 'E', 'NE', 'N', 'NW', 'W'}
            anchor position

          =====  ============
          value  description
          =====  ============
          'C'    Center
          'SW'   bottom left
          'S'    bottom
          'SE'   bottom right
          'E'    right
          'NE'   top right
          'N'    top
          'NW'   top left
          'W'    left
          =====  ============

        """
        if anchor in mtransforms.Bbox.coefs or len(anchor) == 2:
            self._anchor = anchor
        else:
            raise ValueError('argument must be among %s' %
                             ', '.join(mtransforms.BBox.coefs))

    def get_anchor(self):
        "return the anchor"
        return self._anchor

    def set_horizontal(self, h):
        """
        Parameters
        ----------
        h : list of :mod:`~mpl_toolkits.axes_grid.axes_size`
            sizes for horizontal division
        """
        self._horizontal = h

    def get_horizontal(self):
        "return horizontal sizes"
        return self._horizontal

    def set_vertical(self, v):
        """
        Parameters
        ----------
        v : list of :mod:`~mpl_toolkits.axes_grid.axes_size`
            sizes for vertical division
        """
        self._vertical = v

    def get_vertical(self):
        "return vertical sizes"
        return self._vertical

    def set_aspect(self, aspect=False):
        """
        Parameters
        ----------
        aspect : bool
        """
        self._aspect = aspect

    def get_aspect(self):
        "return aspect"
        return self._aspect

    def set_locator(self, _locator):
        self._locator = _locator

    def get_locator(self):
        return self._locator

    def get_position_runtime(self, ax, renderer):
        if self._locator is None:
            return self.get_position()
        else:
            return self._locator(ax, renderer).bounds

    def locate(self, nx, ny, nx1=None, ny1=None, axes=None, renderer=None):
        """
        Parameters
        ----------
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        axes
        renderer
        """

        figW, figH = self._fig.get_size_inches()
        x, y, w, h = self.get_position_runtime(axes, renderer)

        hsizes = self.get_horizontal_sizes(renderer)
        vsizes = self.get_vertical_sizes(renderer)
        k_h = self._calc_k(hsizes, figW*w)
        k_v = self._calc_k(vsizes, figH*h)

        if self.get_aspect():
            k = min(k_h, k_v)
            ox = self._calc_offsets(hsizes, k)
            oy = self._calc_offsets(vsizes, k)

            ww = (ox[-1] - ox[0])/figW
            hh = (oy[-1] - oy[0])/figH
            pb = mtransforms.Bbox.from_bounds(x, y, w, h)
            pb1 = mtransforms.Bbox.from_bounds(x, y, ww, hh)
            pb1_anchored = pb1.anchored(self.get_anchor(), pb)
            x0, y0 = pb1_anchored.x0, pb1_anchored.y0

        else:
            ox = self._calc_offsets(hsizes, k_h)
            oy = self._calc_offsets(vsizes, k_v)
            x0, y0 = x, y

        if nx1 is None:
            nx1 = nx+1
        if ny1 is None:
            ny1 = ny+1

        x1, w1 = x0 + ox[nx]/figW, (ox[nx1] - ox[nx])/figW
        y1, h1 = y0 + oy[ny]/figH, (oy[ny1] - oy[ny])/figH

        return mtransforms.Bbox.from_bounds(x1, y1, w1, h1)

    def new_locator(self, nx, ny, nx1=None, ny1=None):
        """
        Returns a new locator
        (:class:`mpl_toolkits.axes_grid.axes_divider.AxesLocator`) for
        specified cell.

        Parameters
        ----------
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        """
        return AxesLocator(self, nx, ny, nx1, ny1)

    def append_size(self, position, size):

        if position == "left":
            self._horizontal.insert(0, size)
            self._xrefindex += 1
        elif position == "right":
            self._horizontal.append(size)
        elif position == "bottom":
            self._vertical.insert(0, size)
            self._yrefindex += 1
        elif position == "top":
            self._vertical.append(size)
        else:
            raise ValueError("the position must be one of left," +
                             " right, bottom, or top")

    def add_auto_adjustable_area(self,
                                 use_axes, pad=0.1,
                                 adjust_dirs=None,
                                 ):
        if adjust_dirs is None:
            adjust_dirs = ["left", "right", "bottom", "top"]
        from .axes_size import Padded, SizeFromFunc, GetExtentHelper
        for d in adjust_dirs:
            helper = GetExtentHelper(use_axes, d)
            size = SizeFromFunc(helper)
            padded_size = Padded(size, pad)  # pad in inch
            self.append_size(d, padded_size)


class AxesLocator(object):
    """
    A simple callable object, initialized with AxesDivider class,
    returns the position and size of the given cell.
    """
    def __init__(self, axes_divider, nx, ny, nx1=None, ny1=None):
        """
        Parameters
        ----------
        axes_divider : AxesDivider
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        """
        self._axes_divider = axes_divider

        _xrefindex = axes_divider._xrefindex
        _yrefindex = axes_divider._yrefindex

        self._nx, self._ny = nx - _xrefindex, ny - _yrefindex

        if nx1 is None:
            nx1 = nx+1
        if ny1 is None:
            ny1 = ny+1

        self._nx1 = nx1 - _xrefindex
        self._ny1 = ny1 - _yrefindex

    def __call__(self, axes, renderer):

        _xrefindex = self._axes_divider._xrefindex
        _yrefindex = self._axes_divider._yrefindex

        return self._axes_divider.locate(self._nx + _xrefindex,
                                         self._ny + _yrefindex,
                                         self._nx1 + _xrefindex,
                                         self._ny1 + _yrefindex,
                                         axes,
                                         renderer)

    def get_subplotspec(self):
        if hasattr(self._axes_divider, "get_subplotspec"):
            return self._axes_divider.get_subplotspec()
        else:
            return None


from matplotlib.gridspec import SubplotSpec, GridSpec


class SubplotDivider(Divider):
    """
    The Divider class whose rectangle area is specified as a subplot geometry.
    """

    def __init__(self, fig, *args, **kwargs):
        """
        Parameters
        ----------
        fig : :class:`matplotlib.figure.Figure`
        args : tuple (*numRows*, *numCols*, *plotNum*)
            The array of subplots in the figure has dimensions *numRows*,
            *numCols*, and *plotNum* is the number of the subplot
            being created.  *plotNum* starts at 1 in the upper left
            corner and increases to the right.

            If *numRows* <= *numCols* <= *plotNum* < 10, *args* can be the
            decimal integer *numRows* * 100 + *numCols* * 10 + *plotNum*.
        """

        self.figure = fig

        if len(args) == 1:
            if isinstance(args[0], SubplotSpec):
                self._subplotspec = args[0]
            else:
                try:
                    s = str(int(args[0]))
                    rows, cols, num = map(int, s)
                except ValueError:
                    raise ValueError(
                        'Single argument to subplot must be a 3-digit integer')
                self._subplotspec = GridSpec(rows, cols)[num-1]
                # num - 1 for converting from MATLAB to python indexing
        elif len(args) == 3:
            rows, cols, num = args
            rows = int(rows)
            cols = int(cols)
            if isinstance(num, tuple) and len(num) == 2:
                num = [int(n) for n in num]
                self._subplotspec = GridSpec(rows, cols)[num[0]-1:num[1]]
            else:
                self._subplotspec = GridSpec(rows, cols)[int(num)-1]
                # num - 1 for converting from MATLAB to python indexing
        else:
            raise ValueError('Illegal argument(s) to subplot: %s' % (args,))

        # total = rows*cols
        # num -= 1    # convert from matlab to python indexing
        #             # i.e., num in range(0,total)
        # if num >= total:
        #     raise ValueError( 'Subplot number exceeds total subplots')
        # self._rows = rows
        # self._cols = cols
        # self._num = num

        # self.update_params()

        # sets self.fixbox
        self.update_params()

        pos = self.figbox.bounds

        horizontal = kwargs.pop("horizontal", [])
        vertical = kwargs.pop("vertical", [])
        aspect = kwargs.pop("aspect", None)
        anchor = kwargs.pop("anchor", "C")

        if kwargs:
            raise Exception("")

        Divider.__init__(self, fig, pos, horizontal, vertical,
                         aspect=aspect, anchor=anchor)

    def get_position(self):
        "return the bounds of the subplot box"

        self.update_params()  # update self.figbox
        return self.figbox.bounds

    # def update_params(self):
    #     'update the subplot position from fig.subplotpars'

    #     rows = self._rows
    #     cols = self._cols
    #     num = self._num

    #     pars = self.figure.subplotpars
    #     left = pars.left
    #     right = pars.right
    #     bottom = pars.bottom
    #     top = pars.top
    #     wspace = pars.wspace
    #     hspace = pars.hspace
    #     totWidth = right-left
    #     totHeight = top-bottom

    #     figH = totHeight/(rows + hspace*(rows-1))
    #     sepH = hspace*figH

    #     figW = totWidth/(cols + wspace*(cols-1))
    #     sepW = wspace*figW

    #     rowNum, colNum =  divmod(num, cols)

    #     figBottom = top - (rowNum+1)*figH - rowNum*sepH
    #     figLeft = left + colNum*(figW + sepW)

    #     self.figbox = mtransforms.Bbox.from_bounds(figLeft, figBottom,
    #                                                figW, figH)

    def update_params(self):
        'update the subplot position from fig.subplotpars'

        self.figbox = self.get_subplotspec().get_position(self.figure)

    def get_geometry(self):
        'get the subplot geometry, e.g., 2,2,3'
        rows, cols, num1, num2 = self.get_subplotspec().get_geometry()
        return rows, cols, num1+1  # for compatibility

    # COVERAGE NOTE: Never used internally or from examples
    def change_geometry(self, numrows, numcols, num):
        'change subplot geometry, e.g., from 1,1,1 to 2,2,3'
        self._subplotspec = GridSpec(numrows, numcols)[num-1]
        self.update_params()
        self.set_position(self.figbox)

    def get_subplotspec(self):
        'get the SubplotSpec instance'
        return self._subplotspec

    def set_subplotspec(self, subplotspec):
        'set the SubplotSpec instance'
        self._subplotspec = subplotspec


class AxesDivider(Divider):
    """
    Divider based on the pre-existing axes.
    """

    def __init__(self, axes, xref=None, yref=None):
        """
        Parameters
        ----------
        axes : :class:`~matplotlib.axes.Axes`
        xref
        yref
        """
        self._axes = axes
        if xref is None:
            self._xref = Size.AxesX(axes)
        else:
            self._xref = xref
        if yref is None:
            self._yref = Size.AxesY(axes)
        else:
            self._yref = yref

        Divider.__init__(self, fig=axes.get_figure(), pos=None,
                         horizontal=[self._xref], vertical=[self._yref],
                         aspect=None, anchor="C")

    def _get_new_axes(self, **kwargs):
        axes = self._axes

        axes_class = kwargs.pop("axes_class", None)

        if axes_class is None:
            if isinstance(axes, SubplotBase):
                axes_class = axes._axes_class
            else:
                axes_class = type(axes)

        ax = axes_class(axes.get_figure(),
                        axes.get_position(original=True), **kwargs)

        return ax

    def new_horizontal(self, size, pad=None, pack_start=False, **kwargs):
        """
        Add a new axes on the right (or left) side of the main axes.

        Parameters
        ----------
        size : :mod:`~mpl_toolkits.axes_grid.axes_size` or float or string
            A width of the axes. If float or string is given, *from_any*
            function is used to create the size, with *ref_size* set to AxesX
            instance of the current axes.
        pad : :mod:`~mpl_toolkits.axes_grid.axes_size` or float or string
            Pad between the axes. It takes same argument as *size*.
        pack_start : bool
            If False, the new axes is appended at the end
            of the list, i.e., it became the right-most axes. If True, it is
            inserted at the start of the list, and becomes the left-most axes.
        kwargs
            All extra keywords arguments are passed to the created axes.
            If *axes_class* is given, the new axes will be created as an
            instance of the given class. Otherwise, the same class of the
            main axes will be used.
        """

        if pad:
            if not isinstance(pad, Size._Base):
                pad = Size.from_any(pad,
                                    fraction_ref=self._xref)
            if pack_start:
                self._horizontal.insert(0, pad)
                self._xrefindex += 1
            else:
                self._horizontal.append(pad)

        if not isinstance(size, Size._Base):
            size = Size.from_any(size,
                                 fraction_ref=self._xref)

        if pack_start:
            self._horizontal.insert(0, size)
            self._xrefindex += 1
            locator = self.new_locator(nx=0, ny=self._yrefindex)
        else:
            self._horizontal.append(size)
            locator = self.new_locator(nx=len(self._horizontal)-1, ny=self._yrefindex)

        ax = self._get_new_axes(**kwargs)
        ax.set_axes_locator(locator)

        return ax

    def new_vertical(self, size, pad=None, pack_start=False, **kwargs):
        """
        Add a new axes on the top (or bottom) side of the main axes.

        Parameters
        ----------
        size : :mod:`~mpl_toolkits.axes_grid.axes_size` or float or string
            A height of the axes. If float or string is given, *from_any*
            function is used to create the size, with *ref_size* set to AxesX
            instance of the current axes.
        pad : :mod:`~mpl_toolkits.axes_grid.axes_size` or float or string
            Pad between the axes. It takes same argument as *size*.
        pack_start : bool
            If False, the new axes is appended at the end
            of the list, i.e., it became the right-most axes. If True, it is
            inserted at the start of the list, and becomes the left-most axes.
        kwargs
            All extra keywords arguments are passed to the created axes.
            If *axes_class* is given, the new axes will be created as an
            instance of the given class. Otherwise, the same class of the
            main axes will be used.
        """

        if pad:
            if not isinstance(pad, Size._Base):
                pad = Size.from_any(pad,
                                    fraction_ref=self._yref)
            if pack_start:
                self._vertical.insert(0, pad)
                self._yrefindex += 1
            else:
                self._vertical.append(pad)

        if not isinstance(size, Size._Base):
            size = Size.from_any(size,
                                 fraction_ref=self._yref)

        if pack_start:
            self._vertical.insert(0, size)
            self._yrefindex += 1
            locator = self.new_locator(nx=self._xrefindex, ny=0)
        else:
            self._vertical.append(size)
            locator = self.new_locator(nx=self._xrefindex, ny=len(self._vertical)-1)

        ax = self._get_new_axes(**kwargs)
        ax.set_axes_locator(locator)

        return ax

    def append_axes(self, position, size, pad=None, add_to_figure=True,
                    **kwargs):
        """
        create an axes at the given *position* with the same height
        (or width) of the main axes.

         *position*
           ["left"|"right"|"bottom"|"top"]

         *size* and *pad* should be axes_grid.axes_size compatible.
        """

        if position == "left":
            ax = self.new_horizontal(size, pad, pack_start=True, **kwargs)
        elif position == "right":
            ax = self.new_horizontal(size, pad, pack_start=False, **kwargs)
        elif position == "bottom":
            ax = self.new_vertical(size, pad, pack_start=True, **kwargs)
        elif position == "top":
            ax = self.new_vertical(size, pad, pack_start=False, **kwargs)
        else:
            raise ValueError("the position must be one of left," +
                             " right, bottom, or top")

        if add_to_figure:
            self._fig.add_axes(ax)
        return ax

    def get_aspect(self):
        if self._aspect is None:
            aspect = self._axes.get_aspect()
            if aspect == "auto":
                return False
            else:
                return True
        else:
            return self._aspect

    def get_position(self):
        if self._pos is None:
            bbox = self._axes.get_position(original=True)
            return bbox.bounds
        else:
            return self._pos

    def get_anchor(self):
        if self._anchor is None:
            return self._axes.get_anchor()
        else:
            return self._anchor

    def get_subplotspec(self):
        if hasattr(self._axes, "get_subplotspec"):
            return self._axes.get_subplotspec()
        else:
            return None


class HBoxDivider(SubplotDivider):

    def __init__(self, fig, *args, **kwargs):
        SubplotDivider.__init__(self, fig, *args, **kwargs)

    @staticmethod
    def _determine_karray(equivalent_sizes, appended_sizes,
                          max_equivalent_size,
                          total_appended_size):

        n = len(equivalent_sizes)
        import numpy as np
        A = np.mat(np.zeros((n+1, n+1), dtype="d"))
        B = np.zeros((n+1), dtype="d")
        # AxK = B

        # populated A
        for i, (r, a) in enumerate(equivalent_sizes):
            A[i, i] = r
            A[i, -1] = -1
            B[i] = -a
        A[-1, :-1] = [r for r, a in appended_sizes]
        B[-1] = total_appended_size - sum([a for rs, a in appended_sizes])

        karray_H = (A.I*np.mat(B).T).A1
        karray = karray_H[:-1]
        H = karray_H[-1]

        if H > max_equivalent_size:
            karray = ((max_equivalent_size -
                      np.array([a for r, a in equivalent_sizes]))
                      / np.array([r for r, a in equivalent_sizes]))
        return karray

    @staticmethod
    def _calc_offsets(appended_sizes, karray):
        offsets = [0.]

        #for s in l:
        for (r, a), k in zip(appended_sizes, karray):
            offsets.append(offsets[-1] + r*k + a)

        return offsets

    def new_locator(self, nx, nx1=None):
        """
        returns a new locator
        (:class:`mpl_toolkits.axes_grid.axes_divider.AxesLocator`) for
        specified cell.

        Parameters
        ----------
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        """
        return AxesLocator(self, nx, 0, nx1, None)

    def _locate(self, x, y, w, h,
                y_equivalent_sizes, x_appended_sizes,
                figW, figH):
        """
        Parameters
        ----------
        x
        y
        w
        h
        y_equivalent_sizes
        x_appended_sizes
        figW
        figH
        """

        equivalent_sizes = y_equivalent_sizes
        appended_sizes = x_appended_sizes

        max_equivalent_size = figH*h
        total_appended_size = figW*w
        karray = self._determine_karray(equivalent_sizes, appended_sizes,
                                        max_equivalent_size,
                                        total_appended_size)

        ox = self._calc_offsets(appended_sizes, karray)

        ww = (ox[-1] - ox[0])/figW
        ref_h = equivalent_sizes[0]
        hh = (karray[0]*ref_h[0] + ref_h[1])/figH
        pb = mtransforms.Bbox.from_bounds(x, y, w, h)
        pb1 = mtransforms.Bbox.from_bounds(x, y, ww, hh)
        pb1_anchored = pb1.anchored(self.get_anchor(), pb)
        x0, y0 = pb1_anchored.x0, pb1_anchored.y0

        return x0, y0, ox, hh

    def locate(self, nx, ny, nx1=None, ny1=None, axes=None, renderer=None):
        """
        Parameters
        ----------
        axes_divider : AxesDivider
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        axes
        renderer
        """

        figW, figH = self._fig.get_size_inches()
        x, y, w, h = self.get_position_runtime(axes, renderer)

        y_equivalent_sizes = self.get_vertical_sizes(renderer)
        x_appended_sizes = self.get_horizontal_sizes(renderer)
        x0, y0, ox, hh = self._locate(x, y, w, h,
                                      y_equivalent_sizes, x_appended_sizes,
                                      figW, figH)
        if nx1 is None:
            nx1 = nx+1

        x1, w1 = x0 + ox[nx]/figW, (ox[nx1] - ox[nx])/figW
        y1, h1 = y0, hh

        return mtransforms.Bbox.from_bounds(x1, y1, w1, h1)


class VBoxDivider(HBoxDivider):
    """
    The Divider class whose rectangle area is specified as a subplot geometry.
    """

    def new_locator(self, ny, ny1=None):
        """
        returns a new locator
        (:class:`mpl_toolkits.axes_grid.axes_divider.AxesLocator`) for
        specified cell.

        Parameters
        ----------
        ny, ny1 : int
            Integers specifying the row-position of the
            cell. When *ny1* is None, a single *ny*-th row is
            specified. Otherwise location of rows spanning between *ny*
            to *ny1* (but excluding *ny1*-th row) is specified.
        """
        return AxesLocator(self, 0, ny, None, ny1)

    def locate(self, nx, ny, nx1=None, ny1=None, axes=None, renderer=None):
        """
        Parameters
        ----------
        axes_divider : AxesDivider
        nx, nx1 : int
            Integers specifying the column-position of the
            cell. When *nx1* is None, a single *nx*-th column is
            specified. Otherwise location of columns spanning between *nx*
            to *nx1* (but excluding *nx1*-th column) is specified.
        ny, ny1 : int
            Same as *nx* and *nx1*, but for row positions.
        axes
        renderer
        """

        figW, figH = self._fig.get_size_inches()
        x, y, w, h = self.get_position_runtime(axes, renderer)

        x_equivalent_sizes = self.get_horizontal_sizes(renderer)
        y_appended_sizes = self.get_vertical_sizes(renderer)

        y0, x0, oy, ww = self._locate(y, x, h, w,
                                      x_equivalent_sizes, y_appended_sizes,
                                      figH, figW)
        if ny1 is None:
            ny1 = ny+1

        x1, w1 = x0, ww
        y1, h1 = y0 + oy[ny]/figH, (oy[ny1] - oy[ny])/figH

        return mtransforms.Bbox.from_bounds(x1, y1, w1, h1)


class LocatableAxesBase(object):
    def __init__(self, *kl, **kw):

        self._axes_class.__init__(self, *kl, **kw)

        self._locator = None
        self._locator_renderer = None

    def set_axes_locator(self, locator):
        self._locator = locator

    def get_axes_locator(self):
        return self._locator

    def apply_aspect(self, position=None):

        if self.get_axes_locator() is None:
            self._axes_class.apply_aspect(self, position)
        else:
            pos = self.get_axes_locator()(self, self._locator_renderer)
            self._axes_class.apply_aspect(self, position=pos)

    def draw(self, renderer=None, inframe=False):

        self._locator_renderer = renderer

        self._axes_class.draw(self, renderer, inframe)

    def _make_twin_axes(self, *kl, **kwargs):
        """
        Need to overload so that twinx/twiny will work with
        these axes.
        """
        if 'sharex' in kwargs and 'sharey' in kwargs:
            raise ValueError("Twinned Axes may share only one axis.")
        ax2 = type(self)(self.figure, self.get_position(True), *kl, **kwargs)
        ax2.set_axes_locator(self.get_axes_locator())
        self.figure.add_axes(ax2)
        self.set_adjustable('datalim')
        ax2.set_adjustable('datalim')
        self._twinned_axes.join(self, ax2)
        return ax2

_locatableaxes_classes = {}


def locatable_axes_factory(axes_class):

    new_class = _locatableaxes_classes.get(axes_class)
    if new_class is None:
        new_class = type(str("Locatable%s" % (axes_class.__name__)),
                         (LocatableAxesBase, axes_class),
                         {'_axes_class': axes_class})

        _locatableaxes_classes[axes_class] = new_class

    return new_class

#if hasattr(maxes.Axes, "get_axes_locator"):
#    LocatableAxes = maxes.Axes
#else:


def make_axes_locatable(axes):
    if not hasattr(axes, "set_axes_locator"):
        new_class = locatable_axes_factory(type(axes))
        axes.__class__ = new_class

    divider = AxesDivider(axes)
    locator = divider.new_locator(nx=0, ny=0)
    axes.set_axes_locator(locator)

    return divider


def make_axes_area_auto_adjustable(ax,
                                   use_axes=None, pad=0.1,
                                   adjust_dirs=None):
    if adjust_dirs is None:
        adjust_dirs = ["left", "right", "bottom", "top"]
    divider = make_axes_locatable(ax)

    if use_axes is None:
        use_axes = ax

    divider.add_auto_adjustable_area(use_axes=use_axes, pad=pad,
                                     adjust_dirs=adjust_dirs)

#from matplotlib.axes import Axes
from .mpl_axes import Axes
LocatableAxes = locatable_axes_factory(Axes)
