"""
This module provides the routine to adjust subplot layouts so that there are
no overlapping axes or axes decorations.  All axes decorations are dealt with
(labels, ticks, titles, ticklabels) and some dependent artists are also dealt
with (colorbar, suptitle, legend).

Layout is done via :meth:`~matplotlib.gridspec`, with one constraint per
gridspec, so it is possible to have overlapping axes if the gridspecs
overlap (i.e. using :meth:`~matplotlib.gridspec.GridSpecFromSubplotSpec`).
Axes placed using ``figure.subplots()`` or ``figure.add_subplots()`` will
participate in the layout.  Axes manually placed via ``figure.add_axes()``
will not.

See Tutorial: :doc:`/tutorials/intermediate/constrainedlayout_guide`

"""

# Development Notes:

# What gets a layoutbox:
#  - figure
#    - gridspec
#      - subplotspec
#        EITHER:
#         - axes + pos for the axes (i.e. the total area taken by axis and
#            the actual "position" argument that needs to be sent to
#             ax.set_position.)
#           - The axes layout box will also encomapss the legend, and that is
#             how legends get included (axes legeneds, not figure legends)
#         - colorbars are sibblings of the axes if they are single-axes
#           colorbars
#        OR:
#         - a gridspec can be inside a subplotspec.
#           - subplotspec
#           EITHER:
#            - axes...
#           OR:
#            - gridspec... with arbitrary nesting...
#      - colorbars are siblings of the subplotspecs if they are multi-axes
#        colorbars.
#   - suptitle:
#      - right now suptitles are just stacked atop everything else in figure.
#        Could imagine suptitles being gridspec suptitles, but not implimented
#
#   Todo:    AnchoredOffsetbox connected to gridspecs or axes.  This would
#        be more general way to add extra-axes annotations.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import logging
import warnings

from matplotlib.legend import Legend
import matplotlib.transforms as transforms
import matplotlib._layoutbox as layoutbox

_log = logging.getLogger(__name__)


def get_axall_tightbbox(ax, renderer):
    '''
    Get the tight_bbox of the axis ax, and any dependent decorations, like
    a `Legend` instance.
    '''

    # main bbox of the axis....
    bbox = ax.get_tightbbox(renderer=renderer)
    # now add the possibility of the legend...
    for child in ax.get_children():
        if isinstance(child, Legend):
            bboxn = child._legend_box.get_window_extent(renderer)
            bbox = transforms.Bbox.union([bbox, bboxn])
        # add other children here....
    return bbox


def in_same_column(colnum0min, colnum0max, colnumCmin, colnumCmax):
    if colnum0min >= colnumCmin and colnum0min <= colnumCmax:
        return True
    if colnum0max >= colnumCmin and colnum0max <= colnumCmax:
        return True
    return False


def in_same_row(rownum0min, rownum0max, rownumCmin, rownumCmax):
    if rownum0min >= rownumCmin and rownum0min <= rownumCmax:
        return True
    if rownum0max >= rownumCmin and rownum0max <= rownumCmax:
        return True
    return False


######################################################
def do_constrained_layout(fig, renderer, h_pad, w_pad,
        hspace=None, wspace=None):

    """
    Do the constrained_layout.  Called at draw time in
     ``figure.constrained_layout()``

    Parameters
    ----------


    fig: Figure
      is the ``figure`` instance to do the layout in.

    renderer: Renderer
      the renderer to use.

     h_pad, w_pad : float
       are in figure-normalized units, and are a padding around the axes
       elements.

     hspace, wspace : float
        are in fractions of the subplot sizes.

    """

    '''  Steps:

    1. get a list of unique gridspecs in this figure.  Each gridspec will be
    constrained separately.
    2. Check for gaps in the gridspecs.  i.e. if not every axes slot in the
    gridspec has been filled.  If empty, add a ghost axis that is made so
    that it cannot be seen (though visible=True).  This is needed to make
    a blank spot in the layout.
    3. Compare the tight_bbox of each axes to its `position`, and assume that
    the difference is the space needed by the elements around the edge of
    the axes (decorations) like the title, ticklabels, x-labels, etc.  This
    can include legends who overspill the axes boundaries.
    4. Constrain gridspec elements to line up:
        a) if colnum0 neq colnumC, the two subplotspecs are stacked next to
        each other, with the appropriate order.
        b) if colnum0 == columnC line up the left or right side of the
        _poslayoutbox (depending if it is the min or max num that is equal).
        c) do the same for rows...
    5. The above doesn't constrain relative sizes of the _poslayoutboxes at
    all, and indeed zero-size is a solution that the solver often finds more
    convenient than expanding the sizes.  Right now the solution is to compare
    subplotspec sizes (i.e. drowsC and drows0) and constrain the larger
    _poslayoutbox to be larger than the ratio of the sizes.  i.e. if drows0 >
    drowsC,  then ax._poslayoutbox > axc._poslayoutbox * drowsC / drows0. This
    works fine *if* the decorations are similar between the axes.  If the
    larger subplotspec has much larger axes decorations, then the constraint
    above is incorrect.

    We need the greater than in the above, in general, rather than an equals
    sign.  Consider the case of the left column having 2 rows, and the right
    column having 1 row.  We want the top and bottom of the _poslayoutboxes to
    line up. So that means if there are decorations on the left column axes
    they will be smaller than half as large as the right hand axis.

    This can break down if the decoration size for the right hand axis (the
    margins) is very large.  There must be a math way to check for this case.

    '''

    invTransFig = fig.transFigure.inverted().transform_bbox

    # list of unique gridspecs that contain child axes:
    gss = set([])
    for ax in fig.axes:
        if hasattr(ax, 'get_subplotspec'):
            gs = ax.get_subplotspec().get_gridspec()
            if gs._layoutbox is not None:
                gss.add(gs)
    if len(gss) == 0:
        warnings.warn('There are no gridspecs with layoutboxes. '
                      'Possibly did not call parent GridSpec with the figure= '
                      'keyword')

    # check for unoccupied gridspec slots and make ghost axes for these
    # slots...  Do for each gs separately.  This is a pretty big kludge
    # but shoudn't have too much ill effect.  The worst is that
    # someone querrying the figure will wonder why there are more
    # axes than they thought.
    if fig._layoutbox.constrained_layout_called < 1:
        for gs in gss:
            nrows, ncols = gs.get_geometry()
            hassubplotspec = np.zeros(nrows * ncols, dtype=bool)
            axs = []
            for ax in fig.axes:
                if (hasattr(ax, 'get_subplotspec')
                        and ax._layoutbox is not None
                        and ax.get_subplotspec().get_gridspec() == gs):
                    axs += [ax]
            for ax in axs:
                ss0 = ax.get_subplotspec()
                if ss0.num2 is None:
                    ss0.num2 = ss0.num1
                hassubplotspec[ss0.num1:(ss0.num2 + 1)] = True
            for nn, hss in enumerate(hassubplotspec):
                if not hss:
                    # this gridspec slot doesn't have an axis so we
                    # make a "ghost".
                    ax = fig.add_subplot(gs[nn])
                    ax.set_frame_on(False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_facecolor((1, 0, 0, 0))

    # for each axes, make a margin between the *pos* layoutbox and the
    # *axes* layoutbox be a minimum size that can accomodate the
    # decorations on the axis.
    for ax in fig.axes:
        _log.debug(ax._layoutbox)
        if ax._layoutbox is not None:
            pos = ax.get_position(original=True)
            tightbbox = get_axall_tightbbox(ax, renderer)
            bbox = invTransFig(tightbbox)
            # use stored h_pad if it exists
            h_padt = ax._poslayoutbox.h_pad
            if h_padt is None:
                h_padt = h_pad
            w_padt = ax._poslayoutbox.w_pad
            if w_padt is None:
                w_padt = w_pad
            ax._poslayoutbox.edit_left_margin_min(-bbox.x0 +
                    pos.x0 + w_padt)
            ax._poslayoutbox.edit_right_margin_min(bbox.x1 -
                    pos.x1 + w_padt)
            ax._poslayoutbox.edit_bottom_margin_min(
                    -bbox.y0 + pos.y0 + h_padt)
            ax._poslayoutbox.edit_top_margin_min(bbox.y1-pos.y1+h_padt)
            _log.debug('left %f', (-bbox.x0 + pos.x0 + w_pad))
            _log.debug('right %f', (bbox.x1 - pos.x1 + w_pad))
            _log.debug('bottom %f', (-bbox.y0 + pos.y0 + h_padt))
            # Sometimes its possible for the solver to collapse
            # rather than expand axes, so they all have zero height
            # or width.  This stops that...  It *should* have been
            # taken into account w/ pref_width...
            if fig._layoutbox.constrained_layout_called < 1:
                ax._poslayoutbox.constrain_height_min(20, strength='weak')
                ax._poslayoutbox.constrain_width_min(20, strength='weak')
                ax._layoutbox.constrain_height_min(20, strength='weak')
                ax._layoutbox.constrain_width_min(20, strength='weak')
                ax._poslayoutbox.constrain_top_margin(0, strength='weak')
                ax._poslayoutbox.constrain_bottom_margin(0,
                        strength='weak')
                ax._poslayoutbox.constrain_right_margin(0, strength='weak')
                ax._poslayoutbox.constrain_left_margin(0, strength='weak')

    # do layout for suptitle.
    if fig._suptitle is not None:
        sup = fig._suptitle
        bbox = invTransFig(sup.get_window_extent(renderer=renderer))
        height = bbox.y1 - bbox.y0
        sup._layoutbox.edit_height(height+h_pad)

    # OK, the above lines up ax._poslayoutbox with ax._layoutbox
    # now we need to
    #   1) arrange the subplotspecs.  We do it at this level because
    #      the subplotspecs are meant to contain other dependent axes
    #      like colorbars or legends.
    #   2) line up the right and left side of the ax._poslayoutbox
    #      that have the same subplotspec maxes.

    if fig._layoutbox.constrained_layout_called < 1:

        # arrange the subplotspecs...  This is all done relative to each
        # other.  Some subplotspecs conatain axes, and others contain gridspecs
        # the ones that contain gridspecs are a set proportion of their
        # parent gridspec.  The ones that contain axes are not so constrained.
        figlb = fig._layoutbox
        for child in figlb.children:
            if child._is_gridspec_layoutbox():
                # farm the gridspec layout out.
                #
                # This routine makes all the subplot spec containers
                # have the correct arrangement.  It just stacks the
                # subplot layoutboxes in the correct order...
                arange_subplotspecs(child, hspace=hspace, wspace=wspace)

        # - Align right/left and bottom/top spines of appropriate subplots.
        # - Compare size of subplotspec including height and width ratios
        #   and make sure that the axes spines are at least as large
        #   as they should be.
        for gs in gss:
            # for each gridspec...
            nrows, ncols = gs.get_geometry()
            width_ratios = gs.get_width_ratios()
            height_ratios = gs.get_height_ratios()
            if width_ratios is None:
                width_ratios = np.ones(ncols)
            if height_ratios is None:
                height_ratios = np.ones(nrows)

            # get axes in this gridspec....
            axs = []
            for ax in fig.axes:
                if (hasattr(ax, 'get_subplotspec')
                        and ax._layoutbox is not None):
                    if ax.get_subplotspec().get_gridspec() == gs:
                        axs += [ax]
            rownummin = np.zeros(len(axs), dtype=np.int8)
            rownummax = np.zeros(len(axs), dtype=np.int8)
            colnummin = np.zeros(len(axs), dtype=np.int8)
            colnummax = np.zeros(len(axs), dtype=np.int8)
            width = np.zeros(len(axs))
            height = np.zeros(len(axs))

            for n, ax in enumerate(axs):
                ss0 = ax.get_subplotspec()
                if ss0.num2 is None:
                    ss0.num2 = ss0.num1
                rownummin[n], colnummin[n] = divmod(ss0.num1, ncols)
                rownummax[n], colnummax[n] = divmod(ss0.num2, ncols)
                width[n] = np.sum(
                        width_ratios[colnummin[n]:(colnummax[n] + 1)])
                height[n] = np.sum(
                        height_ratios[rownummin[n]:(rownummax[n] + 1)])

            for nn, ax in enumerate(axs[:-1]):
                ss0 = ax.get_subplotspec()

                # now compare ax to all the axs:
                #
                # If the subplotspecs have the same colnumXmax, then line
                # up their right sides.  If they have the same min, then
                # line up their left sides (and vertical equivalents).
                rownum0min, colnum0min = rownummin[nn], colnummin[nn]
                rownum0max, colnum0max = rownummax[nn], colnummax[nn]
                width0, height0 = width[nn], height[nn]
                alignleft = False
                alignright = False
                alignbot = False
                aligntop = False
                alignheight = False
                alignwidth = False
                for mm in range(nn+1, len(axs)):
                    axc = axs[mm]
                    rownumCmin, colnumCmin = rownummin[mm], colnummin[mm]
                    rownumCmax, colnumCmax = rownummax[mm], colnummax[mm]
                    widthC, heightC = width[mm], height[mm]
                    # Horizontally align axes spines if they have the
                    # same min or max:
                    if not alignleft and colnum0min == colnumCmin:
                        # we want the _poslayoutboxes to line up on left
                        # side of the axes spines...
                        layoutbox.align([ax._poslayoutbox,
                                         axc._poslayoutbox],
                                        'left')
                        alignleft = True

                    if not alignright and colnum0max == colnumCmax:
                        # line up right sides of _poslayoutbox
                        layoutbox.align([ax._poslayoutbox,
                                         axc._poslayoutbox],
                                        'right')
                        alignright = True
                    # Vertically align axes spines if they have the
                    # same min or max:
                    if not aligntop and rownum0min == rownumCmin:
                        # line up top of _poslayoutbox
                        _log.debug('rownum0min == rownumCmin')
                        layoutbox.align([ax._poslayoutbox, axc._poslayoutbox],
                                        'top')
                        aligntop = True

                    if not alignbot and rownum0max == rownumCmax:
                        # line up bottom of _poslayoutbox
                        _log.debug('rownum0max == rownumCmax')
                        layoutbox.align([ax._poslayoutbox, axc._poslayoutbox],
                                        'bottom')
                        alignbot = True
                    ###########
                    # Now we make the widths and heights of position boxes
                    # similar. (i.e the spine locations)
                    # This allows vertically stacked subplots to have
                    # different sizes if they occupy different amounts
                    # of the gridspec:  i.e.
                    # gs = gridspec.GridSpec(3,1)
                    # ax1 = gs[0,:]
                    # ax2 = gs[1:,:]
                    # then drows0 = 1, and drowsC = 2, and ax2
                    # should be at least twice as large as ax1.
                    # But it can be more than twice as large because
                    # it needs less room for the labeling.
                    #
                    # For height, this only needs to be done if the
                    # subplots share a column.  For width if they
                    # share a row.

                    drowsC = (rownumCmax - rownumCmin + 1)
                    drows0 = (rownum0max - rownum0min + 1)
                    dcolsC = (colnumCmax - colnumCmin + 1)
                    dcols0 = (colnum0max - colnum0min + 1)

                    if not alignheight and drows0 == drowsC:
                        ax._poslayoutbox.constrain_height(
                                axc._poslayoutbox.height * height0 / heightC)
                        alignheight = True
                    elif in_same_column(colnum0min, colnum0max,
                            colnumCmin, colnumCmax):
                        if height0 > heightC:
                            ax._poslayoutbox.constrain_height_min(
                                axc._poslayoutbox.height * height0 / heightC)
                            # these constraints stop the smaller axes from
                            # being allowed to go to zero height...
                            axc._poslayoutbox.constrain_height_min(
                                ax._poslayoutbox.height * heightC /
                                (height0*1.8))
                        elif height0 < heightC:
                            axc._poslayoutbox.constrain_height_min(
                                ax._poslayoutbox.height * heightC / height0)
                            ax._poslayoutbox.constrain_height_min(
                                ax._poslayoutbox.height * height0 /
                                (heightC*1.8))
                    # widths...
                    if not alignwidth and dcols0 == dcolsC:
                        ax._poslayoutbox.constrain_width(
                                axc._poslayoutbox.width * width0 / widthC)
                        alignwidth = True
                    elif in_same_row(rownum0min, rownum0max,
                            rownumCmin, rownumCmax):
                        if width0 > widthC:
                            ax._poslayoutbox.constrain_width_min(
                                    axc._poslayoutbox.width * width0 / widthC)
                            axc._poslayoutbox.constrain_width_min(
                                    ax._poslayoutbox.width * widthC /
                                    (width0*1.8))
                        elif width0 < widthC:
                            axc._poslayoutbox.constrain_width_min(
                                    ax._poslayoutbox.width * widthC / width0)
                            ax._poslayoutbox.constrain_width_min(
                                    axc._poslayoutbox.width * width0 /
                                    (widthC*1.8))

    fig._layoutbox.constrained_layout_called += 1
    fig._layoutbox.update_variables()
    # Now set the position of the axes...
    for ax in fig.axes:
        if ax._layoutbox is not None:
            newpos = ax._poslayoutbox.get_rect()
            _log.debug('newpos %r', newpos)
            # Now set the new position.
            # ax.set_position will zero out the layout for
            # this axis, allowing users to hard-code the position,
            # so this does the same w/o zeroing layout.
            ax._set_position(newpos, which='original')


def arange_subplotspecs(gs, hspace=0, wspace=0):
    """
    arange the subplotspec children of this gridspec, and then recursively
    do the same of any gridspec children of those gridspecs...
    """
    sschildren = []
    for child in gs.children:
        if child._is_subplotspec_layoutbox():
            for child2 in child.children:
                # check for gridspec children...
                name = (child2.name).split('.')[-1][:-3]
                if name == 'gridspec':
                    arange_subplotspecs(child2, hspace=hspace, wspace=wspace)
            sschildren += [child]
    # now arrange the subplots...
    for child0 in sschildren:
        ss0 = child0.artist
        nrows, ncols = ss0.get_gridspec().get_geometry()
        if ss0.num2 is None:
            ss0.num2 = ss0.num1
        rowNum0min, colNum0min = divmod(ss0.num1, ncols)
        rowNum0max, colNum0max = divmod(ss0.num2, ncols)
        sschildren = sschildren[1:]
        for childc in sschildren:
            ssc = childc.artist
            rowNumCmin, colNumCmin = divmod(ssc.num1, ncols)
            if ssc.num2 is None:
                ssc.num2 = ssc.num1
            rowNumCmax, colNumCmax = divmod(ssc.num2, ncols)
            # OK, this tells us the relative layout of ax
            # with axc
            thepad = wspace / ncols
            if colNum0max < colNumCmin:
                layoutbox.hstack([ss0._layoutbox, ssc._layoutbox],
                        padding=thepad)
            if colNumCmax < colNum0min:
                layoutbox.hstack([ssc._layoutbox, ss0._layoutbox],
                        padding=thepad)

            ####
            # vertical alignment
            thepad = hspace / nrows
            if rowNum0max < rowNumCmin:
                layoutbox.vstack([ss0._layoutbox,
                                 ssc._layoutbox],
                                 padding=thepad)
            if rowNumCmax < rowNum0min:
                layoutbox.vstack([ssc._layoutbox,
                                  ss0._layoutbox],
                                  padding=thepad)


def layoutcolorbarsingle(ax, cax, shrink, aspect, location, pad=0.05):
    """
    Do the layout for a colorbar, to not oeverly pollute colorbar.py

    `pad` is in fraction of the original axis size.
    """
    axlb = ax._layoutbox
    axpos = ax._poslayoutbox
    axsslb = ax.get_subplotspec()._layoutbox
    lb = layoutbox.LayoutBox(
            parent=axsslb,
            name=axsslb.name + '.cbar',
            artist=cax)

    if location in ('left', 'right'):
        lbpos = layoutbox.LayoutBox(
                parent=lb,
                name=lb.name + '.pos',
                tightwidth=False,
                pos=True,
                subplot=False,
                artist=cax)

        if location == 'right':
            # arrange to right of parent axis
            layoutbox.hstack([axlb, lb], padding=pad * axlb.width,
                             strength='strong')
        else:
            layoutbox.hstack([lb, axlb], padding=pad * axlb.width)
        # constrain the height and center...
        layoutbox.match_heights([axpos, lbpos], [1, shrink])
        layoutbox.align([axpos, lbpos], 'v_center')
        # set the width of the pos box
        lbpos.constrain_width(shrink * axpos.height * (1/aspect),
                              strength='strong')
    elif location in ('bottom', 'top'):
        lbpos = layoutbox.LayoutBox(
                parent=lb,
                name=lb.name + '.pos',
                tightheight=True,
                pos=True,
                subplot=False,
                artist=cax)

        if location == 'bottom':
            layoutbox.vstack([axlb, lb], padding=pad * axlb.height)
        else:
            layoutbox.vstack([lb, axlb], padding=pad * axlb.height)
        # constrain the height and center...
        layoutbox.match_widths([axpos, lbpos],
                               [1, shrink], strength='strong')
        layoutbox.align([axpos, lbpos], 'h_center')
        # set the height of the pos box
        lbpos.constrain_height(axpos.width * aspect * shrink,
                                strength='medium')

    return lb, lbpos


def layoutcolorbargridspec(parents, cax, shrink, aspect, location, pad=0.05):
    """
    Do the layout for a colorbar, to not oeverly pollute colorbar.py

    `pad` is in fraction of the original axis size.
    """

    gs = parents[0].get_subplotspec().get_gridspec()
    # parent layout box....
    gslb = gs._layoutbox

    lb = layoutbox.LayoutBox(parent=gslb.parent,
                             name=gslb.parent.name + '.cbar',
                             artist=cax)
    if location in ('left', 'right'):
        lbpos = layoutbox.LayoutBox(
                parent=lb,
                name=lb.name + '.pos',
                tightwidth=False,
                pos=True,
                subplot=False,
                artist=cax)

        if location == 'right':
            # arrange to right of the gridpec sibbling
            layoutbox.hstack([gslb, lb], padding=pad * gslb.width,
                             strength='strong')
        else:
            layoutbox.hstack([lb, gslb], padding=pad * gslb.width)
        # constrain the height and center...
        # This isn't quite right.  We'd like the colorbar
        # pos to line up w/ the axes poss, not the size of the
        # gs.
        maxrow = -100000
        minrow = 1000000
        maxax = None
        minax = None

        for ax in parents:
            subspec = ax.get_subplotspec()
            nrows, ncols = subspec.get_gridspec().get_geometry()
            for num in [subspec.num1, subspec.num2]:
                rownum1, colnum1 = divmod(subspec.num1, ncols)
                if rownum1 > maxrow:
                    maxrow = rownum1
                    maxax = ax
                if rownum1 < minrow:
                    minrow = rownum1
                    minax = ax
        # invert the order so these are bottom to top:
        maxposlb = minax._poslayoutbox
        minposlb = maxax._poslayoutbox
        # now we want the height of the colorbar pos to be
        # set by the top and bottom of these poss
        # bottom              top
        #     b             t
        # h = (top-bottom)*shrink
        # b = bottom + (top-bottom - h) / 2.
        lbpos.constrain_height(
                (maxposlb.top - minposlb.bottom) *
                shrink, strength='strong')
        lbpos.constrain_bottom(
                (maxposlb.top - minposlb.bottom) *
                (1 - shrink)/2 + minposlb.bottom,
                strength='strong')

        # set the width of the pos box
        lbpos.constrain_width(lbpos.height * (shrink / aspect),
                              strength='strong')
    elif location in ('bottom', 'top'):
        lbpos = layoutbox.LayoutBox(
                parent=lb,
                name=lb.name + '.pos',
                tightheight=True,
                pos=True,
                subplot=False,
                artist=cax)

        if location == 'bottom':
            layoutbox.vstack([gslb, lb], padding=pad * gslb.width)
        else:
            layoutbox.vstack([lb, gslb], padding=pad * gslb.width)

        maxcol = -100000
        mincol = 1000000
        maxax = None
        minax = None

        for ax in parents:
            subspec = ax.get_subplotspec()
            nrows, ncols = subspec.get_gridspec().get_geometry()
            for num in [subspec.num1, subspec.num2]:
                rownum1, colnum1 = divmod(subspec.num1, ncols)
                if colnum1 > maxcol:
                    maxcol = colnum1
                    maxax = ax
                if rownum1 < mincol:
                    mincol = colnum1
                    minax = ax
        maxposlb = maxax._poslayoutbox
        minposlb = minax._poslayoutbox
        lbpos.constrain_width((maxposlb.right - minposlb.left) *
                              shrink)
        lbpos.constrain_left(
                (maxposlb.right - minposlb.left) *
                (1-shrink)/2 + minposlb.left)
        # set the height of the pos box
        lbpos.constrain_height(lbpos.width * shrink * aspect,
                               strength='medium')

    return lb, lbpos
