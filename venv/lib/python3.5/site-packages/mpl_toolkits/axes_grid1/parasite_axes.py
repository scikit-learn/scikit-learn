from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib import (
    artist as martist, collections as mcoll, transforms as mtransforms,
    rcParams)
from matplotlib.axes import subplot_class_factory
from matplotlib.transforms import Bbox
from .mpl_axes import Axes

import numpy as np


class ParasiteAxesBase(object):

    def get_images_artists(self):
        artists = {a for a in self.get_children() if a.get_visible()}
        images = {a for a in self.images if a.get_visible()}

        return list(images), list(artists - images)

    def __init__(self, parent_axes, **kargs):

        self._parent_axes = parent_axes
        kargs.update(dict(frameon=False))
        self._get_base_axes_attr("__init__")(self, parent_axes.figure,
                                        parent_axes._position, **kargs)

    def cla(self):
        self._get_base_axes_attr("cla")(self)

        martist.setp(self.get_children(), visible=False)
        self._get_lines = self._parent_axes._get_lines

        # In mpl's Axes, zorders of x- and y-axis are originally set
        # within Axes.draw().
        if self._axisbelow:
            self.xaxis.set_zorder(0.5)
            self.yaxis.set_zorder(0.5)
        else:
            self.xaxis.set_zorder(2.5)
            self.yaxis.set_zorder(2.5)


_parasite_axes_classes = {}
def parasite_axes_class_factory(axes_class=None):
    if axes_class is None:
        axes_class = Axes

    new_class = _parasite_axes_classes.get(axes_class)
    if new_class is None:
        def _get_base_axes_attr(self, attrname):
            return getattr(axes_class, attrname)

        new_class = type(str("%sParasite" % (axes_class.__name__)),
                         (ParasiteAxesBase, axes_class),
                         {'_get_base_axes_attr': _get_base_axes_attr})
        _parasite_axes_classes[axes_class] = new_class

    return new_class

ParasiteAxes = parasite_axes_class_factory()

# #class ParasiteAxes(ParasiteAxesBase, Axes):

#     @classmethod
#     def _get_base_axes_attr(cls, attrname):
#         return getattr(Axes, attrname)



class ParasiteAxesAuxTransBase(object):
    def __init__(self, parent_axes, aux_transform, viewlim_mode=None,
                 **kwargs):

        self.transAux = aux_transform
        self.set_viewlim_mode(viewlim_mode)

        self._parasite_axes_class.__init__(self, parent_axes, **kwargs)

    def _set_lim_and_transforms(self):

        self.transAxes = self._parent_axes.transAxes

        self.transData = \
            self.transAux + \
            self._parent_axes.transData

        self._xaxis_transform = mtransforms.blended_transform_factory(
                self.transData, self.transAxes)
        self._yaxis_transform = mtransforms.blended_transform_factory(
                self.transAxes, self.transData)

    def set_viewlim_mode(self, mode):
        if mode not in [None, "equal", "transform"]:
            raise ValueError("Unknown mode : %s" % (mode,))
        else:
            self._viewlim_mode = mode

    def get_viewlim_mode(self):
        return self._viewlim_mode


    def update_viewlim(self):
        viewlim = self._parent_axes.viewLim.frozen()
        mode = self.get_viewlim_mode()
        if mode is None:
            pass
        elif mode == "equal":
            self.axes.viewLim.set(viewlim)
        elif mode == "transform":
            self.axes.viewLim.set(viewlim.transformed(self.transAux.inverted()))
        else:
            raise ValueError("Unknown mode : %s" % (self._viewlim_mode,))


    def _pcolor(self, method_name, *XYC, **kwargs):
        if len(XYC) == 1:
            C = XYC[0]
            ny, nx = C.shape

            gx = np.arange(-0.5, nx, 1.)
            gy = np.arange(-0.5, ny, 1.)

            X, Y = np.meshgrid(gx, gy)
        else:
            X, Y, C = XYC

        pcolor_routine = self._get_base_axes_attr(method_name)

        if "transform" in kwargs:
            mesh = pcolor_routine(self, X, Y, C, **kwargs)
        else:
            orig_shape = X.shape
            xy = np.vstack([X.flat, Y.flat])
            xyt=xy.transpose()
            wxy = self.transAux.transform(xyt)
            gx, gy = wxy[:,0].reshape(orig_shape), wxy[:,1].reshape(orig_shape)
            mesh = pcolor_routine(self, gx, gy, C, **kwargs)
            mesh.set_transform(self._parent_axes.transData)

        return mesh

    def pcolormesh(self, *XYC, **kwargs):
        return self._pcolor("pcolormesh", *XYC, **kwargs)

    def pcolor(self, *XYC, **kwargs):
        return self._pcolor("pcolor", *XYC, **kwargs)


    def _contour(self, method_name, *XYCL, **kwargs):

        if len(XYCL) <= 2:
            C = XYCL[0]
            ny, nx = C.shape

            gx = np.arange(0., nx, 1.)
            gy = np.arange(0., ny, 1.)

            X,Y = np.meshgrid(gx, gy)
            CL = XYCL
        else:
            X, Y = XYCL[:2]
            CL = XYCL[2:]

        contour_routine = self._get_base_axes_attr(method_name)

        if "transform" in kwargs:
            cont = contour_routine(self, X, Y, *CL, **kwargs)
        else:
            orig_shape = X.shape
            xy = np.vstack([X.flat, Y.flat])
            xyt=xy.transpose()
            wxy = self.transAux.transform(xyt)
            gx, gy = wxy[:,0].reshape(orig_shape), wxy[:,1].reshape(orig_shape)
            cont = contour_routine(self, gx, gy, *CL, **kwargs)
            for c in cont.collections:
                c.set_transform(self._parent_axes.transData)

        return cont

    def contour(self, *XYCL, **kwargs):
        return self._contour("contour", *XYCL, **kwargs)

    def contourf(self, *XYCL, **kwargs):
        return self._contour("contourf", *XYCL, **kwargs)

    def apply_aspect(self, position=None):
        self.update_viewlim()
        self._get_base_axes_attr("apply_aspect")(self)
        #ParasiteAxes.apply_aspect()



_parasite_axes_auxtrans_classes = {}
def parasite_axes_auxtrans_class_factory(axes_class=None):
    if axes_class is None:
        parasite_axes_class = ParasiteAxes
    elif not issubclass(axes_class, ParasiteAxesBase):
        parasite_axes_class = parasite_axes_class_factory(axes_class)
    else:
        parasite_axes_class = axes_class

    new_class = _parasite_axes_auxtrans_classes.get(parasite_axes_class)
    if new_class is None:
        new_class = type(str("%sParasiteAuxTrans" % (parasite_axes_class.__name__)),
                         (ParasiteAxesAuxTransBase, parasite_axes_class),
                         {'_parasite_axes_class': parasite_axes_class,
                         'name': 'parasite_axes'})
        _parasite_axes_auxtrans_classes[parasite_axes_class] = new_class

    return new_class


ParasiteAxesAuxTrans = parasite_axes_auxtrans_class_factory(axes_class=ParasiteAxes)




def _get_handles(ax):
    handles = ax.lines[:]
    handles.extend(ax.patches)
    handles.extend([c for c in ax.collections
                    if isinstance(c, mcoll.LineCollection)])
    handles.extend([c for c in ax.collections
                    if isinstance(c, mcoll.RegularPolyCollection)])
    handles.extend([c for c in ax.collections
                    if isinstance(c, mcoll.CircleCollection)])

    return handles


class HostAxesBase(object):
    def __init__(self, *args, **kwargs):

        self.parasites = []
        self._get_base_axes_attr("__init__")(self, *args, **kwargs)


    def get_aux_axes(self, tr, viewlim_mode="equal", axes_class=None):
        parasite_axes_class = parasite_axes_auxtrans_class_factory(axes_class)
        ax2 = parasite_axes_class(self, tr, viewlim_mode)
        # note that ax2.transData == tr + ax1.transData
        # Anthing you draw in ax2 will match the ticks and grids of ax1.
        self.parasites.append(ax2)
        ax2._remove_method = lambda h: self.parasites.remove(h)
        return ax2

    def _get_legend_handles(self, legend_handler_map=None):
        # don't use this!
        Axes_get_legend_handles = self._get_base_axes_attr("_get_legend_handles")
        all_handles = list(Axes_get_legend_handles(self, legend_handler_map))

        for ax in self.parasites:
            all_handles.extend(ax._get_legend_handles(legend_handler_map))

        return all_handles


    def draw(self, renderer):

        orig_artists = list(self.artists)
        orig_images = list(self.images)

        if hasattr(self, "get_axes_locator"):
            locator = self.get_axes_locator()
            if locator:
                pos = locator(self, renderer)
                self.set_position(pos, which="active")
                self.apply_aspect(pos)
            else:
                self.apply_aspect()
        else:
            self.apply_aspect()

        rect = self.get_position()

        for ax in self.parasites:
            ax.apply_aspect(rect)
            images, artists = ax.get_images_artists()
            self.images.extend(images)
            self.artists.extend(artists)

        self._get_base_axes_attr("draw")(self, renderer)
        self.artists = orig_artists
        self.images = orig_images


    def cla(self):

        for ax in self.parasites:
            ax.cla()

        self._get_base_axes_attr("cla")(self)
        #super(HostAxes, self).cla()


    def twinx(self, axes_class=None):
        """
        create a twin of Axes for generating a plot with a sharex
        x-axis but independent y axis.  The y-axis of self will have
        ticks on left and the returned axes will have ticks on the
        right
        """

        if axes_class is None:
            axes_class = self._get_base_axes()

        parasite_axes_class = parasite_axes_class_factory(axes_class)

        ax2 = parasite_axes_class(self, sharex=self, frameon=False)
        self.parasites.append(ax2)

        self.axis["right"].set_visible(False)

        ax2.axis["right"].set_visible(True)
        ax2.axis["left", "top", "bottom"].set_visible(False)

        def _remove_method(h):
            self.parasites.remove(h)
            self.axis["right"].set_visible(True)
            self.axis["right"].toggle(ticklabels=False, label=False)
        ax2._remove_method = _remove_method

        return ax2

    def twiny(self, axes_class=None):
        """
        create a twin of Axes for generating a plot with a shared
        y-axis but independent x axis.  The x-axis of self will have
        ticks on bottom and the returned axes will have ticks on the
        top
        """

        if axes_class is None:
            axes_class = self._get_base_axes()

        parasite_axes_class = parasite_axes_class_factory(axes_class)

        ax2 = parasite_axes_class(self, sharey=self, frameon=False)
        self.parasites.append(ax2)

        self.axis["top"].set_visible(False)

        ax2.axis["top"].set_visible(True)
        ax2.axis["left", "right", "bottom"].set_visible(False)

        def _remove_method(h):
            self.parasites.remove(h)
            self.axis["top"].set_visible(True)
            self.axis["top"].toggle(ticklabels=False, label=False)
        ax2._remove_method = _remove_method

        return ax2


    def twin(self, aux_trans=None, axes_class=None):
        """
        create a twin of Axes for generating a plot with a sharex
        x-axis but independent y axis.  The y-axis of self will have
        ticks on left and the returned axes will have ticks on the
        right
        """

        if axes_class is None:
            axes_class = self._get_base_axes()

        parasite_axes_auxtrans_class = parasite_axes_auxtrans_class_factory(axes_class)

        if aux_trans is None:
            ax2 = parasite_axes_auxtrans_class(self, mtransforms.IdentityTransform(),
                                               viewlim_mode="equal",
                                               )
        else:
            ax2 = parasite_axes_auxtrans_class(self, aux_trans,
                                               viewlim_mode="transform",
                                               )
        self.parasites.append(ax2)
        ax2._remove_method = lambda h: self.parasites.remove(h)

        self.axis["top", "right"].set_visible(False)

        ax2.axis["top", "right"].set_visible(True)
        ax2.axis["left", "bottom"].set_visible(False)

        def _remove_method(h):
            self.parasites.remove(h)
            self.axis["top", "right"].set_visible(True)
            self.axis["top", "right"].toggle(ticklabels=False, label=False)
        ax2._remove_method = _remove_method

        return ax2

    def get_tightbbox(self, renderer, call_axes_locator=True):

        bbs = [ax.get_tightbbox(renderer, call_axes_locator)
               for ax in self.parasites]
        get_tightbbox = self._get_base_axes_attr("get_tightbbox")
        bbs.append(get_tightbbox(self, renderer, call_axes_locator))

        _bbox = Bbox.union([b for b in bbs if b.width!=0 or b.height!=0])

        return _bbox



_host_axes_classes = {}
def host_axes_class_factory(axes_class=None):
    if axes_class is None:
        axes_class = Axes

    new_class = _host_axes_classes.get(axes_class)
    if new_class is None:
        def _get_base_axes(self):
            return axes_class

        def _get_base_axes_attr(self, attrname):
            return getattr(axes_class, attrname)

        new_class = type(str("%sHostAxes" % (axes_class.__name__)),
                         (HostAxesBase, axes_class),
                         {'_get_base_axes_attr': _get_base_axes_attr,
                          '_get_base_axes': _get_base_axes})

        _host_axes_classes[axes_class] = new_class

    return new_class

def host_subplot_class_factory(axes_class):
    host_axes_class = host_axes_class_factory(axes_class=axes_class)
    subplot_host_class = subplot_class_factory(host_axes_class)
    return subplot_host_class

HostAxes = host_axes_class_factory(axes_class=Axes)
SubplotHost = subplot_class_factory(HostAxes)


def host_axes(*args, **kwargs):
    """
    Create axes that can act as a hosts to parasitic axes.

    Parameters
    ----------
    figure : `matplotlib.figure.Figure`
        Figure to which the axes will be added. Defaults to the current figure
        `pyplot.gcf()`.

    *args, **kwargs :
        Will be passed on to the underlying ``Axes`` object creation.
    """
    import matplotlib.pyplot as plt
    axes_class = kwargs.pop("axes_class", None)
    host_axes_class = host_axes_class_factory(axes_class)
    fig = kwargs.get("figure", None)
    if fig is None:
        fig = plt.gcf()
    ax = host_axes_class(fig, *args, **kwargs)
    fig.add_axes(ax)
    plt.draw_if_interactive()
    return ax

def host_subplot(*args, **kwargs):
    """
    Create a subplot that can act as a host to parasitic axes.

    Parameters
    ----------
    figure : `matplotlib.figure.Figure`
        Figure to which the subplot will be added. Defaults to the current
        figure `pyplot.gcf()`.

    *args, **kwargs :
        Will be passed on to the underlying ``Axes`` object creation.
    """
    import matplotlib.pyplot as plt
    axes_class = kwargs.pop("axes_class", None)
    host_subplot_class = host_subplot_class_factory(axes_class)
    fig = kwargs.get("figure", None)
    if fig is None:
        fig = plt.gcf()
    ax = host_subplot_class(fig, *args, **kwargs)
    fig.add_subplot(ax)
    plt.draw_if_interactive()
    return ax
