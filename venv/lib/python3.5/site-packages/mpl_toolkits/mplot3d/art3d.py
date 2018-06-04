# art3d.py, original mplot3d version by John Porter
# Parts rewritten by Reinier Heeres <reinier@heeres.eu>
# Minor additions by Ben Axelrod <baxelrod@coroware.com>

'''
Module containing 3D artist code and functions to convert 2D
artists into 3D versions which can be added to an Axes3D.
'''
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip

import math

import numpy as np

from matplotlib import (
    artist, cbook, colors as mcolors, lines, text as mtext, path as mpath)
from matplotlib.cbook import _backports
from matplotlib.collections import (
    Collection, LineCollection, PolyCollection, PatchCollection,
    PathCollection)
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from . import proj3d


def norm_angle(a):
    """Return angle between -180 and +180"""
    a = (a + 360) % 360
    if a > 180:
        a = a - 360
    return a


def norm_text_angle(a):
    """Return angle between -90 and +90"""
    a = (a + 180) % 180
    if a > 90:
        a = a - 180
    return a


def get_dir_vector(zdir):
    if zdir == 'x':
        return np.array((1, 0, 0))
    elif zdir == 'y':
        return np.array((0, 1, 0))
    elif zdir == 'z':
        return np.array((0, 0, 1))
    elif zdir is None:
        return np.array((0, 0, 0))
    elif cbook.iterable(zdir) and len(zdir) == 3:
        return zdir
    else:
        raise ValueError("'x', 'y', 'z', None or vector of length 3 expected")


class Text3D(mtext.Text):
    '''
    Text object with 3D position and (in the future) direction.
    '''

    def __init__(self, x=0, y=0, z=0, text='', zdir='z', **kwargs):
        '''
        *x*, *y*, *z*  Position of text
        *text*         Text string to display
        *zdir*         Direction of text

        Keyword arguments are passed onto :func:`~matplotlib.text.Text`.
        '''
        mtext.Text.__init__(self, x, y, text, **kwargs)
        self.set_3d_properties(z, zdir)

    def set_3d_properties(self, z=0, zdir='z'):
        x, y = self.get_position()
        self._position3d = np.array((x, y, z))
        self._dir_vec = get_dir_vector(zdir)
        self.stale = True

    def draw(self, renderer):
        proj = proj3d.proj_trans_points(
            [self._position3d, self._position3d + self._dir_vec], renderer.M)
        dx = proj[0][1] - proj[0][0]
        dy = proj[1][1] - proj[1][0]
        if dx==0. and dy==0.:
            # atan2 raises ValueError: math domain error on 0,0
            angle = 0.
        else:
            angle = math.degrees(math.atan2(dy, dx))
        self.set_position((proj[0][0], proj[1][0]))
        self.set_rotation(norm_text_angle(angle))
        mtext.Text.draw(self, renderer)
        self.stale = False


def text_2d_to_3d(obj, z=0, zdir='z'):
    """Convert a Text to a Text3D object."""
    obj.__class__ = Text3D
    obj.set_3d_properties(z, zdir)


class Line3D(lines.Line2D):
    '''
    3D line object.
    '''

    def __init__(self, xs, ys, zs, *args, **kwargs):
        '''
        Keyword arguments are passed onto :func:`~matplotlib.lines.Line2D`.
        '''
        lines.Line2D.__init__(self, [], [], *args, **kwargs)
        self._verts3d = xs, ys, zs

    def set_3d_properties(self, zs=0, zdir='z'):
        xs = self.get_xdata()
        ys = self.get_ydata()

        try:
            # If *zs* is a list or array, then this will fail and
            # just proceed to juggle_axes().
            zs = float(zs)
            zs = [zs for x in xs]
        except TypeError:
            pass
        self._verts3d = juggle_axes(xs, ys, zs, zdir)
        self.stale = True

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_data(xs, ys)
        lines.Line2D.draw(self, renderer)
        self.stale = False


def line_2d_to_3d(line, zs=0, zdir='z'):
    '''
    Convert a 2D line to 3D.
    '''
    line.__class__ = Line3D
    line.set_3d_properties(zs, zdir)


def path_to_3d_segment(path, zs=0, zdir='z'):
    '''Convert a path to a 3D segment.'''

    zs = _backports.broadcast_to(zs, len(path))
    pathsegs = path.iter_segments(simplify=False, curves=False)
    seg = [(x, y, z) for (((x, y), code), z) in zip(pathsegs, zs)]
    seg3d = [juggle_axes(x, y, z, zdir) for (x, y, z) in seg]
    return seg3d


def paths_to_3d_segments(paths, zs=0, zdir='z'):
    '''
    Convert paths from a collection object to 3D segments.
    '''

    zs = _backports.broadcast_to(zs, len(paths))
    segs = [path_to_3d_segment(path, pathz, zdir)
            for path, pathz in zip(paths, zs)]
    return segs


def path_to_3d_segment_with_codes(path, zs=0, zdir='z'):
    '''Convert a path to a 3D segment with path codes.'''

    zs = _backports.broadcast_to(zs, len(path))
    seg = []
    codes = []
    pathsegs = path.iter_segments(simplify=False, curves=False)
    for (((x, y), code), z) in zip(pathsegs, zs):
        seg.append((x, y, z))
        codes.append(code)
    seg3d = [juggle_axes(x, y, z, zdir) for (x, y, z) in seg]
    return seg3d, codes


def paths_to_3d_segments_with_codes(paths, zs=0, zdir='z'):
    '''
    Convert paths from a collection object to 3D segments with path codes.
    '''

    zs = _backports.broadcast_to(zs, len(paths))
    segments = []
    codes_list = []
    for path, pathz in zip(paths, zs):
        segs, codes = path_to_3d_segment_with_codes(path, pathz, zdir)
        segments.append(segs)
        codes_list.append(codes)
    return segments, codes_list


class Line3DCollection(LineCollection):
    '''
    A collection of 3D lines.
    '''

    def __init__(self, segments, *args, **kwargs):
        '''
        Keyword arguments are passed onto :func:`~matplotlib.collections.LineCollection`.
        '''
        LineCollection.__init__(self, segments, *args, **kwargs)

    def set_sort_zpos(self, val):
        '''Set the position to use for z-sorting.'''
        self._sort_zpos = val
        self.stale = True

    def set_segments(self, segments):
        '''
        Set 3D segments
        '''
        self._segments3d = np.asanyarray(segments)
        LineCollection.set_segments(self, [])

    def do_3d_projection(self, renderer):
        '''
        Project the points according to renderer matrix.
        '''
        xyslist = [
            proj3d.proj_trans_points(points, renderer.M) for points in
            self._segments3d]
        segments_2d = [np.column_stack([xs, ys]) for xs, ys, zs in xyslist]
        LineCollection.set_segments(self, segments_2d)

        # FIXME
        minz = 1e9
        for xs, ys, zs in xyslist:
            minz = min(minz, min(zs))
        return minz

    def draw(self, renderer, project=False):
        if project:
            self.do_3d_projection(renderer)
        LineCollection.draw(self, renderer)


def line_collection_2d_to_3d(col, zs=0, zdir='z'):
    """Convert a LineCollection to a Line3DCollection object."""
    segments3d = paths_to_3d_segments(col.get_paths(), zs, zdir)
    col.__class__ = Line3DCollection
    col.set_segments(segments3d)


class Patch3D(Patch):
    '''
    3D patch object.
    '''

    def __init__(self, *args, **kwargs):
        zs = kwargs.pop('zs', [])
        zdir = kwargs.pop('zdir', 'z')
        Patch.__init__(self, *args, **kwargs)
        self.set_3d_properties(zs, zdir)

    def set_3d_properties(self, verts, zs=0, zdir='z'):
        zs = _backports.broadcast_to(zs, len(verts))
        self._segment3d = [juggle_axes(x, y, z, zdir)
                           for ((x, y), z) in zip(verts, zs)]
        self._facecolor3d = Patch.get_facecolor(self)

    def get_path(self):
        return self._path2d

    def get_facecolor(self):
        return self._facecolor2d

    def do_3d_projection(self, renderer):
        s = self._segment3d
        xs, ys, zs = zip(*s)
        vxs, vys, vzs, vis = proj3d.proj_transform_clip(xs, ys, zs, renderer.M)
        self._path2d = mpath.Path(np.column_stack([vxs, vys]))
        # FIXME: coloring
        self._facecolor2d = self._facecolor3d
        return min(vzs)

    def draw(self, renderer):
        Patch.draw(self, renderer)


class PathPatch3D(Patch3D):
    '''
    3D PathPatch object.
    '''

    def __init__(self, path, **kwargs):
        zs = kwargs.pop('zs', [])
        zdir = kwargs.pop('zdir', 'z')
        Patch.__init__(self, **kwargs)
        self.set_3d_properties(path, zs, zdir)

    def set_3d_properties(self, path, zs=0, zdir='z'):
        Patch3D.set_3d_properties(self, path.vertices, zs=zs, zdir=zdir)
        self._code3d = path.codes

    def do_3d_projection(self, renderer):
        s = self._segment3d
        xs, ys, zs = zip(*s)
        vxs, vys, vzs, vis = proj3d.proj_transform_clip(xs, ys, zs, renderer.M)
        self._path2d = mpath.Path(np.column_stack([vxs, vys]), self._code3d)
        # FIXME: coloring
        self._facecolor2d = self._facecolor3d
        return min(vzs)


def get_patch_verts(patch):
    """Return a list of vertices for the path of a patch."""
    trans = patch.get_patch_transform()
    path =  patch.get_path()
    polygons = path.to_polygons(trans)
    if len(polygons):
        return polygons[0]
    else:
        return []


def patch_2d_to_3d(patch, z=0, zdir='z'):
    """Convert a Patch to a Patch3D object."""
    verts = get_patch_verts(patch)
    patch.__class__ = Patch3D
    patch.set_3d_properties(verts, z, zdir)


def pathpatch_2d_to_3d(pathpatch, z=0, zdir='z'):
    """Convert a PathPatch to a PathPatch3D object."""
    path = pathpatch.get_path()
    trans = pathpatch.get_patch_transform()

    mpath = trans.transform_path(path)
    pathpatch.__class__ = PathPatch3D
    pathpatch.set_3d_properties(mpath, z, zdir)


class Patch3DCollection(PatchCollection):
    '''
    A collection of 3D patches.
    '''

    def __init__(self, *args, **kwargs):
        """
        Create a collection of flat 3D patches with its normal vector
        pointed in *zdir* direction, and located at *zs* on the *zdir*
        axis. 'zs' can be a scalar or an array-like of the same length as
        the number of patches in the collection.

        Constructor arguments are the same as for
        :class:`~matplotlib.collections.PatchCollection`. In addition,
        keywords *zs=0* and *zdir='z'* are available.

        Also, the keyword argument "depthshade" is available to
        indicate whether or not to shade the patches in order to
        give the appearance of depth (default is *True*).
        This is typically desired in scatter plots.
        """
        zs = kwargs.pop('zs', 0)
        zdir = kwargs.pop('zdir', 'z')
        self._depthshade = kwargs.pop('depthshade', True)
        PatchCollection.__init__(self, *args, **kwargs)
        self.set_3d_properties(zs, zdir)

    def set_sort_zpos(self, val):
        '''Set the position to use for z-sorting.'''
        self._sort_zpos = val
        self.stale = True

    def set_3d_properties(self, zs, zdir):
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        offsets = self.get_offsets()
        if len(offsets) > 0:
            xs, ys = zip(*offsets)
        else:
            xs = []
            ys = []
        self._offsets3d = juggle_axes(xs, ys, np.atleast_1d(zs), zdir)
        self._facecolor3d = self.get_facecolor()
        self._edgecolor3d = self.get_edgecolor()
        self.stale = True

    def do_3d_projection(self, renderer):
        xs, ys, zs = self._offsets3d
        vxs, vys, vzs, vis = proj3d.proj_transform_clip(xs, ys, zs, renderer.M)

        fcs = (zalpha(self._facecolor3d, vzs) if self._depthshade else
               self._facecolor3d)
        fcs = mcolors.to_rgba_array(fcs, self._alpha)
        self.set_facecolors(fcs)

        ecs = (zalpha(self._edgecolor3d, vzs) if self._depthshade else
               self._edgecolor3d)
        ecs = mcolors.to_rgba_array(ecs, self._alpha)
        self.set_edgecolors(ecs)
        PatchCollection.set_offsets(self, np.column_stack([vxs, vys]))

        if vzs.size > 0:
            return min(vzs)
        else:
            return np.nan


class Path3DCollection(PathCollection):
    '''
    A collection of 3D paths.
    '''

    def __init__(self, *args, **kwargs):
        """
        Create a collection of flat 3D paths with its normal vector
        pointed in *zdir* direction, and located at *zs* on the *zdir*
        axis. 'zs' can be a scalar or an array-like of the same length as
        the number of paths in the collection.

        Constructor arguments are the same as for
        :class:`~matplotlib.collections.PathCollection`. In addition,
        keywords *zs=0* and *zdir='z'* are available.

        Also, the keyword argument "depthshade" is available to
        indicate whether or not to shade the patches in order to
        give the appearance of depth (default is *True*).
        This is typically desired in scatter plots.
        """
        zs = kwargs.pop('zs', 0)
        zdir = kwargs.pop('zdir', 'z')
        self._depthshade = kwargs.pop('depthshade', True)
        PathCollection.__init__(self, *args, **kwargs)
        self.set_3d_properties(zs, zdir)

    def set_sort_zpos(self, val):
        '''Set the position to use for z-sorting.'''
        self._sort_zpos = val
        self.stale = True

    def set_3d_properties(self, zs, zdir):
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        offsets = self.get_offsets()
        if len(offsets) > 0:
            xs, ys = zip(*offsets)
        else:
            xs = []
            ys = []
        self._offsets3d = juggle_axes(xs, ys, np.atleast_1d(zs), zdir)
        self._facecolor3d = self.get_facecolor()
        self._edgecolor3d = self.get_edgecolor()
        self.stale = True

    def do_3d_projection(self, renderer):
        xs, ys, zs = self._offsets3d
        vxs, vys, vzs, vis = proj3d.proj_transform_clip(xs, ys, zs, renderer.M)

        fcs = (zalpha(self._facecolor3d, vzs) if self._depthshade else
               self._facecolor3d)
        fcs = mcolors.to_rgba_array(fcs, self._alpha)
        self.set_facecolors(fcs)

        ecs = (zalpha(self._edgecolor3d, vzs) if self._depthshade else
               self._edgecolor3d)
        ecs = mcolors.to_rgba_array(ecs, self._alpha)
        self.set_edgecolors(ecs)
        PathCollection.set_offsets(self, np.column_stack([vxs, vys]))

        if vzs.size > 0 :
            return min(vzs)
        else :
            return np.nan


def patch_collection_2d_to_3d(col, zs=0, zdir='z', depthshade=True):
    """
    Convert a :class:`~matplotlib.collections.PatchCollection` into a
    :class:`Patch3DCollection` object
    (or a :class:`~matplotlib.collections.PathCollection` into a
    :class:`Path3DCollection` object).

    Keywords:

    *za*            The location or locations to place the patches in the
                    collection along the *zdir* axis. Defaults to 0.

    *zdir*          The axis in which to place the patches. Default is "z".

    *depthshade*    Whether to shade the patches to give a sense of depth.
                    Defaults to *True*.

    """
    if isinstance(col, PathCollection):
        col.__class__ = Path3DCollection
    elif isinstance(col, PatchCollection):
        col.__class__ = Patch3DCollection
    col._depthshade = depthshade
    col.set_3d_properties(zs, zdir)


class Poly3DCollection(PolyCollection):
    '''
    A collection of 3D polygons.
    '''

    def __init__(self, verts, *args, **kwargs):
        '''
        Create a Poly3DCollection.

        *verts* should contain 3D coordinates.

        Keyword arguments:
        zsort, see set_zsort for options.

        Note that this class does a bit of magic with the _facecolors
        and _edgecolors properties.
        '''
        zsort = kwargs.pop('zsort', True)
        PolyCollection.__init__(self, verts, *args, **kwargs)
        self.set_zsort(zsort)
        self._codes3d = None

    _zsort_functions = {
        'average': np.average,
        'min': np.min,
        'max': np.max,
    }

    def set_zsort(self, zsort):
        '''
        Set z-sorting behaviour:
            boolean: if True use default 'average'
            string: 'average', 'min' or 'max'
        '''

        if zsort is True:
            zsort = 'average'

        if zsort is not False:
            if zsort in self._zsort_functions:
                zsortfunc = self._zsort_functions[zsort]
            else:
                return False
        else:
            zsortfunc = None

        self._zsort = zsort
        self._sort_zpos = None
        self._zsortfunc = zsortfunc
        self.stale = True

    def get_vector(self, segments3d):
        """Optimize points for projection"""
        si = 0
        ei = 0
        segis = []
        points = []
        for p in segments3d:
            points.extend(p)
            ei = si + len(p)
            segis.append((si, ei))
            si = ei

        if len(segments3d):
            xs, ys, zs = zip(*points)
        else :
            # We need this so that we can skip the bad unpacking from zip()
            xs, ys, zs = [], [], []

        ones = np.ones(len(xs))
        self._vec = np.array([xs, ys, zs, ones])
        self._segis = segis

    def set_verts(self, verts, closed=True):
        '''Set 3D vertices.'''
        self.get_vector(verts)
        # 2D verts will be updated at draw time
        PolyCollection.set_verts(self, [], False)
        self._closed = closed

    def set_verts_and_codes(self, verts, codes):
        '''Sets 3D vertices with path codes'''
        # set vertices with closed=False to prevent PolyCollection from
        # setting path codes
        self.set_verts(verts, closed=False)
        # and set our own codes instead.
        self._codes3d = codes

    def set_3d_properties(self):
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        self._sort_zpos = None
        self.set_zsort(True)
        self._facecolors3d = PolyCollection.get_facecolors(self)
        self._edgecolors3d = PolyCollection.get_edgecolors(self)
        self._alpha3d = PolyCollection.get_alpha(self)
        self.stale = True

    def set_sort_zpos(self,val):
        '''Set the position to use for z-sorting.'''
        self._sort_zpos = val
        self.stale = True

    def do_3d_projection(self, renderer):
        '''
        Perform the 3D projection for this object.
        '''
        # FIXME: This may no longer be needed?
        if self._A is not None:
            self.update_scalarmappable()
            self._facecolors3d = self._facecolors

        txs, tys, tzs = proj3d.proj_transform_vec(self._vec, renderer.M)
        xyzlist = [(txs[si:ei], tys[si:ei], tzs[si:ei])
                   for si, ei in self._segis]

        # This extra fuss is to re-order face / edge colors
        cface = self._facecolors3d
        cedge = self._edgecolors3d
        if len(cface) != len(xyzlist):
            cface = cface.repeat(len(xyzlist), axis=0)
        if len(cedge) != len(xyzlist):
            if len(cedge) == 0:
                cedge = cface
            else:
                cedge = cedge.repeat(len(xyzlist), axis=0)

        # if required sort by depth (furthest drawn first)
        if self._zsort:
            z_segments_2d = sorted(
                ((self._zsortfunc(zs), np.column_stack([xs, ys]), fc, ec, idx)
                 for idx, ((xs, ys, zs), fc, ec)
                 in enumerate(zip(xyzlist, cface, cedge))),
                key=lambda x: x[0], reverse=True)
        else:
            raise ValueError("whoops")

        segments_2d = [s for z, s, fc, ec, idx in z_segments_2d]
        if self._codes3d is not None:
            codes = [self._codes3d[idx] for z, s, fc, ec, idx in z_segments_2d]
            PolyCollection.set_verts_and_codes(self, segments_2d, codes)
        else:
            PolyCollection.set_verts(self, segments_2d, self._closed)

        self._facecolors2d = [fc for z, s, fc, ec, idx in z_segments_2d]
        if len(self._edgecolors3d) == len(cface):
            self._edgecolors2d = [ec for z, s, fc, ec, idx in z_segments_2d]
        else:
            self._edgecolors2d = self._edgecolors3d

        # Return zorder value
        if self._sort_zpos is not None:
            zvec = np.array([[0], [0], [self._sort_zpos], [1]])
            ztrans = proj3d.proj_transform_vec(zvec, renderer.M)
            return ztrans[2][0]
        elif tzs.size > 0 :
            # FIXME: Some results still don't look quite right.
            #        In particular, examine contourf3d_demo2.py
            #        with az = -54 and elev = -45.
            return np.min(tzs)
        else :
            return np.nan

    def set_facecolor(self, colors):
        PolyCollection.set_facecolor(self, colors)
        self._facecolors3d = PolyCollection.get_facecolor(self)
    set_facecolors = set_facecolor

    def set_edgecolor(self, colors):
        PolyCollection.set_edgecolor(self, colors)
        self._edgecolors3d = PolyCollection.get_edgecolor(self)
    set_edgecolors = set_edgecolor

    def set_alpha(self, alpha):
        """
        Set the alpha tranparencies of the collection.  *alpha* must be
        a float or *None*.

        ACCEPTS: float or None
        """
        if alpha is not None:
            try:
                float(alpha)
            except TypeError:
                raise TypeError('alpha must be a float or None')
        artist.Artist.set_alpha(self, alpha)
        try:
            self._facecolors = mcolors.to_rgba_array(
                self._facecolors3d, self._alpha)
        except (AttributeError, TypeError, IndexError):
            pass
        try:
            self._edgecolors = mcolors.to_rgba_array(
                    self._edgecolors3d, self._alpha)
        except (AttributeError, TypeError, IndexError):
            pass
        self.stale = True

    def get_facecolors(self):
        return self._facecolors2d
    get_facecolor = get_facecolors

    def get_edgecolors(self):
        return self._edgecolors2d
    get_edgecolor = get_edgecolors

    def draw(self, renderer):
        return Collection.draw(self, renderer)


def poly_collection_2d_to_3d(col, zs=0, zdir='z'):
    """Convert a PolyCollection to a Poly3DCollection object."""
    segments_3d, codes = paths_to_3d_segments_with_codes(col.get_paths(),
                                                         zs, zdir)
    col.__class__ = Poly3DCollection
    col.set_verts_and_codes(segments_3d, codes)
    col.set_3d_properties()


def juggle_axes(xs, ys, zs, zdir):
    """
    Reorder coordinates so that 2D xs, ys can be plotted in the plane
    orthogonal to zdir. zdir is normally x, y or z. However, if zdir
    starts with a '-' it is interpreted as a compensation for rotate_axes.
    """
    if zdir == 'x':
        return zs, xs, ys
    elif zdir == 'y':
        return xs, zs, ys
    elif zdir[0] == '-':
        return rotate_axes(xs, ys, zs, zdir)
    else:
        return xs, ys, zs


def rotate_axes(xs, ys, zs, zdir):
    """
    Reorder coordinates so that the axes are rotated with zdir along
    the original z axis. Prepending the axis with a '-' does the
    inverse transform, so zdir can be x, -x, y, -y, z or -z
    """
    if zdir == 'x':
        return ys, zs, xs
    elif zdir == '-x':
        return zs, xs, ys

    elif zdir == 'y':
        return zs, xs, ys
    elif zdir == '-y':
        return ys, zs, xs

    else:
        return xs, ys, zs


def get_colors(c, num):
    """Stretch the color argument to provide the required number num"""
    return _backports.broadcast_to(
        mcolors.to_rgba_array(c) if len(c) else [0, 0, 0, 0],
        (num, 4))


def zalpha(colors, zs):
    """Modify the alphas of the color list according to depth"""
    # FIXME: This only works well if the points for *zs* are well-spaced
    #        in all three dimensions. Otherwise, at certain orientations,
    #        the min and max zs are very close together.
    #        Should really normalize against the viewing depth.
    colors = get_colors(colors, len(zs))
    if len(zs):
        norm = Normalize(min(zs), max(zs))
        sats = 1 - norm(zs) * 0.7
        colors = [(c[0], c[1], c[2], c[3] * s) for c, s in zip(colors, sats)]
    return colors
