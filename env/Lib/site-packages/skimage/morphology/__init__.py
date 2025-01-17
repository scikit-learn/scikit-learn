"""Morphological algorithms, e.g., closing, opening, skeletonization."""

from .binary import binary_closing, binary_dilation, binary_erosion, binary_opening
from .gray import black_tophat, closing, dilation, erosion, opening, white_tophat
from .isotropic import (
    isotropic_erosion,
    isotropic_dilation,
    isotropic_opening,
    isotropic_closing,
)
from .footprints import (
    ball,
    cube,
    diamond,
    disk,
    ellipse,
    footprint_from_sequence,
    footprint_rectangle,
    mirror_footprint,
    octagon,
    octahedron,
    pad_footprint,
    rectangle,
    square,
    star,
)
from ..measure._label import label
from ._skeletonize import medial_axis, skeletonize, thin
from .convex_hull import convex_hull_image, convex_hull_object
from .grayreconstruct import reconstruction
from .misc import remove_small_holes, remove_small_objects, remove_objects_by_distance
from .extrema import h_maxima, h_minima, local_minima, local_maxima
from ._flood_fill import flood, flood_fill
from .max_tree import (
    area_opening,
    area_closing,
    diameter_closing,
    diameter_opening,
    max_tree,
    max_tree_local_maxima,
)

__all__ = [
    'area_closing',
    'area_opening',
    'ball',
    'binary_closing',
    'binary_dilation',
    'binary_erosion',
    'binary_opening',
    'black_tophat',
    'closing',
    'convex_hull_image',
    'convex_hull_object',
    'diameter_closing',
    'diameter_opening',
    'diamond',
    'dilation',
    'disk',
    'ellipse',
    'erosion',
    'flood',
    'flood_fill',
    'footprint_from_sequence',
    'footprint_rectangle',
    'h_maxima',
    'h_minima',
    'isotropic_closing',
    'isotropic_dilation',
    'isotropic_erosion',
    'isotropic_opening',
    'label',
    'local_maxima',
    'local_minima',
    'max_tree',
    'max_tree_local_maxima',
    'medial_axis',
    'mirror_footprint',
    'octagon',
    'octahedron',
    'opening',
    'pad_footprint',
    'reconstruction',
    'remove_small_holes',
    'remove_small_objects',
    'remove_objects_by_distance',
    'skeletonize',
    'star',
    'thin',
    'white_tophat',
]
