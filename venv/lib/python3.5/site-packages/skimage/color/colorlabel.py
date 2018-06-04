import itertools

import numpy as np

from .._shared.utils import warn
from ..util import img_as_float
from . import rgb_colors
from .colorconv import rgb2gray, gray2rgb

import six
from six.moves import zip


__all__ = ['color_dict', 'label2rgb', 'DEFAULT_COLORS']


DEFAULT_COLORS = ('red', 'blue', 'yellow', 'magenta', 'green',
                  'indigo', 'darkorange', 'cyan', 'pink', 'yellowgreen')


color_dict = dict((k, v) for k, v in six.iteritems(rgb_colors.__dict__)
                  if isinstance(v, tuple))


def _rgb_vector(color):
    """Return RGB color as (1, 3) array.

    This RGB array gets multiplied by masked regions of an RGB image, which are
    partially flattened by masking (i.e. dimensions 2D + RGB -> 1D + RGB).

    Parameters
    ----------
    color : str or array
        Color name in `color_dict` or RGB float values between [0, 1].
    """
    if isinstance(color, six.string_types):
        color = color_dict[color]
    # Slice to handle RGBA colors.
    return np.array(color[:3])


def _match_label_with_color(label, colors, bg_label, bg_color):
    """Return `unique_labels` and `color_cycle` for label array and color list.

    Colors are cycled for normal labels, but the background color should only
    be used for the background.
    """
    # Temporarily set background color; it will be removed later.
    if bg_color is None:
        bg_color = (0, 0, 0)
    bg_color = _rgb_vector([bg_color])

    # map labels to their ranks among all labels from small to large
    unique_labels, mapped_labels = np.unique(label, return_inverse=True)

    # get rank of bg_label
    bg_label_rank_list = mapped_labels[label.flat == bg_label]

    # The rank of each label is the index of the color it is matched to in
    # color cycle. bg_label should always be mapped to the first color, so
    # its rank must be 0. Other labels should be ranked from small to large
    # from 1.
    if len(bg_label_rank_list) > 0:
        bg_label_rank = bg_label_rank_list[0]
        mapped_labels[mapped_labels < bg_label_rank] += 1
        mapped_labels[label.flat == bg_label] = 0
    else:
        mapped_labels += 1

    # Modify labels and color cycle so background color is used only once.
    color_cycle = itertools.cycle(colors)
    color_cycle = itertools.chain(bg_color, color_cycle)

    return mapped_labels, color_cycle


def label2rgb(label, image=None, colors=None, alpha=0.3,
              bg_label=-1, bg_color=(0, 0, 0), image_alpha=1, kind='overlay'):
    """Return an RGB image where color-coded labels are painted over the image.

    Parameters
    ----------
    label : array, shape (M, N)
        Integer array of labels with the same shape as `image`.
    image : array, shape (M, N, 3), optional
        Image used as underlay for labels. If the input is an RGB image, it's
        converted to grayscale before coloring.
    colors : list, optional
        List of colors. If the number of labels exceeds the number of colors,
        then the colors are cycled.
    alpha : float [0, 1], optional
        Opacity of colorized labels. Ignored if image is `None`.
    bg_label : int, optional
        Label that's treated as the background.
    bg_color : str or array, optional
        Background color. Must be a name in `color_dict` or RGB float values
        between [0, 1].
    image_alpha : float [0, 1], optional
        Opacity of the image.
    kind : string, one of {'overlay', 'avg'}
        The kind of color image desired. 'overlay' cycles over defined colors
        and overlays the colored labels over the original image. 'avg' replaces
        each labeled segment with its average color, for a stained-class or
        pastel painting appearance.

    Returns
    -------
    result : array of float, shape (M, N, 3)
        The result of blending a cycling colormap (`colors`) for each distinct
        value in `label` with the image, at a certain alpha value.
    """
    if kind == 'overlay':
        return _label2rgb_overlay(label, image, colors, alpha, bg_label,
                                  bg_color, image_alpha)
    else:
        return _label2rgb_avg(label, image, bg_label, bg_color)


def _label2rgb_overlay(label, image=None, colors=None, alpha=0.3,
                       bg_label=-1, bg_color=None, image_alpha=1):
    """Return an RGB image where color-coded labels are painted over the image.

    Parameters
    ----------
    label : array, shape (M, N)
        Integer array of labels with the same shape as `image`.
    image : array, shape (M, N, 3), optional
        Image used as underlay for labels. If the input is an RGB image, it's
        converted to grayscale before coloring.
    colors : list, optional
        List of colors. If the number of labels exceeds the number of colors,
        then the colors are cycled.
    alpha : float [0, 1], optional
        Opacity of colorized labels. Ignored if image is `None`.
    bg_label : int, optional
        Label that's treated as the background.
    bg_color : str or array, optional
        Background color. Must be a name in `color_dict` or RGB float values
        between [0, 1].
    image_alpha : float [0, 1], optional
        Opacity of the image.

    Returns
    -------
    result : array of float, shape (M, N, 3)
        The result of blending a cycling colormap (`colors`) for each distinct
        value in `label` with the image, at a certain alpha value.
    """
    if colors is None:
        colors = DEFAULT_COLORS
    colors = [_rgb_vector(c) for c in colors]

    if image is None:
        image = np.zeros(label.shape + (3,), dtype=np.float64)
        # Opacity doesn't make sense if no image exists.
        alpha = 1
    else:
        if not image.shape[:2] == label.shape:
            raise ValueError("`image` and `label` must be the same shape")

        if image.min() < 0:
            warn("Negative intensities in `image` are not supported")

        image = img_as_float(rgb2gray(image))
        image = gray2rgb(image) * image_alpha + (1 - image_alpha)

    # Ensure that all labels are non-negative so we can index into
    # `label_to_color` correctly.
    offset = min(label.min(), bg_label)
    if offset != 0:
        label = label - offset  # Make sure you don't modify the input array.
        bg_label -= offset

    new_type = np.min_scalar_type(int(label.max()))
    if new_type == np.bool:
        new_type = np.uint8
    label = label.astype(new_type)

    mapped_labels_flat, color_cycle = _match_label_with_color(label, colors,
                                                              bg_label, bg_color)

    if len(mapped_labels_flat) == 0:
        return image

    dense_labels = range(max(mapped_labels_flat) + 1)

    label_to_color = np.array([c for i, c in zip(dense_labels, color_cycle)])

    mapped_labels = label
    mapped_labels.flat = mapped_labels_flat
    result = label_to_color[mapped_labels] * alpha + image * (1 - alpha)

    # Remove background label if its color was not specified.
    remove_background = 0 in mapped_labels_flat and bg_color is None
    if remove_background:
        result[label == bg_label] = image[label == bg_label]

    return result


def _label2rgb_avg(label_field, image, bg_label=0, bg_color=(0, 0, 0)):
    """Visualise each segment in `label_field` with its mean color in `image`.

    Parameters
    ----------
    label_field : array of int
        A segmentation of an image.
    image : array, shape ``label_field.shape + (3,)``
        A color image of the same spatial shape as `label_field`.
    bg_label : int, optional
        A value in `label_field` to be treated as background.
    bg_color : 3-tuple of int, optional
        The color for the background label

    Returns
    -------
    out : array, same shape and type as `image`
        The output visualization.
    """
    out = np.zeros_like(image)
    labels = np.unique(label_field)
    bg = (labels == bg_label)
    if bg.any():
        labels = labels[labels != bg_label]
        out[bg] = bg_color
    for label in labels:
        mask = (label_field == label).nonzero()
        color = image[mask].mean(axis=0)
        out[mask] = color
    return out
