import numpy as np
from .._shared.utils import deprecated


def join_segmentations(s1, s2):
    """Return the join of the two input segmentations.

    The join J of S1 and S2 is defined as the segmentation in which two
    voxels are in the same segment if and only if they are in the same
    segment in *both* S1 and S2.

    Parameters
    ----------
    s1, s2 : numpy arrays
        s1 and s2 are label fields of the same shape.

    Returns
    -------
    j : numpy array
        The join segmentation of s1 and s2.

    Examples
    --------
    >>> from skimage.segmentation import join_segmentations
    >>> s1 = np.array([[0, 0, 1, 1],
    ...                [0, 2, 1, 1],
    ...                [2, 2, 2, 1]])
    >>> s2 = np.array([[0, 1, 1, 0],
    ...                [0, 1, 1, 0],
    ...                [0, 1, 1, 1]])
    >>> join_segmentations(s1, s2)
    array([[0, 1, 3, 2],
           [0, 5, 3, 2],
           [4, 5, 5, 3]])
    """
    if s1.shape != s2.shape:
        raise ValueError("Cannot join segmentations of different shape. " +
                         "s1.shape: %s, s2.shape: %s" % (s1.shape, s2.shape))
    s1 = relabel_sequential(s1)[0]
    s2 = relabel_sequential(s2)[0]
    j = (s2.max() + 1) * s1 + s2
    j = relabel_sequential(j)[0]
    return j


@deprecated('relabel_sequential')
def relabel_from_one(label_field):
    """Convert labels in an arbitrary label field to {1, ... number_of_labels}.

    This function is deprecated, see ``relabel_sequential`` for more.
    """
    return relabel_sequential(label_field, offset=1)


def relabel_sequential(label_field, offset=1):
    """Relabel arbitrary labels to {`offset`, ... `offset` + number_of_labels}.

    This function also returns the forward map (mapping the original labels to
    the reduced labels) and the inverse map (mapping the reduced labels back
    to the original ones).

    Parameters
    ----------
    label_field : numpy array of int, arbitrary shape
        An array of labels.
    offset : int, optional
        The return labels will start at `offset`, which should be
        strictly positive.

    Returns
    -------
    relabeled : numpy array of int, same shape as `label_field`
        The input label field with labels mapped to
        {offset, ..., number_of_labels + offset - 1}.
    forward_map : numpy array of int, shape ``(label_field.max() + 1,)``
        The map from the original label space to the returned label
        space. Can be used to re-apply the same mapping. See examples
        for usage.
    inverse_map : 1D numpy array of int, of length offset + number of labels
        The map from the new label space to the original space. This
        can be used to reconstruct the original label field from the
        relabeled one.

    Notes
    -----
    The label 0 is assumed to denote the background and is never remapped.

    The forward map can be extremely big for some inputs, since its
    length is given by the maximum of the label field. However, in most
    situations, ``label_field.max()`` is much smaller than
    ``label_field.size``, and in these cases the forward map is
    guaranteed to be smaller than either the input or output images.

    Examples
    --------
    >>> from skimage.segmentation import relabel_sequential
    >>> label_field = np.array([1, 1, 5, 5, 8, 99, 42])
    >>> relab, fw, inv = relabel_sequential(label_field)
    >>> relab
    array([1, 1, 2, 2, 3, 5, 4])
    >>> fw
    array([0, 1, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 5])
    >>> inv
    array([ 0,  1,  5,  8, 42, 99])
    >>> (fw[label_field] == relab).all()
    True
    >>> (inv[relab] == label_field).all()
    True
    >>> relab, fw, inv = relabel_sequential(label_field, offset=5)
    >>> relab
    array([5, 5, 6, 6, 7, 9, 8])
    """
    m = label_field.max()
    if not np.issubdtype(label_field.dtype, np.signedinteger):
        new_type = np.min_scalar_type(int(m))
        label_field = label_field.astype(new_type)
        m = m.astype(new_type)  # Ensures m is an integer
    labels = np.unique(label_field)
    labels0 = labels[labels != 0]
    if m == len(labels0):  # nothing to do, already 1...n labels
        return label_field, labels, labels
    forward_map = np.zeros(m + 1, int)
    forward_map[labels0] = np.arange(offset, offset + len(labels0))
    if not (labels == 0).any():
        labels = np.concatenate(([0], labels))
    inverse_map = np.zeros(offset - 1 + len(labels), dtype=np.intp)
    inverse_map[(offset - 1):] = labels
    relabeled = forward_map[label_field]
    return relabeled, forward_map, inverse_map
