import scipy as sp
import numpy as np
from nose.tools import assert_equal
from ..kmeans_coder import KMeansCoder
from ...feature_extraction.image import PatchExtractor


def test_kmeans_coder_shape():
    lena = sp.lena()[np.newaxis, :, :]
    patches = PatchExtractor(patch_size=(4, 4),
                             max_patches=int(1e2)).transform(lena)
    patches = patches.reshape(len(patches), -1)

    encoder = KMeansCoder(n_atoms=12, max_iter=3,
                          transform_algorithm='threshold', transform_alpha=1.)
    encoder.fit(patches)
    assert_equal(encoder.components_.shape, (12, 4 * 4))
    code = encoder.transform(patches)
    assert_equal(code.shape, (len(patches), 12))
