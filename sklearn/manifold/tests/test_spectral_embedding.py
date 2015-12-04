from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns_message
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.manifold import SpectralEmbedding, spectral_embedding
import numpy as np

@ignore_warnings
def test_spectral_embeding_import():
    random_state = np.random.RandomState(36)
    data = random_state.randn(10, 30)
    sims = rbf_kernel(data)

    assert_warns_message(DeprecationWarning, "spectral_embedding is deprecated",
                         spectral_embedding, sims)
    assert_warns_message(DeprecationWarning, "SpectralEmbedding is deprecated",
                         SpectralEmbedding)
