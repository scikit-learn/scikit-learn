import textwrap

import pytest

from sklearn.utils.testing import assert_run_python_script


# We are deprecating importing anything that isn't in an __init__ file and
# remaming most file.py into _file.py.
# This test makes sure imports are still possible but deprecated, with the
# appropriate error message.

@pytest.mark.parametrize('deprecated_path, importee', (
    ('sklearn.neural_network.rbm', 'BernoulliRBM'),
    ('sklearn.neural_network.multilayer_perceptron', 'MLPClassifier'),

    ('sklearn.utils.mocking', 'MockDataFrame'),
    ('sklearn.utils.weight_vector', 'WeightVector'),
    ('sklearn.utils.seq_dataset', 'ArrayDataset32'),
    ('sklearn.utils.fast_dict', 'IntFloatDict'),

    ('sklearn.cluster.affinity_propagation_', 'AffinityPropagation'),
    ('sklearn.cluster.bicluster', 'SpectralBiclustering'),
    ('sklearn.cluster.birch', 'Birch'),
    ('sklearn.cluster.dbscan_', 'DBSCAN'),
    ('sklearn.cluster.hierarchical', 'FeatureAgglomeration'),
    ('sklearn.cluster.k_means_', 'KMeans'),
    ('sklearn.cluster.mean_shift_', 'MeanShift'),
    ('sklearn.cluster.optics_', 'OPTICS'),
    ('sklearn.cluster.spectral', 'SpectralClustering'),
))
def test_import_is_deprecated(deprecated_path, importee):
    # Make sure that "from deprecated_path import importee" is still possible
    # but raises a warning
    # We only need one entry per file, no need to check multiple imports from
    # the same file.

    # TODO: remove in 0.24

    expected_message = (
        "The {deprecated_path} module is  deprecated in version "
        "0.22 and will be removed in version 0.24. "
        "The corresponding classes / functions "
        "should instead be imported from .*. "
        "Anything that cannot be imported from .* is now "
        "part of the private API."
    ).format(deprecated_path=deprecated_path)

    script = """
    import pytest
    from sklearn.exceptions import SklearnDeprecationWarning

    with pytest.warns(SklearnDeprecationWarning,
                      match="{expected_message}"):
        from {deprecated_path} import {importee}
    """.format(
        expected_message=expected_message,
        deprecated_path=deprecated_path,
        importee=importee
    )
    assert_run_python_script(textwrap.dedent(script))
