import textwrap

import pytest

from sklearn.utils.testing import assert_run_python_script


# We are deprecating importing anything that isn't in an __init__ file and
# remaming most file.py into _file.py.
# This test makes sure imports are still possible but deprecated, with the
# appropriate error message.

@pytest.mark.parametrize('deprecated_path, importee', (
    ('sklearn.tree.tree', 'DecisionTreeClassifier'),
    ('sklearn.tree.tree', 'DecisionTreeRegressor'),
    ('sklearn.tree.tree', 'ExtraTreeClassifier'),
    ('sklearn.tree.tree', 'ExtraTreeRegressor'),
    ('sklearn.tree.export', 'export_graphviz'),
    ('sklearn.tree.export', 'plot_tree'),
    ('sklearn.tree.export', 'export_text'),

    ('sklearn.ensemble.base', 'BaseEnsemble'),
    ('sklearn.ensemble.forest', 'RandomForestClassifier'),
    ('sklearn.ensemble.forest', 'RandomForestRegressor'),
    ('sklearn.ensemble.forest', 'RandomTreesEmbedding'),
    ('sklearn.ensemble.forest', 'ExtraTreesClassifier'),
    ('sklearn.ensemble.forest', 'ExtraTreesRegressor'),
    ('sklearn.ensemble.bagging', 'BaggingClassifier'),
    ('sklearn.ensemble.bagging', 'BaggingRegressor'),
    ('sklearn.ensemble.iforest', 'IsolationForest'),
    ('sklearn.ensemble.weight_boosting', 'AdaBoostClassifier'),
    ('sklearn.ensemble.weight_boosting', 'AdaBoostRegressor'),
    ('sklearn.ensemble.gradient_boosting', 'GradientBoostingClassifier'),
    ('sklearn.ensemble.gradient_boosting', 'GradientBoostingRegressor'),
    ('sklearn.ensemble.voting', 'VotingClassifier'),
    ('sklearn.ensemble.voting', 'VotingRegressor'),

    ('sklearn.neural_network.rbm', 'BernoulliRBM'),
    ('sklearn.neural_network.multilayer_perceptron', 'MLPClassifier'),

    ('sklearn.utils.mocking', 'MockDataFrame'),
))
def test_import_is_deprecated(deprecated_path, importee):
    # Make sure that "from deprecated_path import importee" is still possible
    # but raises a warning

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

    with pytest.warns(DeprecationWarning,
                      match="{expected_message}"):
        from {deprecated_path} import {importee}
    """.format(
        expected_message=expected_message,
        deprecated_path=deprecated_path,
        importee=importee
    )
    assert_run_python_script(textwrap.dedent(script))
