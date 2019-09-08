import textwrap
from sklearn.utils.testing import assert_run_python_script


# This file makes sure importing from sklearn.neural_net.rbm or
# sklearn.neural_net.multilayer_perceptron raises a deprecation warning.


def test_rbm():
    script = """
    import pytest

    with pytest.warns(DeprecationWarning,
                      match="not work anymore in version 0.24"):
        from sklearn.neural_network.rbm import BernoulliRBM
    """
    assert_run_python_script(textwrap.dedent(script))


def test_multilayer_perceptron():
    script = """
    import pytest

    with pytest.warns(DeprecationWarning,
                      match="not work anymore in version 0.24"):
        from sklearn.neural_network.multilayer_perceptron import MLPClassifier
    """
    assert_run_python_script(textwrap.dedent(script))
