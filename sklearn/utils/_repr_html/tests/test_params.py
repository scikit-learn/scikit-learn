import pytest

from sklearn import config_context
from sklearn.utils._repr_html.params import (
    ParamsDict,
    _generate_link_to_param_doc,
    _params_html_repr,
    _read_params,
)


def test_params_dict_content():
    """Check the behavior of the ParamsDict class."""
    params = ParamsDict({"a": 1, "b": 2})
    assert params["a"] == 1
    assert params["b"] == 2
    assert params.non_default == ()

    params = ParamsDict({"a": 1, "b": 2}, non_default=("a",))
    assert params["a"] == 1
    assert params["b"] == 2
    assert params.non_default == ("a",)


def test_params_dict_repr_html_():
    params = ParamsDict({"a": 1, "b": 2}, non_default=("a",), estimator_class="")
    out = params._repr_html_()
    assert "<summary>Parameters</summary>" in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            params._repr_html_()


def test_params_dict_repr_mimebundle():
    params = ParamsDict({"a": 1, "b": 2}, non_default=("a",), estimator_class="")
    out = params._repr_mimebundle_()

    assert "text/plain" in out
    assert "text/html" in out
    assert "<summary>Parameters</summary>" in out["text/html"]
    assert out["text/plain"] == "{'a': 1, 'b': 2}"

    with config_context(display="text"):
        out = params._repr_mimebundle_()
        assert "text/plain" in out
        assert "text/html" not in out


def test_read_params():
    """Check the behavior of the `_read_params` function."""
    out = _read_params("a", 1, tuple())
    assert out["param_type"] == "default"
    assert out["param_name"] == "a"
    assert out["param_value"] == "1"

    # check non-default parameters
    out = _read_params("a", 1, ("a",))
    assert out["param_type"] == "user-set"
    assert out["param_name"] == "a"
    assert out["param_value"] == "1"

    # check that we escape html tags
    tag_injection = "<script>alert('xss')</script>"
    out = _read_params("a", tag_injection, tuple())
    assert (
        out["param_value"]
        == "&quot;&lt;script&gt;alert(&#x27;xss&#x27;)&lt;/script&gt;&quot;"
    )
    assert out["param_name"] == "a"
    assert out["param_type"] == "default"


def test_params_html_repr():
    """Check returned HTML template"""
    params = ParamsDict({"a": 1, "b": 2}, estimator_class="")
    assert "parameters-table" in _params_html_repr(params)
    assert "estimator-table" in _params_html_repr(params)


def test_params_html_repr_with_doc_links():
    """Test `_params_html_repr` with valid and invalid doc links."""

    class MockEstimator:
        """A fake estimator class with a docstring used for testing.

        a : int
            Description of a.
        b : str
        """

        __module__ = "sklearn.mock_module"
        __qualname__ = "MockEstimator"

    params = ParamsDict(
        {"a": 1, "b": "value"},
        non_default=("a",),
        estimator_class=MockEstimator,
        doc_link="mock_module.MockEstimator.html",
    )
    html_output = _params_html_repr(params)

    # Check that the doc links are correctly generated
    assert "Documentation for a" in html_output
    assert "Documentation for b" in html_output
    assert "sk-estimator-doc-link" in html_output
    # Check that the doc links contain the correct URL fragments
    assert "mock_module.MockEstimator.html#:~:text=a,-int" in html_output
    assert "mock_module.MockEstimator.html#:~:text=b,-str" in html_output


def test_params_html_repr_without_doc_links():
    """Test `_params_html_repr` when `link_to_param_doc` returns None."""

    class MockEstimatorWithoutDoc:
        __module__ = "sklearn.mock_module"
        __qualname__ = "MockEstimatorWithoutDoc"
        # No docstring defined on this test class.

    params = ParamsDict(
        {"a": 1, "b": "value"},
        non_default=("a",),
        estimator_class=MockEstimatorWithoutDoc,
    )
    html_output = _params_html_repr(params)
    # Check that no doc links are generated
    assert "?" not in html_output
    assert "Documentation for a" not in html_output
    assert "Documentation for b" not in html_output


def test_generate_link_to_param_doc_basic():
    class MockEstimator:
        __doc__ = """
        alpha : float
            Regularization strength.
        beta : int
            Some integer parameter.
        """

    doc_link = "mock_module.MockEstimator.html"
    # Test for 'alpha'
    url = _generate_link_to_param_doc(MockEstimator, "alpha", doc_link)
    assert url == "mock_module.MockEstimator.html#:~:text=alpha,-float"
    # Test for 'beta'
    url = _generate_link_to_param_doc(MockEstimator, "beta", doc_link)
    assert url == "mock_module.MockEstimator.html#:~:text=beta,-int"


def test_generate_link_to_param_doc_param_not_found():
    class MockEstimator:
        __doc__ = """
        alpha : float
            Regularization strength.
        """

    doc_link = "mock_module.MockEstimator.html"
    url = _generate_link_to_param_doc(MockEstimator, "gamma", doc_link)
    assert url is None


def test_generate_link_to_param_doc_special_characters():
    class MockEstimator:
        __doc__ = """
        param with space : {'stay', 'leave'} default='don't know'
            Parameter with space in name.
        """

    doc_link = "mock_module.MockEstimator.html"
    url = _generate_link_to_param_doc(MockEstimator, "param with space", doc_link)
    # Spaces should be percent-encoded
    assert "param%20with%20space" in url
    assert "-%7B%27stay%27%2C%20%27leave%27%7D" in url


def test_generate_link_to_param_doc_empty_docstring():
    class MockEstimator:
        __doc__ = ""

    doc_link = "mock_module.MockEstimator.html"
    url = _generate_link_to_param_doc(MockEstimator, "alpha", doc_link)
    assert url is None
