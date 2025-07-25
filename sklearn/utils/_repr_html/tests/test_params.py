import pytest

from sklearn import config_context
from sklearn.utils._repr_html.params import ParamsDict, _params_html_repr, _read_params


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
    params = ParamsDict({"a": 1, "b": 2}, non_default=("a",))
    out = params._repr_html_()
    assert "<summary>Parameters</summary>" in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            params._repr_html_()


def test_params_dict_repr_mimebundle():
    params = ParamsDict({"a": 1, "b": 2}, non_default=("a",))
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
    params = ParamsDict({"a": 1, "b": 2})
    assert "parameters-table" in _params_html_repr(params)
    assert "estimator-table" in _params_html_repr(params)
