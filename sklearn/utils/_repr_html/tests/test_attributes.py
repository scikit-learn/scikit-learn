import pytest

from sklearn import config_context
from sklearn.utils._repr_html.fitted_attributes import AttrsDict, _fitted_attr_html_repr


def test_fitted_attrs_dict_content():
    fitted_attrs = AttrsDict({"a": int, "b": bool})
    assert fitted_attrs["a"] == int
    assert fitted_attrs["b"] == bool


def test_fitted_attrs_dict_repr_html_():
    fitted_attrs = AttrsDict({"a": int, "b": bool})
    out = fitted_attrs._repr_html_()
    assert "<summary>Fitted attributes</summary>" in out
    assert "<td><class 'int'></td>" in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            fitted_attrs._repr_html_()


def test_fitted_attrs_dict_repr_mimebundle():
    fitted_attrs = AttrsDict({"a": int, "b": float})
    out = fitted_attrs._repr_mimebundle_()

    assert "text/plain" in out
    assert "text/html" in out
    assert "<summary>Fitted attributes</summary>" in out["text/html"]
    assert out["text/plain"] == "{'a': <class 'int'>, 'b': <class 'float'>}"

    with config_context(display="text"):
        out = fitted_attrs._repr_mimebundle_()
        assert "text/plain" in out
        assert "text/html" not in out


def test_fitted_attr_html_repr():
    out = _fitted_attr_html_repr({"a": int, "b": float})
    assert "<summary>Fitted attributes</summary>" in out
    assert '<table class="body-table">' in out
    assert '<tr class="default">' in out
    assert "<class 'float'></td>" in out
