import numpy as np
import pytest

from sklearn import config_context
from sklearn.utils._repr_html.fitted_attributes import AttrsDict, _fitted_attr_html_repr

fitted_attrs = AttrsDict(
    fitted_attrs={"a": ("int", 6), "b": ("ndarray", (1,), np.dtype("float64"), 8)}
)


def test_fitted_attrs_dict_repr_html_():
    out = fitted_attrs._repr_html_()
    assert "<summary>Fitted attributes</summary>" in out
    assert "<td>int</td>" in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            fitted_attrs._repr_html_()


def test_fitted_attrs_dict_repr_mimebundle():
    out = fitted_attrs._repr_mimebundle_()

    assert "text/plain" in out
    assert "text/html" in out
    assert "<summary>Fitted attributes</summary>" in out["text/html"]
    plain_text = "{'a': ('int', 6), 'b': ('ndarray', (1,), dtype('float64'), 8)}"
    assert out["text/plain"] == plain_text
    with config_context(display="text"):
        out = fitted_attrs._repr_mimebundle_()
        assert "text/plain" in out
        assert "text/html" not in out


def test_fitted_attr_html_repr():
    out = _fitted_attr_html_repr(fitted_attrs)
    assert "<summary>Fitted attributes</summary>" in out
    assert '<table class="parameters-table">' in out
