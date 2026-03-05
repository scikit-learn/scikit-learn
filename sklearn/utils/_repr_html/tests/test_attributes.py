import re

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn import config_context
from sklearn.linear_model import LogisticRegression
from sklearn.utils._repr_html.fitted_attributes import AttrsDict, _fitted_attr_html_repr

fitted_attrs = AttrsDict(
    fitted_attrs={
        "a": {"type_name": "int", "value": 6},
        "b": {
            "type_name": "ndarray",
            "shape": (1,),
            "dtype": np.dtype("float64"),
            "value": 8,
        },
    }
)


def test_fitted_attrs_dict_repr_html_():
    out = fitted_attrs._repr_html_()
    assert "<summary>Fitted attributes</summary>" in out
    assert '<td class="fitted-att-type">int</td>' in out
    assert '<td class="fitted-att-type">ndarray[float64](1,)</td>' in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            fitted_attrs._repr_html_()


def test_fitted_attrs_dict_repr_mimebundle():
    out = fitted_attrs._repr_mimebundle_()

    assert "text/plain" in out
    assert "text/html" in out
    assert "<summary>Fitted attributes</summary>" in out["text/html"]
    plain_text = (
        "{'a': {'type_name': 'int', 'value': 6}, "
        "'b': {'type_name': 'ndarray', 'shape': (1,), "
        "'dtype': dtype('float64'), 'value': 8}}"
    )
    assert out["text/plain"] == plain_text
    with config_context(display="text"):
        out = fitted_attrs._repr_mimebundle_()
        assert "text/plain" in out
        assert "text/html" not in out


def test_fitted_attr_html_repr():
    out = _fitted_attr_html_repr(fitted_attrs)
    assert "<summary>Fitted attributes</summary>" in out
    assert '<table class="parameters-table">' in out


def test_pandas_column_names():
    pd = pytest.importorskip("pandas")
    X = pd.DataFrame({"A": [0, 2, 4], "B": [3, 4, 5], "C": [3, 4, 4]})
    y = pd.DataFrame({"y": [1, 2, 3]})
    model = LogisticRegression()
    model.fit(X, y)

    assert_array_equal(model.feature_names_in_, np.array(["A", "B", "C"], dtype=object))


def test_pandas_series_fitted_attr():
    pd = pytest.importorskip("pandas")
    fitted_attrs_with_series = AttrsDict(
        fitted_attrs={
            "new_": {
                "type_name": "Series",
                "shape": (3,),
                "dtype": np.dtype("int64"),
                "value": pd.Series({"a": 1, "b": 2, "c": 3}),
            }
        }
    )
    html_fitted_attr_out = _fitted_attr_html_repr(fitted_attrs_with_series)
    expected_html_fitted_attr = (
        r'<div class="estimator-table">'
        r"\s*<details>"
        r"\s*<summary>Fitted attributes</summary>"
        r'\s*<table class="parameters-table">'
        r"\s*<tbody>"
        r"\s*<tr>"
        r"\s*<th>Name</th>"
        r"\s*<th>Type</th>"
        r"\s*<th>Value</th>"
        r"\s*</tr>"
        r'\s*<tr class="default">'
        r'\s*<td class="param">'
        r'\s*<a class="param-doc-link" style="text-decoration:none;">new_</a>'
        r"\s*</td>"
        r'\s*<td class="fitted-att-type">Series\[int64\]\(3,\)</td>'
        r"\s*<td>a\s*1\s*b\s*2...3\s*dtype:\s*int64</td>"
        r"\s*</tr>"
        r"\s*</tbody>"
        r"\s*</table>"
        r"\s*</details>"
        r"\s*</div>"
    )
    assert re.search(expected_html_fitted_attr, html_fitted_attr_out, flags=re.DOTALL)


def test_numpy_fitted_attr():
    html_fitted_attr_out = _fitted_attr_html_repr(fitted_attrs)

    expected_html_fitted_attr = (
        r'<div class="estimator-table">'
        r"\s*<details>"
        r"\s*<summary>Fitted attributes</summary>"
        r'\s*<table class="parameters-table">'
        r"\s*<tbody>"
        r"\s*<tr>"
        r"\s*<th>Name</th>"
        r"\s*<th>Type</th>"
        r"\s*<th>Value</th>"
        r"\s*</tr>"
        r'\s*<tr class="default">'
        r'\s*<td class="param">'
        r'\s*<a class="param-doc-link" style="text-decoration:none;">a</a>'
        r"\s*</td>"
        r'\s*<td class="fitted-att-type">int</td>'
        r"\s*<td>6</td>"
        r"\s*</tr>"
        r'\s*<tr class="default">'
        r'\s*<td class="param">'
        r'\s*<a class="param-doc-link" style="text-decoration:none;">b</a>'
        r"\s*</td>"
        r'\s*<td class="fitted-att-type">ndarray\[float64\]\(1,\)</td>'
        r"\s*<td>8</td>"
        r"\s*</tr>"
        r"\s*</tbody>"
        r"\s*</table>"
        r"\s*</details>"
        r"\s*</div>"
    )
    assert re.search(expected_html_fitted_attr, html_fitted_attr_out, flags=re.DOTALL)
