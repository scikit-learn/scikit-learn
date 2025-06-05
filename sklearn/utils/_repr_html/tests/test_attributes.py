from sklearn.utils._repr_html.attributes import AttrsDict


def test_fitted_attrs_dict_content():
    fitted_attrs = AttrsDict({"a": int, "b": bool})
    assert fitted_attrs["a"] == int
    assert fitted_attrs["b"] == bool


def test_fitted_attrs_dict_repr_html_():
    fitted_attrs = AttrsDict({"a": int, "b": bool})
    out = fitted_attrs._repr_html_()
    assert "<summary>Fitted attributes</summary>" in out
    assert "<td><class 'int'></td>" in out
