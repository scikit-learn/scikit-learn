import pytest

from sklearn import config_context
from sklearn.base import BaseEstimator


class MyEstimator(BaseEstimator):
    def __init__(self, l1=0, empty=None):
        self.l1 = l1
        self.empty = empty


def test_repr_html_():
    est = MyEstimator()
    out = est._get_params_html(deep=False)._repr_html_()
    assert "<summary>Parameters</summary>" in out

    with config_context(display="text"):
        msg = "_repr_html_ is only defined when"
        with pytest.raises(AttributeError, match=msg):
            est._get_params_html(deep=False)._repr_html_()


def test_ReprHTMLMixin_repr_mimebundle():
    est = MyEstimator()
    out = est._get_params_html(deep=False)._repr_mimebundle_()

    assert "text/plain" in out
    assert "text/html" in out
    assert "<summary>Parameters</summary>" in out["text/html"]
    assert out["text/plain"] == "{'l1': 0, 'empty': None}"

    with config_context(display="text"):
        out = est._get_params_html(deep=False)._repr_mimebundle_()
        assert "text/plain" in out
        assert "text/html" not in out
