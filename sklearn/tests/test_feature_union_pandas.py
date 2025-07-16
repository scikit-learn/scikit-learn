import pytest

from sklearn.pipeline import FeatureUnion, FunctionTransformer

pd = pytest.importorskip("pandas")


def _double(df):
    return pd.DataFrame({"double_age": df["age"] * 2})


def _triple(df):
    return pd.DataFrame({"triple_age": df["age"] * 3})


@pytest.mark.filterwarnings("ignore:With transform")
def test_feature_union_series_output():
    X = pd.DataFrame({"id": [1, 2, 1], "age": [10, 20, 30]})

    fu = FeatureUnion(
        [
            ("double", FunctionTransformer(_double)),
            ("triple", FunctionTransformer(_triple)),
        ]
    ).set_output(transform="pandas")

    Xt = fu.fit_transform(X)
    assert list(Xt.columns) == ["double_age", "triple_age"]
    assert Xt.iloc[0, 0] == 20 and Xt.iloc[0, 1] == 30
