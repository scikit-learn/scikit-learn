import numpy as np
import pandas as pd
from sklearn.pipeline import FeatureUnion, FunctionTransformer


def _double(x): return (x["age"] * 2).rename("double")
def _triple(x): return (x["age"] * 3).rename("triple")


def test_feature_union_series_output():
    X = pd.DataFrame({"id": [1, 2, 1], "age": [10, 20, 30]})

    fu = FeatureUnion([
        ("double",  FunctionTransformer(_double)),
        ("triple",  FunctionTransformer(_triple)),
    ]).set_output(transform="pandas")

    Xt = fu.fit_transform(X)
    assert list(Xt.columns) == ["double", "triple"]
    assert Xt.iloc[0, 0] == 20 and Xt.iloc[0, 1] == 30
