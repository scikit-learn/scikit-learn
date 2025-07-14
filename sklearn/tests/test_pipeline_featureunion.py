import pandas as pd
from sklearn.pipeline import FeatureUnion
from sklearn.preprocessing import FunctionTransformer

def test_featureunion_series_output_set_output():
    X = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

    def pick_series(x):
        return x["a"]  # returns a Series

    union = FeatureUnion(
        [("id", FunctionTransformer(pick_series, feature_names_out="one-to-one"))]
    ).set_output(transform="pandas")

    Xt = union.fit_transform(X)
    assert isinstance(Xt, pd.DataFrame)
    assert list(Xt.columns) == ["id__a"] 