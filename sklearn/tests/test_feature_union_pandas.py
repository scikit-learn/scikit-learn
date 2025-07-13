import pandas as pd
import numpy as np
from sklearn.pipeline import FunctionTransformer, FeatureUnion

def test_feature_union_series_output():
    df = pd.DataFrame({
        'id': [1,2,1,2],
        'x': [10, 20, 30, 40]
    })
    def double_x(df):
        return df['x'] * 2

    fu = FeatureUnion([('d', FunctionTransformer(double_x))])\
        .set_output(transform='pandas')
    res = fu.fit_transform(df)

    assert isinstance(res, pd.DataFrame)
    # Series 返回时保持原名 'x'
    assert list(res.columns) == ['x']
    assert len(res) == len(df)
    assert (res['x'] == df['x'] * 2).all()
