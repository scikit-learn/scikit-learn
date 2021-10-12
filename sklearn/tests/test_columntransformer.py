import pandas as pd

data = {
    "x1": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    "x2": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "y": [0, 1, 2, 1, 2, 0, 0, 2, 1, 0],
    "equal_sample_weight": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "sample_weight": [.5, 1, .5, 1, .5, 1, .5, 1, .5, 1]
}

df = pd.DataFrame.from_dict(data)

# check outside of columntransformer
from sklearn.preprocessing import StandardScaler

sc1 = StandardScaler()
sc1_xWeight = sc1.fit_transform(X=df[['x1']], y=df['y'])
sc1_wWeight = sc1.fit_transform(X=df[['x1']], y=df['y'], sample_weight=df['sample_weight'])
sc1_wEqualWeight = sc1.fit_transform(X=df[['x1']], y=df['y'], sample_weight=df['equal_sample_weight'])

print(f"xWeight: {sc1_xWeight}")
print(f"wWeight: {sc1_wWeight}")
print(f"wEqualWeight: {sc1_wEqualWeight}")
from sklearn.utils._testing import assert_array_equal
assert_array_equal(sc1_xWeight, sc1_wEqualWeight, err_msg= "These should be equal")
assert_array_equal(sc1_xWeight, sc1_wWeight, err_msg= "You should see this message because we shouldn't be equal")

# test as part of a columntransformer
from sklearn.compose import ColumnTransformer
ct = ColumnTransformer(
    [
        ("standard_scaler", StandardScaler(), ["x1"]),
    ]
)

ct_xWeight = ct.fit_transform(X=df[['x1']], y=df['y'])
fit_params = {'standard_scaler__sample_weight': df['sample_weight']}
ct_wWeight = ct.fit_transform(X=df[['x1']], y=df['y'], **fit_params)
fit_params = {'standard_scaler__sample_weight': df['equal_sample_weight']}
ct_wEqualWeight = ct.fit_transform(X=df[['x1']], y=df['y'], **fit_params)

assert_array_equal(sc1_xWeight, ct_xWeight, err_msg= "These should be equal")
assert_array_equal(sc1_wWeight, ct_wWeight, err_msg= "These should be equal")
assert_array_equal(sc1_wEqualWeight, ct_wEqualWeight, err_msg= "These should be equal")

# Test with Pipeline
from sklearn.pipeline import Pipeline
ct2 = ColumnTransformer(
    [
        ("standard_scaler", StandardScaler(), ["x1"]),
    ]
)
pt = Pipeline([('ct_step', ct2), ('passthrough_test',"passthrough" )])
pt_xWeight = pt.fit_transform(X=df[['x1']], y=df['y'])
fit_params = {'ct_step__standard_scaler__sample_weight': df['sample_weight']}
pt_wWeight = pt.fit_transform(X=df[['x1']], y=df['y'], **fit_params)
fit_params = {'ct_step__standard_scaler__sample_weight': df['equal_sample_weight']}
pt_wEqualWeight = pt.fit_transform(X=df[['x1']], y=df['y'], **fit_params)

assert_array_equal(sc1_xWeight, pt_xWeight, err_msg= "These should be equal")
assert_array_equal(sc1_wWeight, pt_wWeight, err_msg= "These should be equal")
assert_array_equal(sc1_wEqualWeight, pt_wEqualWeight, err_msg= "These should be equal")
