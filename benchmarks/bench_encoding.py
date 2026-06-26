import argparse
from timeit import timeit

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, TargetEncoder

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scale", type=int, default=100)
parser.add_argument("-e", "--encoders", default="onehot,ordinal,target")
parser.add_argument(
    "-f", "--format", default="string", choices=["string", "categorical", "integer"]
)
parser.add_argument("-w", "--with-numericals", action="store_true")
parser.add_argument("--polars", action="store_true")
args = parser.parse_args()
scale = args.scale

# fetch & preprocess data:
X, y = fetch_openml(data_id=42165, as_frame=True, return_X_y=True)

string_cols = X.select_dtypes(include=["object", "string"]).columns

if not args.with_numericals:
    X = X.loc[:, string_cols]

if args.format in ("categorical", "integer"):
    X[string_cols] = X[string_cols].astype("category")
if args.format == "integer":
    for col in string_cols:
        X[col] = X[col].cat.codes

X = pd.concat([X] * scale)
y = np.concat([y.astype(np.float32).values] * scale)


def _as_polars(X, string_cols, format):
    import polars as pl

    string_cols = set(string_cols)
    data = {}
    for col in X.columns:
        if col in string_cols and format != "integer":
            data[col] = [
                None if missing else str(value)
                for value, missing in zip(X[col], X[col].isna())
            ]
        else:
            data[col] = X[col].to_list()

    X_polars = pl.DataFrame(data)
    if format == "categorical":
        categorical_cols = [col for col in X_polars.columns if col in string_cols]
        X_polars = X_polars.with_columns(pl.col(categorical_cols).cast(pl.Categorical))
    return X_polars


if args.polars:
    X = _as_polars(X, string_cols, args.format)

# define encoders:
encoders = {
    "onehot": OneHotEncoder(handle_unknown="ignore", max_categories=10),
    "ordinal": OrdinalEncoder(max_categories=10),
    "target": TargetEncoder(target_type="continuous"),
}


# small measurement helper:
def measure_ms(func):
    n_repeats = max(300 // scale, 1)
    return round(timeit(func, number=n_repeats) * 1000 / n_repeats, 1)


# actual benchmark:
for name in args.encoders.split(","):
    encoder = encoders[name]
    if args.with_numericals:
        encoder = ColumnTransformer(
            transformers=[("categorical", encoder, list(string_cols))],
            remainder="passthrough",
            sparse_threshold=0.9,
        )
    print(name, "fit", measure_ms(lambda: encoder.fit(X, y)), "ms")
    encoder.fit(X, y)
    print(name, "transform", measure_ms(lambda: encoder.transform(X)), "ms")
    print(name, "fit_transform", measure_ms(lambda: encoder.fit_transform(X, y)), "ms")
