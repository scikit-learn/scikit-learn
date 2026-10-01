import argparse
from math import ceil
from statistics import median
from timeit import Timer

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, TargetEncoder

DATASETS = {
    # ordered by increasing number of rows:
    "house_prices": 42165,
    "bank_marketing": 45211,
    "amazon_employee_access": 4135,
    "kddcup09_churn": 1112,
    "kick": 41162,
}

parser = argparse.ArgumentParser()
parser.add_argument(
    "-n",
    "--n-samples",
    type=int,
    default=None,
    help="tile/truncate each dataset to this many rows; default keeps each "
    "dataset's native row count",
)
parser.add_argument("-f", "--formats", default="object,string,categorical")
parser.add_argument("-l", "--libraries", default="pandas,polars")
parser.add_argument("-e", "--encoders", default="ordinal,ordinal-max10")
parser.add_argument("-w", "--with-numericals", action="store_true")

parser.add_argument("-o", "--output", default="bench_encoders_results.csv")
args = parser.parse_args()

# column widths for the table-like progress output below
LIB_W = max(len(lib) for lib in args.libraries.split(","))
FMT_W = max(len(fmt) for fmt in args.formats.split(","))
ENC_W = max(len(name) for name in args.encoders.split(","))

# define encoders:
encoders = {
    "onehot": OneHotEncoder(handle_unknown="ignore", max_categories=20),
    "ordinal-max10": OrdinalEncoder(max_categories=10),
    "ordinal": OrdinalEncoder(encoded_missing_value=-1),
    "target": TargetEncoder(target_type="continuous"),
}


def measure_ms(func, prec=1):
    times = Timer(func).repeat(repeat=3, number=1)
    if sum(times) < 0.1:
        times += Timer(func).repeat(repeat=10, number=1)
    return round(median(times) * 1000, prec)


results = []


for dataset_name, data_id in DATASETS.items():
    # fetch & preprocess data:
    X, y = fetch_openml(data_id=data_id, as_frame=True, return_X_y=True)
    # print(len(X))
    if y is None:
        # Some high-cardinality OpenML datasets have no default target. A synthetic
        # continuous target keeps TargetEncoder timings available for those datasets.
        y = np.linspace(0.0, 1.0, X.shape[0], endpoint=False)

    categorical_cols = X.select_dtypes(include=["object", "string", "category"]).columns

    if not args.with_numericals:
        X = X.loc[:, categorical_cols]

    if args.n_samples is None:
        n_samples = X.shape[0]
        X_base = X
        y_arr = np.asarray(y)
    else:
        n_samples = args.n_samples
        repeats = ceil(n_samples / X.shape[0])
        X_base = pd.concat([X] * repeats, ignore_index=True).iloc[:n_samples]
        y_arr = np.concat([np.asarray(y)] * repeats)[:n_samples]

    print(
        f"---- {dataset_name} {n_samples} rows with {len(categorical_cols)}"
        " categorical columns ----"
    )

    for fmt in args.formats.split(","):
        X_fmt = X_base.copy()
        if fmt == "object":
            X_fmt[categorical_cols] = (
                X_fmt[categorical_cols]
                .astype("string")
                .astype(object)
                .replace({pd.NA: np.nan})
            )
        elif fmt == "string":
            X_fmt[categorical_cols] = X_fmt[categorical_cols].astype("string")
        else:  # categorical
            X_fmt[categorical_cols] = X_fmt[categorical_cols].astype("category")

        for lib in args.libraries.split(","):
            if lib == "polars":
                import polars as pl

                # needs pyarrow to work:
                X_run = pl.from_pandas(X_fmt, include_index=False)
                if fmt == "categorical":
                    X_run = X_run.with_columns(
                        pl.col(categorical_cols).cast(pl.Categorical)
                    )
            else:
                X_run = X_fmt

            # actual benchmark:
            for name in args.encoders.split(","):
                label = f"{lib:<{LIB_W}} | {fmt:<{FMT_W}} | {name:<{ENC_W}} || "
                encoder = encoders[name]

                if args.with_numericals:
                    encoder = ColumnTransformer(
                        transformers=[("categorical", encoder, list(categorical_cols))],
                        remainder="passthrough",
                        sparse_threshold=0,
                    )

                row = dict(
                    dataset=dataset_name,
                    encoder=name,
                    lib=lib,
                    format=fmt,
                    n_samples=X_run.shape[0],
                    n_features=X_run.shape[1],
                )

                try:
                    row["fit_ms"] = measure_ms(lambda: encoder.fit(X_run, y_arr))
                    encoder.fit(X_run, y_arr)
                    row["transform_ms"] = measure_ms(lambda: encoder.transform(X_run))
                    if name == "target":
                        # for other encoders, fit_transform is fit+transform
                        row["fit_transform_ms"] = measure_ms(
                            lambda: encoder.fit_transform(X_run, y_arr)
                        )
                    encoder.transform(X_run[:1])
                    row["transform_single_ms"] = measure_ms(
                        lambda: encoder.transform(X_run[:1])
                    )
                    print(
                        f"{label}  f={row['fit_ms']:>6.1f}"
                        f"  t={row['transform_ms']:>6.1f}"
                        f"  t(1)={row['transform_single_ms']:>6.1f}"
                    )
                except Exception as e:
                    row["error"] = f"{type(e).__name__}: {e}"
                    print(f"{label}  FAILED: {row['error']}")
                results.append(row)

results_df = pd.DataFrame(results)
results_df.to_csv(args.output, index=False)
print(f"\nWrote {len(results_df)} rows to {args.output}")
