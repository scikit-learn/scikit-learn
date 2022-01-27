# %%
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn import metrics

import numpy as np
import time
import functools
from joblib import Memory
import pandas as pd
import argparse
import streamlit as st
import altair as alt

parser = argparse.ArgumentParser(
    description="Determine whether to save/load the dataframe."
)
parser.add_argument(
    "--save",
    type=str,
    default="",
    metavar="s",
    help="Saves the dataframe to the path if provided",
)
parser.add_argument(
    "--load",
    type=str,
    default="",
    metavar="l",
    help="Loads the dataframe at the path if provided",
)
parser.add_argument(
    "--val",
    action="store_true",
    help="Validates convergence of trust-ncg against lbfgs",
)

args = parser.parse_args()


def main() -> pd.DataFrame:
    # %%
    m = Memory(cachedir=".")

    # %%
    """
    TOY_DSETS = {
        'iris':functools.partial(datasets.load_iris,return_X_y=True),
        'wine':functools.partial(datasets.load_wine,return_X_y=True),
        'breast_cancer':functools.partial(datasets.load_breast_cancer,return_X_y=True),
        'news':functools.partial(datasets.fetch_20newsgroups_vectorized,return_X_y=True)
        }
    """
    TOY_DSETS = {}
    BASE_SYNTH_DSETS = {
        "base": datasets.make_classification,
        "dense": functools.partial(datasets.make_classification, n_informative=15),
        "redundant": functools.partial(datasets.make_classification, n_redundant=15),
        "sparse": functools.partial(datasets.make_classification, n_features=200),
        "repeated": functools.partial(datasets.make_classification, n_repeated=15),
    }
    N_SAMPLES = (10, 50, 250, 1000, 5000, 15000)
    N_CLASSES = (2, 5, 15)
    SYNTH_DSETS = {
        f"{base}_{n_samples}_{n_classes}": functools.partial(
            BASE_SYNTH_DSETS[base], n_samples=n_samples
        )
        for base in BASE_SYNTH_DSETS
        for n_samples in N_SAMPLES
        for n_classes in N_CLASSES
    }
    DSETS = TOY_DSETS | SYNTH_DSETS

    # %%
    @m.cache()
    def make_dset(name: str):
        return DSETS[name]()

    # %%
    MULTI_CLASS = ("ovr", "multinomial")
    SOLVER = ("lbfgs", "trust-ncg", "sag", "saga", "newton-cg")
    PENALTY = ("l2", "none")
    CONFIG = tuple(
        {"multi_class": mc, "solver": s, "penalty": p}
        for mc in MULTI_CLASS
        for s in SOLVER
        for p in PENALTY
    )
    _LogisticRegression = functools.partial(LogisticRegression, max_iter=1000)
    scaler = StandardScaler()

    # %%
    def build_row(X, y, conf, dset: str):
        lr = _LogisticRegression(**conf)
        t_start = time.perf_counter()
        lr.fit(X, y)
        t_total = time.perf_counter() - t_start
        base, n_samples, n_classes = dset.split("_")
        return {
            "cpu_time": t_total,
            "n_iter_": np.mean(lr.n_iter_),
            "NLL": metrics.log_loss(lr.predict(X), y),
            "dset": base,
            "n_samples": n_samples,
            "n_classes": n_classes,
            **conf,
        }

    # %%
    df = None
    if args.load:
        df = on_load(args.load)
    else:
        data = []
        for i, dset in enumerate(DSETS):
            print(f"Progress: {i*100./len(DSETS):.2f}")
            X, y = make_dset(dset)
            X = scaler.fit_transform(X)
            data.extend([build_row(X, y, conf, dset) for conf in CONFIG])
        df = pd.DataFrame(data)

    if args.save and not args.load:
        df.to_pickle(args.save)
        print(f"Dataframe saved to {args.save}")

    dset_option = st.selectbox(
        "Which dataset would you like to view?", ("base", "dense", "sparse")
    )
    y_option = st.selectbox(
        "Which metric would you like charted?", ("NLL", "cpu_time", "n_iter_")
    )
    n_classes_option = st.selectbox(
        "Choose the number of classes in the dataset.", ("2", "5", "15")
    )
    multi_class_option = st.selectbox(
        "Which multi_class strategy should the solvers use?", ("ovr", "multinomial")
    )
    penalty_option = st.selectbox(
        "Which penalty should the solvers use?", ("none", "l2")
    )

    conditions = (
        (df["n_classes"] == n_classes_option)
        & (df["multi_class"] == multi_class_option)
        & (df["penalty"] == penalty_option)
        & (df["dset"] == dset_option)
    )
    chart = (
        alt.Chart(df[conditions], height=500, width=900)
        .mark_line()
        .encode(
            x="n_samples",
            y=y_option,
            color="solver",
        )
    )

    st.altair_chart(chart)

    def validate():
        count_match = 0
        count = 0
        for dset in DSETS:
            X, y = make_dset(dset)
            for conf in CONFIG:
                count += 1
                conf_ = dict(conf)
                del conf_["solver"]
                lr_trust = LogisticRegression(solver="trust-ncg", **conf_)
                lr_lbfgs = LogisticRegression(solver="lbfgs", **conf_)
                lr_trust.fit(X, y)
                lr_lbfgs.fit(X, y)
                if np.isclose(lr_trust.coef_, lr_lbfgs.coef_).min():
                    count_match += 1
        return count_match, count

    if args.val:
        count_match, count = validate()
        print(count_match * 100.0 / count)


def on_load(pth: str) -> pd.DataFrame:
    print(f"Dataframe loaded from {args.load}")
    return pd.read_pickle(pth)


if __name__ == "__main__":
    main()
