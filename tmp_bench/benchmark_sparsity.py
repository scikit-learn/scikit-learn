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
    def make_classification_scaled():
        X, y = datasets.make_classification(
            n_samples=15000, n_informative=15, n_classes=15
        )
        return StandardScaler().fit_transform(X), y

    DSETS = {
        "dense": make_classification_scaled,
        "sparse": functools.partial(
            datasets.fetch_20newsgroups_vectorized, return_X_y=True
        ),
    }

    # %%
    @m.cache()
    def make_dset(name: str):
        return DSETS[name]()

    # %%
    MULTI_CLASS = ("ovr", "multinomial")
    SOLVER = ("lbfgs", "trust-ncg")
    PENALTY = ("l2", "none")
    CONFIG = tuple(
        {"multi_class": mc, "solver": s, "penalty": p}
        for mc in MULTI_CLASS
        for s in SOLVER
        for p in PENALTY
    )
    _LogisticRegression = functools.partial(LogisticRegression, max_iter=1000)

    # %%
    def build_row(X, y, conf, dset: str):
        lr = _LogisticRegression(**conf)
        t_start = time.perf_counter()
        lr.fit(X, y)
        t_total = time.perf_counter() - t_start
        return {
            "cpu_time": t_total,
            "n_iter_": np.mean(lr.n_iter_),
            "NLL": metrics.log_loss(lr.predict(X), y),
            "dset": dset,
            **conf,
        }

    # %%
    df = None
    if args.load:
        df = on_load(args.load)
    else:
        data = []
        for dset in DSETS:
            X, y = make_dset(dset)
            data.extend([build_row(X, y, conf, dset) for conf in CONFIG])
        df = pd.DataFrame(data)

    if args.save and not args.load:
        df.to_pickle(args.save)
        print(f"Dataframe saved to {args.save}")

    dset_option = st.selectbox(
        "Which dataset would you like to view?", ("dense", "sparse")
    )
    y_option = st.selectbox(
        "Which metric would you like charted?", ("NLL", "cpu_time", "n_iter_")
    )
    multi_class_option = st.selectbox(
        "Which multi_class strategy should the solvers use?", ("ovr", "multinomial")
    )
    penalty_option = st.selectbox(
        "Which penalty should the solvers use?", ("none", "l2")
    )

    conditions = (
        (df["multi_class"] == multi_class_option)
        & (df["penalty"] == penalty_option)
        & (df["dset"] == dset_option)
    )
    chart = (
        alt.Chart(df[conditions], height=500, width=900)
        .mark_bar()
        .encode(x="solver", y=y_option, column="multi_class", color="penalty")
        .properties(title=f"{y_option} using the {dset_option}-style dataset")
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

# %%
