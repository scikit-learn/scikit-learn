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
import hvplot.pandas

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
SOLVER = ("trust-ncg", "lbfgs", "sag", "saga", "newton-cg")
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
data = []
for i, dset in enumerate(DSETS):
    print(f"Progress: {i*100./len(DSETS):.2f}")
    X, y = make_dset(dset)
    X = scaler.fit_transform(X)
    data.extend([build_row(X, y, conf, dset) for conf in CONFIG])
# %%
df = pd.DataFrame(data)
plot = df.hvplot.line(
    x="n_samples",
    y="NLL",
    groupby=["solver", "n_classes", "multi_class", "penalty", "dset"],
)
hvplot.show(plot)
