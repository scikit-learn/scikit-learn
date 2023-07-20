from itertools import product
from timeit import timeit

import numpy as np
import pandas as pd

from sklearn.datasets import fetch_openml
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import OneHotEncoder


def get_data(trunc_ncat):
    # the data is located here: https://www.openml.org/d/4135
    X, y = fetch_openml(data_id=4135, return_X_y=True)
    X = pd.DataFrame(X)

    Xdicts = []
    for trunc in trunc_ncat:
        X_trunc = X % trunc if trunc > 0 else X
        keep_idx = np.array(
            [idx[0] for idx in X_trunc.groupby(list(X.columns)).groups.values()]
        )
        X_trunc = X_trunc.values[keep_idx]
        y_trunc = y[keep_idx]

        X_ohe = OneHotEncoder(categories="auto").fit_transform(X_trunc)

        Xdicts.append({"X": X_trunc, "y": y_trunc, "ohe": False, "trunc": trunc})
        Xdicts.append({"X": X_ohe, "y": y_trunc, "ohe": True, "trunc": trunc})

    return Xdicts


# Training dataset
trunc_factor = [2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 64, 0]
data = get_data(trunc_factor)
results = []
# Loop over classifiers and datasets
for Xydict, clf_type in product(data, [RandomForestClassifier, ExtraTreesClassifier]):
    # Can't use non-truncated categorical data with RandomForest
    # and it becomes intractable with too many categories
    if (
        clf_type is RandomForestClassifier
        and not Xydict["ohe"]
        and (not Xydict["trunc"] or Xydict["trunc"] > 16)
    ):
        continue

    X, y = Xydict["X"], Xydict["y"]
    tech = "One-hot" if Xydict["ohe"] else "NOCATS"
    trunc = "truncated({})".format(Xydict["trunc"]) if Xydict["trunc"] > 0 else "full"
    cat = "none" if Xydict["ohe"] else "all"
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=17).split(X, y)

    traintimes = []
    testtimes = []
    aucs = []
    name = "({}, {}, {})".format(clf_type.__name__, trunc, tech)

    for train, test in cv:
        # Train
        clf = clf_type(
            n_estimators=10,
            max_features=None,
            min_samples_leaf=1,
            random_state=23,
            bootstrap=False,
            max_depth=None,
            categorical=cat,
        )

        traintimes.append(
            timeit(
                "clf.fit(X[train], y[train])".format(),
                "from __main__ import clf, X, y, train",
                number=1,
            )
        )

        """
        # Check that all leaf nodes are pure
        for est in clf.estimators_:
            leaves = est.tree_.children_left < 0
            print(np.max(est.tree_.impurity[leaves]))
            #assert(np.all(est.tree_.impurity[leaves] == 0))
        """

        # Test
        probs = []
        testtimes.append(
            timeit(
                "probs.append(clf.predict_proba(X[test]))",
                "from __main__ import probs, clf, X, test",
                number=1,
            )
        )

        aucs.append(roc_auc_score(y[test], probs[0][:, 1]))

    traintimes = np.array(traintimes)
    testtimes = np.array(testtimes)
    aucs = np.array(aucs)
    results.append(
        [
            name,
            traintimes.mean(),
            traintimes.std(),
            testtimes.mean(),
            testtimes.std(),
            aucs.mean(),
            aucs.std(),
        ]
    )

    results_df = pd.DataFrame(results)
    results_df.columns = [
        "name",
        "train time mean",
        "train time std",
        "test time mean",
        "test time std",
        "auc mean",
        "auc std",
    ]
    results_df = results_df.set_index("name")
    print(results_df)
