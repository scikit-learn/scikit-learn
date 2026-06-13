import argparse
import warnings
from time import perf_counter

import numpy as np

from sklearn._loss import HalfBinomialLoss, HalfMultinomialLoss
from sklearn.datasets import fetch_20newsgroups_vectorized, make_classification
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.model_selection import train_test_split


def obj(clf, X, y):
    """Objective function for Elastic-Net Logistic Regression."""
    n_classes = len(clf.classes_)
    fit_intercept = clf.fit_intercept
    lml = LinearModelLoss(
        base_loss=HalfBinomialLoss()
        if n_classes < 3
        else HalfMultinomialLoss(n_classes=n_classes),
        fit_intercept=fit_intercept,
    )
    reg = 1 / (clf.C * X.shape[0])
    l1_reg = clf.l1_ratio * reg
    l2_reg = (1 - clf.l1_ratio) * reg
    if fit_intercept:
        if n_classes <= 2:
            coef = np.r_[clf.coef_[0], clf.intercept_]
        else:
            coef = np.c_[clf.coef_, clf.intercept_]
    else:
        coef = clf.coef_ if n_classes >= 3 else clf.coef_[0]
    return lml.loss(
        coef=coef,
        X=X,
        y=y.astype(X.dtype),
        l1_reg_strength=l1_reg,
        l2_reg_strength=l2_reg,
    )


def fit_solvers(solvers, exclude_solvers, X, y, **params):
    for solver in solvers:
        if solver in exclude_solvers:
            print(f"{solver:>15} was excluded")
            continue

        with warnings.catch_warnings(record=True) as record:
            tic = perf_counter()
            clf = LogisticRegression(solver=solver, **params).fit(X, y)
            toc = perf_counter()
        if len(record) > 0 and any(
            issubclass(rec.category, ConvergenceWarning) for rec in record
        ):
            converged = False
        else:
            converged = True

        objective = obj(clf, X, y)
        msg = (
            f"{solver:>15} iterations={clf.n_iter_[0]:6_d} in {toc - tic:6.3f} seconds"
            f" at {objective=: 13.12f}"
        )
        if clf.l1_ratio > 0:
            n_zeros = np.sum(clf.coef_ == 0)
            msg += f" {n_zeros} of {clf.coef_.size} coefficients are zero."
        msg += f" {converged=}"
        print(msg)


def fit_one_config(
    n_classes=2,
    n_samples=100,
    n_features=5,
    C=1,
    l1_ratio=0,
    fit_intercept=False,
    tol=1e-5,
    exclude_solvers=[],
):
    X, y = make_classification(
        n_samples=n_samples,
        n_classes=n_classes,
        n_features=n_features,
        n_informative=max(1, n_features // 2),
        n_redundant=1,
        random_state=42,
    )
    print(
        f"\n{n_samples=:_d} {n_features=:_d} {n_classes=:d} {l1_ratio=} "
        f"{tol=} {fit_intercept=}"
    )
    params = dict(
        C=C,
        fit_intercept=fit_intercept,
        l1_ratio=l1_ratio,
        tol=tol,
        max_iter=10_000,
        random_state=444,
    )

    # Make X hot in memory
    X.max()

    solvers = ["saga", "newton-cd-gram", "newton-cd"]
    if n_classes == 2:
        solvers += ["liblinear"]
    if l1_ratio == 0:
        solvers += ["newton-cg", "newton-cholesky", "lbfgs"]

    fit_solvers(solvers, exclude_solvers, X, y, **params)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--penalty",
        choices=["L2", "L1"],
        default="L1",
        type=str,
        help="Penalty type. Solvers are selected automatically.",
    )
    parser.add_argument(
        "--data",
        choices=["make_classification", "20newsgroup"],
        default="make_classification",
        type=str,
        help=(
            "The data to use. Solvers are selected automatically.\n"
            "Option 'make_classification' creates dense data with different number "
            "of classes.\n"
            "Option '20newsgroup' uses a down-sampled variant of the vectorized "
            "20-news-group data, it is sparse and has 20 classes."
        ),
    )

    args = parser.parse_args()

    print("Start simple benchmarks for logistic regression.")

    if args.data == "make_classification":
        print("dataset from make_classification, tol=1e-5")
        if args.penalty == "L2":
            print("\nL2 logistic regression")
            params = dict(n_samples=100_000, n_features=10, fit_intercept=False)
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params, exclude_solvers=["newton-cd"])

            print()
            params = dict(n_samples=100_000, n_features=10, fit_intercept=True)
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params, exclude_solvers=["newton-cd"])

            print()
            params = dict(n_samples=100, n_features=1000, fit_intercept=False)
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params, exclude_solvers=["newton-cd-gram"])

            print()
            params = dict(n_samples=100, n_features=1000, fit_intercept=True)
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params, exclude_solvers=["newton-cd-gram"])

        else:  # args.penalty == "L1":
            print("\nL1 logistic regression")
            params = dict(
                n_samples=10_000, n_features=10, l1_ratio=1, fit_intercept=False
            )
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params)

            print()
            params = dict(
                n_samples=10_000, n_features=10, l1_ratio=1, fit_intercept=True
            )
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params)

            print()
            params = dict(
                n_samples=100, n_features=1000, l1_ratio=1, fit_intercept=False
            )
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params)

            print()
            params = dict(
                n_samples=100, n_features=1000, l1_ratio=1, fit_intercept=True
            )
            fit_one_config(n_classes=2, **params)
            fit_one_config(n_classes=3, **params)
            fit_one_config(n_classes=10, **params)

    else:
        warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn")
        # Turn down for faster run time
        n_samples = 5000
        X, y = fetch_20newsgroups_vectorized(subset="all", return_X_y=True)
        X = X[:n_samples]
        y = y[:n_samples]

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, random_state=42, stratify=y, test_size=0.1
        )
        train_samples, n_features = X_train.shape
        n_classes = np.unique(y).shape[0]

        print(
            f"dataset 20newsgroup, {train_samples=:_}, {n_features=:_}, "
            f"{n_classes=:_}, tol=1e-4"
        )

        exclude_solvers = ["newton-cd-gram"]  # too many features!!!
        if args.penalty == "L2":
            exclude_solvers += ["newton-cd"]  # takes long time
            exclude_solvers += ["newton-cholesky"]  # too many features!!!
            solvers = [
                "saga",
                "lbfgs",
                "newton-cd",
                "newton-cd-gram",
                "newton-cg",
                "newton-cholesky",
                "lbfgs",
            ]
            l1_ratio = 0
        else:
            solvers = ["saga", "newton-cd"]
            l1_ratio = 1
        params = dict(
            C=1,
            fit_intercept=True,
            l1_ratio=l1_ratio,
            tol=1e-4,
            max_iter=1_000,
            random_state=444,
        )

        fit_solvers(solvers, exclude_solvers, X, y, **params)
