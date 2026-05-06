import numpy as np

from sklearn.metrics import balanced_accuracy_score, r2_score


def neg_mean_inertia(X, labels, centers):
    return -(np.asarray(X - centers[labels]) ** 2).sum(axis=1).mean()


def make_gen_classif_scorers(caller):
    caller.train_scorer = balanced_accuracy_score
    caller.test_scorer = balanced_accuracy_score


def make_gen_reg_scorers(caller):
    caller.test_scorer = r2_score
    caller.train_scorer = r2_score


def neg_mean_data_error(X, U, V):
    return -np.sqrt(((X - U.dot(V)) ** 2).mean())


def make_dict_learning_scorers(caller):
    caller.train_scorer = lambda _, __: (
        neg_mean_data_error(
            caller.X, caller.estimator.transform(caller.X), caller.estimator.components_
        )
    )
    caller.test_scorer = lambda _, __: (
        neg_mean_data_error(
            caller.X_val,
            caller.estimator.transform(caller.X_val),
            caller.estimator.components_,
        )
    )


def explained_variance_ratio(Xt, X):
    return np.var(Xt, axis=0).sum() / np.var(X, axis=0).sum()


def make_pca_scorers(caller):
    caller.train_scorer = lambda _, __: caller.estimator.explained_variance_ratio_.sum()
    caller.test_scorer = lambda _, __: (
        explained_variance_ratio(caller.estimator.transform(caller.X_val), caller.X_val)
    )
