import numpy as np
import pandas as pd
from datetime import datetime
import time
from sklearn.ensemble import RandomForestClassifier, ObliqueRandomForestClassifier
from sklearn.model_selection import GridSearchCV
import _pickle as cPickle


def get_forest_sizes(X, y, df, random_state):
    N = df.shape[0]
    tree_size = []
    max_n_leaves = []

    df.reset_index(drop=True, inplace=True)

    t_i = time.time()

    # for every parameter combination in GridSearch, refit the clf
    # with the parameters and estimate the size of the forest.
    for i in range(N):
        row = df.iloc[i, :]

        if row["clf"] == "RF":
            clf = RandomForestClassifier(random_state=random_state, **row["params"])
        elif row["clf"] == "OF":
            clf = ObliqueRandomForestClassifier(
                random_state=random_state, **row["params"]
            )
        else:
            print("Cannot identify estimator")

        clf.fit(X, y)
        tree_size.append(cPickle.dumps(clf).__sizeof__())

        max_leaves = max([est.tree_.n_leaves for est in clf.estimators_])
        max_n_leaves.append(max_leaves)

    # create the new dataframe
    new_df = pd.concat([df, pd.DataFrame(tree_size, columns=["clf_size"])], axis=1)
    new_df = pd.concat(
        [new_df, pd.DataFrame(max_n_leaves, columns=["max_leaves"])], axis=1
    )
    print(f"Refitting trees took {int(time.time()-t_i)} seconds")
    return new_df


#########################
random_state = 123456
t0 = datetime.now()
params = {
    "max_features": ["sqrt", "log2", None],
    "n_estimators": [i for i in range(100, 300, 500)],
    "max_depth": [5, 10, 15, 20, None],
}

df_cv = pd.DataFrame()
clfs = {
    "RF": RandomForestClassifier(random_state=random_state),
    "OF": ObliqueRandomForestClassifier(random_state=random_state),
}


def sparse_parity(n_samples, p=20, p_star=3, random_seed=None, **kwarg):
    if random_seed:
        np.random.seed(random_seed)

    X = np.random.uniform(-1, 1, (n_samples, p))
    y = np.zeros(n_samples)

    for i in range(0, n_samples):
        y[i] = sum(X[i, :p_star] > 0) % 2

    return X, y


X, y = sparse_parity(n_samples=10000, random_seed=random_state)
for clf_lab, clf in clfs.items():
    t_i = time.time()

    if clf_lab == "OF":
        of_params = params.copy()
        of_params["max_features"] = ["sqrt", "log2", None, 2 * X.shape[1]]
        param_grid = of_params
    else:
        param_grid = params.copy()

    search = GridSearchCV(
        estimator=clf,
        param_grid=param_grid,
        return_train_score=True,
        random_state=random_state,
        cv=None,  # default 5-fold
        n_jobs=1,
    )
    search.fit(X, y)
    print(f"Grid Searching took {int(time.time()-t_i)} seconds for {clf_lab} algorithm")

    df_tmp = pd.DataFrame(search.cv_results_)
    df_tmp.columns = [i.replace("param_", "") for i in df_tmp.columns]
    df_tmp["clf"] = clf_lab
    # df_tmp.fillna('None', inplace=True)
    df_tmp["mean_test_score"] = df_tmp.apply(
        lambda x: round(x["mean_test_score"], 3), axis=1
    )

    # refit trees to get the forest sizes
    df_tmp = get_forest_sizes(X, y, df_tmp, random_state)
    df_cv = pd.concat([df_cv, df_tmp])

df_cv["clf_size_MB"] = df_cv.apply(lambda x: x["clf_size"] / 10e6, axis=1)
df_cv.to_csv("OFvsRF_grid_search_cv_results_sparse_parity.csv")
