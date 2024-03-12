"""
==========================================
IsolationForest prediction benchmark
==========================================
A test of IsolationForest on classical anomaly detection datasets.

The benchmark is run as follows:
1. The dataset is randomly split into a training set and a test set, both
assumed to contain outliers.
2. Isolation Forest is trained on the training set fixed at 1000 samples.
3. The test samples are scored using the trained model at:
    - 1000, 10000, 100000 samples
    - 10, 100, 1000 features
    - 0.01, 0.1, 0.5 contamination
    - 1, 2, 3, 4 n_jobs

We compare the prediction time at the very end.
"""

from time import time
from pathlib import Path
from collections import defaultdict

import numpy as np
from joblib import parallel_backend
import pandas as pd

from sklearn.ensemble import IsolationForest

print(__doc__)


def get_data(
    n_samples_train, n_samples_test, n_features, contamination=0.1, random_state=0
):
    """Function based on code from: https://scikit-learn.org/stable/
    auto_examples/ensemble/plot_isolation_forest.html#sphx-glr-auto-
    examples-ensemble-plot-isolation-forest-py
    """
    rng = np.random.RandomState(random_state)

    X = 0.3 * rng.randn(n_samples_train, n_features)
    X_train = np.r_[X + 2, X - 2]

    X = 0.3 * rng.randn(n_samples_test, n_features)
    X_test = np.r_[X + 2, X - 2]

    n_outliers = int(np.floor(contamination * n_samples_test))
    X_outliers = rng.uniform(low=-4, high=4, size=(n_outliers, n_features))

    outlier_idx = rng.choice(np.arange(0, n_samples_test), n_outliers, replace=False)
    X_test[outlier_idx, :] = X_outliers

    return X_train, X_test


random_state = 1

results = defaultdict(list)

# Loop over all datasets for fitting and scoring the estimator:
n_samples_train = 1000
for n_samples_test in [1000, 10000, 100000]:
    for n_features in [10, 100, 1000]:
        for contamination in [0.01, 0.1, 0.5]:
            for n_jobs in [1, 2, 3, 4]:
                X_train, X_test = get_data(
                    n_samples_train,
                    n_samples_test,
                    n_features,
                    contamination,
                    random_state,
                )

                print("--- Fitting the IsolationForest estimator...")
                model = IsolationForest(n_jobs=-1, random_state=random_state)
                tstart = time()
                model.fit(X_train)
                fit_time = time() - tstart

                tstart = time()
                scoring = -model.decision_function(
                    X_test
                )  # the lower, the more abnormal
                predict_time = time() - tstart

                # clearcache
                for _ in range(1000):
                    1 + 1

                with parallel_backend("loky", n_jobs=n_jobs):
                    tstart = time()
                    scoring = -model.decision_function(
                        X_test
                    )  # the lower, the more abnormal
                    predict_time = time() - tstart

                results["predict_time"].append(predict_time)
                results["fit_time"].append(fit_time)
                results["n_samples_train"].append(n_samples_train)
                results["n_samples_test"].append(n_samples_test)
                results["n_features"].append(n_features)
                results["contamination"].append(contamination)
                results["n_jobs"].append(n_jobs)

df = pd.DataFrame(results)
df.to_csv("bench_results/isolation_forest_predict.csv", index=False)


def plot(bench_results, pr_name, main_name, image_path):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    results_path = Path(bench_results)
    pr_path = results_path / f"{pr_name}.csv"
    main_path = results_path / f"{main_name}.csv"
    image_path = results_path / image_path

    df_pr = pd.read_csv(pr_path).assign(branch=pr_name)
    df_main = pd.read_csv(main_path).assign(branch=main_name)
    df_all = pd.concat((df_pr, df_main), ignore_index=True)

    gb = df_all.groupby(["n_jobs", "n_samples_test"])
    groups = gb.groups

    n_rows, n_cols = 2, 4
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 8), constrained_layout=True)
    axes_flat = axes.ravel()
    for i, (keys, idx) in enumerate(groups.items()):
        ax = axes_flat[i]
        ax.set_title(" | ".join(keys))
        sns.boxplot(data=df_all.loc[idx], y="predict_time", x="branch", ax=ax)
        if i % n_cols != 0:
            ax.set_ylabel("")

    axes_flat[-1].set_visible(False)

    fig.savefig(image_path)
    print(f"Saved image to {image_path}")


bench_results = Path("~/bench_results_forest")
pr_name = "pr"
main_name = "main"
image_path = "results_image.png"
plot(
    bench_results=bench_results,
    pr_name=pr_name,
    main_name=main_name,
    image_path=image_path,
)
