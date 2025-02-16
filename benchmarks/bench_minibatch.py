import os
import time
from collections import defaultdict
from datetime import datetime
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn import metrics
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.datasets import fetch_20newsgroups, fetch_openml
from sklearn.decomposition._truncated_svd import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import StandardScaler


# MARK: -Load and preprocess data
def load_and_preprocess_data(dataset, n_sample=None):
    if dataset == "20newsgroups":
        newsgroups = fetch_20newsgroups(subset="train")
        vectorizer = TfidfVectorizer()
        X = vectorizer.fit_transform(newsgroups.data)
        y = newsgroups.target
        n_components = 100
        # Apply Truncated SVD to the sparse matrix X
        svd = TruncatedSVD(n_components=n_components)
        X = svd.fit_transform(X)
        return X, y
    else:
        try:
            X, y = fetch_openml(
                name=dataset,
                version=1,
                as_frame=False,
                return_X_y=True,
                data_home=None,
                cache=True,
                parser="auto",
            )
        except Exception as e:
            raise Exception(f"Could not load dataset {dataset}. Error: {e}")
    normalize = False
    if dataset == "pendigits":
        normalize = True
    if dataset == "letter":
        normalize = True
    if n_sample is not None:
        shuffle = np.random.permutation(X.shape[0])
        X = X[shuffle]
        y = y[shuffle]
        X = X[: min(X.shape[0], n_sample)]
        y = y[: min(X.shape[0], n_sample)]
    if normalize:
        X = StandardScaler().fit_transform(X)
    return X, y


# MARK: -Evaluation
def evaluate(kms, X, labels, num_iters, n_clusters, batch_size, n_runs=50):

    evaluations = []

    for name, km in kms.items():
        train_times = []
        print(f"Evaluating {name}")
        scores = defaultdict(list)
        for seed in range(n_runs):
            km.random_state = seed
            t0 = time.time()
            km.fit(X)
            # include the time it took to construct the kernel matrix in
            # the training time
            train_times.append(time.time() - t0)
            scores["NMI"].append(
                metrics.normalized_mutual_info_score(labels, km.labels_)
            )
            scores["ARI"].append(metrics.adjusted_rand_score(labels, km.labels_))
        train_times = np.asarray(train_times)

        evaluation = {
            "estimator": name,
            "num_iters": num_iters,
            "n_clusters": n_clusters,
            "batch_size": batch_size,
            "train_time_mean": train_times.mean(),
            "train_time_std": train_times.std(),
        }
        for score_name, score_values in scores.items():
            mean_score, std_score = np.mean(score_values), np.std(score_values)
            evaluation[score_name + "_mean"] = mean_score
            evaluation[score_name + "_std"] = std_score

        evaluations.append(evaluation)

        print(
            f"\n {name}, num_iters: {num_iters}, n_clusters: {n_clusters},\
                 batch size: {batch_size}"
        )
        for score_name, score_values in scores.items():
            mean_score, std_score = np.mean(score_values), np.std(score_values)
            print(f"{score_name}: {mean_score:.3f} ± {std_score:.3f}")
    return evaluations


colors = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # safety orange
    "#ff7f0e",  # safety orange
    "#9467bd",  # muted purple
    "#9467bd",  # muted purple
    "#8c564b",  # chestnut brown
    "#e377c2",  # raspberry yogurt pink
    "#7f7f7f",  # middle gray
    "#bcbd22",  # curry yellow-green
    "#17becf",  # blue-teal
]
hatches = ["", "//", "", "//", ""]


def plot_results(to_plot):
    plt.rcParams.update({"font.size": 24})
    plt.figure(figsize=(10, 6))
    num_res = len(to_plot)  # Number of rows in the grid
    # assume all DFs have the same batch sizes
    batch_sizes = to_plot[0]["batch_size"].unique()
    num_batches = len(batch_sizes)
    print(num_batches, num_res)
    fig, axes = plt.subplots(
        num_batches, num_res, figsize=(7 * num_res, 6 * num_batches)
    )
    for j in range(num_batches):
        for i, df1 in enumerate(to_plot):
            b = batch_sizes[j]
            name = df1["dataset"].iloc[0]
            df = df1[df1["batch_size"] == batch_sizes[j]]

            if num_batches == 1 and num_res == 1:
                ax = axes
            elif num_batches == 1:
                ax = axes[i]
            else:
                ax = axes[j][i]
            ax1, ax2 = plot_results_bars(df, ax, i == 0 and j == 0)
            if i == 0:
                ax1.set_ylabel("Score")
            if i == num_res - 1:
                ax2.set_ylabel("Time (s)")
            ax.set_title(f"{name} (batch size: {b})")

    fig.legend(loc="lower center", bbox_to_anchor=(0.5, 1.04), ncol=5, fontsize=34)

    plt.tight_layout()
    # write to results directory
    plt.savefig("minibatch_results/results.png", bbox_inches="tight")


def plot_results_bars(df, ax1, set_labels=True):
    metric_names = ["ARI", "NMI"]
    time_metric = "train_time"
    sorted(df["estimator"].unique())
    ax2 = ax1.twinx()
    n_metrics = len(metric_names) + 1  # Including train_time
    bar_width = 0.4
    positions = np.arange(n_metrics) * (len(df["estimator"].unique()) * bar_width + 0.5)
    df_comb = df

    for i, metric in enumerate(metric_names + [time_metric]):
        metric_mean = metric + "_mean"
        metric_std = metric + "_std"
        for j, name in enumerate(sorted(df["estimator"].unique())):
            position = positions[i] + j * bar_width - 0.5
            ax = ax1
            if metric == time_metric:
                ax = ax2
            alg_name = name[2:]
            ax.bar(
                position,
                df_comb[df_comb["estimator"] == name][metric_mean].iloc[0],
                bar_width,
                color=colors[j],
                label=(alg_name) if i == 0 and set_labels else "",
                yerr=df_comb[df_comb["estimator"] == name][metric_std].iloc[0],
                capsize=5,
                hatch=hatches[j],
                edgecolor="black",
                linewidth=1,
            )
    ax1.set_xticks(positions + bar_width / 2)
    ax1.set_xticklabels(metric_names + ["runtime"])
    return ax1, ax2


result_files = []


n_runs = 10
n_iters = [100]
batch_size_values = [1024]
tol = 1e-4

to_plot = []

dataset_names = [
    "pendigits",
    "har",
    "mnist_784",
    "letter",
    "20newsgroups",
]
print("Running on datasets:", dataset_names)
for dataset_name in dataset_names:
    n = None
    X, Y = load_and_preprocess_data(dataset_name, n_sample=n)
    if n is None:
        n = X.shape[0]

    num_clusters = np.unique(Y).shape[0]
    print(f"num clusters: {num_clusters}")
    n_clusters = np.unique(Y).shape[0]
    n_clusters_values = [n_clusters]
    evaluations = []
    current_datetime = datetime.now()
    print(f"dataset: {dataset_name}")
    evaluations = []

    for num_iters, n_clusters, batch_size in product(
        n_iters, n_clusters_values, batch_size_values
    ):
        print("#" * 20)
        mbk_newlr = MiniBatchKMeans(
            n_clusters=n_clusters,
            batch_size=batch_size,
            max_iter=num_iters,
            adaptive_lr=True,
            tol=tol,
        )
        mbk_oldlr = MiniBatchKMeans(
            n_clusters=n_clusters,
            batch_size=batch_size,
            max_iter=num_iters,
            adaptive_lr=False,
            tol=tol,
        )
        km_full = KMeans(n_clusters=n_clusters, max_iter=num_iters, tol=tol)
        mbks = {
            "1.new lr MiniBatch": mbk_newlr,
            "2.MiniBatch": mbk_oldlr,
            "3.KMeans": km_full,
        }

        evaluations += evaluate(mbks, X, Y, num_iters, n_clusters, batch_size, n_runs)

    # Convert evaluations to DataFrame
    df = pd.DataFrame(evaluations)
    metric_names = [
        "Homogeneity",
        "Completeness",
        "V-measure",
        "ARI",
        "Silhouette Coefficient",
        "NMI",
    ]
    param_vals = {
        "num_iters": n_iters,
        "n_clusters": n_clusters_values,
        "batch_size": batch_size_values,
        "n_runs": n_runs,
        "n": n,
    }

    if not os.path.exists("minibatch_results"):
        os.makedirs("minibatch_results")

    df["dataset"] = dataset_name
    to_plot.append(df)

plot_results(to_plot)
