#%%
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import make_blobs, make_moons
import matplotlib.pyplot as plt
#from sklearn.cluster import MiniBatchKMeans
from sklearn.datasets import fetch_openml
from collections import defaultdict
from time import time
from tqdm import tqdm
from sklearn import metrics
from datetime import datetime
from sklearn.preprocessing import StandardScaler
#import mini batch k-means from sklrean
from sklearn.cluster import MiniBatchKMeans
import pandas as pd
import time
import os
from itertools import product
from copy import deepcopy
import argparse



# MARK: -Load and preprocess data
def load_and_preprocess_data(dataset, n_sample=None):
    try:
        X,y = fetch_openml(name=dataset, version=1, as_frame=False, return_X_y=True,data_home="data/",cache=True,parser="auto")
    except Exception as e:
        print(e)
        X= None
    
    normalize = False
    knn = 2000
    heat_knn = 250
    if dataset == "pendigits":
        kernel_sigma = 2.23606797749979#3.5#5.6496
        normalize = True
        knn = 1000
        heat_knn = 500

    if dataset == "har":
        #kernel_sigma = 5.5819
        kernel_sigma = 10.414010843544517 #unnormalized
        kernel_sigma = 4.25
        knn = 2000
        heat_knn = 500
        
    if dataset == "moons":
        kernel_sigma = 1.3369
        #generate data set
        X, y = make_moons(n_samples=5000, noise=0.1, random_state=42)

    if dataset == "mnist_784":
        kernel_sigma = 5.8301
        X = X/255
        normalize = False
        knn = 2000
        heat_knn = 200

    if dataset == "EMNIST_Balanced":
        #kernel_sigma = 11.510648727416992 #wang
        kernel_sigma = 3.5 #good for n=100000
        X = X/255
        normalize = False

    if dataset == "letter":
        kernel_sigma = 2.23606797749979#2.3351
        kernel_sigma = 3.0
        normalize = True
        knn = 1000
        heat_knn = 500
        
    
    if n_sample is not None:
        shuffle = np.random.permutation(X.shape[0])
        X = X[shuffle]
        y = y[shuffle]
        X = X[:min(X.shape[0],n_sample)]
        y = y[:min(X.shape[0],n_sample)]


    if normalize:
        X = StandardScaler().fit_transform(X)
    return X,y, kernel_sigma, knn, heat_knn


def initialize_kmeans_plusplus(X, k, rng):
    n_samples, n_features = X.shape
    centroids = np.empty((k, n_features), dtype=X.dtype)

    # Step 1: Randomly choose the first centroid using the provided rng
    centroids[0] = X[rng.integers(n_samples)]

    # Step 2 and 3: Choose the remaining centroids
    for c in range(1, k):
        squared_diff = np.square(X[:, np.newaxis, :] - centroids[np.newaxis, :c, :])
        distances = np.min(squared_diff, axis=1)
        total_distances = np.sum(distances, axis=1)
        probabilities = total_distances / np.sum(total_distances)
        centroids[c] = X[rng.choice(n_samples, p=probabilities)]

    return centroids


# MARK: -Evaluation
def evaluate(kms, X, labels, num_iters, n_clusters, batch_size, n_runs=50):
    
    evaluations = []
                    
    for name, km in kms.items():
        train_times = []
        print(f"Evaluating {name}")
        scores = defaultdict(list)
        for seed in tqdm(range(n_runs)):
            km.random_state=seed
            t0 = time.time()
            km.fit(X)
            
            # include the time it took to construct the kernel matrix in the training time
            train_times.append(time.time() - t0)
            scores["NMI"].append(metrics.normalized_mutual_info_score(labels, km.labels_))
            #scores["Homogeneity"].append(metrics.homogeneity_score(labels, km.labels_))
            #scores["Completeness"].append(metrics.completeness_score(labels, km.labels_))
            #scores["V-measure"].append(metrics.v_measure_score(labels, km.labels_))
            scores["ARI"].append(metrics.adjusted_rand_score(labels, km.labels_))
            #scores["Silhouette Coefficient"].append(metrics.silhouette_score(X, km.labels_, sample_size=2000))
    
        train_times = np.asarray(train_times)

        evaluation = {
            "estimator": name,
            "num_iters": num_iters,
            "n_clusters": n_clusters,
            "batch_size": batch_size,
            "tau": 0, #temp
            "train_time_mean": train_times.mean(),
            "train_time_std": train_times.std()
        }

        for score_name, score_values in scores.items():
            mean_score, std_score = np.mean(score_values), np.std(score_values)
            evaluation[score_name + "_mean"] = mean_score
            evaluation[score_name + "_std"] = std_score
    
        evaluations.append(evaluation)

        print(f"\n {name}, num_iters: {num_iters}, n_clusters: {n_clusters}, batch size: {batch_size}")
        for score_name, score_values in scores.items():
            mean_score, std_score = np.mean(score_values), np.std(score_values)
            print(f"{score_name}: {mean_score:.3f} ± {std_score:.3f}")
        

    return evaluations




colors = [
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#ff7f0e',  # safety orange
    '#9467bd',  # muted purple
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
    # Add more as needed
]

hatches = [
    '','//','','//',''  
]

def plot_results(to_plot):
    plt.rcParams.update({'font.size': 24}) 
    plt.figure(figsize=(10, 6))

    num_res = len(to_plot)  # Number of rows in the grid
    #assume all DFs have the same batch sizes
    batch_sizes = to_plot[0]['batch_size'].unique()
    num_batches = len(batch_sizes)
    print(num_batches,num_res)
    fig, axes = plt.subplots(num_batches, num_res, figsize=(7*num_res, 6*num_batches))
    #fig, axes = plt.subplots(num_batches, num_res, figsize=(7*num_res, 6))
    
    for j in range(num_batches):
        for i, df1 in enumerate(to_plot):
            b = batch_sizes[j]
            name = df1['dataset'].iloc[0]
            df = df1[df1['batch_size'] == batch_sizes[j]]

            if num_batches == 1 and num_res == 1:
                ax = axes
            elif num_batches == 1:
                ax = axes[i]
            else:
                ax = axes[j][i]
            ax1,ax2 = plot_results_bars(df, ax,i==0 and j==0)
            if i == 0:
                ax1.set_ylabel('Score')
            if i == num_res - 1:
                ax2.set_ylabel('Time (s) log scale')
            ax.set_title(f"{name} (batch size: {b})")
    
    fig.legend(loc='lower center', bbox_to_anchor=(0.5, 1.04),ncol=5, fontsize=34)
       
    plt.tight_layout()
    # write to results directory
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y-%m-%d-%H-%M-%S')
    plt.savefig(f"results/{formatted_datetime}_results.png", bbox_inches='tight')



def plot_results_bars(df, ax1, set_labels=True):
    metric_names = ["ARI", "NMI"]
    time_metric = "train_time"
    sorted(df['estimator'].unique())
    ax2 = ax1.twinx()  # Create a second y-axis to plot the train_time
    n_metrics = len(metric_names) + 1  # Including train_time
    
    bar_width = 0.4
    positions = np.arange(n_metrics) * (len(df['estimator'].unique()) * bar_width + 0.5)

    df_comb = df
    #ax2.set_ylabel('Time (s) log scale')
    #set ax2 to log scale
    ax2.set_yscale('log')
    #ax1.set_ylabel('Score')

    
    for i, metric in enumerate(metric_names + [time_metric]):
        metric_mean = metric + "_mean"
        metric_std = metric + "_std"
        for j, name in enumerate(sorted(df['estimator'].unique())):
            position = positions[i] + j * bar_width-0.5
            ax = ax1
            if metric == time_metric:
                ax = ax2
            alg_name = name[2:]
            ax.bar(position, df_comb[df_comb['estimator'] == name][metric_mean].iloc[0], bar_width,
                    color=colors[j], label=(alg_name) if i == 0 and set_labels else "", yerr=df_comb[df_comb['estimator'] == name][metric_std].iloc[0],
                    capsize=5, hatch=hatches[j], edgecolor='black', linewidth=1)
            '''
            ax.bar(position,df["kernel_construction_time"], bottom=df_comb[df_comb['estimator'] == name][metric_mean].iloc[0], width=bar_width,
                    color="black", label=(alg_name) if i == 0 and set_labels else "",
                    capsize=5, hatch=hatches[j], edgecolor='black', linewidth=1)
            '''
            
    ax1.set_xticks(positions + bar_width / 2)
    ax1.set_xticklabels(metric_names + ["runtime"])
    return ax1,ax2

    
#%%

result_files = []

# make parameters global

n_runs = 10
n_iters = [5]
# Define parameter ranges
batch_size_values = [1024]
#taus = [50, 100, 200, 300]

to_plot = []

dataset_names = [
        "pendigits",
        "har",
        "mnist_784",
        "letter"
    ]
print("Running on datasets:", dataset_names)
for dataset_name in dataset_names:
    n = None
    X, Y, kernel_sigma, knn, heat_knn = load_and_preprocess_data(dataset_name, n_sample = n)
    if n is None:
        n = X.shape[0]
    
    print(dataset_name)
    print("dataset size", X.shape, "kernel sigma", kernel_sigma)
    num_clusters = np.unique(Y).shape[0]
    print(f"num clusters: {num_clusters}")

    n_clusters = np.unique(Y).shape[0]
    n_clusters_values = [n_clusters]
    
    evaluations = []
    current_datetime = datetime.now()

    # Format the date and time
    formatted_datetime = current_datetime.strftime('%Y_%m_%d__%H_%M_%S')
    print(formatted_datetime)
    print(f"dataset: {dataset_name}")
    evaluations = []
    
    for num_iters, n_clusters, batch_size in product(n_iters, n_clusters_values, batch_size_values):
        print("#"*20)
        
        #mbk_newlr = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters,max_no_improvement=None, reassignment_ratio = 0, new_lr=True)
        #mbk_oldlr = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters,max_no_improvement=None, reassignment_ratio = 0, new_lr=False)
        mbk_newlr_default = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters, new_lr=True)
        mbk_oldlr_default = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters, new_lr=False)

        mbks = {
                #"1.new lr MiniBatch": mbk_newlr,
                #"2.MiniBatch": mbk_oldlr,
                "3.default new lr MiniBatch": mbk_newlr_default,
                "4.default MiniBatch": mbk_oldlr_default,
                }
        
        evaluations += evaluate(mbks,X, Y, num_iters, n_clusters, batch_size, n_runs)
        #evaluations.append(temp_evals)
        
    

    # Convert evaluations to DataFrame
    df = pd.DataFrame(evaluations)
    metric_names = ["Homogeneity", "Completeness", "V-measure", "ARI", "Silhouette Coefficient", "NMI"]
    param_vals = {'num_iters': n_iters, 'n_clusters': n_clusters_values, 'batch_size': batch_size_values, 'n_runs': n_runs, 'n': n}
    

    if not os.path.exists("results"):
        os.makedirs("results")

    
    df['dataset'] = dataset_name
    to_plot.append(df)

plot_results(to_plot)


# %%
