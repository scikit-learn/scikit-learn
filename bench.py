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
            rng = np.random.default_rng(seed) #set random seed
            
        
            km.random_state=seed
            t0 = time.time()
            init_centroids = None
            km.fit(X, init_centroids)
            
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
    '#9467bd',  # muted purple
    '#17becf',   # blue-teal
    # Add more as needed
]

hatches = [ 
    "x","","//","//","","//"
]

hatches_map = {
    "Kernel": "x",
    "MiniBatch Kernel": "",
    "$\\beta$-MiniBatch": "//",
    "Truncated $\\beta$-MiniBatch Kernel": "//",
    "MiniBatch": "",
    "$\\beta$-MiniBatch Kernel": "//",
}

colors_map = {
    "Kernel": colors[0],
    "Truncated $\\beta$-MiniBatch Kernel": colors[3],
    "$\\beta$-MiniBatch Kernel": colors[1],
    "MiniBatch Kernel": colors[1],
    "$\\beta$-MiniBatch": colors[2],
    "MiniBatch": colors[2],
}

def plot_results(to_plot):
    plt.rcParams.update({'font.size': 24}) 
    
    # Extract unique batch sizes and taus for the grid
    batch_sizes = sorted(to_plot[0]['batch_size'].unique())
    num_batches = len(batch_sizes)
    

    # Now plot the figure that will go in the body of the paper. This will
    # be one plot for batch size 1024 and tau 200 with one column per dataset.
    # It will not include the _knn or _heat datasets.

    # Create a figure with a single row and as many columns as datasets
    datasets = [df['dataset'].iloc[0] for df in to_plot]
    datasets = [d for d in datasets if not d.endswith("_knn") and not d.endswith("_heat")]
    num_datasets = len(datasets)

    #greg: changed sharey to False
    fig, axes = plt.subplots(1, num_datasets, figsize=(7*num_datasets, 6),sharey=False)

    # Adjust the title size and space below the title
    fig.suptitle("Results for Batch size 1024 and $\\tau$ 200", fontsize=40, y=0.99)

    # Iterate over datasets that are not _knn or _heat
    #greg: changed tau to 200
    for i, data_set_group in enumerate([d for d in to_plot if (not d['dataset'].iloc[0].endswith("_knn")) and (not d['dataset'].iloc[0].endswith("_heat"))]):
        dataset_name = data_set_group['dataset'].iloc[0]

        df_tau_bs = data_set_group[(data_set_group['tau'] == 200) & (data_set_group['batch_size'] == 1024)]


        ax1, ax2 = plot_results_bars(df_tau_bs, axes[i], i == 0)

        # Label y-axes
        if i == 0:
            ax1.set_ylabel('Score')
        if i == num_datasets - 1:  # Apply time label to all rows
            ax2.set_ylabel('Time (s) log scale')

        ax1.set_title(f"{dataset_name}")

    # Position the legend outside of the plotting area, reducing space between the plots and legend
    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=5, fontsize=28, frameon=False)

    # Adjust layout: Reduce space between title and plots, and between plots and the legend
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust rect to control the space left for title/legend

    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y_%m_%d__%H_%M_%S')

    # Save the figure for each dataset
    plt.savefig(f"results/paper_body_results.png", bbox_inches='tight')
    plt.close()


    #temp fix
    num_taus=1
    taus = [200]
    print(f"num batches: {num_batches}, num taus: {num_taus}")
    
    # Iterate over datasets
    for dataset_group in to_plot:
        dataset_name = dataset_group['dataset'].iloc[0]
        
        # Create a figure for each dataset
        fig, axes = plt.subplots(num_taus, num_batches, figsize=(7*num_batches, 6*num_taus),sharey=True)

        # Adjust the title size and space below the title
        fig.suptitle(f"Dataset: {dataset_name}", fontsize=40, y=0.99)


        maximum_time = dataset_group['train_time_mean'].max()*1.08

        # Iterate over taus and batch sizes
        for i, tau in enumerate(taus):
            for j, batch_size in enumerate(batch_sizes):
                df_tau_bs = dataset_group[(dataset_group['tau'] == tau) & (dataset_group['batch_size'] == batch_size)]

                if num_taus == 1 and num_batches == 1:
                    ax = axes
                elif num_taus == 1:
                    ax = axes[j]
                elif num_batches == 1:
                    ax = axes[i]
                else:
                    ax = axes[i, j]

                ax1, ax2 = plot_results_bars(df_tau_bs, ax, i == 0 and j == 0)

                # Add batch size label on the top of each column
                if i == 0:
                    ax.set_title(f"Batch size: {batch_size}")
                
                # Add tau label to first column of each row with increased font size and vertically aligned
                if j == 0:
                    # Increase font size of Tau labels
                    ax1.text(-0.33, 0.5, f"$\\tau$: {tau}", va='center', ha='center', 
                             rotation=0, transform=ax1.transAxes, fontsize=26)

                # Label y-axes
                if j == 0:
                    ax1.set_ylabel('Score')
                if j == num_batches - 1:  # Apply time label to all rows
                    ax2.set_ylabel('Time (s) log scale')

                ax2.set_ylim(0.1, maximum_time)

        # Position the legend outside of the plotting area, reducing space between the plots and legend
        fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=5, fontsize=28, frameon=False)

        # Adjust layout: Reduce space between title and plots, and between plots and the legend
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust rect to control the space left for title/legend

        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime('%Y_%m_%d__%H_%M_%S')

        # Save the figure for each dataset
        plt.savefig(f"results/{formatted_datetime}_{dataset_name}_results.png", bbox_inches='tight')
        plt.close()








def plot_results_bars(df, ax1, set_labels=True):
    metric_names = ["ARI", "NMI"]
    time_metric = "train_time"
    
    # Sort estimators explicitly
    estimators_ordered = df['estimator'][::-1]
    estimators_ordered = [
        "1.Kernel",
        "9.Truncated $\\beta$-MiniBatch Kernel",
        "2.$\\beta$-MiniBatch Kernel",
        "3.MiniBatch Kernel",
        "4.$\\beta$-MiniBatch",
        "5.MiniBatch",
    ]
    
    ax2 = ax1.twinx()  # Create a second y-axis to plot the train_time
    n_metrics = len(metric_names) + 1  # Including train_time
    
    bar_width = 0.4
    positions = np.arange(n_metrics) * (len(estimators_ordered) * bar_width + 0.5)

    df_comb = df
    ax2.set_yscale('log')

    # Iterate over metrics and estimators in the same order
    for i, metric in enumerate(metric_names + [time_metric]):
        metric_mean = metric + "_mean"
        metric_std = metric + "_std"
        
        for j, name in enumerate(estimators_ordered):
            position = positions[i] + j * bar_width - 0.5
            ax = ax1
            if metric == time_metric:
                ax = ax2
            
            alg_name = name[2:]
            kernel_construction_time = 0

            # Add kernel construction time only for relevant methods
            if alg_name in ["Truncated $\\beta$-MiniBatch Kernel",
                            "$\\beta$-MiniBatch Kernel", 
                            "Kernel",
                            "MiniBatch Kernel"] and metric == "train_time":
                kernel_construction_time = df["kernel_construction_time"].iloc[0]

            # Plot the kernel construction time (black part of the bars)
            if kernel_construction_time > 0:
                ax.bar(position, kernel_construction_time, bar_width, color="black", edgecolor='black', linewidth=1)

            # Plot the main bar with color and hatch
            hatch = hatches_map[alg_name]
            color = colors_map[alg_name]
            bar_height = df_comb[df_comb['estimator'] == name][metric_mean].iloc[0]
            ax.bar(position, bar_height, bar_width, color=color, 
                   label=(alg_name) if i == 0 and set_labels else None, 
                   yerr=df_comb[df_comb['estimator'] == name][metric_std].iloc[0],
                   capsize=5, hatch=hatch, bottom=kernel_construction_time, edgecolor='black', linewidth=1)

            # # Add text annotation above bars
            # total_height = bar_height + kernel_construction_time
            # ax.text(position, total_height + df_comb[df_comb['estimator'] == name][metric_std].iloc[0] + 0.01,  
            #         f'{total_height:.2f}', ha='center', va='bottom', fontsize=10)
            

    # Set x-ticks and labels
    ax1.set_xticks(positions + bar_width / 2)
    ax1.set_xticklabels(metric_names + ["runtime"])
    
    return ax1, ax2

    
#%%
# MARK: -Main
# read arguments using argparse

#parser = argparse.ArgumentParser()
#parser.add_argument("mode", type=str, help="Mode of operation: 'run' or 'plot'", choices=["run","plot"])

# add an optional list of result files to plot
#parser.add_argument("-f", "--files", nargs="+", help="List of result files to plot", default=[])

#args = parser.parse_args()

#mode = args.mode
#result_files = args.files
print("Mode:!")

# mode ="run"
mode = "run"
result_files = []

# make parameters global

n_runs = 10
n_iters = [200]
# Define parameter ranges
batch_size_values = [256,512,1024,2048]
#taus = [50, 100, 200, 300]


if mode == "run":

    skip_full = False
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
            
            mbk_newlr = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters)#, new_lr=True)
            mbk_oldlr = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch_size, max_iter=num_iters)#, new_lr=False)

            mbks = {
                    "4.$\\beta$-MiniBatch": mbk_newlr,
                    "5.MiniBatch": mbk_oldlr,
                    #"2.$\\beta$-MiniBatch Kernel": mbkk_newlr,
                    #"3.MiniBatch Kernel": mbkk_oldlr,
                    }
            
            temp_evals = evaluate(mbks,X, Y, num_iters, n_clusters, batch_size, n_runs)
            


        

        # Convert evaluations to DataFrame
        df = pd.DataFrame(evaluations)
        metric_names = ["Homogeneity", "Completeness", "V-measure", "ARI", "Silhouette Coefficient", "NMI"]
        param_vals = {'num_iters': n_iters, 'n_clusters': n_clusters_values, 'batch_size': batch_size_values, 'n_runs': n_runs, 'n': n}
        

        if not os.path.exists("results"):
            os.makedirs("results")

        #  clear the directory
        #for filename in os.listdir("results"):
        #    if filename.endswith(".csv"):
        #        os.remove(f"results/{filename}")

        #add dataset name to df
        df['dataset'] = dataset_name
        #add kernel construction time to df
        
        if result_files != []:
            path = os.path.join("results",f"{result_files[0]}")
        else:
            path = os.path.join("results",f"{dataset_name}_{formatted_datetime}_results.csv")
        df.to_csv(path, index=False)

elif mode == "plot":
    # Set the directory for results

    filepaths = []
    if result_files != []:
        filepaths = result_files
    else:
        directory = "results/"
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                filepaths.append(directory + filename)


    to_plot = []
    # Assume num_epochs_values contains a single value
    num_iters = n_iters[0]  # Extract the single value of num_epochs

    # Iterate over CSV files in the directory
    for filename in filepaths:
        df = pd.read_csv(filename)
        to_plot.append(df)

    # Prepare a list of dataframes filtered for the single num_epochs value
    to_plot_filtered = []

    for df in to_plot:
        # Filter dataframe by num_iters only
        df_filtered = df[df['num_iters'] == num_iters]

        # If there is data for this num_epochs value, add it to the list
        if not df_filtered.empty:
            to_plot_filtered.append(df_filtered)

    # Only plot if there is data to plot
    if to_plot_filtered:
        plot_results(to_plot_filtered)


# %%
# %%
