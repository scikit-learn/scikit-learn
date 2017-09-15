# som.py
# version 1.0
# (c) 2017-2017 Lutz Hamel, Li Yuan, University of Rhode Island
#
# This file constitues a set of routines which are useful in constructing
# and evaluating self-organizing maps (SOMs).
# The main utilities available in this file are:
# build --------- constructs a map
# convergence --- reports the map convergence index
# embed --------- reports the embedding of the map in terms of modeling the
#                 underlying data distribution (100% if all feature
#                 distributions are modeled correctly, 0% if none are)
# embed_ks ------ reports the embedding of the map using the
#                 Kolmogorov-Smirnov Test
# topo ---------- reports the estimated topographic accuracy
# significance -- graphically reports the significance of each feature with
#                 respect to the self-organizing map model
# starburst ----- displays the starburst representation of the SOM model,
#                 the centers of starbursts are the centers of clusters
# projection ---- print a table with the associations of labels with map
#                 elements
# neuron -------- returns the contents of a neuron at (x,y) on the map as a
#                 vector
# marginal ------ displays a density plot of a training dataframe dimension
#                 overlayed with the neuron density for that same dimension
#                 or index.
#
#
# ## License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# ##

import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns                   # Plot density by marginal
import statsmodels.stats.api as sms     # t-test
import statistics as stat               # F-test
from random import randint
from sklearn.metrics.pairwise import euclidean_distances
from scipy import stats                 # KS Test
from scipy.stats import f               # F-test
from itertools import combinations


def build(data, labels, xdim=10, ydim=5, alpha=.3, train=1000,
          algorithm="som"):
    """ build -- construct a SOM, returns an object of class 'map'

     parameters:
     - data   - a dataframe where each row contains an unlabeled training
                instance
     - labels - a vector or dataframe with one label for each observation
                in data
     - xdim, ydim - the dimensions of the map
     - alpha  - the learning rate, should be a positive non-zero real number
     - train  - number of training iterations
     - algorithm - selection switch

     retuns:
     - an object of type 'map' -- see below

     Example:
     - m = som.build(data,labels,algorithm=algorithm)

     NOTE: default algorithm: "som" also available: "som_f"

    """
    algorithms = ["som"]

    # check if the dims are reasonable
    if (xdim < 3 or ydim < 3):
        sys.exit("build: map is too small.")

    try:
        index_algorithm = algorithms.index(algorithm)
    except ValueError:
        sys.exit("build only supports 'som','som_f'")

    if index_algorithm == 0:        # som by python
        neurons = vsom_p(data,
                         xdim=xdim,
                         ydim=ydim,
                         alpha=alpha,
                         train=train)

    else:
        sys.exit("build only supports 'som'")

    map = {
            'data': data,
            'labels': labels,
            'xdim': xdim,
            'ydim': ydim,
            'alpha': alpha,
            'train': train,
            'algorithm': algorithm,
            'neurons': neurons}

    visual = []

    for i in range(data.shape[0]):
        b = best_match(map, data.iloc[[i]])
        visual.extend([b])

    map['visual'] = visual

    return(map)


def convergence(map, conf_int=.95, k=50, verb=False, ks=False):
    """ map.convergence - the convergence index of a map

     parameters:
     - map is an object if type 'map'
     - conf_int is the confidence interval of the quality assessment
                (default 95%)
     - k is the number of samples used for the estimated topographic accuracy
         computation
     - verb if true reports the two convergence components separately,
            otherwise it will report the linear combination of the two
     - ks is a switch, true for ks-test, false for standard var and means test

     - return value is the convergence index

     Example:
     - som.convergence(m)
    """
    if ks:
        embed = embed_ks(map, conf_int, verb=False)
    else:
        embed = embed_vm(map, conf_int, verb=False)

    topo_ = topo(map, k, conf_int, verb=False, interval=False)

    if verb:
        return {"embed": embed, "topo": topo_}
    else:
        return (0.5*embed + 0.5*topo_)


def embed(map, conf_int=.95, verb=False, ks=False):
    """ map.embed - evaluate the embedding of a map using the F-test and
                    a Bayesian estimate of the variance in the training data.
    parameters:
    - map is an object if type 'map'
    - conf_int is the confidence interval of the convergence test (default 95%)
    - verb is switch that governs the return value false: single convergence
      value is returned, true: a vector of individual feature congences is
      returned.

    - return value is the cembedding of the map (variance captured by the map
      so far)

    Example:
    - som.embed(m)
    - som.embed(m,verb=True)

    Hint: the embedding index is the variance of the training data captured by
          the map; maps with convergence of less than 90% are typically not
          trustworthy. Of course, the precise cut-off depends on the noise
          level in your training data.
    """

    if ks:
        return embed_ks(map, conf_int, verb)
    else:
        return embed_vm(map, conf_int, verb)


def topo(map, k=50, conf_int=.95, verb=False, interval=True):
    """ topo - measure the topographic accuracy of the map using sampling

    parameters:
    - map is an object if type 'map'
    - k is the number of samples used for the accuracy computation
    - conf.int is the confidence interval of the accuracy test (default 95%)
    - verb is switch that governs the return value, false: single accuracy
      value is returned, true: a vector of individual feature accuracies is
      returned.
    - interval is a switch that controls whether the confidence interval is
      computed.

    return value is the estimated topographic accuracy

    Example:
    - som.topo(m)
    """

    # data.df is a matrix that contains the training data
    data_df = map['data']

    if (k > data_df.shape[0]):
        sys.exit("topo: sample larger than training data.")

    data_sample_ix = [randint(1, data_df.shape[0]) for _ in range(k)]

    # compute the sum topographic accuracy - the accuracy of a single sample
    # is 1 if the best matching unit is a neighbor otherwise it is 0
    acc_v = []
    for i in range(k):
        acc_v.append(accuracy(map,
                              data_df.iloc[data_sample_ix[i]-1],
                              data_sample_ix[i]))

    # compute the confidence interval values using the bootstrap
    if interval:
        bval = bootstrap(map, conf_int, data_df, k, acc_v)

    # the sum topographic accuracy is scaled by the number of samples -
    # estimated
    # topographic accuracy
    if verb:
        return acc_v
    else:
        val = np.sum(acc_v)/k
        if interval:
            return {'val': val, 'lo': bval['lo'], 'hi': bval['hi']}
        else:
            return val


def starburst(map, explicit=False, smoothing=2, merge_clusters=True,
              merge_range=.25):
    """ starburst - compute and display the starburst representation of clusters

    parameters:
    - map is an object if type 'map'
    - explicit controls the shape of the connected components
    - smoothing controls the smoothing level of the umat (NULL,0,>0)
    - merge_clusters is a switch that controls if the starburst clusters are
                     merged together
    - merge_range - a range that is used as a percentage of a certain distance
                    in the code to determine whether components are closer to
                    their centroids or centroids closer to each other.

    Example:
    - som.starburst(m)

    """

    umat = compute_umat(map, smoothing=smoothing)
    plot_heat(map,
              umat,
              explicit=explicit,
              comp=True,
              merge=merge_clusters,
              merge_range=merge_range)


def projection(map):
    """ projection - print the association of labels with map elements

    parameters:
    - map is an object if type 'map'

    return values:
    - a dataframe containing the projection onto the map for each observation

    Example:
    - som.projection(m)

    """

    labels_v = map['labels']
    x_v = []
    y_v = []

    for i in range(len(labels_v)):

        ix = map['visual'][i]
        coord = coordinate(map, ix)
        x_v.append(coord[0])
        y_v.append(coord[1])

    return pd.DataFrame({'labels': labels_v, 'x': x_v, 'y': y_v})


def significance(map, graphics=True, feature_labels=False):
    """ significance - compute the relative significance of each feature and
                       plot it

    parameters:
    - map is an object if type 'map'
    - graphics is a switch that controls whether a plot is generated or not
    - feature.labels is a switch to allow the plotting of feature names vs
      feature indices

    return value:
    - a vector containing the significance for each feature

    Example:
    - som.significance(m)

    """

    data_df = map['data']
    nfeatures = data_df.shape[1]

    # Compute the variance of each feature on the map
    var_v = [randint(1, 1) for _ in range(nfeatures)]

    for i in range(nfeatures):
        var_v[i] = np.var(np.array(data_df)[:, i])

    # we use the variance of a feature as likelihood of
    # being an important feature, compute the Bayesian
    # probability of significance using uniform priors

    var_sum = np.sum(var_v)
    prob_v = var_v/var_sum

    # plot the significance
    if graphics:
        y = max(prob_v)

        plt.axis([0, nfeatures+1, 0, y])

        x = np.arange(1, nfeatures+1)
        tag = list(data_df)

        plt.xticks(x, tag)
        plt.yticks = np.linspace(0, y, 5)

        i = 1
        for xc in prob_v:
            plt.axvline(x=i, ymin=0, ymax=xc)
            i += 1

        plt.xlabel('Features')
        plt.ylabel('Significance')
        plt.show()
    else:
        return prob_v


def neuron(map, x, y):
    """ neuron - returns the contents of a neuron at (x, y) on the map as
                 a vector

    parameters:
      map - the neuron map
      x - map x-coordinate of neuron
      y - map y-coordinate of neuron

    return value:
      a vector representing the neuron

    Example:
    - som.neuron(m,3,3)

    """

    ix = rowix(map, x, y)
    return map['neurons'][ix]


def marginal(map, marginal):
    """ marginal - plot that shows the marginal probability distribution
                   of the neurons and data

    parameters:
    - map is an object of type 'map'
    - marginal is the name of a training data frame dimension or index

    Example:
    - som.marginal(m,2)

    """

    # check if the second argument is of type character
    if type(marginal) == str and marginal in list(map['data']):

        f_ind = list(map['data']).index(marginal)
        f_name = marginal
        train = np.matrix(map['data'])[:, f_ind]
        neurons = map['neurons'][:, f_ind]
        plt.ylabel('Density')
        plt.xlabel(f_name)
        sns.kdeplot(np.ravel(train),
                    label="training data",
                    shade=True,
                    color="b")
        sns.kdeplot(neurons, label="neurons", shade=True, color="r")
        plt.legend(fontsize=15)
        plt.show()

    elif (type(marginal) == int and marginal < map['ydim'] - 1 and
          marginal >= 0):

        f_ind = marginal
        f_name = list(map['data'])[marginal]
        train = np.matrix(map['data'])[:, f_ind]
        neurons = map['neurons'][:, f_ind]
        plt.ylabel('Density')
        plt.xlabel(f_name)
        sns.kdeplot(np.ravel(train),
                    label="training data",
                    shade=True,
                    color="b")
        sns.kdeplot(neurons, label="neurons", shade=True, color="r")
        plt.legend(fontsize=15)
        plt.show()

    else:
        sys.exit("marginal: second argument is not the name of a training \
                 data frame dimension or index")

# --------------------- local functions ---------------------------#


def embed_vm(map, conf_int=.95, verb=False):
    """ embed_vm -- using variance test and mean test to evaluate the
                    map quality

    """

    # map_df is a dataframe that contains the neurons
    map_df = map['neurons']

    # data_df is a dataframe that contain the training data
    data_df = np.array(map['data'])

    # do the F-test on a pair of datasets
    vl = df_var_test1(map_df, data_df, conf_int)

    # do the t-test on a pair of datasets
    ml = df_mean_test(map_df, data_df, conf=conf_int)

    # compute the variance captured by the map --
    # but only if the means have converged as well.
    nfeatures = map_df.shape[1]
    prob_v = significance(map, graphics=False)
    var_sum = 0

    for i in range(nfeatures):

        if (vl['conf_int_lo'][i] <= 1.0 and vl['conf_int_hi'][i] >= 1.0 and
           ml['conf_int_lo'][i] <= 0.0 and ml['conf_int_hi'][i] >= 0.0):

            var_sum = var_sum + prob_v[i]
        else:
            # not converged - zero out the probability
            prob_v[i] = 0

    # return the variance captured by converged features
    if verb:
        return prob_v
    else:
        return var_sum


def embed_ks(map, conf_int=0.95, verb=False):
    """ embed_ks -- using the kolgomorov-smirnov test to evaluate the map
                    quality

    """

    # map_df is a dataframe that contains the neurons
    map_df = map['neurons']

    # data_df is a dataframe that contain the training data
    data_df = np.array(map['data'])

    nfeatures = map_df.shape[1]

    # use the Kolmogorov-Smirnov Test to test whether the neurons and training
    # data appear
    # to come from the same distribution
    ks_vector = []
    for i in range(nfeatures):
        ks_vector.append(stats.mstats.ks_2samp(map_df[:, i], data_df[:, i]))

    prob_v = significance(map, graphics=False)
    var_sum = 0

    # compute the variance captured by the map
    for i in range(nfeatures):

        # the second entry contains the p-value
        if ks_vector[i][1] > (1 - conf_int):
            var_sum = var_sum + prob_v[i]
        else:
            # not converged - zero out the probability
            prob_v[i] = 0

    # return the variance captured by converged features
    if verb:
        return prob_v
    else:
        return var_sum


def bootstrap(map, conf_int, data_df, k, sample_acc_v):
    """ bootstrap -- compute the topographic accuracies for the given confide

    """

    ix = int(100 - conf_int*100)
    bn = 200

    bootstrap_acc_v = [np.sum(sample_acc_v)/k]

    for i in range(2, bn+1):

        bs_v = np.array([randint(1, k) for _ in range(k)])-1
        a = np.sum(list(np.array(sample_acc_v)[list(bs_v)]))/k
        bootstrap_acc_v.append(a)

    bootstrap_acc_sort_v = np.sort(bootstrap_acc_v)

    lo_val = bootstrap_acc_sort_v[ix-1]
    hi_val = bootstrap_acc_sort_v[bn-ix-1]

    return {'lo': lo_val, 'hi': hi_val}


def best_match(map, obs, full=False):
    """ best_match -- given observation obs, return the best matching neuron

    """

    # NOTE: replicate obs so that there are nr rows of obs
    obs_m = np.tile(obs, (map['neurons'].shape[0], 1))
    diff = map['neurons'] - obs_m
    squ = diff * diff
    s = np.sum(squ, axis=1)
    d = np.sqrt(s)
    o = np.argsort(d)

    if full:
        return o
    else:
        return o[0]


def accuracy(map, sample, data_ix):
    """ accuracy -- the topographic accuracy of a single sample is 1 is the
                    best matching unit and the second best matching unit are
                    neighbors otherwise it is 0

    """

    o = best_match(map, sample, full=True)
    best_ix = o[0]
    second_best_ix = o[1]

    # sanity check
    coord = coordinate(map, best_ix)
    coord_x = coord[0]
    coord_y = coord[1]

    map_ix = map['visual'][data_ix-1]
    coord = coordinate(map, map_ix)
    map_x = coord[0]
    map_y = coord[1]

    if (coord_x != map_x or coord_y != map_y or best_ix != map_ix):
        print("Error: best_ix: ", best_ix, " map_ix: ", map_ix, "\n")

    # determine if the best and second best are neighbors on the map
    best_xy = coordinate(map, best_ix)
    second_best_xy = coordinate(map, second_best_ix)
    diff_map = np.array(best_xy) - np.array(second_best_xy)
    diff_map_sq = diff_map * diff_map
    sum_map = np.sum(diff_map_sq)
    dist_map = np.sqrt(sum_map)

    # it is a neighbor if the distance on the map
    # between the bmu and 2bmu is less than 2,   should be 1 or 1.414
    if dist_map < 2:
        return 1
    else:
        return 0


def coordinate(map, rowix):
    """  coordinate -- convert from a row index to a map xy-coordinate  """
    x = (rowix) % map['xdim']
    y = (rowix) // map['xdim']
    return [x, y]


def rowix(map, x, y):
    """  rowix -- convert from a map xy-coordinate to a row index  """
    rix = x + y*map['xdim']
    return rix


def plot_heat(map, heat, explicit=False, comp=True,
              merge=False, merge_range=0.25):
    """ plot_heat - plot a heat map based on a 'map', this plot also contains
                    the connected components of the map based on the landscape
                    of the heat map

    parameters:
    - map is an object if type 'map'
    - heat is a 2D heat map of the map returned by 'map'
    - labels is a vector with labels of the original training data set
    - explicit controls the shape of the connected components
    - comp controls whether we plot the connected components on the heat map
    - merge controls whether we merge the starbursts together.
    - merge_range - a range that is used as a percentage of a certain distance
                    in the code to determine whether components are closer to
                    their centroids or centroids closer to each other.

    """

    # keep an unaltered copy of the unified distance matrix,
    # required for merging the starburst clusters
    umat = heat

    x = map['xdim']
    y = map['ydim']
    nobs = map['data'].shape[0]
    count = np.matrix([[0]*y]*x)

    # need to make sure the map doesn't have a dimension of 1
    if (x <= 1 or y <= 1):
        sys.exit("plot_heat: map dimensions too small")

    tmp = pd.cut(heat, bins=100, labels=False)

    tmp_1 = np.array(np.matrix.transpose(tmp))
    fig, ax = plt.subplots()
    ax.pcolor(tmp_1, cmap=plt.cm.YlOrRd)
    ax.set_xticks(np.arange(x)+0.5, minor=False)
    ax.set_yticks(np.arange(y)+0.5, minor=False)
    plt.xlabel("x")
    plt.ylabel("y")
    ax.set_xticklabels(np.arange(x), minor=False)
    ax.set_yticklabels(np.arange(y), minor=False)
    ax.xaxis.set_tick_params(labeltop='on')
    ax.yaxis.set_tick_params(labelright='on')

    # put the connected component lines on the map
    if comp:
        if not merge:
            # find the centroid for each neuron on the map
            centroids = compute_centroids(map,
                                          heat,
                                          explicit)
        else:
            # find the unique centroids for the neurons on the map
            centroids = compute_combined_clusters(map,
                                                  umat,
                                                  explicit,
                                                  merge_range)

        # connect each neuron to its centroid
        for ix in range(x):
            for iy in range(y):
                cx = centroids['centroid_x'][ix, iy]
                cy = centroids['centroid_y'][ix, iy]
                plt.plot([ix+0.5, cx+0.5],
                         [iy+0.5, cy+0.5],
                         color='grey',
                         linestyle='-',
                         linewidth=1.0)

    # put the labels on the map if available
    if not (map['labels'] is None) and (len(map['labels']) != 0):

        # count the labels in each map cell
        for i in range(nobs):

            nix = map['visual'][i]
            c = coordinate(map, nix)
            ix = c[0]
            iy = c[1]

            count[ix-1, iy-1] = count[ix-1, iy-1]+1

        for i in range(nobs):

            c = coordinate(map, map['visual'][i])
            ix = c[0]
            iy = c[1]

            # we only print one label per cell
            if count[ix-1, iy-1] > 0:

                count[ix-1, iy-1] = 0
                ix = ix - .5
                iy = iy - .5
                l = map['labels'][i]
                plt.text(ix+1, iy+1, l)

    plt.show()


def compute_centroids(map, heat, explicit=False):
    """ compute_centroids -- compute the centroid for each point on the map

    parameters:
    - map is an object if type 'map'
    - heat is a matrix representing the heat map representation
    - explicit controls the shape of the connected component

    return value:
    - a list containing the matrices with the same x-y dims as the original
      map containing the centroid x-y coordinates

    """

    xdim = map['xdim']
    ydim = map['ydim']
    centroid_x = np.matrix([[-1] * ydim for _ in range(xdim)])
    centroid_y = np.matrix([[-1] * ydim for _ in range(xdim)])

    heat = np.matrix(heat)

    def compute_centroid(ix, iy):

        if (centroid_x[ix, iy] > -1) and (centroid_y[ix, iy] > -1):
            return {"bestx": centroid_x[ix, iy], "besty": centroid_y[ix, iy]}

        min_val = heat[ix, iy]
        min_x = ix
        min_y = iy

        # (ix, iy) is an inner map element
        if ix > 0 and ix < xdim-1 and iy > 0 and iy < ydim-1:

            if heat[ix-1, iy-1] < min_val:
                min_val = heat[ix-1, iy-1]
                min_x = ix-1
                min_y = iy-1

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix+1, iy-1] < min_val:
                min_val = heat[ix+1, iy-1]
                min_x = ix+1
                min_y = iy-1

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

            if heat[ix+1, iy+1] < min_val:
                min_val = heat[ix+1, iy+1]
                min_x = ix+1
                min_y = iy+1

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

            if heat[ix-1, iy+1] < min_val:
                min_val = heat[ix-1, iy+1]
                min_x = ix-1
                min_y = iy+1

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # (ix, iy) is bottom left corner
        elif ix == 0 and iy == 0:

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

            if heat[ix+1, iy+1] < min_val:
                min_val = heat[ix+1, iy+1]
                min_x = ix+1
                min_y = iy+1

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

        # (ix, iy) is bottom right corner
        elif ix == xdim-1 and iy == 0:

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

            if heat[ix-1, iy+1] < min_val:
                min_val = heat[ix-1, iy+1]
                min_x = ix-1
                min_y = iy+1

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # (ix, iy) is top right corner
        elif ix == xdim-1 and iy == ydim-1:

            if heat[ix-1, iy-1] < min_val:
                min_val = heat[ix-1, iy-1]
                min_x = ix-1
                min_y = iy-1

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # (ix, iy) is top left corner
        elif ix == 0 and iy == ydim-1:

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix+1, iy-1] < min_val:
                min_val = heat[ix+1, iy-1]
                min_x = ix+1
                min_y = iy-1

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

        # (ix, iy) is a left side element
        elif ix == 0 and iy > 0 and iy < ydim-1:

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix+1, iy-1] < min_val:
                min_val = heat[ix+1, iy-1]
                min_x = ix+1
                min_y = iy-1

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

            if heat[ix+1, iy+1] < min_val:
                min_val = heat[ix+1, iy+1]
                min_x = ix+1
                min_y = iy+1

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

        # (ix, iy) is a bottom side element
        elif ix > 0 and ix < xdim-1 and iy == 0:

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

            if heat[ix+1, iy+1] < min_val:
                min_val = heat[ix+1, iy+1]
                min_x = ix+1
                min_y = iy+1

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

            if heat[ix-1, iy+1] < min_val:
                min_val = heat[ix-1, iy+1]
                min_x = ix-1
                min_y = iy+1

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # (ix, iy) is a right side element
        elif ix == xdim-1 and iy > 0 and iy < ydim-1:

            if heat[ix-1, iy-1] < min_val:
                min_val = heat[ix-1, iy-1]
                min_x = ix-1
                min_y = iy-1

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix, iy+1] < min_val:
                min_val = heat[ix, iy+1]
                min_x = ix
                min_y = iy+1

            if heat[ix-1, iy+1] < min_val:
                min_val = heat[ix-1, iy+1]
                min_x = ix-1
                min_y = iy+1

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # (ix, iy) is a top side element
        elif ix > 0 and ix < xdim-1 and iy == ydim-1:

            if heat[ix-1, iy-1] < min_val:
                min_val = heat[ix-1, iy-1]
                min_x = ix-1
                min_y = iy-1

            if heat[ix, iy-1] < min_val:
                min_val = heat[ix, iy-1]
                min_x = ix
                min_y = iy-1

            if heat[ix+1, iy-1] < min_val:
                min_val = heat[ix+1, iy-1]
                min_x = ix+1
                min_y = iy-1

            if heat[ix+1, iy] < min_val:
                min_val = heat[ix+1, iy]
                min_x = ix+1
                min_y = iy

            if heat[ix-1, iy] < min_val:
                min_val = heat[ix-1, iy]
                min_x = ix-1
                min_y = iy

        # if successful
        # move to the square with the smaller value, i_e_, call
        # compute_centroid on this new square
        # note the RETURNED x-y coords in the centroid_x and
        # centroid_y matrix at the current location
        # return the RETURNED x-y coordinates

        if min_x != ix or min_y != iy:
            r_val = compute_centroid(min_x, min_y)

            # if explicit is set show the exact connected component
            # otherwise construct a connected componenent where all
            # nodes are connected to a centrol node
            if explicit:

                centroid_x[ix, iy] = min_x
                centroid_y[ix, iy] = min_y
                return {"bestx": min_x, "besty": min_y}

            else:
                centroid_x[ix, iy] = r_val['bestx']
                centroid_y[ix, iy] = r_val['besty']
                return r_val

        else:
            centroid_x[ix, iy] = ix
            centroid_y[ix, iy] = iy
            return {"bestx": ix, "besty": iy}

    for i in range(xdim):
        for j in range(ydim):
            compute_centroid(i, j)

    return {"centroid_x": centroid_x, "centroid_y": centroid_y}


def compute_umat(map, smoothing=None):
    """ compute_umat -- compute the unified distance matrix

    parameters:
    - map is an object if type 'map'
    - smoothing is either NULL, 0, or a positive floating point value
      controlling the smoothing of the umat representation
    return value:
    - a matrix with the same x-y dims as the original map containing
      the umat values

    """

    d = euclidean_distances(map['neurons'], map['neurons'])
    umat = compute_heat(map, d, smoothing)

    return umat


def compute_heat(map, d, smoothing=None):
    """ compute_heat -- compute a heat value map representation of the
                        given distance matrix

    parameters:
    - map is an object if type 'map'
    - d is a distance matrix computed via the 'dist' function
    - smoothing is either NULL, 0, or a positive floating point value
      controlling the smoothing of the umat representation
    return value:
    - a matrix with the same x-y dims as the original map containing the heat

    """

    x = map['xdim']
    y = map['ydim']
    heat = np.matrix([[0.0] * y for _ in range(x)])

    if x == 1 or y == 1:
        sys.exit("compute_heat: heat map can not be computed for a map \
                 with a dimension of 1")

    # this function translates our 2-dim map coordinates
    # into the 1-dim coordinates of the neurons
    def xl(ix, iy):

        return ix + iy * x

    # check if the map is larger than 2 x 2 (otherwise it is only corners)
    if x > 2 and y > 2:
        # iterate over the inner nodes and compute their umat values
        for ix in range(1, x-1):
            for iy in range(1, y-1):
                sum = (d[xl(ix, iy), xl(ix-1, iy-1)] +
                       d[xl(ix, iy), xl(ix, iy-1)] +
                       d[xl(ix, iy), xl(ix+1, iy-1)] +
                       d[xl(ix, iy), xl(ix+1, iy)] +
                       d[xl(ix, iy), xl(ix+1, iy+1)] +
                       d[xl(ix, iy), xl(ix, iy+1)] +
                       d[xl(ix, iy), xl(ix-1, iy+1)] +
                       d[xl(ix, iy), xl(ix-1, iy)])

                heat[ix, iy] = sum/8

        # iterate over bottom x axis
        for ix in range(1, x-1):
            iy = 0
            sum = (d[xl(ix, iy), xl(ix+1, iy)] +
                   d[xl(ix, iy), xl(ix+1, iy+1)] +
                   d[xl(ix, iy), xl(ix, iy+1)] +
                   d[xl(ix, iy), xl(ix-1, iy+1)] +
                   d[xl(ix, iy), xl(ix-1, iy)])

            heat[ix, iy] = sum/5

        # iterate over top x axis
        for ix in range(1, x-1):
            iy = y-1
            sum = (d[xl(ix, iy), xl(ix-1, iy-1)] +
                   d[xl(ix, iy), xl(ix, iy-1)] +
                   d[xl(ix, iy), xl(ix+1, iy-1)] +
                   d[xl(ix, iy), xl(ix+1, iy)] +
                   d[xl(ix, iy), xl(ix-1, iy)])

            heat[ix, iy] = sum/5

        # iterate over the left y-axis
        for iy in range(1, y-1):
            ix = 0
            sum = (d[xl(ix, iy), xl(ix, iy-1)] +
                   d[xl(ix, iy), xl(ix+1, iy-1)] +
                   d[xl(ix, iy), xl(ix+1, iy)] +
                   d[xl(ix, iy), xl(ix+1, iy+1)] +
                   d[xl(ix, iy), xl(ix, iy+1)])

            heat[ix, iy] = sum/5

        # iterate over the right y-axis
        for iy in range(1, y-1):
            ix = x-1
            sum = (d[xl(ix, iy), xl(ix-1, iy-1)] +
                   d[xl(ix, iy), xl(ix, iy-1)] +
                   d[xl(ix, iy), xl(ix, iy+1)] +
                   d[xl(ix, iy), xl(ix-1, iy+1)] +
                   d[xl(ix, iy), xl(ix-1, iy)])

            heat[ix, iy] = sum/5

    # compute umat values for corners
    if x >= 2 and y >= 2:
        # bottom left corner
        ix = 0
        iy = 0
        sum = (d[xl(ix, iy), xl(ix+1, iy)] +
               d[xl(ix, iy), xl(ix+1, iy+1)] +
               d[xl(ix, iy), xl(ix, iy+1)])

        heat[ix, iy] = sum/3

        # bottom right corner
        ix = x-1
        iy = 0
        sum = (d[xl(ix, iy), xl(ix, iy+1)] +
               d[xl(ix, iy), xl(ix-1, iy+1)] +
               d[xl(ix, iy), xl(ix-1, iy)])
        heat[ix, iy] = sum/3

        # top left corner
        ix = 0
        iy = y-1
        sum = (d[xl(ix, iy), xl(ix, iy-1)] +
               d[xl(ix, iy), xl(ix+1, iy-1)] +
               d[xl(ix, iy), xl(ix+1, iy)])
        heat[ix, iy] = sum/3

        # top right corner
        ix = x-1
        iy = y-1
        sum = (d[xl(ix, iy), xl(ix-1, iy-1)] +
               d[xl(ix, iy), xl(ix, iy-1)] +
               d[xl(ix, iy), xl(ix-1, iy)])
        heat[ix, iy] = sum/3

    # smooth the heat map
    pts = []

    for i in range(y):
        for j in range(x):
            pts.extend([[j, i]])

    if smoothing is not None:
        if smoothing == 0:
            heat = smooth_2d(heat,
                             nrow=x,
                             ncol=y,
                             surface=False)
        elif smoothing > 0:
            heat = smooth_2d(heat,
                             nrow=x,
                             ncol=y,
                             surface=False,
                             theta=smoothing)
        else:
            sys.exit("compute_heat: bad value for smoothing parameter")

    return heat


def df_var_test1(df1, df2, conf=.95):
    """ df_var_test -- a function that applies the F-test testing the ratio
                       of the variances of the two data frames

    parameters:
    - df1,df2 - data frames with the same number of columns
    - conf - confidence level for the F-test (default .95)

    """

    if df1.shape[1] != df2.shape[1]:
        sys.exit("df_var_test: cannot compare variances of data frames")

    # init our working arrays
    var_ratio_v = [randint(1, 1) for _ in range(df1.shape[1])]
    var_confintlo_v = [randint(1, 1) for _ in range(df1.shape[1])]
    var_confinthi_v = [randint(1, 1) for _ in range(df1.shape[1])]

    def var_test(x, y, ratio=1, conf_level=0.95):

        DF_x = len(x) - 1
        DF_y = len(y) - 1
        V_x = stat.variance(x.tolist())
        V_y = stat.variance(y.tolist())

        ESTIMATE = V_x / V_y

        BETA = (1 - conf_level) / 2
        CINT = [ESTIMATE / f.ppf(1 - BETA, DF_x, DF_y),
                ESTIMATE / f.ppf(BETA, DF_x, DF_y)]

        return {"estimate": ESTIMATE, "conf_int": CINT}

    # compute the F-test on each feature in our populations
    for i in range(df1.shape[1]):

        t = var_test(df1[:, i], df2[:, i], conf_level=conf)
        var_ratio_v[i] = t['estimate']
        var_confintlo_v[i] = t['conf_int'][0]
        var_confinthi_v[i] = t['conf_int'][1]

    # return a list with the ratios and conf intervals for each feature
    return {"ratio": var_ratio_v,
            "conf_int_lo": var_confintlo_v,
            "conf_int_hi": var_confinthi_v}


def df_mean_test(df1, df2, conf=0.95):

    """ df_mean_test -- a function that applies the t-test testing the difference
                        of the means of the two data frames

    parameters:
    - df1,df2 - data frames with the same number of columns
    - conf - confidence level for the t-test (default .95)

    """

    if df1.shape[1] != df2.shape[1]:
        sys.exit("df_mean_test: cannot compare means of data frames")

    # init our working arrays
    mean_diff_v = [randint(1, 1) for _ in range(df1.shape[1])]
    mean_confintlo_v = [randint(1, 1) for _ in range(df1.shape[1])]
    mean_confinthi_v = [randint(1, 1) for _ in range(df1.shape[1])]

    def t_test(x, y, conf_level=0.95):
        estimate_x = np.mean(x)
        estimate_y = np.mean(y)
        cm = sms.CompareMeans(sms.DescrStatsW(x), sms.DescrStatsW(y))
        conf_int_lo = cm.tconfint_diff(alpha=1-conf_level, usevar='unequal')[0]
        conf_int_hi = cm.tconfint_diff(alpha=1-conf_level, usevar='unequal')[1]

        return {"estimate": [estimate_x, estimate_y],
                "conf_int": [conf_int_lo, conf_int_hi]}

    # compute the F-test on each feature in our populations
    for i in range(df1.shape[1]):
        t = t_test(x=df1[:, i], y=df2[:, i], conf_level=conf)
        mean_diff_v[i] = t['estimate'][0] - t['estimate'][1]
        mean_confintlo_v[i] = t['conf_int'][0]
        mean_confinthi_v[i] = t['conf_int'][1]

    # return a list with the ratios and conf intervals for each feature
    return {"diff": mean_diff_v,
            "conf_int_lo": mean_confintlo_v,
            "conf_int_hi": mean_confinthi_v}


def vsom_p(data, xdim, ydim, alpha, train):
    """ vsom_p - vectorized, unoptimized version of the stochastic SOM
                 training algorithm written entirely in python

    """
    # some constants
    dr = data.shape[0]
    dc = data.shape[1]
    nr = xdim*ydim
    nc = dc  # dim of data and neurons is the same

    # build and initialize the matrix holding the neurons
    cells = nr * nc  # No. of neurons times number of data dimensions

    # vector with small init values for all neurons
    v = np.random.uniform(-1, 1, cells)

    # NOTE: each row represents a neuron, each column represents a dimension.
    neurons = np.reshape(v, (nr, nc))  # rearrange the vector as matrix

    # compute the initial neighborhood size and step
    nsize = max(xdim, ydim) + 1
    nsize_step = np.ceil(train/nsize)
    step_counter = 0  # counts the number of epochs per nsize_step

    # convert a 1D rowindex into a 2D map coordinate
    def coord2D(rowix):

        x = (np.array(rowix) - 1) % xdim
        y = (np.array(rowix) - 1) // xdim

        return np.concatenate((x, y))

    # constants for the Gamma function
    m = [i for i in range(1, nr+1)]  # a vector with all neuron 1D addresses

    # x-y coordinate of ith neuron: m2Ds[i,] = c(xi, yi)
    m2Ds = np.matrix.transpose(coord2D(m).reshape(2, nr))

    # neighborhood function
    def Gamma(c):

        # lookup the 2D map coordinate for c
        c2D = m2Ds[c, ]
        # a matrix with each row equal to c2D
        c2Ds = np.outer(np.linspace(1, 1, nr), c2D)
        # distance vector of each neuron from c in terms of map coords!
        d = np.sqrt(np.dot((c2Ds - m2Ds)**2, [1, 1]))
        # if m on the grid is in neigh then alpha else 0.0
        hood = np.where(d < nsize*1.5, alpha, 0.0)

        return hood

    # training #
    # the epochs loop
    for epoch in range(1, train):

        # hood size decreases in disrete nsize.steps
        step_counter = step_counter + 1
        if step_counter == nsize_step:

            step_counter = 0
            nsize = nsize - 1

        # create a sample training vector
        ix = randint(0, dr-1)
        xk = data.iloc[[ix]]

        # competitive step
        xk_m = np.outer(np.linspace(1, 1, nr), xk)

        diff = neurons - xk_m
        squ = diff * diff
        s = np.dot(squ, np.linspace(1, 1, nc))
        o = np.argsort(s)
        c = o[0]

        # update step
        gamma_m = np.outer(Gamma(c), np.linspace(1, 1, nc))
        neurons = neurons - diff * gamma_m

    return neurons


def compute_combined_clusters(map, heat, explicit, rang):
    """ compute_combined_clusters -- to Combine connected components that
                                     represent the same cluster

        TOP LEVEL FOR DISTANCE MATRIX ORIENTED CLUSTER COMBINE
    """

    # compute the connected components
    centroids = compute_centroids(map, heat, explicit)
    # Get unique centroids
    unique_centroids = get_unique_centroids(map, centroids)
    # Get distance from centroid to cluster elements for all centroids
    within_cluster_dist = distance_from_centroids(map,
                                                  centroids,
                                                  unique_centroids,
                                                  heat)
    # Get average pairwise distance between clusters
    between_cluster_dist = distance_between_clusters(map,
                                                     centroids,
                                                     unique_centroids,
                                                     heat)
    # Get a boolean matrix of whether two components should be combined
    combine_cluster_bools = combine_decision(within_cluster_dist,
                                             between_cluster_dist,
                                             rang)
    # Create the modified connected components grid
    ne_centroid = new_centroid(combine_cluster_bools,
                               centroids,
                               unique_centroids,
                               map)

    return ne_centroid


def get_unique_centroids(map, centroids):
    """ get_unique_centroids -- a function that computes a list of unique
                                centroids from a matrix of centroid
                                locations.

    parameters:
    - map is an object of type 'map'
    - centroids - a matrix of the centroid locations in the map

    """

    # get the dimensions of the map
    xdim = map['xdim']
    ydim = map['ydim']
    xlist = []
    ylist = []
    x_centroid = centroids['centroid_x']
    y_centroid = centroids['centroid_y']

    for ix in range(xdim):
        for iy in range(ydim):
            cx = x_centroid[ix, iy]
            cy = y_centroid[ix, iy]

    # Check if the x or y of the current centroid is not in the list
    # and if not
    # append both the x and y coordinates to the respective lists
            if not(cx in xlist) or not(cy in ylist):
                xlist.append(cx)
                ylist.append(cy)

    # return a list of unique centroid positions
    return {"position_x": xlist, "position_y": ylist}


def distance_from_centroids(map, centroids, unique_centroids, heat):
    """ distance_from_centroids -- A function to get the average distance
                                   from centroid by cluster.

    parameters:
    - map is an object of type 'map'
    - centroids - a matrix of the centroid locations in the map
    - unique.centroids - a list of unique centroid locations
    - heat is a unified distance matrix

    """

    centroids_x_positions = unique_centroids['position_x']
    centroids_y_positions = unique_centroids['position_y']
    within = []

    for i in range(len(centroids_x_positions)):
        cx = centroids_x_positions[i]
        cy = centroids_y_positions[i]

        # compute the average distance
        distance = cluster_spread(cx, cy, np.matrix(heat), centroids, map)

        # append the computed distance to the list of distances
        within.append(distance)

    return within


def cluster_spread(x, y, umat, centroids, map):
    """ cluster_spread -- Function to calculate the average distance in
                          one cluster given one centroid.

    parameters:
    - x - x position of a unique centroid
    - y - y position of a unique centroid
    - umat - a unified distance matrix
    - centroids - a matrix of the centroid locations in the map
    - map is an object of type 'map'

    """

    centroid_x = x
    centroid_y = y
    sum = 0
    elements = 0
    xdim = map['xdim']
    ydim = map['ydim']
    centroid_weight = umat[centroid_x, centroid_y]

    for xi in range(xdim):
        for yi in range(ydim):
            cx = centroids['centroid_x'][xi, yi]
            cy = centroids['centroid_y'][xi, yi]

            if(cx == centroid_x and cy == centroid_y):
                cweight = umat[xi, yi]
                sum = sum + abs(cweight - centroid_weight)
                elements = elements + 1

    average = sum / elements

    return average


def distance_between_clusters(map, centroids, unique_centroids, umat):
    """ distance_between_clusters -- A function to compute the average
                                     pairwise distance between clusters.

    parameters:
    - map is an object of type 'map'
    - centroids - a matrix of the centroid locations in the map
    - unique_centroids - a list of unique centroid locations
    - umat - a unified distance matrix

    """

    cluster_elements = list_clusters(map, centroids, unique_centroids, umat)

    tmp_1 = np.zeros(shape=(max([len(cluster_elements[i]) for i in range(
            len(cluster_elements))]), len(cluster_elements)))

    for i in range(len(cluster_elements)):
        for j in range(len(cluster_elements[i])):
            tmp_1[j, i] = cluster_elements[i][j]

    columns = tmp_1.shape[1]

    tmp = np.matrix.transpose(np.array(list(combinations(
          [i for i in range(columns)], 2))))

    tmp_3 = np.zeros(shape=(tmp_1.shape[0], tmp.shape[1]))

    for i in range(tmp.shape[1]):
        tmp_3[:, i] = np.where(tmp_1[:, tmp[1, i]]*tmp_1[:, tmp[0, i]] != 0,
                               abs(tmp_1[:, tmp[0, i]]-tmp_1[:, tmp[1, i]]), 0)
        # both are not equals 0

    mean = np.true_divide(tmp_3.sum(0), (tmp_3 != 0).sum(0))
    index = 0
    mat = np.zeros((columns, columns))

    for xi in range(columns-1):
        for yi in range(xi, columns-1):
            mat[xi, yi + 1] = mean[index]
            mat[yi + 1, xi] = mean[index]
            index = index + 1

    return mat


def list_clusters(map, centroids, unique_centroids, umat):
    """ list_clusters -- A function to get the clusters as a list of lists.

    parameters:
    - map is an object of type 'map'
    - centroids - a matrix of the centroid locations in the map
    - unique.centroids - a list of unique centroid locations
    - umat - a unified distance matrix

    """

    centroids_x_positions = unique_centroids['position_x']
    centroids_y_positions = unique_centroids['position_y']
    cluster_list = []

    for i in range(len(centroids_x_positions)):
        cx = centroids_x_positions[i]
        cy = centroids_y_positions[i]

    # get the clusters associated with a unique centroid and store it in a list
        cluster_list.append(list_from_centroid(map, cx, cy, centroids, umat))

    return cluster_list


def list_from_centroid(map, x, y, centroids, umat):
    """ list.from.centroid -- A funtion to get all cluster elements
                              associated to one centroid.

    parameters:
    - map is an object of type 'map'
    - x - the x position of a centroid
    - y - the y position of a centroid
    - centroids - a matrix of the centroid locations in the map
    - umat - a unified distance matrix

    """

    centroid_x = x
    centroid_y = y
    xdim = map['xdim']
    ydim = map['ydim']

    cluster_list = []
    for xi in range(xdim):
        for yi in range(ydim):
            cx = centroids['centroid_x'][xi, yi]
            cy = centroids['centroid_y'][xi, yi]

            if(cx == centroid_x and cy == centroid_y):
                cweight = np.matrix(umat)[xi, yi]
                cluster_list.append(cweight)

    return cluster_list


def combine_decision(within_cluster_dist, distance_between_clusters, rang):
    """ combine_decision -- A function that produces a boolean matrix
                            representing which clusters should be combined.

    parameters:
    - within_cluster_dist - A list of the distances from centroid to
                            cluster elements for all centroids
    - distance_between_clusters - A list of the average pairwise distance
                                  between clusters
    - range is the distance where the clusters are merged together.

    """

    inter_cluster = distance_between_clusters
    centroid_dist = within_cluster_dist
    dim = inter_cluster.shape[1]
    to_combine = np.matrix([[False]*dim]*dim)

    for xi in range(dim):
        for yi in range(dim):
            cdist = inter_cluster[xi, yi]
            if cdist != 0:
                rx = centroid_dist[xi] * rang
                ry = centroid_dist[yi] * rang
                if (cdist < centroid_dist[xi] + rx or
                   cdist < centroid_dist[yi] + ry):
                    to_combine[xi, yi] = True

    return to_combine


def swap_centroids(map, x1, y1, x2, y2, centroids):
    """ swap_centroids -- A function that changes every instance of a
                          centroid to one that it should be combined with.
    parameters:
    - map is an object of type 'map'
    - x1 -
    - y1 -
    - x2 -
    - y2 -
    - centroids - a matrix of the centroid locations in the map

    """

    xdim = map['xdim']
    ydim = map['ydim']
    compn_x = centroids['centroid_x']
    compn_y = centroids['centroid_y']
    for xi in range(xdim):
        for yi in range(ydim):
            if compn_x[xi, 0] == x1 and compn_y[yi, 0] == y1:
                compn_x[xi, 0] = x2
                compn_y[yi, 0] = y2

    return {"centroid_x": compn_x, "centroid_y": compn_y}


def new_centroid(bmat, centroids, unique_centroids, map):
    """ new_centroid -- A function to combine centroids based on matrix of
                        booleans.

    parameters:
    - bmat - a boolean matrix containing the centroids to merge
    - centroids - a matrix of the centroid locations in the map
    - unique.centroids - a list of unique centroid locations
    - map is an object of type 'map'

    """

    bmat_rows = bmat.shape[0]
    bmat_columns = bmat.shape[1]
    centroids_x = unique_centroids['position_x']
    centroids_y = unique_centroids['position_y']
    components = centroids

    for xi in range(bmat_rows):
        for yi in range(bmat_columns):
            if bmat[xi, yi]:
                x1 = centroids_x[xi]
                y1 = centroids_y[xi]
                x2 = centroids_x[yi]
                y2 = centroids_y[yi]
                components = swap_centroids(map, x1, y1, x2, y2, components)

    return components


def smooth_2d(Y, ind=None, weight_obj=None, grid=None, nrow=64, ncol=64,
              surface=True, theta=None):
    """ smooth_2d -- Kernel Smoother For Irregular 2-D Data """

    def exp_cov(x1, x2, theta=2, p=2, distMat=0):

        x1 = x1*(1/theta)
        x2 = x2*(1/theta)

        distMat = euclidean_distances(x1, x2)

        distMat = distMat**p

        return np.exp(-distMat)

    NN = [[1]*ncol] * nrow
    grid = {'x': [i for i in range(nrow)], "y": [i for i in range(ncol)]}

    if weight_obj is None:
        dx = grid['x'][1] - grid['x'][0]
        dy = grid['y'][1] - grid['y'][0]
        m = len(grid['x'])
        n = len(grid['y'])
        M = 2 * m
        N = 2 * n

        xg = []

        for i in range(N):
            for j in range(M):
                xg.extend([[j, i]])

        xg = np.matrix(xg)

        center = []
        center.append([int(dx * M/2-1), int((dy * N)/2-1)])

        out = exp_cov(xg, np.matrix(center))
        out = np.matrix.transpose(np.reshape(out, (N, M)))
        temp = np.zeros((M, N))
        temp[int(M/2-1)][int(N/2-1)] = 1

        wght = np.fft.fft2(out)/(np.fft.fft2(temp) * M * N)
        weight_obj = {"m": m, "n": n, "N": N, "M": M, "wght": wght}

    temp = np.zeros((weight_obj['M'], weight_obj['N']))
    temp[0:m, 0:n] = Y
    temp2 = np.fft.ifft2(np.fft.fft2(temp) *
                         weight_obj['wght']).real[0:weight_obj['m'],
                                                  0:weight_obj['n']]

    temp = np.zeros((weight_obj['M'], weight_obj['N']))
    temp[0:m, 0:n] = NN
    temp3 = np.fft.ifft2(np.fft.fft2(temp) *
                         weight_obj['wght']).real[0:weight_obj['m'],
                                                  0:weight_obj['n']]

    return temp2/temp3
