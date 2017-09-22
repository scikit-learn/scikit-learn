"""
This module implement a regression tree specifically design for the
gradient-boosting regression trees.
"""
from __future__ import division, print_function

import numbers

from collections import defaultdict
from math import ceil

import numpy as np

from ..tree.tree import BaseDecisionTree
from ..base import RegressorMixin
from ..utils.validation import check_array, check_random_state, check_X_y
from ..externals import six

from ._splitter import Splitter
from ._split_record import SplitRecord
from ._criterion import impurity_mse

from ._tree import Tree
from . import _tree

TREE_UNDEFINED, TREE_LEAF, FEAT_UNKNOWN = -2, -1, -3
DTYPE = _tree.DTYPE
DOUBLE = _tree.DOUBLE


class RegressionTree(BaseDecisionTree, RegressorMixin):
    """A regression tree specifically designed for gradient-boosting.

    Parameters
    ----------
    criterion : string, optional (default="mse")
        The function to measure the quality of a split. Supported criteria
        are "mse" for the mean squared error, which is equal to variance
        reduction as feature selection criterion, and "mae" for the mean
        absolute error.

        .. versionadded:: 0.18
           Mean Absolute Error (MAE) criterion.

    splitter : string, optional (default="best")
        The strategy used to choose the split at each node. Supported
        strategies are "best" to choose the best split and "random" to choose
        the best random split.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a percentage and
          `int(max_features * n_features)` features are considered at each
          split.
        - If "auto", then `max_features=n_features`.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_depth : int or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int, float, optional (default=2)
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a percentage and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for percentages.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node:

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a percentage and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for percentages.

    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    min_impurity_split : float, optional (default=1e-7)
        Threshold for early stopping in tree growth. If the impurity
        of a node is below the threshold, the node is a leaf.

        .. versionadded:: 0.18

    Attributes
    ----------
    feature_importances_ : array of shape = [n_features]
        The feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the
        (normalized) total reduction of the criterion brought
        by that feature. It is also known as the Gini importance [4]_.

    max_features_ : int,
        The inferred value of max_features.

    n_features_ : int
        The number of features when ``fit`` is performed.

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    tree_ : Tree object
        The underlying Tree object.

    See also
    --------
    DecisionTreeClassifier, DecisionTreeRegressor

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    .. [4] L. Breiman, and A. Cutler, "Random Forests",
           http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
    """
    def __init__(self,
                 criterion="mse",
                 splitter="best",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_weight_fraction_leaf=0.,
                 max_features=None,
                 max_leaf_nodes=None,  # this parameter is not active
                 random_state=None,
                 min_impurity_split=1e-7):
        super(RegressionTree, self).__init__(
            criterion=criterion,
            splitter=splitter,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_features=max_features,
            max_leaf_nodes=None,
            random_state=random_state,
            min_impurity_split=min_impurity_split,
            presort=True)

    def fit(self, X, y, sample_weight=None, check_input=True,
            X_idx_sorted=None):
        """Build a decision tree regressor from the training set (X, y).

        Parameters
        ----------
        X : array-like or sparse matrix, shape = [n_samples, n_features]
            The training input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csc_matrix``.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (real numbers). Use ``dtype=np.float64`` and
            ``order='C'`` for maximum efficiency.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node.

        check_input : boolean, (default=True)
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.

        X_idx_sorted : array-like, shape = [n_samples, n_features], optional
            The indexes of the sorted training input samples. If many tree
            are grown on the same dataset, this allows the ordering to be
            cached between trees. If None, the data will be sorted here.
            Don't use this parameter unless you know what to do.

        Returns
        -------
        self : object
            Returns self.
        """
        random_state = check_random_state(self.random_state)
        if check_input:
            # FIXME do not accept sparse data for the moment
            X, y = check_X_y(X, y)
            X = check_array(X, dtype=DTYPE)
            y = check_array(y, ensure_2d=False, dtype=None)

        # Determine output settings
        n_samples, self.n_features_ = X.shape

        y = np.atleast_1d(y)
        expanded_class_weight = None

        if y.ndim == 1:
            # reshape is necessary to preserve the data contiguity against vs
            # [:, np.newaxis] that does not.
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]

        self.classes_ = [None] * self.n_outputs_
        self.n_classes_ = [1] * self.n_outputs_

        self.n_classes_ = np.array(self.n_classes_, dtype=np.intp)

        if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        # Check parameters
        max_depth = ((2 ** 31) - 1 if self.max_depth is None
                     else self.max_depth)
        max_leaf_nodes = (-1 if self.max_leaf_nodes is None
                          else self.max_leaf_nodes)

        if isinstance(self.min_samples_leaf, (numbers.Integral, np.integer)):
            if not 1 <= self.min_samples_leaf:
                raise ValueError("min_samples_leaf must be at least 1 "
                                 "or in (0, 0.5], got %s"
                                 % self.min_samples_leaf)
            min_samples_leaf = self.min_samples_leaf
        else:  # float
            if not 0. < self.min_samples_leaf <= 0.5:
                raise ValueError("min_samples_leaf must be at least 1 "
                                 "or in (0, 0.5], got %s"
                                 % self.min_samples_leaf)
            min_samples_leaf = int(ceil(self.min_samples_leaf * n_samples))

        if isinstance(self.min_samples_split, (numbers.Integral, np.integer)):
            if not 2 <= self.min_samples_split:
                raise ValueError("min_samples_split must be an integer "
                                 "greater than 1 or a float in (0.0, 1.0]; "
                                 "got the integer %s"
                                 % self.min_samples_split)
            min_samples_split = self.min_samples_split
        else:  # float
            if not 0. < self.min_samples_split <= 1.:
                raise ValueError("min_samples_split must be an integer "
                                 "greater than 1 or a float in (0.0, 1.0]; "
                                 "got the float %s"
                                 % self.min_samples_split)
            min_samples_split = int(ceil(self.min_samples_split * n_samples))
            min_samples_split = max(2, min_samples_split)

        min_samples_split = max(min_samples_split, 2 * min_samples_leaf)

        if isinstance(self.max_features, six.string_types):
            if self.max_features == "auto":
                max_features = self.n_features_
            elif self.max_features == "sqrt":
                max_features = max(1, int(np.sqrt(self.n_features_)))
            elif self.max_features == "log2":
                max_features = max(1, int(np.log2(self.n_features_)))
            else:
                raise ValueError(
                    'Invalid value for max_features. Allowed string '
                    'values are "auto", "sqrt" or "log2".')
        elif self.max_features is None:
            max_features = self.n_features_
        elif isinstance(self.max_features, (numbers.Integral, np.integer)):
            max_features = self.max_features
        else:  # float
            if self.max_features > 0.0:
                max_features = max(1,
                                   int(self.max_features * self.n_features_))
            else:
                max_features = 0

        self.max_features_ = max_features

        if len(y) != n_samples:
            raise ValueError("Number of labels=%d does not match "
                             "number of samples=%d" % (len(y), n_samples))
        if not 0 <= self.min_weight_fraction_leaf <= 0.5:
            raise ValueError("min_weight_fraction_leaf must in [0, 0.5]")
        if max_depth <= 0:
            raise ValueError("max_depth must be greater than zero. ")
        if not (0 < max_features <= self.n_features_):
            raise ValueError("max_features must be in (0, n_features]")
        if not isinstance(max_leaf_nodes, (numbers.Integral, np.integer)):
            raise ValueError("max_leaf_nodes must be integral number but was "
                             "%r" % max_leaf_nodes)
        if -1 < max_leaf_nodes < 2:
            raise ValueError(("max_leaf_nodes {0} must be either smaller than "
                              "0 or larger than 1").format(max_leaf_nodes))

        if sample_weight is not None:
            if (getattr(sample_weight, "dtype", None) != DOUBLE or
                    not sample_weight.flags.contiguous):
                sample_weight = np.ascontiguousarray(
                    sample_weight, dtype=DOUBLE)
            if len(sample_weight.shape) > 1:
                raise ValueError("Sample weights array has more "
                                 "than one dimension: %d" %
                                 len(sample_weight.shape))
            if len(sample_weight) != n_samples:
                raise ValueError("Number of weights=%d does not match "
                                 "number of samples=%d" %
                                 (len(sample_weight), n_samples))

        if expanded_class_weight is not None:
            if sample_weight is not None:
                sample_weight = sample_weight * expanded_class_weight
            else:
                sample_weight = expanded_class_weight

        # Set min_weight_leaf from min_weight_fraction_leaf
        if sample_weight is None:
            sample_weight = np.ones(y.size)
            min_weight_leaf = (self.min_weight_fraction_leaf *
                               n_samples)
        else:
            min_weight_leaf = (self.min_weight_fraction_leaf *
                               np.sum(sample_weight))

        if self.min_impurity_split < 0.:
            raise ValueError("min_impurity_split must be greater than "
                             "or equal to 0")

        # FIXME: to have cython buffer compatibility
        y = np.squeeze(y)

        # If multiple trees are built on the same dataset, we only want to
        # presort once. Splitters now can accept presorted indices if desired,
        # but do not handle any presorting themselves. Ensemble algorithms
        # which desire presorting must do presorting themselves and pass that
        # matrix into each tree.
        if X_idx_sorted is None and self.presort:
            X_idx_sorted = np.asfortranarray(np.argsort(X, axis=0),
                                             dtype=np.int32)

        if self.presort and X_idx_sorted.shape != X.shape:
            raise ValueError("The shape of X (X.shape = {}) doesn't match "
                             "the shape of X_idx_sorted (X_idx_sorted"
                             ".shape = {})".format(X.shape,
                                                   X_idx_sorted.shape))

        self.tree_ = Tree(self.n_features_, self.n_classes_, self.n_outputs_)

        if max_depth <= 10:
            init_capacity = (2 ** (max_depth + 1)) - 1
        else:
            init_capacity = 2047
        self.tree_._resize_py(init_capacity)

        weighted_n_samples = np.sum(sample_weight)
        # initialize the number of splitter
        n_splitters = 0
        # the array to map the samples to the correct splitter
        X_nid = np.zeros(y.size, dtype=int)
        # the list of splitter at each round
        splitter_list = []
        # the output split record
        split_record_map = defaultdict(lambda: None)

        # create the root node statistics
        root_sum_y = np.sum(np.ravel(y) * sample_weight)
        root_sum_sq_y = np.sum(np.ravel(y ** 2) * sample_weight)
        root_n_samples = n_samples
        root_sum_weighted_samples = weighted_n_samples
        # create the parent split record
        parent_split_record = SplitRecord()
        # affect the stats to the record
        parent_split_record.init_stats(root_sum_y, root_sum_sq_y,
                                       root_n_samples,
                                       root_sum_weighted_samples,
                                       0., 0., 0, 0.,
                                       root_sum_y, root_sum_sq_y,
                                       root_n_samples,
                                       root_sum_weighted_samples,)
        # compute the impurity for the parent node
        # FIXME only MSE impurity for the moment
        parent_split_record.impurity = impurity_mse(root_sum_y,
                                                    root_sum_sq_y,
                                                    root_n_samples,
                                                    root_sum_weighted_samples)

        parent_split_record.nid = self.tree_._add_node_py(
            parent=TREE_UNDEFINED,
            is_left=1, is_leaf=TREE_LEAF,
            feature=FEAT_UNKNOWN,
            threshold=TREE_UNDEFINED,
            impurity=parent_split_record.impurity,
            n_node_samples=n_samples,
            weighted_n_node_samples=weighted_n_samples,
            node_value=(root_sum_y /
                        root_sum_weighted_samples))

        # create a list to keep track of the constant features
        # constant_features = []
        constant_features = defaultdict(list)

        # Create a dictionary to store the parents split overtime
        parent_split_map = {parent_split_record.nid: parent_split_record}
        # find the node to be extended
        expandable_nids = list(parent_split_map.keys())

        current_depth = 0
        while current_depth < max_depth:
            # see if we should add or remove splitter
            n_splitters = len(expandable_nids)
            curr_n_splitters = len(splitter_list)

            # add splitters
            if n_splitters - curr_n_splitters > 0:
                splitter_list += [Splitter(X, y, sample_weight,
                                           weighted_n_samples,
                                           FEAT_UNKNOWN, TREE_UNDEFINED,
                                           parent_split_map[nid],
                                           min_samples_leaf,
                                           min_weight_leaf)
                                  for nid in expandable_nids[
                                          curr_n_splitters:]]

            # drop splitters
            else:
                splitter_list = splitter_list[:n_splitters]

            # create a dictionary from the list of splitter
            splitter_map = {nid: splitter_list[i]
                            for i, nid in enumerate(expandable_nids)}

            # Create an array from where to select randomly the feature
            shuffled_feature_idx = random_state.permutation(
                np.arange(X.shape[1]))

            # get the feature
            n_visited_feature = 0
            for feat_i in shuffled_feature_idx:

                # break the loop when enough features have been seen
                if n_visited_feature >= self.max_features_:
                    break

                # Get the sorted index
                X_col = X_idx_sorted[:, feat_i]

                # reset the splitter
                for i, nid in enumerate(expandable_nids):
                    splitter_map[nid].reset(feat_i, X_col[0],
                                            parent_split_map[nid])

                # scans all samples and evaluate all possible splits for all
                # the different splitters
                for sample_idx_sorted in X_col:
                    # Samples which are not in a leaf
                    if X_nid[sample_idx_sorted] != -1:
                        # check that the feature is not consider as constant
                        if feat_i not in constant_features[
                                X_nid[sample_idx_sorted]]:
                            # check that the sample value are different enough
                            splitter_map[X_nid[
                                sample_idx_sorted]].node_evaluate_split(
                                    sample_idx_sorted)

                b_constant = []
                # copy the split_record if the improvement is better
                for nid in expandable_nids:
                    if ((split_record_map[nid] is None) or
                            (splitter_map[
                                nid].best_split_record.impurity_improvement >
                             split_record_map[nid].impurity_improvement)):
                        split_record_map[nid] = SplitRecord()
                        splitter_map[nid].best_split_record.copy_to(
                            split_record_map[nid])
                    # declare the feature as constant if no best split was
                    # found
                    if np.isnan(
                            splitter_map[nid].best_split_record.threshold):
                        constant_features[nid].append(feat_i)
                        b_constant.append(True)
                    else:
                        b_constant.append(False)

                if not np.all(b_constant):
                    n_visited_feature += 1

            # all features have been marked as constant and we need to clean
            # the list of nodes which should have grown
            if not n_visited_feature:
                for nid in expandable_nids:
                    del parent_split_map[nid]
                expandable_nids = []

            feature_update_X_nid = []
            b_grow = False
            for nid in expandable_nids:
                # store the feature to visit for the update of X_nid
                feature_update_X_nid.append(split_record_map[nid].feature)

                # expand the tree structure
                if not np.isnan(split_record_map[nid].threshold):
                    best_split = split_record_map[nid]

                    # create the left and right which have been found
                    # from the parent splits
                    left_sr, right_sr = best_split.expand_record()

                    # the statistics for the children are not computed yet
                    # add a node for left child
                    # find out if the next node will be a lead or not
                    l_stats = left_sr.c_stats
                    left_nid = self.tree_._add_node_py(
                        parent=nid,
                        is_left=1,
                        is_leaf=TREE_LEAF,
                        feature=FEAT_UNKNOWN,
                        threshold=TREE_UNDEFINED,
                        impurity=left_sr.impurity,
                        n_node_samples=l_stats.n_samples,
                        weighted_n_node_samples=l_stats.sum_weighted_samples,
                        node_value=(l_stats.sum_y /
                                    l_stats.sum_weighted_samples))

                    # add a node for the right child
                    r_stats = right_sr.c_stats
                    right_nid = self.tree_._add_node_py(
                        parent=nid,
                        is_left=0,
                        is_leaf=TREE_LEAF,
                        feature=FEAT_UNKNOWN,
                        threshold=TREE_UNDEFINED,
                        impurity=right_sr.impurity,
                        n_node_samples=r_stats.n_samples,
                        weighted_n_node_samples=r_stats.sum_weighted_samples,
                        node_value=(r_stats.sum_y /
                                    r_stats.sum_weighted_samples))

                    # the depth has changed
                    b_grow = True

                    # Update the parent node with the found best split
                    c_stats = best_split.c_stats
                    self.tree_._update_node_py(
                        node_id=nid,
                        left_child=left_nid,
                        right_child=right_nid,
                        threshold=best_split.threshold,
                        impurity=best_split.impurity,
                        feature=best_split.feature,
                        n_node_samples=c_stats.n_samples,
                        weighted_n_node_samples=c_stats.sum_weighted_samples)

                    # update the dictionary with the new record
                    # add only the record if the impurity at the node is large
                    # enough and that the number of samples in the leaf and for
                    # the split is large enough.

                    # left child
                    b_impurity = left_sr.impurity > self.min_impurity_split
                    b_samples_split = (left_sr.c_stats.n_samples >=
                                       min_samples_split)
                    b_samples_lead = (left_sr.c_stats.n_samples >=
                                      min_samples_leaf)
                    if (b_impurity and b_samples_split and b_samples_lead):
                        parent_split_map.update({left_nid: left_sr})
                        # propagate the info about constant feature
                        constant_features[left_nid] = constant_features[nid]

                    # right child
                    b_impurity = right_sr.impurity > self.min_impurity_split
                    b_samples_split = (right_sr.c_stats.n_samples >=
                                       min_samples_split)
                    b_samples_lead = (right_sr.c_stats.n_samples >=
                                      min_samples_leaf)
                    if (b_impurity and b_samples_split and b_samples_lead):
                        parent_split_map.update({right_nid: right_sr})
                        # propagate the info about constant feature
                        constant_features[left_nid] = constant_features[nid]

                    self.counter_X_nid_labels_ = np.zeros(
                        max(parent_split_map.keys()), dtype=int)

                # we can flush the data from the parent_split_map for the
                # current node
                del parent_split_map[nid]
                del constant_features[nid]

            # the depth increased
            if b_grow:
                current_depth += 1

            # update of the expandable nodes
            expandable_nids = list(parent_split_map.keys())

            # remove redundant index of feature to visit when updating X_nid
            feature_update_X_nid = np.unique(feature_update_X_nid)

            # check that some node need to be extended before to update
            # the node index
            # make a copy of X_nid
            X_nid_tmp = X_nid.copy()
            if not expandable_nids:
                # break if we cannot grow anymore
                break
            else:
                for sample_idx in range(X.shape[0]):
                    for feat_i in feature_update_X_nid:
                        # get the index of samples to update
                        X_idx = X_idx_sorted[sample_idx, feat_i]
                        parent_nid = X_nid[X_idx]
                        # only if the sample was not a leaf
                        if parent_nid != -1:
                            if split_record_map[parent_nid].feature == feat_i:
                                # if the feature correspond, we can update the
                                # feature
                                parent_n_left_samples = split_record_map[
                                    parent_nid].l_stats.n_samples

                                # no threshold found -> this is a leaf
                                if np.isnan(split_record_map[
                                        parent_nid].threshold):
                                    X_nid_tmp[X_idx] = -1
                                else:
                                    # counter to know how many samples we
                                    # checked per splitter
                                    self.counter_X_nid_labels_[parent_nid] += 1
                                    # track how many samples we send to the
                                    # left child handle the time that several
                                    # samples are equal
                                    if (self.counter_X_nid_labels_[
                                            parent_nid] <=
                                            parent_n_left_samples):
                                        # is it a leaf
                                        if (self.tree_.children_left[
                                                parent_nid] in
                                                expandable_nids):
                                            # FIXME -> PEP8
                                            x = X_idx
                                            X_nid_tmp[
                                                x] = self.tree_.children_left[
                                                parent_nid]
                                        else:
                                            X_nid_tmp[X_idx] = -1
                                    else:
                                        # is it a leaf
                                        if (self.tree_.children_right[
                                                parent_nid] in
                                                expandable_nids):
                                            # FIXME -> PEP8
                                            x = X_idx
                                            X_nid_tmp[
                                                x] = self.tree_.children_right[
                                                parent_nid]
                                        else:
                                            X_nid_tmp[X_idx] = -1
                X_nid = X_nid_tmp

        # shrink the tree to the node count only
        rc = self.tree_._resize_c_py(self.tree_.node_count)

        if rc >= 0:
            self.tree_.max_depth = current_depth

        return self
