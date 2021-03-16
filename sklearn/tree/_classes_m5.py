"""
This module gathers methods built arount the M5 (model trees) methodology.
"""

# Authors: Sylvain Marie <sylvain.marie@schneider-electric.com>
#
# License: BSD 3 clause
from collections import namedtuple
from copy import copy
from warnings import warn

import numpy as np
from scipy.sparse import issparse

from ..metrics import mean_squared_error
from ..base import RegressorMixin, clone
from ..base import is_classifier
from ..preprocessing import StandardScaler
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..linear_model import LinearRegression

from . import _tree
from ._classes import BaseDecisionTree, DTYPE, DOUBLE
from ._classes_m5_linreg_utils import DeNormalizableMixIn, \
    linreg_model_to_text, DeNormalizableLinearRegression

__all__ = ["M5Base",
           "M5Prime"]


_SmoothingDetails = namedtuple('_SmoothingDetails', ('A', 'B', 'C'))
"""Internal structure to contain the recursive smoothed coefficients details 
in the smoothing algorithm"""


class M5Base(BaseDecisionTree):
    """
    M5Base. Implements base routines for generating M5 PredictionModel trees and rules.

    The original algorithm M5 was invented by Quinlan:

    Quinlan J. R. (1992). Learning with continuous classes. Proceedings of the
    Australian Joint Conference on Artificial Intelligence. 343--348. World
    Scientific, Singapore.


    Yong Wang made improvements and created M5':

    Wang, Y and Witten, I. H. (1997). Induction of model trees for predicting
    continuous classes. Proceedings of the poster papers of the European
    Conference on Machine Learning. University of Economics, Faculty of
    Informatics and Statistics, Prague.

    Pruning and Smoothing can be activated and deactivated on top of the base
    model. TODO check if 'rules' should be supported too
    
    Inspired by Weka (https://github.com/bnjmn/weka) M5Base class, from Mark Hall
    """

    def __init__(self,
                 criterion='friedman_mse',
                 # According to M5' paper, mse should lead to similar results than (std/rmse), that is not implemented in sklearn.
                 splitter='best',  # take the feature with the best split
                 max_depth=None,  # no maximum depth limit
                 min_samples_split=4,  # in weka this parameter is named "M"
                 min_samples_leaf=2,  # from the M5' article : otherwise (n+v)/(n-v) is infinite
                 min_weight_fraction_leaf=0.,  # TODO this would actually maybe be better than min_sample_leaf ?
                 max_features=None,  # no feature reduction: take all features
                 max_leaf_nodes=None,  # no limitation in number of leaves
                 min_impurity_decrease=0.,
                 # min_impurity_split_as_initial_ratio = 0.05,  # TODO The paper suggests to use 5% of all STDR
                 min_impurity_split=None,
                 random_state=None,
                 class_weight=None,
                 leaf_model=None,  # the regression model used in the leaves (it will be cloned for each leaf)
                 use_pruning=True,
                 use_smoothing=None,  # False, True, or 'installed', or 'on_prediction'. None and True means 'installed' by default except if smoothing constant or smoothing constant ratio is 0.0
                 smoothing_constant=None,  # The number of samples (default=15) that a parent model represents when mixed with a child model. 0 means disable smoothing
                 smoothing_constant_ratio=None,  # if provided, smoothing_constant will not be used. 0.0 means disable smoothing
                 debug_prints=False,
                 ccp_alpha=0.0  # TODO is this relevant ? does that conflict with "use_pruning" ?
                 ):

        # ------ WEKA VERSION FOR REFERENCE ---
        # From https://github.com/bnjmn/weka/blob/master/weka/src/main/java/weka/classifiers/trees/m5/RuleNode.java
        # it is like the M5' paper: a minimum number + a rule on standard deviation ratio
        # if ((m_numInstances < m_splitNum)
        #         | | (Rule.stdDev(m_classIndex, m_instances) < (m_globalDeviation * m_devFraction))) {
        # m_isLeaf = true;
        # Note: their impurity criterion for splitting is a standard deviation reduction (SDR) as stated in the paper:
        # see https://github.com/Waikato/weka-trunk/blob/cd79b7f15ffeb078392aae08f5a4e9ae38620a98/weka/src/main/java/weka/classifiers/trees/m5/Impurity.java
        # it should be equivalent to the mse one

        # ------ THIS (SKLEARN BASE) VERSION ----
        # From https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/tree.py
        # max_leaf_nodes = -1 if max_leaf_nodes is None
        # if max_leaf_nodes < 0:
        #     builder = DepthFirstTreeBuilder(splitter, min_samples_split,
        #                                     min_samples_leaf,
        #                                     min_weight_leaf,
        #                                     max_depth,
        #                                     self.min_impurity_decrease,
        #                                     min_impurity_split)

        # From https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx
        # is_leaf = (depth >= max_depth or
        #            n_node_samples < min_samples_split or
        #            n_node_samples < 2 * min_samples_leaf or
        #            weighted_n_node_samples < 2 * min_weight_leaf or
        #            impurity <= min_impurity_split)  (this one is deprecated)
        #
        # if not is_leaf:
        #     DO THE SPLIT :)
        #     splitter.node_split(impurity, & split, & n_constant_features)
        #
        #     is_leaf = (is_leaf or split.pos >= end or
        #                (split.improvement + EPSILON < min_impurity_decrease))

        # TODO the paper suggests to do this with 5%
        # min_impurity_split = min_impurity_split_as_initial_ratio * initial_impurity

        super(M5Base, self).__init__(
            criterion=criterion,
            splitter=splitter,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_features=max_features,
            max_leaf_nodes=max_leaf_nodes,
            min_impurity_decrease=min_impurity_decrease,
            min_impurity_split=min_impurity_split,
            random_state=random_state,
            class_weight=class_weight,
            ccp_alpha=ccp_alpha
        )

        # warning : if the field names are different from constructor params,
        # then clone(self) will not work.
        if leaf_model is None:
            # to handle case when the model is learnt on normalized data and
            # we wish to be able to read the model equations.
            leaf_model = DeNormalizableLinearRegression()
        self.leaf_model = leaf_model

        # smoothing related
        if smoothing_constant_ratio is not None and smoothing_constant is not None:
            raise ValueError("Only one of `smoothing_constant` and "
                             "`smoothing_constant_ratio` should be provided")
        elif (smoothing_constant_ratio == 0.0 or smoothing_constant == 0) \
                and (use_smoothing is True or use_smoothing == 'installed'):
            raise ValueError("`use_smoothing` was explicitly enabled with "
                             "pre-installed models, while smoothing "
                             "constant/ratio are explicitly set to zero")

        self.use_pruning = use_pruning
        self.use_smoothing = use_smoothing
        self.smoothing_constant = smoothing_constant
        self.smoothing_constant_ratio = smoothing_constant_ratio

        self.debug_prints = debug_prints

    def as_pretty_text(self, **kwargs):
        """
        Returns a multi-line representation of this decision tree, using 
        `tree_to_text`.

        :return: a multi-line string representing this decision tree
        """
        from ._export_m5p import export_text_m5
        return export_text_m5(self, out_file=None, **kwargs)

    def fit(self, X, y, sample_weight=None, check_input=True,
            X_idx_sorted="deprecated"):

        # (0) smoothing default values behaviour
        if (self.smoothing_constant_ratio == 0.0 
                or self.smoothing_constant == 0):
            self.use_smoothing = False
        elif self.use_smoothing is None or self.use_smoothing is True:
            # default behaviour
            if isinstance(self.leaf_model, LinearRegression):
                self.use_smoothing = 'installed'
            else:
                self.use_smoothing = 'on_prediction'
        # finally make sure we are ok now
        assert self.use_smoothing in (False, np.bool_(False), 'installed',
                                      'on_prediction')

        # (1) Build the initial tree as usual
        super(M5Base, self).fit(
            X, y,
            sample_weight=sample_weight,
            check_input=check_input,
            X_idx_sorted=X_idx_sorted)

        if self.debug_prints:
            print("(debug_prints) Initial tree:")
            print(self.as_pretty_text(node_ids=True))

        # (1b) prune initial tree to take into account min impurity in splits
        # TODO since we're not on sklearn 0.17 anymore this can maybe be
        #  replaced with proper usage of self.min_impurity_split
        prune_on_min_impurity(self.tree_)

        if self.debug_prints:
            print("(debug_prints) Postprocessed tree:")
            print(self.as_pretty_text(node_ids=True))

        # (2) Now prune tree and replace pruned branches with linear models
        # -- unfortunately we have to re-do this input validation
        # step, because it also converts the input to float32.
        if check_input:
            X = check_array(X, dtype=DTYPE, accept_sparse="csc")
            if issparse(X):
                X.sort_indices()

                if X.indices.dtype != np.intc or X.indptr.dtype != np.intc:
                    raise ValueError("No support for np.int64 index based "
                                     "sparse matrices")

        # -- initialise the structure to contain the leaves and node models
        self.node_models = np.empty((self.tree_.node_count,), dtype=object)

        # -- Pruning requires to know the global deviation of the target
        global_std_dev = np.nanstd(y)
        # global_abs_dev = np.nanmean(np.abs(y))

        # -- Pruning requires to know the samples that reached each node.
        # From http://scikit-learn.org/stable/auto_examples/tree/plot_unveil_tree_structure.html
        # Retrieve the decision path of each sample.
        samples_to_nodes = self.decision_path(X)
        # * row i is to see the nodes (non-empty j) in which sample i appears.
        #   sparse row-first (CSR) format is OK
        # * column j is to see the samples (non-empty i) that fall into that
        #   node. To do that, we need to make it CSC
        nodes_to_samples = samples_to_nodes.tocsc()

        # -- execute the pruning
        # TODO if "self.use_pruning"
        # self.features_usage is a dict feature_idx -> nb times used, only
        #   for used features.
        self.features_usage = build_models_and_get_pruning_info(
            self.tree_, X, y, nodes_to_samples, self.leaf_model,
            self.node_models, global_std_dev, use_pruning=self.use_pruning
        )

        # -- cleanup to compress inner structures: only keep non-pruned ones
        self._cleanup_tree()

        if self.debug_prints:
            print("(debug_prints) Pruned tree:")
            print(self.as_pretty_text(node_ids=True))

        if self.use_smoothing == 'installed':
            # Retrieve the NEW decision path of each sample.
            samples_to_nodes = self.decision_path(X)
            # * row i is to see the nodes (non-empty j) in which sample i
            #   appears. sparse row-first (CSR) format is OK
            # * column j is to see the samples (non-empty i) that fall into
            #   that node. To do that, we need to make it CSC
            nodes_to_samples = samples_to_nodes.tocsc()

            # default behaviour for smoothing constant and ratio
            smoothing_constant = self._get_smoothing_constant_to_use(X)

            self.install_smoothing(X, y, nodes_to_samples,
                                   smoothing_constant=smoothing_constant)

            if self.debug_prints:
                print("(debug_prints) Pruned and smoothed tree:")
                print(self.as_pretty_text(node_ids=True))

        return self

    def _get_smoothing_constant_to_use(self, X):
        """
        Returns the smoothing_constant to use for smoothing, based on current 
        settings and X data
        """
        if self.smoothing_constant_ratio is not None:
            nb_training_samples = X.shape[0]
            smoothing_cstt = self.smoothing_constant_ratio * nb_training_samples
            if smoothing_cstt < 15:
                warn("smoothing constant ratio %s is leading to an extremely "
                     "small smoothing constant %s because nb training samples"
                     " is %s. Clipping to 15."
                     % (self.smoothing_constant_ratio, smoothing_cstt, 
                        nb_training_samples))
                smoothing_cstt = 15
        else:
            smoothing_cstt = self.smoothing_constant
            if smoothing_cstt is None:
                smoothing_cstt = 15  # default from the original Quinlan paper
                
        return smoothing_cstt

    def _cleanup_tree(self):
        """
        Reduces the size of this object by removing from internal structures 
        all items that are not used any more (all leaves that have been pruned)
        """
        old_tree = self.tree_
        old_node_modls = self.node_models

        # Get all information to create a copy of the inner tree. 
        # Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = _tree.Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # cleanup: only keep the nodes that are not undefined.
        # note: this is identical to
        # valid_nodes_indices = dct["nodes"]['left_child'] != TREE_UNDEFINED
        valid_nodes_indices = old_tree.children_left != _tree.TREE_UNDEFINED
        # valid_nodes_indices2 = old_tree.children_right != TREE_UNDEFINED
        # assert np.all(valid_nodes_indices == valid_nodes_indices2)
        new_node_count = sum(valid_nodes_indices)

        # create empty new structures
        n_shape = (new_node_count, *dct["nodes"].shape[1:])
        new_nodes = _empty_contig_ar(n_shape, dtype=dct["nodes"].dtype)
        v_shape = (new_node_count, *dct["values"].shape[1:])
        new_values = _empty_contig_ar(v_shape, dtype=dct["values"].dtype)
        m_shape = (new_node_count, *old_node_modls.shape[1:])
        new_node_models = _empty_contig_ar(m_shape, dtype=old_node_modls.dtype)

        # Fill structures while reindexing the tree and remembering the depth
        global next_free_id
        next_free_id = 0

        def _compress(old_node_id):
            """

            :param old_node_id:
            :return: the depth and new indices of left and right children
            """
            global next_free_id
            new_node_id = next_free_id
            next_free_id += 1

            # use the old tree to walk
            old_node = dct["nodes"][old_node_id]
            left_id = old_node['left_child']
            right_id = old_node['right_child']

            # Create the new node with a copy of the old
            new_nodes[new_node_id] = old_node  # this is an entire row so it is probably copied already by doing so.
            new_values[new_node_id] = dct["values"][old_node_id]
            new_node_models[new_node_id] = copy(old_node_modls[old_node_id])

            if left_id == _tree.TREE_LEAF:
                # ... and right_id == _tree.TREE_LEAF
                # => node_id is a leaf. Nothing to do
                return 1, new_node_id
            else:

                left_depth, new_id_left = _compress(left_id)
                right_depth, new_id_right = _compress(right_id)

                # store the new indices
                new_nodes[new_node_id]['left_child'] = new_id_left
                new_nodes[new_node_id]['right_child'] = new_id_right

                return 1 + max(left_depth, right_depth), new_node_id

        # complete definition of the new tree
        dct["max_depth"] = _compress(0)[0] - 1  # root node has depth 0, not 1
        dct["node_count"] = new_node_count  # new_nodes.shape[0]
        dct["nodes"] = new_nodes
        dct["values"] = new_values
        new_tree.__setstate__(dct)

        # Fix an issue on sklearn 0.17.1: setstate was not updating max_depth
        # see https://github.com/scikit-learn/scikit-learn/blob/0.17.1/sklearn/tree/_tree.pyx#L623
        new_tree.max_depth = dct["max_depth"]

        # update self fields
        self.tree_ = new_tree
        self.node_models = new_node_models

    def install_smoothing(
            self, 
            X_train_all, y_train_all, nodes_to_samples, smoothing_constant
    ):
        """
        Executes the smoothing procedure described in the M5 and M5P paper, 
        "once for all". This means that all models are modified so that after 
        this method has completed, each model in the tree is already a smoothed
        model.

        This has pros (prediction speed) and cons (the model equations are 
        harder to read - lots of redundancy)

        It can only be done if leaf models are instances of `LinearRegression`
        """
        # check validity: leaf models have to support pre-computed smoothing
        if not isinstance(self.leaf_model, LinearRegression):
            raise TypeError("`install_smoothing` is only available if leaf "
                            "models are instances of `LinearRegression` or a "
                            "subclass")

        # Select the Error metric to compute model errors (used to compare 
        # node with subtree for pruning)
        # TODO in Weka they use RMSE, but in papers they use MAE. this could be a parameter.
        # TODO Shall we go further and store the residuals, or a bunch of metrics? not sure
        err_metric = root_mean_squared_error  # mean_absolute_error

        # --- Node info
        # TODO we should probably not do this once in each method, but once or give access directly (no copy)
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = self.tree_.__getstate__().copy()
        old_node_models = self.node_models

        m_shape = old_node_models.shape
        new_node_models = _empty_contig_ar(m_shape, dtype=old_node_models.dtype)

        def smooth__(coefs_,
                     n_samples,
                     features=None,
                     parent=None,  # type: _SmoothingDetails
                     parent_features=None,
                     k=smoothing_constant  #
                     ):
            # type: (...) -> _SmoothingDetails
            """
            Smoothes the model coefficients or intercept `coefs_` (an array or
            a scalar), using the smoothing results at parent node.
            At each node we keep in memory three results A, B, C.

             - the new coef at each node is A + B
             - the recursion equations are
                - B(n) = (k / (n_samples(n) + k)) * A(n-1) + B(n-1)
                - C(n) = (n_samples(n) / (n_samples(n) + k)) * C(n-1)
                - A(n) = coef(n)*C(n)
            """
            if parent is None:
                # A0 = coef(0), C0 = 1, B0 = 0
                if features is None:
                    # single scalar value
                    return _SmoothingDetails(A=coefs_, B=0, C=1)
                else:
                    # vector
                    if len(coefs_) == 0:
                        coefs_ = np.asarray(coefs_)
                    if coefs_.shape[0] != len(features):
                        raise ValueError("nb features does not match the nb "
                                         "of coefficients")
                    return _SmoothingDetails(
                        A=coefs_, B=np.zeros(coefs_.shape, dtype=coefs_.dtype),
                        C=np.ones(coefs_.shape, dtype=coefs_.dtype)
                    )
            else:
                # B(n) = k/(n_samples(n)+k)*A(n-1) + B(n-1)
                # C(n) = (n_samples(n) / (n_samples(n) + k)) * C(n-1)
                Bn = (k / (n_samples + k)) * parent.A + parent.B
                Cn = (n_samples / (n_samples + k)) * parent.C

                # A(n) = coef(n) * C(n)
                if features is None:
                    # single scalar value: easy
                    An = coefs_ * Cn
                    return _SmoothingDetails(A=An, B=Bn, C=Cn)
                else:
                    # vector of coefs: we have to 'expand' the coefs array
                    # because the coefs at this node apply for (features)
                    # while coefs at parent node apply for (parent_features)
                    An = np.zeros(Cn.shape, dtype=Cn.dtype)
                    parent_features = np.array(parent_features)
                    features = np.array(features)

                    # thanks https://stackoverflow.com/a/8251757/7262247 !
                    index = np.argsort(parent_features)
                    sorted_parents = parent_features[index]
                    sorted_index = np.searchsorted(sorted_parents, features)

                    features_index = np.take(index, sorted_index, mode="clip")
                    if np.any(parent_features[features_index] != features):
                        # result = np.ma.array(features_index, mask=mask)
                        raise ValueError("Internal error - please report this."
                                         "One feature was found in the child "
                                         "node, that was not in the parent "
                                         "node.")

                    if len(features_index) > 0:
                        An[features_index] = coefs_ * Cn[features_index]

                    # old inefficient version
                    # for i, feature in enumerate(parent_features):
                    #     try:
                    #         j = features.index(feature)
                    #     except IndexError:
                    #         pass
                    #     else:
                    #         An[i] = coefs_[features[j]] * Cn[i]

                    return _SmoothingDetails(A=An, B=Bn, C=Cn)

        def _smooth(node_id, parent_features=None,
                    parent_coefs=None,     # type: _SmoothingDetails
                    parent_intercept=None  # type: _SmoothingDetails
                    ):
            # gather all info on this node
            # --base regression tree
            node_info = dct["nodes"][node_id]
            left_id = node_info['left_child']
            right_id = node_info['right_child']
            # --additional model
            node_model = old_node_models[node_id]
            # --samples
            samples_at_this_node = get_samples_at_node(node_id,
                                                       nodes_to_samples)
            n_samples_at_this_node = samples_at_this_node.shape[0]
            # note: should be equal to tree.n_node_samples[node_id]
            assert n_samples_at_this_node == self.tree_.n_node_samples[node_id]
            y_true_this_node = y_train_all[samples_at_this_node]
            X_this_node = X_train_all[samples_at_this_node, :]

            # (1) smooth current node
            # parent_model = old_node_modls[parent_id] if parent_id is not None
            parent_features = parent_features if parent_features is not None else None
            is_constant_leaf = False
            if (left_id == _tree.TREE_LEAF
                    and isinstance(node_model, ConstantLeafModel)):
                is_constant_leaf = True
                node_features = ()
                smoothed_features = (parent_features
                                     if parent_features is not None
                                     else node_features)
                node_coefs = ()
                node_intercept = dct["values"][node_id]
                # extract the unique scalar value
                assert len(node_intercept) == 1
                node_intercept = node_intercept.item()
            else:
                # a leaf LinRegLeafModel or a split SplitNodeModel
                node_features = node_model.features
                smoothed_features = (parent_features
                                     if parent_features is not None
                                     else node_features)
                node_coefs = node_model.model.coef_
                node_intercept = node_model.model.intercept_

            # create a new linear regression model with smoothed coefficients
            smoothed_sklearn_model = clone(self.leaf_model)
            smoothed_coefs = smooth__(node_coefs, features=node_features,
                                      n_samples=n_samples_at_this_node,
                                      parent=parent_coefs,
                                      parent_features=parent_features)
            smoothed_intercept = smooth__(node_intercept,
                                          n_samples=n_samples_at_this_node,
                                          parent=parent_intercept)
            smoothed_sklearn_model.coef_ = smoothed_coefs.A + smoothed_coefs.B
            smoothed_sklearn_model.intercept_ = (smoothed_intercept.A
                                                 + smoothed_intercept.B)

            # finally update the node
            if is_constant_leaf:
                smoothed_node_model = LinRegLeafModel(
                    smoothed_features, smoothed_sklearn_model, None
                )
            else:
                smoothed_node_model = copy(node_model)
                smoothed_node_model.features = smoothed_features
                smoothed_node_model.model = smoothed_sklearn_model

            # remember the new smoothed model
            new_node_models[node_id] = smoothed_node_model

            if left_id == _tree.TREE_LEAF:
                # If this is a leaf, update the prediction error on X
                y_pred_this_node = smoothed_node_model.predict(X_this_node)
                smoothed_node_model.error = err_metric(
                    y_true_this_node, y_pred_this_node
                )

            else:
                # If this is a split node - recurse on each subtree
                _smooth(left_id, parent_features=smoothed_features,
                        parent_coefs=smoothed_coefs,
                        parent_intercept=smoothed_intercept)
                _smooth(right_id, parent_features=smoothed_features,
                        parent_coefs=smoothed_coefs,
                        parent_intercept=smoothed_intercept)

                # Update the error using the same formula than the one we used
                # in build_models_and_get_pruning_info
                y_pred_children = predict_from_leaves_no_smoothing(
                    self.tree_, new_node_models, X_this_node
                )
                err_children = err_metric(y_true_this_node, y_pred_children)

                # the total number of parameters is the sum of params in each
                # branch PLUS 1 for the split
                n_params_splitmodel = (new_node_models[left_id].n_params
                                       + new_node_models[right_id].n_params
                                       + 1)
                smoothed_node_model.n_params = n_params_splitmodel
                # do not adjust the error now, simply store the raw one
                # smoothed_node_model.error = compute_adjusted_error(
                #    err_children, n_samples_at_this_node, n_params_splitmodel)
                smoothed_node_model.error = err_children

            return

        # smooth the whole tree
        _smooth(0)

        # use the new node models now
        self.node_models = new_node_models

        # remember the smoothing constant installed
        self.installed_smoothing_constant = smoothing_constant
        return

    def denormalize(self, x_scaler, y_scaler):
        """
        De-normalizes this model according to the provided x and y normalization scalers.
        Currently only StandardScaler issupported.

        :param x_scaler: a StandardScaler or None
        :param y_scaler: a StandardScaler or None
        :return:
        """
        # perform denormalization
        self._denormalize_tree(x_scaler, y_scaler)

    def _denormalize_tree(self, x_scaler, y_scaler):
        """
        De-normalizes all models in the tree
        :return:
        """
        old_tree = self.tree_
        old_node_models = self.node_models

        # Get all information to create a copy of the inner tree.
        # Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = _tree.Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # create empty new structures for nodes and models
        n_shape = dct["nodes"].shape
        new_nodes = _empty_contig_ar(n_shape, dtype=dct["nodes"].dtype)
        v_shape = dct["values"].shape
        new_values = _empty_contig_ar(v_shape, dtype=dct["values"].dtype)
        m_shape = old_node_models.shape
        new_node_models = _empty_contig_ar(m_shape, dtype=old_node_models.dtype)

        def _denormalize(node_id, x_scaler, y_scaler):
            """
            denormalizes the subtree below node `node_id`.

             - `new_nodes[node_id]` is filled with a copy of the old node
               (see old_node.dtype to see the various fields). If the node is
               a split node, the split threshold
               `new_nodes[node_id]['threshold']` is denormalized.

             - `new_values[node_id]` is filled with a denormalized copy of the
                old constant prediction at this node. Reminder: this constant
                prediction is actually used only when on the leaves now, but
                we keep it for ref.

             - `new_node_models[node_id]` is filled with a copy of the old
               model `old_node_models[node_id]`, that is denormalized if it
               is not a constant leaf

            :param node_id:
            :return: (nothing)
            """
            # use the old tree to walk
            old_node = dct["nodes"][node_id]
            left_id = old_node['left_child']
            right_id = old_node['right_child']

            # Create the new node with a copy of the old
            new_nodes[node_id] = old_node  # this is an entire row so it is probably copied already by doing so.
            new_model = copy(old_node_models[node_id])
            new_node_models[node_id] = new_model

            # Create the new value by de-scaling y
            if y_scaler is not None:
                # Note: if this is a split node with a linear regression model
                # the value will never be used. However to preserve consistency
                # of the whole values structure and for debugging purposes, we
                # choose this safe path of denormalizing ALL.
                # TODO we could also do it at once outside of the recursive
                #  calls, but we should check for side effects
                new_values[node_id] = y_scaler.inverse_transform(
                    dct["values"][node_id]
                )
            else:
                # no denormalization: simply copy
                new_values[node_id] = dct["values"][node_id]

            if left_id == _tree.TREE_LEAF:
                # ... and right_id == _tree.TREE_LEAF
                # => node_id is a leaf
                if isinstance(new_model, ConstantLeafModel):
                    # nothing to do: we already re-scaled the value
                    return
                elif isinstance(new_model, LinRegLeafModel):
                    # denormalize model
                    new_model.denormalize(x_scaler, y_scaler)
                    return
                else:
                    raise TypeError("Internal error - Leafs can only be"
                                    "constant or linear regression")
            else:
                # this is a split node, denormalize each subtree
                _denormalize(left_id, x_scaler, y_scaler)
                _denormalize(right_id, x_scaler, y_scaler)

                # denormalize the split value if needed
                if x_scaler is not None:
                    split_feature = old_node['feature']
                    # The denormalizer requires a vector with all the features,
                    # even if we only want to denormalize one.
                    # -- put split value in a vector where it has pos 'feature'
                    old_threshold_and_zeros = np.zeros((self.n_features_, ), dtype=dct["nodes"]['threshold'].dtype)
                    old_threshold_and_zeros[split_feature] = old_node['threshold']
                    # -- denormalize the vector and retrieve value 'feature'
                    new_nodes[node_id]['threshold'] = x_scaler.inverse_transform(old_threshold_and_zeros)[split_feature]
                else:
                    # no denormalization: simply copy
                    new_nodes[node_id]['threshold'] = old_node['threshold']

                if isinstance(new_model, SplitNodeModel):
                    # denormalize model at split node too, even if it is not
                    # always used (depending on smoothing mode)
                    new_model.denormalize(x_scaler, y_scaler)
                else:
                    raise TypeError("Internal error: all intermediate nodes"
                                    "should be SplitNodeModel")

                return

        # denormalize the whole tree and put the result in (new_nodes,
        # new_values, new_node_models) recursively
        _denormalize(0, x_scaler, y_scaler)

        # complete definition of the new tree
        dct["nodes"] = new_nodes
        dct["values"] = new_values
        new_tree.__setstate__(dct)

        # update the self fields
        # self.features_usage
        self.tree_ = new_tree
        self.node_models = new_node_models

    def compress_features(self):
        """
        Compresses the model and returns a vector of required feature indices.
        This model input will then be X[:, features] instead of X.
        """
        used_features = sorted(self.features_usage.keys())
        new_features_lookup_dct = {old_feature_idx: i
                                   for i, old_feature_idx
                                   in enumerate(used_features)}

        if used_features == list(range(self.n_features_)):
            # NOTHING TO DO: we need all features
            return used_features

        # Otherwise, We can compress. For this we have to create a copy of the
        # tree because we will change its internals
        old_tree = self.tree_
        old_node_modls = self.node_models

        # Get all information to create a copy of the inner tree.
        # Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = _tree.Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # create empty new structures for nodes and models
        n_shape = dct["nodes"].shape
        new_nodes = _empty_contig_ar(n_shape, dtype=dct["nodes"].dtype)
        m_shape = old_node_modls.shape
        new_node_models = _empty_contig_ar(m_shape, dtype=old_node_modls.dtype)

        def _compress_features(node_id):
            """

            :param node_id:
            :return: the depth and new indices of left and right children
            """
            # use the old tree to walk
            old_node = dct["nodes"][node_id]
            left_id = old_node['left_child']
            right_id = old_node['right_child']

            # Create the new node with a copy of the old
            new_nodes[node_id] = old_node  # this is an entire row so it is probably copied already by doing so.
            new_model = copy(old_node_modls[node_id])
            new_node_models[node_id] = new_model

            if left_id == _tree.TREE_LEAF:
                # ... and right_id == _tree.TREE_LEAF
                # => node_id is a leaf
                if isinstance(new_model, ConstantLeafModel):
                    # no features used
                    return
                elif isinstance(new_model, LinRegLeafModel):
                    # compress that model
                    new_model.reindex_features(new_features_lookup_dct)
                    return
                else:
                    raise TypeError("Internal error - Leafs can only be"
                                    "constant or linear regression")
            else:

                _compress_features(left_id)
                _compress_features(right_id)

                # store the new split feature index in the node
                new_nodes[node_id]['feature'] = (
                    new_features_lookup_dct[old_node['feature']]
                )

                if isinstance(new_model, SplitNodeModel):
                    # TODO now that split node models have a linear regression
                    #  model too, we should update them here
                    pass
                else:
                    raise TypeError("Internal error: all intermediate nodes"
                                    "should be SplitNodeModel")

                return

        _compress_features(0)

        # complete definition of the new tree
        dct["nodes"] = new_nodes
        new_tree.__setstate__(dct)

        # update the self fields
        self.features_usage = {new_features_lookup_dct[k]: v
                               for k, v in self.features_usage.items()}
        self.tree_ = new_tree
        self.node_models = new_node_models
        self.n_features_ = len(self.features_usage)

        # return the vector of used features
        return used_features

    @property
    def feature_importances_(self):
        """Return the feature importances (the higher, the more important).

        Returns
        -------
        feature_importances_ : array, shape = [n_features]
        """
        check_is_fitted(self, 'tree_')

        # TODO when m5tree inherits from tree,compute_feature_importances below
        features = np.array([self.features_usage[k]
                             for k in sorted(self.features_usage.keys())],
                            dtype=int)

        return features

    # def compute_feature_importances(self, normalize=True):
    #     """Computes the importance of each feature (aka variable)."""
    #     nodes = self.nodes
    #     node = nodes
    #     end_node = node + self.node_count
    #
    #     normalizer = 0.
    #
    #     np.ndarray[np.float64_t, ndim=1] importances
    #     importances = np.zeros((self.n_features,))
    #     cdef DOUBLE_t* importance_data = <DOUBLE_t*>importances.data
    #
    #     with nogil:
    #         while node != end_node:
    #             if node.left_child != _TREE_LEAF:
    #                 # ... and node.right_child != _TREE_LEAF:
    #                 left = &nodes[node.left_child]
    #                 right = &nodes[node.right_child]
    #
    #                 importance_data[node.feature] += (
    #                     node.weighted_n_node_samples * node.impurity -
    #                     left.weighted_n_node_samples * left.impurity -
    #                     right.weighted_n_node_samples * right.impurity)
    #             node += 1
    #
    #     importances /= nodes[0].weighted_n_node_samples
    #
    #     if normalize:
    #         normalizer = np.sum(importances)
    #
    #         if normalizer > 0.0:
    #             # Avoid dividing by zero (e.g., when root is pure)
    #             importances /= normalizer
    #
    #     return importances

    def predict(self, X, check_input=True, smooth_predictions=None,
                smoothing_constant=None):
        """Predict class or regression value for X.

        For a classification model, the predicted class for each sample in X is
        returned. For a regression model, the predicted value based on X is
        returned.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        check_input : boolean, (default=True)
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.

        smooth_predictions : boolean, (default=None)
            None means "use self config"

        smoothing_constant: int, (default=self.smoothing_constant)
            Smoothing constant used when smooth_predictions is True. During
            smoothing, child node models are recursively mixed with their
            parent node models, and the mix is done using a weighted sum. The
            weight given to the child model is the number of training samples
            that reached its node, while the weight given to the parent model
            is the smoothing constant. Therefore it can be seen as an
            equivalent number of samples that parent models represent when
            injected into the recursive weighted sum.

        Returns
        -------
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes, or the predict values.
        """
        check_is_fitted(self, 'tree_')

        # If this is just a constant node only check the input's shape.
        if self.n_features_ == 0:
            # perform the input checking manually to set ensure_min_features=0
            X = check_array(X, dtype=DTYPE, accept_sparse="csr",
                            ensure_min_features=0)
            if issparse(X) and (X.indices.dtype != np.intc or
                                X.indptr.dtype != np.intc):
                raise ValueError("No support for np.int64 index based "
                                 "sparse matrices")
            # skip it then
            check_input = False

        # validate and convert dtype
        X = self._validate_X_predict(X, check_input)

        # -------- This is the only change wrt parent class. TODO maybe rather replace self.tree_ with a proxy ---------
        if smooth_predictions is None:
            # Default: smooth prediction at prediction time if configured as
            # such in the model. Note: if self.use_smoothing == 'installed',
            # models are already smoothed models, so no need to smooth again
            smooth_predictions = (self.use_smoothing == 'on_prediction')
        else:
            # user provided an explicit value for smooth_predictions
            if not smooth_predictions and self.use_smoothing == 'installed':
                raise ValueError("Smoothing has been pre-installed on this "
                                 "tree, it is not anymore possible to make "
                                 "predictions without smoothing")

        if smooth_predictions and smoothing_constant is None:
            # default parameter for smoothing is the one defined in the model
            # (with possible ratio)
            smoothing_constant = self._get_smoothing_constant_to_use(X)

        # Do not use the embedded tree (like in super): it has been pruned but
        # still has constant nodes
        # proba = self.tree_.predict(X)
        if not smooth_predictions:
            proba = predict_from_leaves_no_smoothing(self.tree_,
                                                     self.node_models, X)
        else:
            proba = predict_from_leaves(self, X, smoothing=True,
                                        smoothing_constant=smoothing_constant)
        if len(proba.shape) < 2:
            proba = proba.reshape(-1, 1)
        # ------------------------------------------

        n_samples = X.shape[0]

        # Classification
        if is_classifier(self):
            if self.n_outputs_ == 1:
                return self.classes_.take(np.argmax(proba, axis=1), axis=0)
            else:
                predictions = np.zeros((n_samples, self.n_outputs_))
                for k in range(self.n_outputs_):
                    predictions[:, k] = self.classes_[k].take(
                        np.argmax(proba[:, k], axis=1),
                        axis=0)
                return predictions

        # Regression
        else:
            if self.n_outputs_ == 1:
                return proba[:, 0]
            else:
                return proba[:, :, 0]


class ConstantLeafModel:
    """
    Represents the additional information about a leaf node that is not pruned.
    It contains the error associated with the training samples at this node.

    Note: the constant value is not stored here, as it is already available in
    the sklearn tree struct. So to know the prediction at this node, use
    `tree.value[node_id]`
    """
    __slots__ = ('error',)

    def __init__(self, error):
        self.error = error

    # @property
    # def n_drivers(self):
    #     """
    #     Returns the number of drivers used by this model. (This does not
    #     include the constant driver)
    #     """
    #     return 0

    @property
    def n_params(self):
        """
        Returns the number of parameters used by this model, including the
        constant driver
        """
        return 1

    @staticmethod
    def predict_cstt(tree, node_id, n_samples):
        """
        This static method is a helper to get an array of constant predictions
        associated with this node.

        :param tree:
        :param node_id:
        :param n_samples:
        :return:
        """
        cstt_prediction = tree.value[node_id]

        # cstt_prediction can be multioutput so it is an array. replicate it
        from numpy import matlib  # note: np.matlib not available in 1.10.x
        return matlib.repmat(cstt_prediction, n_samples, 1)


class LinRegNodeModel(DeNormalizableMixIn):
    """
    Represents the additional information about a tree node that contains a
    linear regression model.

    It contains
     - the features used by the model,
     - the scikit learn model object itself,
     - and the error for the training samples that reached this node

    """
    __slots__ = ('features', 'model', 'error')

    def __init__(self, features, model, error):
        self.features = features
        self.model = model
        self.error = error

    def to_text(self, feature_names=None, target_name=None, precision=3,
                line_breaks=False):
        """ Returns a text representation of the linear regression model """
        return linreg_model_to_text(self.model, feature_names=feature_names,
                                    target_name=target_name,
                                    precision=precision,
                                    line_breaks=line_breaks)

    @property
    def n_params(self):
        """
        Returns the number of parameters used by this model, including the
        constant driver
        """
        return len(self.features) + 1

    def reindex_features(self, new_features_lookup_dct):
        """
        Reindexes the required features using the provided lookup dictionary.

        :param new_features_lookup_dct:
        :return:
        """
        self.features = [new_features_lookup_dct[f] for f in self.features]

    def denormalize(self,
                    x_scaler: StandardScaler=None,
                    y_scaler: StandardScaler=None
                    ):
        """
        De-normalizes the linear model.

        :param x_scaler:
        :param y_scaler:
        :return:
        """
        # create a clone of the x scaler with only the used features
        x_scaler = copy(x_scaler)
        if len(self.features) > 0:
            x_scaler.scale_ = x_scaler.scale_[self.features]
            x_scaler.mean_ = x_scaler.mean_[self.features]
        else:
            # in that particular case, the above expression is not working
            x_scaler.scale_ = x_scaler.scale_[0:0]
            x_scaler.mean_ = x_scaler.mean_[0:0]

        # use the denormalize function on the internal model
        # TODO what if internal model is not able ?
        self.model.denormalize(x_scaler, y_scaler)

    def predict(self, X):
        """ Performs a prediction for X, only using the required features """
        if len(self.features) < 1 and isinstance(self.model, LinearRegression):
            # unfortunately LinearRegression models do not like it when no
            # features are needed: their input validation requires at least 1.
            # so we do it ourselves.
            return self.model.intercept_ * np.ones((X.shape[0], 1))
        else:
            return self.model.predict(X[:, self.features])


class LinRegLeafModel(LinRegNodeModel):
    """
    Represents the additional information about a leaf node with a linear
    regression model
    """
    pass


class SplitNodeModel(LinRegNodeModel):
    """
    Represents the additional information about a split node, with a linear
    regression model.
    """
    __slots__ = ('n_params', )

    def __init__(self, n_params, error, features, model):
        self.n_params = n_params
        super(SplitNodeModel, self).__init__(features=features, model=model,
                                             error=error)


PRUNING_MULTIPLIER = 2
# see https://github.com/bnjmn/weka/blob/master/weka/src/main/java/weka/classifiers/trees/m5/RuleNode.java#L124
# TODO check why they use 2 instead of 1 (from the article) ??


def root_mean_squared_error(*args, **kwargs):
    return np.sqrt(mean_squared_error(*args, **kwargs))


def prune_on_min_impurity(tree):
    """
    Edits the given tree so as to prune subbranches that do not respect the min
    impurity criterion

    The paper suggests to do this with 5% but min_impurity_split is not
    available in 0.17 (and criterion is not std but mse so we have to square)

    :param tree:
    :return:
    """
    left_children = tree.children_left
    right_children = tree.children_right
    impurities = tree.impurity

    # The paper suggests to do this with 5% but min_impurity_split is not
    # available in 0.17 (and criterion is not std but mse so we have to square)
    root_impurity = impurities[0]
    impurity_threshold = root_impurity * (0.05 ** 2)

    def stop_on_min_impurity(node_id):
        # note: in the paper that is 0.05 but criterion is on std. Here
        # impurity is mse so a squared equivalent of std.
        left_id = left_children[node_id]
        right_id = right_children[node_id]
        if left_id != _tree.TREE_LEAF:  # a split node
            if impurities[node_id] < impurity_threshold:
                # stop here, this will be a leaf
                removed_ns = prune_children(node_id, tree)
                # print("removed nodes below [%s]: %s" % (node_id, removed_ns))
            else:
                stop_on_min_impurity(left_id)
                stop_on_min_impurity(right_id)

    stop_on_min_impurity(0)


def build_models_and_get_pruning_info(tree, X_train_all, y_train_all,
                                      nodes_to_samples, leaf_model,
                                      node_models, global_std_dev,
                                      use_pruning, node_id=0):
    """


    :param tree: a tree that will be pruned on the way
    :param X:
    :param y:
    :param nodes_to_samples:
    :param leaf_model:
    :param node_models:
    :param global_std_dev:
    :param use_pruning:
    :param node_id:
    :return: a dictionary where the key is the feature index and the value is
        the number of samples where this feature is used
    """

    # Select the Error metric to compute model errors (used to compare node
    # with subtree for pruning)
    # TODO in Weka they use RMSE, but in papers they use MAE. could be a param
    # TODO Shall we go further and store the residuals, or a bunch of metrics?
    err_metric = root_mean_squared_error  # mean_absolute_error

    # Get the samples associated with this node
    samples_at_this_node = get_samples_at_node(node_id, nodes_to_samples)
    n_samples_at_this_node = samples_at_this_node.shape[0]
    y_true_this_node = y_train_all[samples_at_this_node]

    # Is this a leaf node or a split node ?
    left_node = tree.children_left[node_id]  # this way we do not have to query it again in the else.
    if left_node == _tree.TREE_LEAF:
        # --> Current node is a LEAF. See is_leaf(node_id, tree) if you have doubts <--

        # -- create a linear model for this node
        # leaves do not have the right to use any features since they have no
        # subtree: keep the constant prediction
        y_pred_this_node = ConstantLeafModel.predict_cstt(
            tree, node_id, y_true_this_node.shape[0]
        )
        err_this_node = err_metric(y_true_this_node, y_pred_this_node)

        # TODO when use_pruning = False, should we allow all features to be
        #   used instead of having to stick to the M5 rule of "only use a
        #   feature if the subtree includes a split with this feature" ?
        #   OR alternate proposal: should we transform the boolean use_pruning
        #   into a use_pruning_max_nb integer to say for example "only 2 level
        #   of pruning" ?

        # -- store the model information
        node_models[node_id] = ConstantLeafModel(err_this_node)

        # -- return an empty dict for features used
        return dict()  # np.array([], dtype=int)

    else:
        # --> Current node is a SPLIT <--
        right_node = tree.children_right[node_id]

        # (1) prune left and right subtree and get some information
        features_l = build_models_and_get_pruning_info(
            tree, X_train_all, y_train_all, nodes_to_samples, leaf_model,
            node_models, global_std_dev, use_pruning, node_id=left_node
        )
        features_r = build_models_and_get_pruning_info(
            tree, X_train_all, y_train_all, nodes_to_samples, leaf_model,
            node_models, global_std_dev, use_pruning, node_id=right_node
        )

        # (2) select only the samples that reach this node
        X_this_node = X_train_all[samples_at_this_node, :]

        # (3) Create a model for this node
        # -- fit a linear regression model taking into account only the
        # features used in the subtrees + the split one
        # TODO should we normalize=True (variance scaling) or use a whole
        #  pipeline here ? Not sure..
        # estimators = [('scale', StandardScaler()), ('clf', model_type())];
        # pipe = Pipeline(estimators)
        # skmodel_this_node = leaf_model_type(**leaf_model_params)
        skmodel_this_node = clone(leaf_model)

        # -- old - we used to only store the array of features
        # selected_features = np.union1d(features_l, features_r)
        # selected_features = np.union1d(selected_features, tree.feature[node_id])
        # -- new - we also gather the nb samples where this feature is used
        selected_features_dct = features_l
        for feature_id, n_samples_feat_used in features_r.items():
            if feature_id in selected_features_dct.keys():
                selected_features_dct[feature_id] += n_samples_feat_used
            else:
                selected_features_dct[feature_id] = n_samples_feat_used
        selected_features_dct[tree.feature[node_id]] = n_samples_at_this_node
        # -- use only the selected features, in the natural integer order
        selected_features = sorted(selected_features_dct.keys())

        X_train_this_node = X_this_node[:, selected_features]
        skmodel_this_node.fit(X_train_this_node, y_true_this_node)
        # -- predict and compute error
        y_pred_this_node = skmodel_this_node.predict(X_train_this_node)
        err_this_node = err_metric(y_true_this_node, y_pred_this_node)
        # -- create the object
        # TODO the paper suggest to perform recursive feature elimination in
        #  this model until adjusted_err_model is minimal, is it same in Weka ?
        model_this_node = LinRegLeafModel(selected_features, skmodel_this_node,
                                          err_this_node)

        # (4) compute adj error criterion = ERR * (n+v)/(n-v) for both models
        adjusted_err_model = compute_adjusted_error(err_this_node,
                                                    n_samples_at_this_node,
                                                    model_this_node.n_params)

        # (5) predict and compute adj error for the combination of child models
        # -- Note: this is recursive so the leaves may contain linear models
        # already
        y_pred_children = predict_from_leaves_no_smoothing(tree, node_models,
                                                           X_this_node)
        err_children = err_metric(y_true_this_node, y_pred_children)

        # TODO the Weka implem (below) differs from the paper that suggests a
        #  weigthed sum of adjusted errors. Why?
        # the total number of parameters if we do not prune, is the sum of
        # params in each branch PLUS 1 for the split
        n_params_splitmodel = (node_models[left_node].n_params
                               + node_models[right_node].n_params + 1)
        adjusted_err_children = compute_adjusted_error(
            err_children, n_samples_at_this_node, n_params_splitmodel
        )

        # (6) compare and either prune at this node or keep the subtrees
        std_dev_this_node = np.nanstd(y_true_this_node)
        # note: the first criterion is now already checked before that
        # function call, in `prune_on_min_impurity`
        if use_pruning and (
                std_dev_this_node < (global_std_dev * 0.05)
                or (adjusted_err_model <= adjusted_err_children)
                or (adjusted_err_model < (global_std_dev * 0.00001))  # so this means very good R² already
        ):
            # Choose model for this node rather than subtree model
            # -- prune from this node on
            removed_nodes = prune_children(node_id, tree)
            node_models[removed_nodes] = _tree.TREE_UNDEFINED

            # store the model information
            node_models[node_id] = model_this_node

            # update and return the features used
            selected_features_dct = {k: n_samples_at_this_node
                                     for k in selected_features_dct.keys()}
            return selected_features_dct

        else:
            # The subtrees are good or we do not want pruning: keep them.
            # This node will remain a split, and will only contain the digest
            # about the subtree
            node_models[node_id] = SplitNodeModel(
                n_params_splitmodel, err_children, selected_features,
                skmodel_this_node
            )

            # return the features used
            return selected_features_dct


def compute_adjusted_error(err, n_samples, n_parameters):
    """
    TODO why is this different from the term on the left of eq (4) in the
        paper ? (PRUNING_MULTIPLIER ?)
    TODO why 10 in the special case if n samples < n parameters ?

    :param n_samples: the number of samples
    :param n_parameters: number of parameters in the model
    :return:
    """
    #
    if n_samples <= n_parameters:
        # denominator is zero or negative: use a large factor so as to penalize
        # a lot this overly complex model.
        factor = 10.0  # Caution says Yong in his code
    else:
        factor = ((n_samples + PRUNING_MULTIPLIER * n_parameters)
                  / (n_samples - n_parameters))

    return err * factor


def predict_from_leaves_no_smoothing(tree, node_models, X):
    """
    Returns the prediction for all samples in X, based on using the appropriate
    model leaf.

    This function
     - uses tree.apply(X) to know in which tree leaf each sample in X goes
     - then for each of the leaves that are actually touched, uses
       node_models[leaf_id] to predict, for the X reaching that leaf

    This function assumes that node_models contains non-empty LinRegLeafModel
    entries for all leaf nodes that will be reached by the samples X.

    :param node_models:
    :param X: should contain all features !
    :return:
    """
    # **** This does the job, but we have one execution of model.predict() per
    # sample: probably not efficient
    # sample_ids_to_leaf_node_ids = tree.apply(X)
    # model_and_x = np.concatenate((node_models[leaf_node_ids].reshape(-1, 1), X), axis=1)
    # def pred(m_and_x):
    #     return m_and_x[0].model.predict(m_and_x[1:].reshape(1,-1))[0]
    # y_predicted = np.array(list(map(pred, model_and_x)))

    # **** This should be faster because the number of calls to predict() is
    # equal to the number of leaves touched
    sample_ids_to_leaf_node_ids = tree.apply(X)
    y_predicted = -np.ones(sample_ids_to_leaf_node_ids.shape, dtype=DOUBLE)

    # -- find the unique list of leaves touched
    leaf_node_ids, inverse_idx = np.unique(sample_ids_to_leaf_node_ids,
                                           return_inverse=True)

    # -- for each leaf perform the prediction for samples reaching that leaf
    for leaf_node_id in leaf_node_ids:
        # get the indices of the samples reaching that leaf
        sample_indices = np.nonzero(
            sample_ids_to_leaf_node_ids == leaf_node_id)[0]

        # predict
        node_model = node_models[leaf_node_id]
        if isinstance(node_model, LinRegLeafModel):
            y_predicted[sample_indices] = node_model.predict(
                X[sample_indices, :])
        else:
            # isinstance(node_model, ConstantLeafModel)
            y_predicted[sample_indices] = tree.value[leaf_node_id]

    # **** Recursive strategy: not used anymore
    # left_node = tree.children_left[node_id]  # this way we do not have to query it again in the else.
    # if left_node == _tree.TREE_LEAF:
    #     # --> Current node is a LEAF. See is_leaf(node_id, tree) <--
    #     y_predicted = node_models[node_id].model.predict()
    # else:
    #     # --> Current node is a SPLIT <--
    #     right_node = tree.children_right[node_id]
    #
    #
    #     samples_at_this_node = get_samples_at_node(node_id, nodes_to_samples)
    #     y_true_this_node = y_train_all[samples_at_this_node]
    #     # As surprising as it may seem, in numpy [samples_at_this_node, selected_features] does something else.
    #     X_train_this_node = X_train_all[samples_at_this_node, :][:, selected_features]
    #
    #     X_left

    return y_predicted


def predict_from_leaves(m5p, X, smoothing=True, smoothing_constant=15):
    """
    Predicts using the M5P tree, without using the compiled sklearn
    `tree.apply` subroutine.

    The main purpose of this function is to apply smoothing to a M5P model tree
    where smoothing has not been pre-installed on the models. For examples to
    enable a model to be used both without and with smoothing for comparisons
    purposes, or for models whose leaves are not Linear Models and therefore
    for which no pre-installation method exist.

    Note: this method is slower than `predict_from_leaves_no_smoothing` when
    `smoothing=False`.

    :param m5p:
    :param X: input data for which the predictions should be made.
    :return:
    """
    # validate and converts dtype just in case this was directly called
    # e.g. in unit tests
    X = m5p._validate_X_predict(X, check_input=True)

    tree = m5p.tree_
    node_models = m5p.node_models
    nb_samples = X.shape[0]
    y_predicted = -np.ones((nb_samples, 1), dtype=DOUBLE)

    # sample_ids_to_leaf_node_ids = tree.apply(X)
    def smooth_predictions(ancestor_nodes, X_at_node, y_pred_at_node):
        # note: y_predicted_at_node can be a constant
        current_node_model_id = ancestor_nodes[-1]
        for i, parent_model_id in enumerate(reversed(ancestor_nodes[:-1])):
            # warning: this is the nb of TRAINING samples at this node
            node_nb_train_samples = tree.n_node_samples[current_node_model_id]
            parent_model = node_models[parent_model_id]
            parent_predictions = parent_model.predict(X_at_node)
            y_pred_at_node = ((node_nb_train_samples * y_pred_at_node
                               + smoothing_constant * parent_predictions)
                              / (node_nb_train_samples + smoothing_constant))
            current_node_model_id = parent_model_id

        return y_pred_at_node

    def apply_prediction(node_id, ids=None, ancestor_nodes=None):
        first_call = False
        if ids is None:
            ids = slice(None)
            first_call = True
        if smoothing:
            if ancestor_nodes is None:
                ancestor_nodes = [node_id]
            else:
                ancestor_nodes.append(node_id)

        left_id = tree.children_left[node_id]
        if left_id == _tree.TREE_LEAF:
            # ... and tree.children_right[node_id] == _tree.TREE_LEAF
            # LEAF node: predict
            # predict
            node_model = node_models[node_id]
            # assert (ids == (sample_ids_to_leaf_node_ids == node_id)).all()
            if isinstance(node_model, LinRegLeafModel):
                X_at_node = X[ids, :]
                predictions = node_model.predict(X_at_node)
            else:
                # isinstance(node_model, ConstantLeafModel)
                predictions = tree.value[node_id]
                if smoothing:
                    X_at_node = X[ids, :]

            # finally apply smoothing
            if smoothing:
                y_predicted[ids] = smooth_predictions(
                    ancestor_nodes, X_at_node, predictions)
            else:
                y_predicted[ids] = predictions
        else:
            right_id = tree.children_right[node_id]
            # non-leaf node: split samples and recurse
            left_group = np.zeros(nb_samples, dtype=bool)
            left_group[ids] = (X[ids, tree.feature[node_id]]
                               <= tree.threshold[node_id])
            right_group = ((~left_group)
                           if first_call
                           else (ids & (~left_group)))

            # important: copy ancestor_nodes BEFORE calling anything, otherwise
            # it will be modified
            apply_prediction(left_id, ids=left_group,
                             ancestor_nodes=(ancestor_nodes.copy()
                                             if ancestor_nodes is not None
                                             else None))
            apply_prediction(right_id, ids=right_group,
                             ancestor_nodes=ancestor_nodes)

    # recurse to fill all predictions
    apply_prediction(0)

    return y_predicted


def prune_children(node_id, tree):
    """
    Prunes the children of node_id in the given `tree`.

    Inspired by https://github.com/shenwanxiang/sklearn-post-prune-tree/blob/master/tree_prune.py#L122

    IMPORTANT this relies on the fact that `children_left` and `children_right`
    are modificable (and are the only things we need to modify to fix the
    tree). It seems to be the case as of now.

    :return: a list of removed nodes
    """

    def _prune_below(_node_id):
        left_child = tree.children_left[_node_id]
        right_child = tree.children_right[_node_id]
        if left_child == _tree.TREE_LEAF:
            # _node_id is a leaf: left_ & right_child say "leaf".Nothing to do
            return []
        else:
            # Make sure everything is pruned below: children should be leaves
            removed_l = _prune_below(left_child)
            removed_r = _prune_below(right_child)

            # -- first declare that they are not leaves anymore but they do not exist at all
            assert tree.children_left[left_child] == _tree.TREE_LEAF  # left_child was a leaf
            assert tree.children_right[left_child] == _tree.TREE_LEAF  # left_child was a leaf
            tree.children_left[left_child] = _tree.TREE_UNDEFINED
            tree.children_right[left_child] = _tree.TREE_UNDEFINED
            assert tree.children_left[right_child] == _tree.TREE_LEAF  # right_child was a leaf
            assert tree.children_right[right_child] == _tree.TREE_LEAF  # right_child was a leaf
            tree.children_left[right_child] = _tree.TREE_UNDEFINED
            tree.children_right[right_child] = _tree.TREE_UNDEFINED

            # -- then declare that current node is a leaf
            tree.children_left[_node_id] = _tree.TREE_LEAF
            tree.children_right[_node_id] = _tree.TREE_LEAF

            # return the list of removed nodes
            return removed_l + removed_r + [left_child, right_child]

    # true_node_count = tree.node_count - sum(tree.children_left == _tree.TREE_UNDEFINED)
    # do we need to change node_count here ? not sure, it may have side effects... and we clean all of this later so ok.
    # tree.node_count -= 2*len(nodes_to_remove)

    return _prune_below(node_id)


def is_leaf(node_id, tree):
    """
    Returns true if node with id `node_id` is a leaf in tree `tree`.
    Is is not actually used in this file because we always need the left child
    node id in our code.

    But it is kept here as an easy way to remember how it works.

    :param node_id:
    :param tree:
    :return:
    """
    if node_id == _tree.TREE_LEAF or node_id == _tree.TREE_UNDEFINED:
        raise ValueError("Invalid node_id %s" % node_id)

    return tree.children_left[node_id] == _tree.TREE_LEAF


def get_samples_at_node(node_id, nodes_to_samples):
    """
    Returns an array containing the ids of the samples for node `node_id`.
    This method requires the user to
     - first compute the decision path for the sample matrix X
     - then convert it to a csc

    >>> samples_to_nodes = estimator.decision_path(X)  # returns a Scipy compressed sparse row matrix (CSR)
    >>> nodes_to_samples = samples_to_nodes.tocsc()    # we need the column equivalent (CSC)
    >>> samples = get_samples_at_node(node_id, nodes_to_samples)

    :param node_id:
    :param nodes_to_samples_csc:
    :return:
    """
    return nodes_to_samples.indices[nodes_to_samples.indptr[node_id]:
                                    nodes_to_samples.indptr[node_id + 1]]


# def prune_nodes_in_tree(nodes_to_remove, tree):
#     """
#     Transforms all nodes in nodes_to_remove into leaf nodes, and returns the
#     resulting tree (work is done on a copy)
#
#     :param nodes_to_remove:
#     :param tree:
#     :return:
#     """
#     # Create a copy of the inner tree
#     # tree._copy is gone, but this does the same thing
#     out_tree = _tree.Tree(*tree.__reduce__()[1])
#     out_tree.__setstate__(tree.__getstate__().copy())
#
#     for node in nodes_to_remove:
#         # TO DO: Add a Tree method to remove a branch of a tree
#         out_tree.children_left[out_tree.children_left[node]] = _tree.TREE_UNDEFINED
#         out_tree.children_right[out_tree.children_left[node]] = _tree.TREE_UNDEFINED
#         out_tree.children_left[out_tree.children_right[node]] = _tree.TREE_UNDEFINED
#         out_tree.children_right[out_tree.children_right[node]] = _tree.TREE_UNDEFINED
#         out_tree.children_left[node] = _tree.TREE_LEAF
#         out_tree.children_right[node] = _tree.TREE_LEAF
#
#     # FIX ME: currently should not change node_count, after deletion
#     # this is not number of nodes in the tree
#     # out_tree.node_count -= 2*len(nodes_to_remove)
#
#     return out_tree


class M5Prime(M5Base, RegressorMixin):
    """An M5' (M five prime) model tree regressor.

    The original M5 algorithm was invented by R. Quinlan. Y. Wang made
    improvements and named the resulting algorithm M5 Prime.

    This implementation was inspired by Weka (https://github.com/bnjmn/weka)
    M5Prime class, from Mark Hall

    Read more in the :ref:`User Guide <tree>`.

    See also
    --------
    M5Base

    References
    ----------

    .. [1] Ross J. Quinlan, "Learning with Continuous Classes", 5th Australian
           Joint Conference on Artificial Intelligence, pp343-348, 1992.
    .. [2] Y. Wang and I. H. Witten, "Induction of model trees for predicting
           continuous classes", Poster papers of the 9th European Conference
           on Machine Learning, 1997.

    """
    # def __init__(self, save_instances: bool=False):
    #     super(M5Prime, self).__init__(generate_rule=False,
    #                                   save_instances=save_instances)


def _empty_contig_ar(shape, dtype):
    return np.ascontiguousarray(np.empty(shape, dtype=dtype))
