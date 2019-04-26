"""
This module gathers methods built arount the M5 (model trees) methodology.
"""

# Authors: Sylvain Marie <sylvain.marie@schneider-electric.com>
#
# License: BSD 3 clause
from copy import copy

import numpy as np
from scipy.sparse import issparse

from ..metrics import mean_squared_error
from ..base import RegressorMixin, clone
from ..base import is_classifier
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..linear_model import LinearRegression

from ._tree import Tree, TREE_UNDEFINED, TREE_LEAF
from .tree import BaseDecisionTree, DTYPE, DOUBLE


__all__ = ["M5Base",
           "M5Prime"]


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

    Pruning and Smoothing (TODO) can be activated and deactivated on top of the base model.
    TODO check if 'rules' should be supported too
    
    Inspired by Weka (https://github.com/bnjmn/weka) M5Base class, from Mark Hall
    """

    def __init__(self,
                 criterion='mse',
                 # According to M5' paper, this should lead to similar results than (std/rmse), that is not implemented in sklearn. TODO check FriedmanMSE
                 splitter='best',  # take the feature with the best split
                 max_depth=None,  # no maximum depth limit
                 min_samples_split=4,  # in weka this parameter is named "M"
                 min_samples_leaf=2,  # from the M5' article : otherwise (n+v)/(n-v) is infinite
                 min_weight_fraction_leaf=0.,  # TODO this would actually maybe be better than min_sample_leaf ?
                 max_features=None,  # no feature reduction: take all features
                 max_leaf_nodes=None,  # no limitation in number of leaves
                 min_impurity_decrease=0.,  #
                 # min_impurity_split_as_initial_ratio = 0.05,  # TODO The paper suggests to use 5% of all STDR
                 min_impurity_split=None,  # TODO this is deprecated and will be removed in version 0.21
                 random_state=None,
                 class_weight=None,
                 presort=False,
                 leaf_model=None,  # the regression model used in the leaves (it will be cloned for each leaf)
                 ):

        # ------ WEKA VERSION FOR REFERENCE ---
        # From https://github.com/bnjmn/weka/blob/master/weka/src/main/java/weka/classifiers/trees/m5/RuleNode.java
        # it is like the M5' paper: a minimum number + a rule on standard deviation ratio
        # if ((m_numInstances < m_splitNum)
        #         | | (Rule.stdDev(m_classIndex, m_instances) < (m_globalDeviation * m_devFraction))) {
        # m_isLeaf = true;
        # TODO check their impurity criterion for splitting

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

        # TODO the paper suggests to do this with 5% but min_impurity_split is deprecated in 0.21 ... worth it ?
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
            presort=presort)

        # warning : if the field names are different from constructor params,
        # then clone(self) will not work.
        if leaf_model is None:
            leaf_model = LinearRegression()
        self.leaf_model = leaf_model

    # TODO where shall we put our tree pretty-printing function in sklearn?
    # def as_pretty_text(self, **kwargs):
    #     """
    #     Returns a multi-line representation of this decision tree, using
    #     `tree_to_text`.
    #
    #     :return: a multi-line string representing this decision tree
    #     """
    #     from sepycom.sklearn_utils.tree_to_text import tree_to_text
    #     return tree_to_text(self, out_file=None, **kwargs)

    def fit(self, X, y, sample_weight=None, check_input=True,
            X_idx_sorted=None):

        # (1) Build the initial tree as usual
        super(M5Base, self).fit(
            X, y,
            sample_weight=sample_weight,
            check_input=check_input,
            X_idx_sorted=X_idx_sorted)

        # ** Debug **
        # from skl_mfivep.tree_to_text import tree_to_text
        # print(tree_to_text(self, out_file=None, node_ids=True))

        # (2) Now prune the tree (M5P step) and replace pruned branches with linear models
        # -- unfortunately we have to re-do this input validation step, because it also converts the input to float32.
        if check_input:
            X = check_array(X, dtype=DTYPE, accept_sparse="csc")
            if issparse(X):
                X.sort_indices()

                if X.indices.dtype != np.intc or X.indptr.dtype != np.intc:
                    raise ValueError("No support for np.int64 index based "
                                     "sparse matrices")

        # -- initialise the structure that will contain the leaves and nodes models
        self.node_models = np.empty((self.tree_.node_count,), dtype=object)

        # -- The pruning procedure requires to know the samples that reached each node.
        # From http://scikit-learn.org/stable/auto_examples/tree/plot_unveil_tree_structure.html
        # Retrieve the decision path of each sample.
        samples_to_nodes = self.decision_path(X)
        # * row i is to see the nodes (non-empty j) in which sample i appears. sparse row-first (CSR) format is OK
        # * column j is to see the samples (non-empty i) that fall into that node. To do that, we need to make it CSC
        nodes_to_samples = samples_to_nodes.tocsc()

        # -- execute the pruning
        # self.features_usage is a dict feature_idx -> nb times used, only for used features.
        self.features_usage = build_models_and_get_pruning_info(self.tree_, X, y, nodes_to_samples, self.leaf_model,
                                                                self.node_models)

        # ** Debug **
        # from skl_mfivep.tree_to_text import tree_to_text
        # a = tree_to_text(self, out_file=None, node_ids=False)  # use node_ids=False for comparison
        # print(tree_to_text(self, out_file=None, node_ids=True))

        # -- cleanup to compress inner structures: only keep the non-pruned ones
        self._cleanup_tree()

        # ** Debug **
        # from skl_mfivep.tree_to_text import tree_to_text
        # b = tree_to_text(self, out_file=None, node_ids=False)  # use node_ids=False for comparison
        # print(tree_to_text(self, out_file=None, node_ids=True))
        # assert a == b

        return self

    def _cleanup_tree(self):
        """
        Reduces the size of this object by removing from internal structures all items that are not used any more.
        So, all leaves that have been pruned.
        :return:
        """
        old_tree = self.tree_
        old_node_models = self.node_models

        # Get all information to create a copy of the inner tree. Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # cleanup: only keep the nodes that are not undefined.
        valid_nodes_indices = dct["nodes"]['left_child'] != TREE_UNDEFINED
        new_node_count = sum(valid_nodes_indices)

        # create empty new structures
        n_shape = (new_node_count, *dct["nodes"].shape[1:])
        new_nodes = np.ascontiguousarray(np.empty(n_shape, dtype=dct["nodes"].dtype))
        v_shape = (new_node_count, *dct["values"].shape[1:])
        new_values = np.ascontiguousarray(np.empty(v_shape, dtype=dct["values"].dtype))
        m_shape = (new_node_count, *old_node_models.shape[1:])
        new_node_models = np.ascontiguousarray(np.empty(m_shape, dtype=old_node_models.dtype))

        # Fill the structures while reindexing the tree and remembering the depth
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
            new_node_models[new_node_id] = copy(old_node_models[old_node_id])

            if left_id == TREE_LEAF:
                # ... and right_id == TREE_LEAF
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

        # update self fields
        self.tree_ = new_tree
        self.node_models = new_node_models

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

        # Get all information to create a copy of the inner tree. Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # create empty new structures for nodes and models
        n_shape = dct["nodes"].shape
        new_nodes = np.ascontiguousarray(np.empty(n_shape, dtype=dct["nodes"].dtype))
        v_shape = dct["values"].shape
        new_values = np.ascontiguousarray(np.empty(v_shape, dtype=dct["values"].dtype))
        m_shape = old_node_models.shape
        new_node_models = np.ascontiguousarray(np.empty(m_shape, dtype=old_node_models.dtype))

        def _denormalize(node_id, x_scaler, y_scaler):
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
            new_model = copy(old_node_models[node_id])
            new_node_models[node_id] = new_model

            # Create the new value by de-scaling y
            if y_scaler is not None:
                new_values[node_id] = y_scaler.inverse_transform(dct["values"][node_id])
            else:
                new_values[node_id] = dct["values"][node_id]

            if left_id == TREE_LEAF:
                # ... and right_id == TREE_LEAF
                # => node_id is a leaf
                if isinstance(new_model, ConstantLeafModel):
                    # nothing to do: we already re-scaled the value
                    return
                elif isinstance(new_model, LinRegLeafModel):
                    # denormalize model
                    new_model.denormalize(x_scaler, y_scaler)
                    return
                else:
                    raise TypeError("Internal error - Leafs can only be constant or linear regression")
            else:

                _denormalize(left_id, x_scaler, y_scaler)
                _denormalize(right_id, x_scaler, y_scaler)

                # denormalize the split value if needed
                if x_scaler is not None:
                    split_feature = old_node['feature']
                    # -- put the split value in a vector where it has the position 'feature'
                    old_threshold_and_zeros = np.zeros((self.n_features_, ), dtype=dct["nodes"]['threshold'].dtype)
                    old_threshold_and_zeros[split_feature] = old_node['threshold']
                    # -- denormalize the vector and retrieve the value 'feature'
                    new_nodes[node_id]['threshold'] = x_scaler.inverse_transform(old_threshold_and_zeros)[split_feature]
                else:
                    new_nodes[node_id]['threshold'] = old_node['threshold']

                if isinstance(new_model, SplitNodeModel):
                    # TODO when split node models will have linear regression models too, update them here
                    pass
                else:
                    raise TypeError("Internal error: all intermediate nodes should be SplitNodeModel")

                return

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
        Compresses the model and returns a vector of required feature indices. This model input will then be
        X[:, features] instead of X.

        :return:
        """

        used_features = sorted(self.features_usage.keys())
        new_features_lookup_dct = {old_feature_idx: i for i, old_feature_idx in enumerate(used_features)}

        if used_features == list(range(self.n_features_)):
            # NOTHING TO DO: we need all features
            return used_features

        # Otherwise, We can compress. For this we have to create a copy of the tree because we will change its internals
        old_tree = self.tree_
        old_node_models = self.node_models

        # Get all information to create a copy of the inner tree. Note: _tree.copy() is gone so we use the pickle way
        # --- Base info: nb features, nb outputs, output classes
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L631
        # [1] = (self.n_features, sizet_ptr_to_ndarray(self.n_classes, self.n_outputs), self.n_outputs)
        # So these remain, we are just interested in changing the node-related arrays
        new_tree = Tree(*old_tree.__reduce__()[1])

        # --- Node info
        # see https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/tree/_tree.pyx#L637
        dct = old_tree.__getstate__().copy()

        # create empty new structures for nodes and models
        n_shape = dct["nodes"].shape
        new_nodes = np.ascontiguousarray(np.empty(n_shape, dtype=dct["nodes"].dtype))
        m_shape = old_node_models.shape
        new_node_models = np.ascontiguousarray(np.empty(m_shape, dtype=old_node_models.dtype))

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
            new_model = copy(old_node_models[node_id])
            new_node_models[node_id] = new_model

            if left_id == TREE_LEAF:
                # ... and right_id == TREE_LEAF
                # => node_id is a leaf
                if isinstance(new_model, ConstantLeafModel):
                    # no features used
                    return
                elif isinstance(new_model, LinRegLeafModel):
                    # compress that model
                    new_model.reindex_features(new_features_lookup_dct)
                    return
                else:
                    raise TypeError("Internal error - Leafs can only be constant or linear regression")
            else:

                _compress_features(left_id)
                _compress_features(right_id)

                # store the new split feature index in the node
                new_nodes[node_id]['feature'] = new_features_lookup_dct[old_node['feature']]

                if isinstance(new_model, SplitNodeModel):
                    # TODO when split node models will have linear regression models too, update them here
                    pass
                else:
                    raise TypeError("Internal error: all intermediate nodes should be SplitNodeModel")

                return

        _compress_features(0)

        # complete definition of the new tree
        dct["nodes"] = new_nodes
        new_tree.__setstate__(dct)

        # update the self fields
        self.features_usage = {new_features_lookup_dct[k]: v for k, v in self.features_usage.items()}
        self.tree_ = new_tree
        self.node_models = new_node_models
        self.n_features_ = len(self.features_usage)

        # return the vector of used features
        return used_features

    @property
    def feature_importances_(self):
        """Return the feature importances (the higher, the more important the feature).

        Returns
        -------
        feature_importances_ : array, shape = [n_features]
        """
        check_is_fitted(self, 'tree_')

        # TODO when m5tree inherits from tree, compute_feature_importances below
        features = np.array([self.features_usage[k] for k in sorted(self.features_usage.keys())], dtype=int)

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

    def predict(self, X, check_input=True):
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

        Returns
        -------
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes, or the predict values.
        """
        check_is_fitted(self, 'tree_')
        # If the model is just a constant node, do not check the input, only its shape.
        if self.n_features_ == 0:
            # perform the input checking manually to set ensure_min_features=0
            X = check_array(X, dtype=DTYPE, accept_sparse="csr", ensure_min_features=0)
            if issparse(X) and (X.indices.dtype != np.intc or
                                X.indptr.dtype != np.intc):
                raise ValueError("No support for np.int64 index based "
                                 "sparse matrices")
            # skip it then
            check_input = False

        X = self._validate_X_predict(X, check_input)

        # -------- This is the only change wrt parent class. TODO maybe rather replace self.tree_ with a proxy ---------
        # Do not use the embedded tree (like in super): it has been pruned but still has constant nodes
        # proba = self.tree_.predict(X)
        proba = predict_from_leafs_no_smoothing(self.tree_, self.node_models, X)
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

    To know the prediction at this node, use tree.value[node_id]
    """
    __slots__ = ('error',)

    def __init__(self, error):
        self.error = error

    # @property
    # def n_drivers(self):
    #     """ Returns the number of drivers used by this model. (This does not include the constant driver) """
    #     return 0

    @property
    def n_params(self):
        """ Returns the number of parameters used by this model, including the constant driver """
        return 1

    @staticmethod
    def predict_cstt(tree, node_id, n_samples):
        """
        This static method is a helper to get an array of constant predictions associated with this node.

        :param tree:
        :param node_id:
        :param n_samples:
        :return:
        """
        cstt_prediction = tree.value[node_id]

        # cstt_prediction can be multioutput so it is an array. replicate it
        from numpy import matlib  # note: np.matlib is not available in 1.10.x series
        return matlib.repmat(cstt_prediction, n_samples, 1)


class LinRegLeafModel(object):  # TODO DeNormalizable API
    """
    Represents the additional information about a tree leaf that contains a linear regression model.
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

    # def to_text(self, feature_names=None, target_name=None, precision=3, line_breaks=False):
    #     """ Returns a text representation of the linear regression model """
    #     return linreg_model_to_text(self.model, feature_names=feature_names, target_name=target_name,
    #                                 precision=precision, line_breaks=line_breaks)

    @property
    def n_params(self):
        """ Returns the number of parameters used by this model, including the constant driver """
        return len(self.features) + 1

    def reindex_features(self, new_features_lookup_dct):
        """
        Reindexes the required features using the provided lookup dictionary.

        :param new_features_lookup_dct:
        :return:
        """
        self.features = [new_features_lookup_dct[f] for f in self.features]

    # TODO is there a denormalization API now in scikit learn ?
    # def denormalize(self, x_scaler: StandardScaler=None, y_scaler: StandardScaler=None):
    #     """
    #     De-normalizes the linear model.
    #
    #     :param x_scaler:
    #     :param y_scaler:
    #     :return:
    #     """
    #     # create a clone of the x scaler with only the used features
    #     x_scaler = copy(x_scaler)
    #     x_scaler.scale_ = x_scaler.scale_[self.features]
    #     x_scaler.mean_ = x_scaler.mean_[self.features]
    #
    #     # use the denormalize function on the internal model
    #     # TODO what if internal model is not able ?
    #     self.model.denormalize(x_scaler, y_scaler)

    def predict(self, X):
        """ Performs a prediction for X, only using the required features """
        return self.model.predict(X[:, self.features])


class SplitNodeModel:
    __slots__ = ('n_params', 'error')

    def __init__(self, n_params, error):
        self.n_params = n_params
        self.error = error


PRUNING_MULTIPLIER = 2
# see https://github.com/bnjmn/weka/blob/master/weka/src/main/java/weka/classifiers/trees/m5/RuleNode.java#L124
# TODO check why they use 2 instead of 1 (from the article) ??


def root_mean_squared_error(*args, **kwargs):
    return np.sqrt(mean_squared_error(*args, **kwargs))


def build_models_and_get_pruning_info(tree, X_train_all, y_train_all, nodes_to_samples, leaf_model,
                                      node_models, node_id=0):
    """


    :param tree: a tree that will be pruned on the way
    :param X:
    :param y:
    :param nodes_to_samples:
    :param leaf_model:
    :param node_models:
    :param node_id:
    :return: a dictionary where the key is the feature index and the value is the number of samples where this feature
        is used
    """

    # Select the Error metric to compute model errors (used to compare node with subtree for pruning)
    # TODO in Weka they use RMSE, but in papers they use MAE. this could be a parameter.
    # TODO Shall we go further and store the residuals, or a bunch of metrics? not sure
    err_metric = root_mean_squared_error  # mean_absolute_error

    # Get the samples associated with this node
    samples_at_this_node = get_samples_at_node(node_id, nodes_to_samples)
    n_samples_at_this_node = samples_at_this_node.shape[0]
    y_true_this_node = y_train_all[samples_at_this_node]

    # Is this a leaf node or a split node ?
    left_node = tree.children_left[node_id]  # this way we do not have to query it again in the else.
    if left_node == TREE_LEAF:
        # --> Current node is a LEAF. See is_leaf(node_id, tree) if you have doubts <--

        # -- create a linear model for this node
        # leaves do not have the right to use any features since they have no subtree: keep the constant prediction
        y_pred_this_node = ConstantLeafModel.predict_cstt(tree, node_id, y_true_this_node.shape[0])
        err_this_node = err_metric(y_true_this_node, y_pred_this_node)

        # -- store the model information
        node_models[node_id] = ConstantLeafModel(err_this_node)

        # -- return an empty dict for features used
        return dict()  # np.array([], dtype=int)

    else:
        # --> Current node is a SPLIT <--
        right_node = tree.children_right[node_id]

        # (1) prune left and right subtree and get some information
        features_l = build_models_and_get_pruning_info(tree, X_train_all, y_train_all, nodes_to_samples,
                                                       leaf_model, node_models, node_id=left_node)
        features_r = build_models_and_get_pruning_info(tree, X_train_all, y_train_all, nodes_to_samples,
                                                       leaf_model, node_models, node_id=right_node)

        # (2) select only the samples that reach this node
        X_this_node = X_train_all[samples_at_this_node, :]

        # (3) Create a model for this node
        # -- fit a linear regression model taking into account only the features used in the subtrees + the split one
        # TODO should we normalize=True (variance scaling) or use a whole pipeline here ? Not sure..
        # estimators = [('scale', StandardScaler()), ('clf', model_type())]; pipe = Pipeline(estimators)
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
        model_this_node = LinRegLeafModel(selected_features, skmodel_this_node, err_this_node)

        # (4) predict and compute adjusted error for the combination of child nodes models
        # -- Note: the subtrees are already pruned so the leafs may contain linear models already
        y_pred_children = predict_from_leafs_no_smoothing(tree, node_models, X_this_node)
        err_children = err_metric(y_true_this_node, y_pred_children)

        # (5) compute adjusted error criterion = ERR * (n+v)/(n-v) for both models
        adjusted_err_model = compute_adjusted_error(err_this_node, n_samples_at_this_node, model_this_node.n_params)

        # TODO the Weka implem below differs from the paper that suggests a weigthed sum of adjusted errors. Why?
        # the total number of parameters if we do not prune, is the sum of params in each branch PLUS 1 for the split
        n_params_splitmodel = node_models[left_node].n_params + node_models[right_node].n_params + 1
        adjusted_err_children = compute_adjusted_error(err_children, n_samples_at_this_node, n_params_splitmodel)

        # (6) compare and either prune at this node or keep the subtrees
        m_globalDeviation = 0  # TODO where do we get it from ?
        if adjusted_err_model <= adjusted_err_children \
                or adjusted_err_model < (m_globalDeviation * 0.00001):  # so this means very good R² already
            # Choose model for this node rather than subtree model
            # -- prune from this node on
            removed_nodes = prune_children(node_id, tree)
            node_models[removed_nodes] = TREE_UNDEFINED

            # store the model information
            node_models[node_id] = model_this_node

            # update and return the features used
            selected_features_dct = {k: n_samples_at_this_node for k in selected_features_dct.keys()}
            return selected_features_dct

        else:
            # The subtrees are good, keep them.
            # This node will remain a split, and will only contain the digest about the subtree
            node_models[node_id] = SplitNodeModel(n_params_splitmodel, err_children)

            # return the features used
            return selected_features_dct


def compute_adjusted_error(err, n_samples, n_parameters):
    """
    TODO why is this different from the term on the left of eq (4) in the paper ? (PRUNING_MULTIPLIER ?)
    TODO why 10 in the special case if n samples < n parameters ?

    :param n_samples: the number of samples
    :param n_parameters: number of parameters in the model
    :return:
    """
    #
    if n_samples <= n_parameters:
        # denominator might be zero or negative: use a large factor so as to penalize a lot this overly complex model.
        factor = 10.0  # Caution says Yong in his code
    else:
        factor = (n_samples + PRUNING_MULTIPLIER * n_parameters) / (n_samples - n_parameters)

    return err * factor


def predict_from_leafs_no_smoothing(tree, node_models, X):
    """
    Returns the prediction for all samples in X, based on using the appropriate model leaf.

    This function
     - uses tree.apply(X) to know in which tree leaf each sample in X goes
     - then for each of the leaves that are actually touched, uses node_models[leaf_id] to predict, for the X
     reaching that leaf

    This function assumes that node_models contains non-empty LinRegLeafModel entries for all leaf nodes that will be
    reached by the samples X.

    :param node_models:
    :param X: should contain all features !
    :return:
    """
    # **** This does the job, but we have one execution of model.predict() per sample: probably not efficient
    # sample_ids_to_leaf_node_ids = tree.apply(X)
    # model_and_x = np.concatenate((node_models[leaf_node_ids].reshape(-1, 1), X), axis=1)
    # def pred(m_and_x):
    #     return m_and_x[0].model.predict(m_and_x[1:].reshape(1,-1))[0]
    # y_predicted = np.array(list(map(pred, model_and_x)))

    # **** This should be faster because the number of calls to predict() is equal to the number of leaves touched
    sample_ids_to_leaf_node_ids = tree.apply(X)
    y_predicted = -np.ones(sample_ids_to_leaf_node_ids.shape, dtype=DOUBLE)

    # -- find the unique list of leaves touched
    leaf_node_ids, inverse_idx = np.unique(sample_ids_to_leaf_node_ids, return_inverse=True)

    # -- for each leaf, perform the prediction for the samples reaching that leaf
    for leaf_node_id in leaf_node_ids:
        # get the indices of the samples reaching that leaf
        sample_indices = np.nonzero(sample_ids_to_leaf_node_ids == leaf_node_id)

        # predict
        node_model = node_models[leaf_node_id]
        if isinstance(node_model, LinRegLeafModel):
            y_predicted[sample_indices] = node_model.predict(X[sample_indices])
        else:
            # isinstance(node_model, ConstantLeafModel)
            y_predicted[sample_indices] = tree.value[leaf_node_id]

    # **** Recursive strategy: not used anymore
    # left_node = tree.children_left[node_id]  # this way we do not have to query it again in the else.
    # if left_node == TREE_LEAF:
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


def prune_children(node_id, tree):
    """
    Prunes the children of node_id in the given `tree`.

    Inspired by https://github.com/shenwanxiang/sklearn-post-prune-tree/blob/master/tree_prune.py#L122

    :return: a list of removed nodes
    """

    def _prune_below(_node_id):
        left_child = tree.children_left[_node_id]
        right_child = tree.children_right[_node_id]
        if left_child == TREE_LEAF:
            # _node_id is a leaf: left_child and right_child say "leaf". Nothing to do
            return []
        else:
            # Make sure that everything is pruned below: the children should be leaves
            removed_l = _prune_below(left_child)
            removed_r = _prune_below(right_child)

            # -- first declare that they are not leaves anymore but they do not exist at all
            assert tree.children_left[left_child] == TREE_LEAF  # left_child was a leaf
            assert tree.children_right[left_child] == TREE_LEAF  # left_child was a leaf
            tree.children_left[left_child] = TREE_UNDEFINED
            tree.children_right[left_child] = TREE_UNDEFINED
            assert tree.children_left[right_child] == TREE_LEAF  # right_child was a leaf
            assert tree.children_right[right_child] == TREE_LEAF  # right_child was a leaf
            tree.children_left[right_child] = TREE_UNDEFINED
            tree.children_right[right_child] = TREE_UNDEFINED

            # -- then declare that current node is a leaf
            tree.children_left[_node_id] = TREE_LEAF
            tree.children_right[_node_id] = TREE_LEAF

            # return the list of removed nodes
            return removed_l + removed_r + [left_child, right_child]

    # true_node_count = tree.node_count - sum(tree.children_left == TREE_UNDEFINED)
    # TODO do we need to change node_count here ? not sure, it may have side effects...
    # tree.node_count -= 2*len(nodes_to_remove)

    return _prune_below(node_id)


def is_leaf(node_id, tree):
    """
    Returns true if node with id `node_id` is a leaf in tree `tree`.
    Is is not actually used in this file because we always need the left child node id in our code.

    But it is kept here as an easy way to remember how it works.

    :param node_id:
    :param tree:
    :return:
    """
    if node_id == TREE_LEAF or node_id == TREE_UNDEFINED:
        raise ValueError("Invalid node_id %s" % node_id)

    return tree.children_left[node_id] == TREE_LEAF


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
#     Transforms all nodes in nodes_to_remove into leaf nodes, and returns the resulting tree (work is done on a copy)
#
#     :param nodes_to_remove:
#     :param tree:
#     :return:
#     """
#     # Create a copy of the inner tree
#     # tree._copy is gone, but this does the same thing
#     out_tree = Tree(*tree.__reduce__()[1])
#     out_tree.__setstate__(tree.__getstate__().copy())
#
#     for node in nodes_to_remove:
#         # TODO: Add a Tree method to remove a branch of a tree
#         out_tree.children_left[out_tree.children_left[node]] = TREE_UNDEFINED
#         out_tree.children_right[out_tree.children_left[node]] = TREE_UNDEFINED
#         out_tree.children_left[out_tree.children_right[node]] = TREE_UNDEFINED
#         out_tree.children_right[out_tree.children_right[node]] = TREE_UNDEFINED
#         out_tree.children_left[node] = TREE_LEAF
#         out_tree.children_right[node] = TREE_LEAF
#
#     # FIXME: currently should not change node_count, after deletion
#     # this is not number of nodes in the tree
#     # out_tree.node_count -= 2*len(nodes_to_remove)
#
#     return out_tree


class M5Prime(M5Base, RegressorMixin):
    """An M5' (M five prime) model tree regressor.

    The original M5 algorithm was invented by R. Quinlan. Y. Wang made improvements and named the resulting algorithm
    M5 Prime.

    This implementation was inspired by Weka (https://github.com/bnjmn/weka) M5Prime class, from Mark Hall

    Read more in the :ref:`User Guide <tree>`.

    See also
    --------
    M5Base

    References
    ----------

    .. [1] Ross J. Quinlan, "Learning with Continuous Classes", 5th Australian Joint Conference on Artificial
           Intelligence, pp343-348, 1992.
    .. [2] Y. Wang and I. H. Witten, "Induction of model trees for predicting continuous classes",
           Poster papers of the 9th European Conference on Machine Learning, 1997.

    """
    # def __init__(self, save_instances: bool=False):
    #     super(M5Prime, self).__init__(generate_rule=False,
    #                                   save_instances=save_instances)

