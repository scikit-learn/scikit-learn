#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

from libc.stdlib cimport malloc, free, realloc
from libc.math cimport NAN, INFINITY, isnan

from ._tree cimport Node

from ._splitter cimport Splitter
from ._splitter cimport splitter_init
from ._splitter cimport splitter_reset
from ._splitter cimport splitter_node_evaluate_split
from ._splitter cimport splitter_expand
from ._splitter cimport splitter_set_nid
from ._splitter cimport splitter_copy_to

from ._split_record cimport SplitRecord
from ._split_record cimport split_record_reset

from ._stats_node cimport StatsNode
from ._stats_node cimport stats_node_reset
from ._stats_node cimport stats_node_clear

from ._criterion cimport _impurity_mse

# from sklearn.utils import check_random_state


cdef:
    int TREE_UNDEFINED = -2
    int FEAT_UNKNOWN = -3
    int TREE_LEAF = -1
    bint TREE_NOT_LEAF = 0


cdef void weighted_sum_y(double[::1] y, double[::1] sample_weight,
                         double* p_sum_y, double* p_sum_sq_y):
    cdef int i
    p_sum_y[0] = 0.0
    p_sum_sq_y[0] = 0.0
    for i in range(y.shape[0]):
        p_sum_y[0] += y[i] * sample_weight[i]
        p_sum_sq_y[0] += y[i] ** 2 * sample_weight[i]


cdef class ExactTreeBuilder(TreeBuilder):

    def __cinit__(self, int min_samples_split,
                  int min_samples_leaf, double min_weight_leaf,
                  double min_impurity_split,
                  int max_depth, int max_features):
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf
        self.min_impurity_split = min_impurity_split
        self.max_depth = max_depth
        self.max_features = max_features

    cpdef build(self, Tree tree, float[:, ::1] X, int[::1, :] X_idx_sorted,
                double[::1] y, double[::1] sample_weight,
                double sum_total_weighted_samples):

        # FIXME: don't set the random state here
        # rng = check_random_state(0)

        # we don't do any checking for the moment

        cdef:
            int n_samples = X.shape[0]
            int n_features = X.shape[1]
            int init_capacity
            int n_splitter
            int max_n_splitter
            Splitter* splitters_
            int* expanding_splitters
            int n_expanding_splitter
            int next_n_splitter
            int expanding_n_splitter
            int start_expanding_splitter

        if tree.max_depth <= 10:
            init_capacity = (2 ** (tree.max_depth + 1)) - 1
        else:
            init_capacity = 2047

        tree._resize(init_capacity)

        ##################################################################
        # Initialization
        ##################################################################
        start_expanding_splitter = 0
        n_splitter = 1
        n_expanding_splitter = 1
        max_n_splitter = 1
        next_n_splitter = 1
        splitters_ = <Splitter*> malloc(n_splitter * sizeof(Splitter))
        expanding_splitters = <int*> malloc(n_splitter * sizeof(int))

        # compute the stats for the root node
        cdef:
            double root_sum_y = 0.0
            double root_sum_sq_y = 0.0
            int root_n_samples = n_samples
            double root_sum_weighted_samples = sum_total_weighted_samples
        weighted_sum_y(y, sample_weight, &root_sum_y, &root_sum_sq_y)
        # init a stats_node and split_record which will be used to create
        # the splitter
        cdef:
            SplitRecord split_record
            StatsNode current_stats_node, left_stats_node, right_stats_node
        stats_node_clear(&left_stats_node)
        stats_node_reset(&right_stats_node, root_sum_y,
                         root_sum_sq_y, root_n_samples,
                         root_sum_weighted_samples)
        stats_node_reset(&current_stats_node, root_sum_y,
                         root_sum_sq_y, root_n_samples,
                         root_sum_weighted_samples)
        split_record_reset(&split_record, 0, 0, NAN,
                           _impurity_mse(&right_stats_node),
                           -INFINITY, 0,
                           &current_stats_node, &left_stats_node,
                           &right_stats_node)
        # create the tree node
        split_record.nid = tree._add_node_with_value(
            TREE_UNDEFINED,
            0, TREE_NOT_LEAF,
            FEAT_UNKNOWN,
            TREE_UNDEFINED,
            split_record.impurity,
            n_samples,
            sum_total_weighted_samples,
            root_sum_y / root_sum_weighted_samples)
        # init the splitter
        splitter_init(&splitters_[0],
                      FEAT_UNKNOWN, TREE_UNDEFINED,
                      &split_record, self.min_samples_leaf,
                      self.min_weight_leaf)
        # add the correspondence nid to splitter idx
        expanding_splitters[0] = 0

        cdef:
            int i
            int* X_nid = <int*> malloc(n_samples * sizeof(int))
            int start_reset_count_X = 0
        for i in range(n_samples):
            X_nid[i] = 0

        cdef:
            int j
            int current_depth = 0
            int n_visited_feature = 0
            int feat_idx
            int[::1] shuffled_feature_idx
            int sample_idx_sorted
            int[::1] X_col
            bint b_grow
            int X_idx
            int parent_n_left_samples
            int splitter_idx
            bint b_impurity
            bint b_samples_split
            bint b_samples_leaf
            int X_idx_init
            float X_init
            double threshold

        #######################################################################
        # Tree growing
        #######################################################################
        while current_depth < self.max_depth:
            # shuffled_feature_idx = rng.permutation(range(n_features)).astype(np.int32)

            n_visited_feature = 0

            # ENUMERATE ALL SPLITS
            for i in range(n_features):
                # feat_idx = shuffled_feature_idx[i]
                feat_idx = i

                if n_visited_feature >= self.max_features:
                    break

                # reset the split record for the current feature
                # keeping the best split record untouched
                X_idx_init = X_idx_sorted[0, feat_idx]
                X_init = X[X_idx_init, feat_idx]
                for j in range(n_splitter):
                    splitter_reset(&splitters_[expanding_splitters[
                        j + start_expanding_splitter]],
                                   feat_idx, X_idx_init, X_init)
                # evaluate all possible split for the different samples
                for j in range(n_samples):
                    sample_idx_sorted = X_idx_sorted[j, feat_idx]
                    if X_nid[sample_idx_sorted] != -1:
                        # FIXME: need to incorporate constant feature
                        splitter_node_evaluate_split(
                            &splitters_[X_nid[sample_idx_sorted]],
                            X, y, sample_weight, sum_total_weighted_samples,
                            sample_idx_sorted)

                n_visited_feature += 1

            # EXPAND SPLITTERS
            start_reset_count_X = max_n_splitter
            max_n_splitter += 2 * n_splitter
            n_expanding_splitter = 0
            splitters_ = <Splitter*> realloc(splitters_, max_n_splitter * sizeof(Splitter))
            expanding_splitters = <int*> realloc(expanding_splitters, max_n_splitter * sizeof(int))
            b_grow = 0
            for i in range(n_splitter):
                splitter_idx = expanding_splitters[i + start_expanding_splitter]
                if isnan(splitters_[splitter_idx].best_split_record.threshold):
                    tree._set_node_as_leaf(&splitters_[splitter_idx])
                else:
                    splitter_expand(&splitters_[splitter_idx],
                                    &splitters_[next_n_splitter],
                                    &splitters_[next_n_splitter + 1])
                    parent_nid = splitters_[splitter_idx].split_record.nid
                    left_nid = tree._add_node_set_id(
                        parent_nid, &splitters_[next_n_splitter], 1)
                    right_nid = tree._add_node_set_id(
                        parent_nid, &splitters_[next_n_splitter + 1], 0)
                    tree._update_parent_splitter(&splitters_[splitter_idx],
                                                 left_nid, right_nid)

                    if self._check_minimum_stats(&splitters_[next_n_splitter]):
                        b_grow = 1
                        expanding_splitters[start_expanding_splitter + n_splitter +
                                            n_expanding_splitter] = left_nid
                        n_expanding_splitter += 1
                    else:
                        tree._set_node_as_leaf(&splitters_[next_n_splitter])

                    if self._check_minimum_stats(&splitters_[next_n_splitter + 1]):
                        b_grow = 1
                        expanding_splitters[start_expanding_splitter + n_splitter +
                                            n_expanding_splitter] = right_nid
                        n_expanding_splitter += 1
                    else:
                        tree._set_node_as_leaf(&splitters_[next_n_splitter + 1])

                    next_n_splitter += 2

            # REASSIGN SPLITTERS ID TO SAMPLE
            if b_grow:
                for i in range(n_samples):
                    parent_nid = X_nid[i]
                    if parent_nid != -1:
                        threshold = splitters_[parent_nid].best_split_record.threshold
                        feat_idx = splitters_[parent_nid].best_split_record.feature
                        if X[i, feat_idx] <= threshold:
                            if tree.nodes[parent_nid].left_child == TREE_LEAF:
                                X_nid[i] = -1
                            else:
                                if (tree.nodes[tree.nodes[parent_nid].left_child].left_child == TREE_LEAF and
                                    tree.nodes[tree.nodes[parent_nid].left_child].right_child == TREE_LEAF):
                                    X_nid[i] = -1
                                else:
                                    X_nid[i] = tree.nodes[parent_nid].left_child
                        else:
                            if tree.nodes[parent_nid].right_child == TREE_LEAF:
                                X_nid[i] = -1
                            else:
                                if (tree.nodes[tree.nodes[parent_nid].right_child].left_child == TREE_LEAF and
                                    tree.nodes[tree.nodes[parent_nid].right_child].right_child == TREE_LEAF):
                                    X_nid[i] = -1
                                else:
                                    X_nid[i] = tree.nodes[parent_nid].right_child

            start_expanding_splitter += n_splitter
            n_splitter = n_expanding_splitter
            if b_grow:
                current_depth += 1
            else:
                break

        # Set all remaining nodes as leaf
        for i in range(n_splitter):
            splitter_idx = start_expanding_splitter + i
            tree._set_node_as_leaf(&splitters_[expanding_splitters[splitter_idx]])

        # Deallocate X_nid and splitters_
        free(splitters_)
        free(X_nid)

        rc = tree._resize_c(tree.node_count)

        if rc >= 0:
            tree.max_depth = current_depth
