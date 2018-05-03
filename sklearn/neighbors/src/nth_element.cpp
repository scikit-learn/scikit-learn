#include "nth_element.h"
#include <algorithm>

class IndexComparator {
    double *data;
    Py_intptr_t split_dim, n_features;

public:
    IndexComparator(double *data, Py_intptr_t split_dim, Py_intptr_t n_features):
        data(data), split_dim(split_dim), n_features(n_features) {}

    bool operator()(Py_intptr_t a, Py_intptr_t b) {
        return data[a * n_features + split_dim]
            < data[b * n_features + split_dim];
    }
};

int partition_node_indices(double *data,
                           Py_intptr_t *node_indices,
                           Py_intptr_t split_dim,
                           Py_intptr_t split_index,
                           Py_intptr_t n_features,
                           Py_intptr_t n_points) {
    IndexComparator index_comparator(data, split_dim, n_features);
    std::nth_element(node_indices,
                     node_indices + split_index,
                     node_indices + n_points,
                     index_comparator);
    return 0;
}

