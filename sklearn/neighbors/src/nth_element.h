#include <Python.h>

int partition_node_indices(double *data,
                           Py_intptr_t *node_indices,
                           Py_intptr_t split_dim,
                           Py_intptr_t split_index,
                           Py_intptr_t n_features,
                           Py_intptr_t n_points);

