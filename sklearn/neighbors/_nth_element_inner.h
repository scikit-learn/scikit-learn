#include <algorithm>

template<class D, class I>
class IndexComparator {
private:
    const D *data;
    I split_dim, n_features;
public:
    IndexComparator(const D *data, const I &split_dim, const I &n_features):
        data(data), split_dim(split_dim), n_features(n_features) {}

    bool operator()(const I &a, const I &b) const {
        D a_value = data[a * n_features + split_dim];
        D b_value = data[b * n_features + split_dim];
        return a_value == b_value ? a < b : a_value < b_value;
    }
};

template<class D, class I>
void partition_node_indices_inner(
    const D *data,
    I *node_indices,
    const I &split_dim,
    const I &split_index,
    const I &n_features,
    const I &n_points) {
    IndexComparator<D, I> index_comparator(data, split_dim, n_features);
    std::nth_element(
        node_indices,
        node_indices + split_index,
        node_indices + n_points,
        index_comparator);
}
