from graphviz import Digraph

from sklearn.ensemble._hist_gradient_boosting.gradient_boosting import BaseHistGradientBoosting


def plot_tree(est_or_grower, est_lightgbm=None, tree_index=0, view=True,
              **kwargs):
    """Plot the i'th predictor tree of an estimator, or a grower's tree

    est_or_grower can either be a GradientBoostingMachine instance or a
    TreeGrower. In this latter case tree_index is ignored, and more debugging
    info are displayed. Trees displayed from TreeGrower has additional
    profiling information that are not kept in the predictor trees that
    result from fitting a GradientBoostingMachine.

    tree_index corresponds to the ith built tree. In a multiclass setting, the
    ith tree isn't necessarily the tree built durint the ith iteration because
    there are K trees per iteration. For example with 3 classes,
    tree_index=5 will print the third tree of the second iteration.

    Can also plot a LightGBM estimator (on the left) for comparison.

    Requires matplotlib and graphviz (both python package and binary program).

    kwargs are passed to graphviz.Digraph()

    Example: plotting.plot_tree(est_pygbm, est_lightgbm, view=False,
    filename='output') will silently save output to output.pdf
    """
    def make_pygbm_tree():
        def add_predictor_node(node_idx, parent=None, decision=None):
            iteration = tree_index // est_or_grower.n_trees_per_iteration_
            k = tree_index % est_or_grower.n_trees_per_iteration_
            predictor_tree = est_or_grower._predictors[iteration][k]
            node = predictor_tree.nodes[node_idx]
            name = 'split__{}'.format(node_idx)
            label = 'split_feature_index: {}'.format(
                node['feature_idx'])
            label += r'\nthreshold: {:.3f}'.format(node['threshold'])
            label += r'\ngain: {:.3E}'.format(node['gain'])
            label += r'\nvalue: {:.3f}'.format(node['value'])
            label += r'\ncount: {:,}'.format(node['count'])
            label += r'\nnans_go_left: {:,}'.format(node['missing_go_to_left'])

            graph.node(name, label=label)
            if not node['is_leaf']:
                add_predictor_node(node['left'], name, decision='<=')
                add_predictor_node(node['right'], name, decision='>')

            if parent is not None:
                graph.edge(parent, name, decision)

        def add_grower_node(node, parent=None, decision=None):
            name = 'split__{0}'.format(id(node))
            si = node.split_info
            if si is None:
                feature_idx = 0
                bin_idx = 0
                gain = 0.
                sum_gradients = 0.
                sum_hessians = 0.
            else:
                feature_idx = si.feature_idx
                gain = 0. if si.gain is None else si.gain
                bin_idx = si.bin_idx
                sum_gradients = si.gradient_left + si.gradient_right
                sum_hessians = si.hessian_left + si.hessian_right

            value = 0. if node.value is None else node.value
            label = 'split_feature_index: {}'.format(feature_idx)
            label += r'\nbin threshold: {}'.format(bin_idx)
            label += r'\ngain: {:.3E}'.format(gain)
            label += r'\nvalue: {:.3f}'.format(value)
            label += r'\ncount: {:,}'.format(node.sample_indices.shape[0])
            label += r'\nhist substration: {}'.format(node.hist_subtraction)
            label += r'\nhist speed: {:.3E}'.format(
                node.construction_speed)
            label += r'\nfind split time: {:.4f}'.format(node.find_split_time)
            label += r'\napply split time: {:.4f}'.format(
                node.apply_split_time)
            label += r'\nsum gradients: {:.3E}'.format(sum_gradients)
            label += r'\nsum hessians: {:.3E}'.format(sum_hessians)

            graph.node(name, label=label)
            if node.value is None:  # not a leaf node
                add_grower_node(node.left_child, name, decision='<=')
                add_grower_node(node.right_child, name, decision='>')

            if parent is not None:
                graph.edge(parent, name, decision)

        if isinstance(est_or_grower, BaseHistGradientBoosting):
            add_predictor_node(0)
        # elif isinstance(est_or_grower, pygbm.grower.TreeGrower):
        #     add_grower_node(est_or_grower.root)

    # make lightgbm tree
    if est_lightgbm is not None:
        import lightgbm as lb
        graph = lb.create_tree_digraph(
            est_lightgbm,
            tree_index=tree_index,
            show_info=['split_gain', 'internal_value', 'internal_count',
                       'leaf_count'],
            **kwargs)
    else:
        graph = Digraph(**kwargs)

    # make pygbm tree
    make_pygbm_tree()

    graph.render(view=view)
