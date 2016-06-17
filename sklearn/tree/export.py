"""
This module defines export functions for decision trees.
"""

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Noel Dawe <noel@dawe.me>
#          Satrajit Gosh <satrajit.ghosh@gmail.com>
#          Trevor Stephens <trev.stephens@gmail.com>
<<<<<<< 59bf211789d29d70642fd64d6d4bef95ff518e8d
# License: BSD 3 clause
=======
#          nvdv <modestdev@gmail.com>
# Licence: BSD 3 clause
>>>>>>> Add decision tree plotting.

import numpy as np

from ..externals import six

from . import _criterion
from . import _tree


def _color_brew(n):
    """Generate n colors with equally spaced hues.

    Parameters
    ----------
    n : int
        The number of colors required.

    Returns
    -------
    color_list : list, length n
        List of n tuples of form (R, G, B) being the components of each color.
    """
    color_list = []

    # Initialize saturation & value; calculate chroma & value shift
    s, v = 0.75, 0.9
    c = s * v
    m = v - c

    for h in np.arange(25, 385, 360. / n).astype(int):
        # Calculate some intermediate values
        h_bar = h / 60.
        x = c * (1 - abs((h_bar % 2) - 1))
        # Initialize RGB with same hue & chroma as our color
        rgb = [(c, x, 0),
               (x, c, 0),
               (0, c, x),
               (0, x, c),
               (x, 0, c),
               (c, 0, x),
               (c, x, 0)]
        r, g, b = rgb[int(h_bar)]
        # Shift the initial RGB values to match value and store
        rgb = [(int(255 * (r + m))),
               (int(255 * (g + m))),
               (int(255 * (b + m)))]
        color_list.append(rgb)

    return color_list



def get_graphviz_source(decision_tree, max_depth=None, special_characters=False,
                        feature_names=None, class_names=None, label='all',
                        filled=False, leaves_parallel=False, impurity=True,
                        node_ids=False, proportion=False, rotate=False,
                        rounded=False):
    """Returns decision tree source in GraphViz format."""
    def get_color(value):
        # Find the appropriate color & intensity for a node
        if colors['bounds'] is None:
            # Classification tree
            color = list(colors['rgb'][np.argmax(value)])
            sorted_values = sorted(value, reverse=True)
            if len(sorted_values) == 1:
                alpha = 0
            else:
                alpha = int(np.round(255 * (sorted_values[0] - sorted_values[1]) /
                                           (1 - sorted_values[1]), 0))
        else:
            # Regression tree or multi-output
            color = list(colors['rgb'][0])
            alpha = int(np.round(255 * ((value - colors['bounds'][0]) /
                                        (colors['bounds'][1] -
                                         colors['bounds'][0])), 0))

        # Return html color code in #RRGGBBAA format
        color.append(alpha)
        hex_codes = [str(i) for i in range(10)]
        hex_codes.extend(['a', 'b', 'c', 'd', 'e', 'f'])
        color = [hex_codes[c // 16] + hex_codes[c % 16] for c in color]

        return '#' + ''.join(color)

    def node_to_str(tree, node_id, criterion):
        # Generate the node content string
        if tree.n_outputs == 1:
            value = tree.value[node_id][0, :]
        else:
            value = tree.value[node_id]

        # Should labels be shown?
        labels = (label == 'root' and node_id == 0) or label == 'all'

        # PostScript compatibility for special characters
        if special_characters:
            characters = ['&#35;', '<SUB>', '</SUB>', '&le;', '<br/>', '>']
            node_string = '<'
        else:
            characters = ['#', '[', ']', '<=', '\\n', '"']
            node_string = '"'

        # Write node ID
        if node_ids:
            if labels:
                node_string += 'node '
            node_string += characters[0] + str(node_id) + characters[4]

        # Write decision criteria
        if tree.children_left[node_id] != _tree.TREE_LEAF:
            # Always write node decision criteria, except for leaves
            if feature_names is not None:
                feature = feature_names[tree.feature[node_id]]
            else:
                feature = "X%s%s%s" % (characters[1],
                                       tree.feature[node_id],
                                       characters[2])
            node_string += '%s %s %s%s' % (feature,
                                           characters[3],
                                           round(tree.threshold[node_id], 4),
                                           characters[4])

        # Write impurity
        if impurity:
            if isinstance(criterion, _criterion.FriedmanMSE):
                criterion = "friedman_mse"
            elif not isinstance(criterion, six.string_types):
                criterion = "impurity"
            if labels:
                node_string += '%s = ' % criterion
            node_string += (str(round(tree.impurity[node_id], 4)) +
                            characters[4])

        # Write node sample count
        if labels:
            node_string += 'samples = '
        if proportion:
            percent = (100. * tree.n_node_samples[node_id] /
                       float(tree.n_node_samples[0]))
            node_string += (str(round(percent, 1)) + '%' +
                            characters[4])
        else:
            node_string += (str(tree.n_node_samples[node_id]) +
                            characters[4])

        # Write node class distribution / regression value
        if proportion and tree.n_classes[0] != 1:
            # For classification this will show the proportion of samples
            value = value / tree.weighted_n_node_samples[node_id]
        if labels:
            node_string += 'value = '
        if tree.n_classes[0] == 1:
            # Regression
            value_text = np.around(value, 4)
        elif proportion:
            # Classification
            value_text = np.around(value, 2)
        elif np.all(np.equal(np.mod(value, 1), 0)):
            # Classification without floating-point weights
            value_text = value.astype(int)
        else:
            # Classification with floating-point weights
            value_text = np.around(value, 4)
        # Strip whitespace
        value_text = str(value_text.astype('S32')).replace("b'", "'")
        value_text = value_text.replace("' '", ", ").replace("'", "")
        if tree.n_classes[0] == 1 and tree.n_outputs == 1:
            value_text = value_text.replace("[", "").replace("]", "")
        value_text = value_text.replace("\n ", characters[4])
        node_string += value_text + characters[4]

        # Write node majority class
        if (class_names is not None and
                tree.n_classes[0] != 1 and
                tree.n_outputs == 1):
            # Only done for single-output classification trees
            if labels:
                node_string += 'class = '
            if class_names is not True:
                class_name = class_names[np.argmax(value)]
            else:
                class_name = "y%s%s%s" % (characters[1],
                                          np.argmax(value),
                                          characters[2])
            node_string += class_name

        # Clean up any trailing newlines
        if node_string[-2:] == '\\n':
            node_string = node_string[:-2]
        if node_string[-5:] == '<br/>':
            node_string = node_string[:-5]

        return node_string + characters[5]

    def recurse(tree, node_id, criterion, parent=None, depth=0):
        if node_id == _tree.TREE_LEAF:
            raise ValueError("Invalid node_id %s" % _tree.TREE_LEAF)

        left_child = tree.children_left[node_id]
        right_child = tree.children_right[node_id]

        # Add node with description
        if max_depth is None or depth <= max_depth:

            # Collect ranks for 'leaf' option in plot_options
            if left_child == _tree.TREE_LEAF:
                ranks['leaves'].append(str(node_id))
            elif str(depth) not in ranks:
                ranks[str(depth)] = [str(node_id)]
            else:
                ranks[str(depth)].append(str(node_id))

            source.append('%d [label=%s'
                           % (node_id,
                              node_to_str(tree, node_id, criterion)))

            if filled:
                # Fetch appropriate color for node
                if 'rgb' not in colors:
                    # Initialize colors and bounds if required
                    colors['rgb'] = _color_brew(tree.n_classes[0])
                    if tree.n_outputs != 1:
                        # Find max and min impurities for multi-output
                        colors['bounds'] = (np.min(-tree.impurity),
                                            np.max(-tree.impurity))
                    elif tree.n_classes[0] == 1 and len(np.unique(tree.value)) != 1:
                        # Find max and min values in leaf nodes for regression
                        colors['bounds'] = (np.min(tree.value),
                                            np.max(tree.value))
                if tree.n_outputs == 1:
                    node_val = (tree.value[node_id][0, :] /
                                tree.weighted_n_node_samples[node_id])
                    if tree.n_classes[0] == 1:
                        # Regression
                        node_val = tree.value[node_id][0, :]
                else:
                    # If multi-output color node by impurity
                    node_val = -tree.impurity[node_id]
                source.append(', fillcolor="%s"' % get_color(node_val))
            source.append('] ;\n')

            if parent is not None:
                # Add edge to parent
                source.append('%d -> %d' % (parent, node_id))
                if parent == 0:
                    # Draw True/False labels if parent is root node
                    angles = np.array([45, -45]) * ((rotate - .5) * -2)
                    source.append(' [labeldistance=2.5, labelangle=')
                    if node_id == 1:
                        source.append('%d, headlabel="True"]' % angles[0])
                    else:
                        source.append('%d, headlabel="False"]' % angles[1])
                source.append(' ;\n')

            if left_child != _tree.TREE_LEAF:
                recurse(tree, left_child, criterion=criterion, parent=node_id,
                        depth=depth + 1)
                recurse(tree, right_child, criterion=criterion, parent=node_id,
                        depth=depth + 1)

        else:
            ranks['leaves'].append(str(node_id))

            source.append('%d [label="(...)"' % node_id)
            if filled:
                # color cropped nodes grey
                source.append(', fillcolor="#C0C0C0"')
            source.append('] ;\n' % node_id)

            if parent is not None:
                # Add edge to parent
                source.append('%d -> %d ;\n' % (parent, node_id))

    source = []

    # The depth of each node for plotting with 'leaf' option
    ranks = {'leaves': []}
    # The colors to render each node with
    colors = {'bounds': None}

    source.append('digraph Tree {\n')

    # Specify node aesthetics
    source.append('node [shape=box')
    rounded_filled = []
    if filled:
        rounded_filled.append('filled')
    if rounded:
        rounded_filled.append('rounded')
    if len(rounded_filled) > 0:
        source.append(', style="%s", color="black"' % ", ".join(
            rounded_filled))
    if rounded:
        source.append(', fontname=helvetica')
    source.append('] ;\n')

    # Specify graph & edge aesthetics
    if leaves_parallel:
        source.append('graph [ranksep=equally, splines=polyline] ;\n')
    if rounded:
        source.append('edge [fontname=helvetica] ;\n')
    if rotate:
        source.append('rankdir=LR ;\n')

    # Now recurse the tree and add node & edge attributes
    if isinstance(decision_tree, _tree.Tree):
        recurse(decision_tree, 0, criterion="impurity")
    else:
        recurse(decision_tree.tree_, 0, criterion=decision_tree.criterion)

    # If required, draw leaf nodes at same depth as each other
    if leaves_parallel:
        for rank in sorted(ranks):
            source.append("{rank=same ; " +
                          "; ".join(r for r in ranks[rank]) + "} ;\n")
    source.append("}")
    return ''.join(source)


def export_graphviz(decision_tree, out_file="tree.dot", max_depth=None,
                    feature_names=None, class_names=None, label='all',
                    filled=False, leaves_parallel=False, impurity=True,
                    node_ids=False, proportion=False, rotate=False,
                    rounded=False, special_characters=False):
    """Export a decision tree in DOT format."""
    tree_source = get_graphviz_source(
        tree, max_depth, feature_names, class_names, label, filled,
        leaves_parallel, impurity, node_ids,proportion, rotate, rounded,
        special_characters)
    try:
        if six.PY3:
            source_file = open(out_file, "w", encoding="utf-8")
        else:
            source_file = open(out_file, "wb")
        source_file.write(tree_source)
    finally:
        source_file.close()
