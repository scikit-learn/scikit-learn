from numbers import Integral
from warnings import warn

import numpy as np
import six

from ._export import SENTINEL
from ._criterion import FriedmanMSE
from ._classes_m5 import is_leaf, ConstantLeafModel, M5Base, check_is_fitted
from . import _tree


def export_text_m5(decision_tree, out_file=SENTINEL, max_depth=None,
                   feature_names=None, class_names=None, label='all',
                   target_name=None,
                   # filled=False, leaves_parallel=False,
                   impurity=True,
                   node_ids=False, proportion=False,
                   # rounded=False, rotate=False,
                   special_characters=False, precision=3, **kwargs):
    """Export a decision tree in TXT format.

    Note: this should be merged with ._export.export_text

    Inspired by WEKA and by
    >>> from sklearn.tree import export_graphviz

    This function generates a human-readable, text representation of the
    decision tree, which is then written into `out_file`.

    The sample counts that are shown are weighted with any sample_weights that
    might be present.

    Read more in the :ref:`User Guide <tree>`.

    Parameters
    ----------
    decision_tree : decision tree classifier
        The decision tree to be exported to text.

    out_file : file object or string, optional (default='tree.dot')
        Handle or name of the output file. If ``None``, the result is
        returned as a string. This will the default from version 0.20.

    max_depth : int, optional (default=None)
        The maximum depth of the representation. If None, the tree is fully
        generated.

    feature_names : list of strings, optional (default=None)
        Names of each of the features.

    class_names : list of strings, bool or None, optional (default=None)
        Names of each of the target classes in ascending numerical order.
        Only relevant for classification and not supported for multi-output.
        If ``True``, shows a symbolic representation of the class name.

    label : {'all', 'root', 'none'}, optional (default='all')
        Whether to show informative labels for impurity, etc.
        Options include 'all' to show at every node, 'root' to show only at
        the top root node, or 'none' to not show at any node.

    target_name : optional string with the target name. If not provided, the
        target will not be displayed in the equations

    impurity : bool, optional (default=True)
        When set to ``True``, show the impurity at each node.

    node_ids : bool, optional (default=False)
        When set to ``True``, show the ID number on each node.

    proportion : bool, optional (default=False)
        When set to ``True``, change the display of 'values' and/or 'samples'
        to be proportions and percentages respectively.

    special_characters : bool, optional (default=False)
        When set to ``False``, ignore special characters for PostScript
        compatibility.

    precision : int, optional (default=3)
        Number of digits of precision for floating point in the values of
        impurity, threshold and value attributes of each node.

    kwargs : other keyword arguments for the linear model printer

    Returns
    -------
    dot_data : string
        String representation of the input tree in GraphViz dot format.
        Only returned if ``out_file`` is None.

        .. versionadded:: 0.18

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn import tree

    >>> clf = tree.DecisionTreeClassifier()
    >>> iris = load_iris()

    >>> clf = clf.fit(iris.data, iris.target)
    >>> tree_to_text(clf, out_file='tree.txt')                # doctest: +SKIP

    """

    models = []

    def add_model(node_model):
        models.append(node_model)
        return len(models)

    def node_to_str(tree, node_id, criterion, node_models=None):
        """ Generates the node content string """

        # Should labels be shown?
        labels = (label == 'root' and node_id == 0) or label == 'all'

        # PostScript compatibility for special characters
        if special_characters:
            characters = ['&#35;', '<SUB>', '</SUB>', '&le;', '<br/>', '>']
            node_string = '<'
        else:
            characters = ['#', '[', ']', '<=', '\\n', '']
            node_string = ''

        # -- If this node not a leaf,Write the split decision criteria (x <= y)
        leaf = is_leaf(node_id, tree)
        if not leaf:
            if feature_names is not None:
                feature = feature_names[tree.feature[node_id]]
            else:
                feature = "X%s%s%s" % (characters[1],
                                       tree.feature[node_id],  # feature id for the split
                                       characters[2])
            node_string += '%s %s %s' % (feature,
                                         characters[3],  # <=
                                         round(tree.threshold[node_id],   # threshold for the split
                                               precision))
        else:
            node_string += 'LEAF'

        # Node details - start bracket [
        node_string += ' %s' % characters[1]

        # -- Write impurity
        if impurity:
            if isinstance(criterion, FriedmanMSE):
                criterion = "friedman_mse"
            elif not isinstance(criterion, str):
                criterion = "impurity"
            if labels:
                node_string += '%s=' % criterion
            node_string += str(round(tree.impurity[node_id], precision)) + ', '

        # -- Write node sample count
        if labels:
            node_string += 'samples='
        if proportion:
            percent = (100. * tree.n_node_samples[node_id] /
                       float(tree.n_node_samples[0]))
            node_string += str(round(percent, 1)) + '%'
        else:
            node_string += str(tree.n_node_samples[node_id])

        # Node details - end bracket ]
        node_string += '%s' % characters[2]

        # -- Write node class distribution / regression value
        if tree.n_outputs == 1:
            value = tree.value[node_id][0, :]
        else:
            value = tree.value[node_id]

        if proportion and tree.n_classes[0] != 1:
            # For classification this will show the proportion of samples
            value = value / tree.weighted_n_node_samples[node_id]
        if tree.n_classes[0] == 1:
            # Regression
            value_text = np.around(value, precision)
        elif proportion:
            # Classification
            value_text = np.around(value, precision)
        elif np.all(np.equal(np.mod(value, 1), 0)):
            # Classification without floating-point weights
            value_text = value.astype(int)
        else:
            # Classification with floating-point weights
            value_text = np.around(value, precision)

        # Strip whitespace
        value_text = str(value_text.astype('S32')).replace("b'", "'")
        value_text = value_text.replace("' '", ", ").replace("'", "")
        if tree.n_classes[0] == 1 and tree.n_outputs == 1:
            value_text = value_text.replace("[", "").replace("]", "")
        value_text = value_text.replace("\n ", characters[4])

        if node_models is None:
            node_string += ' : '
            if labels:
                node_string += 'value='
        else:
            nodemodel = node_models[node_id]
            model_err_val = np.around(nodemodel.error, precision)
            if leaf:
                if isinstance(nodemodel, ConstantLeafModel):
                    # the model does not contain the value. rely on the value_text computed from tree
                    value_text = " : %s (err=%s, params=%s)" % (value_text, model_err_val, nodemodel.n_params)
                else:
                    # put the model in the stack, we'll write it later
                    model_id = add_model(nodemodel)
                    value_text = " : LM%s (err=%s, params=%s)" % (model_id, model_err_val, nodemodel.n_params)
            else:
                # replace the value text with error at this node and number of parameters
                value_text = " (err=%s, params=%s)" % (model_err_val, nodemodel.n_params)

        node_string += value_text

        # Write node majority class
        if (class_names is not None and
                tree.n_classes[0] != 1 and
                tree.n_outputs == 1):
            # Only done for single-output classification trees
            node_string += ', '
            if labels:
                node_string += 'class='
            if class_names is not True:
                class_name = class_names[np.argmax(value)]
            else:
                class_name = "y%s%s%s" % (characters[1],
                                          np.argmax(value),
                                          characters[2])
            node_string += class_name

        return node_string + characters[5]

    def recurse(tree, node_id, criterion, parent=None, depth=0, node_models=None):
        if node_id == _tree.TREE_LEAF:
            raise ValueError("Invalid node_id %s" % _tree.TREE_LEAF)

        # Add node with description
        if max_depth is None or depth <= max_depth:
            indent_str = ("|   " * depth)
            if node_ids:
                out_file.write('%d| %s%s\n' % (node_id, indent_str, node_to_str(tree, node_id, criterion,
                                                                                node_models=node_models)))
            else:
                out_file.write('%s%s\n' % (indent_str, node_to_str(tree, node_id, criterion, node_models=node_models)))

            # Recurse on Children if needed
            left_child = tree.children_left[node_id]
            right_child = tree.children_right[node_id]

            if left_child != _tree.TREE_LEAF:
                # that means that node_id is not a leaf (see is_leaf() below.): recurse on children
                recurse(tree, left_child, criterion=criterion, parent=node_id, depth=depth + 1,
                        node_models=node_models)
                recurse(tree, right_child, criterion=criterion, parent=node_id, depth=depth + 1,
                        node_models=node_models)

        else:
            ranks['leaves'].append(str(node_id))
            out_file.write('%d| (...)\n')

    def write_models(models):
        for i, model in enumerate(models):
            out_file.write("LM%s: %s\n" % (i + 1, model.to_text(feature_names=feature_names, precision=precision,
                                                                target_name=target_name, **kwargs)))

    # Main
    check_is_fitted(decision_tree, 'tree_')
    own_file = False
    return_string = False
    try:
        if out_file == SENTINEL:
            warn("out_file can be set to None starting from 0.18. This will be"
                 "the default in 0.20.", DeprecationWarning)
            out_file = "tree.txt"

        if isinstance(out_file, str):
            out_file = open(out_file, "w", encoding="utf-8")
            own_file = True

        if out_file is None:
            return_string = True
            out_file = six.StringIO()

        if isinstance(precision, Integral):
            if precision < 0:
                raise ValueError("'precision' should be greater or equal to 0."
                                 " Got {} instead.".format(precision))
        else:
            raise ValueError("'precision' should be an integer. Got {}"
                             " instead.".format(type(precision)))

        # Check length of feature_names before getting into the tree node
        # Raise error if length of feature_names does not match
        # n_features_ in the decision_tree
        if feature_names is not None:
            if len(feature_names) != decision_tree.n_features_:
                raise ValueError("Length of feature_names, %d "
                                 "does not match number of features, %d"
                                 % (len(feature_names),
                                    decision_tree.n_features_))

        # The depth of each node for plotting with 'leaf' option TODO probably remove
        ranks = {'leaves': []}

        # Tree title
        if isinstance(decision_tree, M5Base):
            if hasattr(decision_tree, 'installed_smoothing_constant'):
                details = "pre-smoothed with constant %s" % decision_tree.installed_smoothing_constant
            else:
                if decision_tree.use_smoothing == 'installed':
                    details = "under construction - not pre-smoothed yet"
                else:
                    details = "unsmoothed - but this can be done at prediction time"

            # add more info or M5P
            out_file.write('%s (%s):\n' % (type(decision_tree).__name__, details))
        else:
            # generic title
            out_file.write('%s :\n' % type(decision_tree).__name__)

        # some space for readability
        out_file.write('\n')

        # Now recurse the tree and add node & edge attributes
        if isinstance(decision_tree, _tree.Tree):
            recurse(decision_tree, 0, criterion="impurity")
        elif isinstance(decision_tree, M5Base) and hasattr(decision_tree, 'node_models'):
            recurse(decision_tree.tree_, 0, criterion=decision_tree.criterion, node_models=decision_tree.node_models)

            # extra step: write all models
            out_file.write("\n")
            write_models(models)
        else:
            recurse(decision_tree.tree_, 0, criterion=decision_tree.criterion)

        # Return the text if needed
        if return_string:
            return out_file.getvalue()

    finally:
        if own_file:
            out_file.close()
