# License: BSD 3 clause
# Authors: the scikit-learn developers


class ComputationNode:
    """A node in a computation tree.

    Parameters
    ----------
    estimator_name : str
        The name of the estimator this computation node belongs to.

    stage : str, default=None
        A description of the stage this computation node belongs to.
        None means it's a leaf.

    n_children : int, default=None
        The number of its children. None means it's a leaf.

    idx : int, default=0
        The index of this node among its siblings.

    parent : ComputationNode instance, default=None
        The parent node. None means this is the root.

    Attributes
    ----------
    children : list
        The list of its children nodes. For a leaf, it's an empty list
    """

    def __init__(
        self,
        estimator_name,
        stage=None,
        n_children=None,
        idx=0,
        parent=None,
    ):
        # estimator_name and description are tuples because an estimator can be
        # a sub-estimator of a meta-estimator. In that case, the root of the computation
        # tree of the sub-estimator and a leaf of the computation tree of the
        # meta-estimator correspond to the same computation step. Therefore, both
        # nodes are merged into a single node, retaining the information of both.
        self.estimator_name = (estimator_name,)
        self.stage = (stage,)

        self.parent = parent
        self.n_children = n_children
        self.idx = idx

        self.children = []

    @property
    def depth(self):
        """The depth of this node in the computation tree."""
        return 0 if self.parent is None else self.parent.depth + 1

    @property
    def path(self):
        """List of all the nodes in the path from the root to this node."""
        return [self] if self.parent is None else self.parent.path + [self]

    def __iter__(self):
        """Pre-order depth-first traversal"""
        yield self
        for node in self.children:
            yield from node


def build_computation_tree(estimator_name, tree_structure, parent=None, idx=0):
    """Build the computation tree from the description of the levels.

    Parameters
    ----------
    estimator_name : str
        The name of the estimator this computation tree belongs to.

    tree_structure : list of dict
        The description of the stages of the computation tree. Each dict must have
        the following keys:
            - stage: str
                A human readable description of the stage.
            - n_children: int or None
                The number of its children. None means it's a leaf.

    parent : ComputationNode instance, default=None
        The parent node. None means this is the root.

    idx : int, default=0
        The index of this node among its siblings.

    Returns
    -------
    computation_tree : ComputationNode instance
        The root of the computation tree.
    """
    this_stage = tree_structure[0]

    node = ComputationNode(
        estimator_name=estimator_name,
        parent=parent,
        n_children=this_stage["n_children"],
        stage=this_stage["stage"],
        idx=idx,
    )

    if parent is not None and parent.n_children is None:
        # parent node is a leaf of the computation tree of an outer estimator. It means
        # that this node is the root of the computation tree of this estimator. They
        # both correspond the same computation step, so we merge both nodes.
        node.stage = parent.stage + node.stage
        node.estimator_name = parent.estimator_name + node.estimator_name
        node.parent = parent.parent
        node.idx = parent.idx
        parent.parent.children[node.idx] = node

    if node.n_children is not None:
        for i in range(node.n_children):
            node.children.append(
                build_computation_tree(
                    estimator_name, tree_structure[1:], parent=node, idx=i
                )
            )

    return node
