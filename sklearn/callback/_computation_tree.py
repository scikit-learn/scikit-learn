# License: BSD 3 clause
# Authors: the scikit-learn developers


class ComputationNode:
    """A node in a computation tree

    Parameters
    ----------
    estimator_name : str
        The name of the estimator this computation node belongs to.

    parent : ComputationNode instance, default=None
        The parent node. None means this is the root.

    max_iter : int, default=None
        The number of its children. None means it's a leaf.

    description : str, default=None
        A description of this computation node. None means it's a leaf.

    Attributes
    ----------
    children : list
        The list of its children nodes. For a leaf, it's an empty list
    """

    def __init__(
        self,
        estimator_name,
        description=None,
        max_iter=None,
        idx=0,
        parent=None,
    ):
        # estimator_name and description are tuples because an estimator can be
        # a sub-estimator of a meta-estimator. In that case, the root of the computation
        # tree of the sub-estimator and a leaf of the computation tree of the
        # meta-estimator correspond to the same computation step. Therefore, both
        # nodes are merged into a single node, retaining the information of both.
        self.estimator_name = (estimator_name,)
        self.description = (description,)

        self.parent = parent
        self.max_iter = max_iter
        self.idx = idx

        self.children = []

    @property
    def depth(self):
        """The depth of this node in the computation tree"""
        return 0 if self.parent is None else self.parent.depth + 1
    
    @property
    def path(self):
        """List of all the nodes in the path from the root to this node"""
        return [self] if self.parent is None else self.parent.path + [self]

    def __iter__(self):
        """Pre-order depth-first traversal"""
        yield self
        for node in self.children:
            yield from node


def build_computation_tree(estimator_name, levels, parent=None, idx=0):
    """Build the computation tree from the description of the levels"""
    this_level = levels[0]

    node = ComputationNode(
        estimator_name=estimator_name,
        parent=parent,
        max_iter=this_level["max_iter"],
        description=this_level["descr"],
        idx=idx,
    )

    if parent is not None and parent.max_iter is None:
        # parent node is a leaf of the computation tree of an outer estimator. It means
        # that this node is the root of the computation tree of this estimator. They
        # both correspond the same computation step, so we merge both nodes.
        node.description = parent.description + node.description
        node.estimator_name = parent.estimator_name + node.estimator_name
        node.parent = parent.parent
        node.idx = parent.idx
        parent.parent.children[node.idx] = node

    if node.max_iter is not None:
        for i in range(node.max_iter):
            node.children.append(
                build_computation_tree(estimator_name, levels[1:], parent=node, idx=i)
            )

    return node
