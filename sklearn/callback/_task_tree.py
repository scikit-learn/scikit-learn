# License: BSD 3 clause
# Authors: the scikit-learn developers


class TaskNode:
    """A node in a task tree.

    Parameters
    ----------
    estimator_name : str
        The name of the estimator this task node belongs to.

    name : str, default=None
        The name of the task this node represents.

    max_subtasks : int or None, default=0
        The maximum number of its children. 0 means it's a leaf.
        None means the number of children is not known in advance.

    idx : int, default=0
        The index of this node among its siblings.

    parent : TaskNode instance, default=None
        The parent node. None means this is the root.

        Note that the root task of an estimator can become an intermediate node
        of a meta-estimator.

    Attributes
    ----------
    children : dict
        A mapping from the index of a child to the child node `{idx: TaskNode}`.
        For a leaf, it's an empty dictionary.
    """

    def __init__(
        self,
        estimator_name,
        name="fit",
        max_subtasks=0,
        idx=0,
        parent=None,
    ):
        # estimator_name and name are tuples because an estimator can be
        # a sub-estimator of a meta-estimator. In that case, the root of the task
        # tree of the sub-estimator and a leaf of the task tree of the
        # meta-estimator correspond to the same computation step. Therefore, both
        # nodes are merged into a single node, retaining the information of both.
        self.estimator_name = (estimator_name,)
        self.name = (name,)

        self.max_subtasks = max_subtasks
        self.idx = idx
        self.parent = parent

        # Children stored in a dict indexed by their idx for easy access because the
        # order in which self.children is populated is not guaranteed to follow the
        # order of the idx du to parallelism.
        self.children = {}

    def _add_child(self, name, max_subtasks, idx):
        child = TaskNode(
            estimator_name=self.estimator_name[-1],
            name=name,
            max_subtasks=max_subtasks,
            idx=idx,
            parent=self,
        )
        self.children[idx] = child

        return child

    def _merge_with(self, task_node):
        self.parent = task_node.parent
        self.idx = task_node.idx
        task_node.parent.children[self.idx] = self

        self.name = task_node.name + self.name
        self.estimator_name = task_node.estimator_name + self.estimator_name

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
        for node in self.children.values():
            yield from node
