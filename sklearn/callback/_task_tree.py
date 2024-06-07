# License: BSD 3 clause
# Authors: the scikit-learn developers


class TaskNode:
    """A node in a task tree.

    Parameters
    ----------
    task_name : str
        The name of the task this node represents.

    task_id : int
        An identifier for this task that distinguishes it from its siblings. Usually
        the index of this node among its siblings.

    max_tasks : int or None
        The maximum number of its siblings. None means the maximum number of siblings
        is not known in advance.

    estimator_name : str
        The name of the estimator this task node belongs to.

    Attributes
    ----------
    parent : TaskNode instance or None
        The parent node. None means this is the root.

        Note that it's dynamic since the root task of an estimator can become an
        intermediate node of a meta-estimator.

    children_map : dict
        A mapping from the task_id of a child to the child node `{task_id: TaskNode}`.
        For a leaf, it's an empty dictionary.

    max_subtasks : int or None
        The maximum number of subtasks of this node. 0 means it's a leaf. None
        means the maximum number of subtasks is not known in advance.

    prev_estimator_name : str or None
        The estimator name of the node this node was merged with. None if it was not
        merged with another node.

    prev_task_name : str
        The task name of the node this node was merged with. None if it was not
        merged with another node.
    """

    def __init__(self, *, task_name, task_id, max_tasks, estimator_name):
        self.task_name = task_name
        self.task_id = task_id
        self.max_tasks = max_tasks
        self.estimator_name = estimator_name

        self.parent = None
        self.children_map = {}
        self.max_subtasks = 0

        # When an estimator is a sub-estimator of a meta-estimator, the root task of
        # the estimator is merged with the corresponding leaf task of the
        # meta-estimator because both correspond to the same computation step.
        # The root task of the estimator takes the place of the leaf task of the
        # meta-estimator in the task tree but we keep the information about the
        # leaf task it was merged with to fully describe the merged node.
        self.prev_estimator_name = None
        self.prev_task_name = None

    def _add_child(self, task_node):
        if task_node.task_id in self.children_map:
            raise ValueError(
                f"Task node {self.task_name} of estimator {self.estimator_name} "
                f"already has a child with task_id={task_node.task_id}."
            )

        if len(self.children_map) == task_node.max_tasks:
            raise ValueError(
                f"Cannot add child to task node {self.task_name} of estimator "
                f"{self.estimator_name} because it already has its maximum "
                f"number of children ({task_node.max_tasks})."
            )

        self.children_map[task_node.task_id] = task_node
        self.max_subtasks = task_node.max_tasks
        task_node.parent = self

    def _merge_with(self, task_node):
        # Set the parent of the sub-estimator's root task node to the parent
        # of the meta-estimator's leaf task node
        self.parent = task_node.parent
        self.task_id = task_node.task_id
        self.max_tasks = task_node.max_tasks
        task_node.parent.children_map[self.task_id] = self

        # Keep information about the node it was merged with
        self.prev_task_name = task_node.task_name
        self.prev_estimator_name = task_node.estimator_name

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
        for task_node in self.children_map.values():
            yield from task_node
