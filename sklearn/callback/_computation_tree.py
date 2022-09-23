# License: BSD 3 clause

import os
from pathlib import Path
import pickle
from tempfile import mkdtemp
from uuid import uuid4

import numpy as np


class ComputationNode:
    """A node in a ComputationTree

    Parameters
    ----------
    computation_tree : ComputationTree instance
        The computation tree it belongs to.

    parent : ComputationNode instance, default=None
        The parent node. None means this is the root.

    max_iter : int, default=None
        The number of its children. None means it's a leaf.

    description : str, default=None
        A description of this computation node. None means it's a leaf.

    tree_status_idx : int, default=0
        The index of the status of this node in the `tree_status` array of its
        computation tree.

    idx : int, default=0
        The index of this node in the children list of its parent.

    Attributes
    ----------
    children : list
        The list of its children nodes. For a leaf, it's an empty list

    depth : int
        The depth of this node in its computation tree. The root has a depth of 0.
    """

    def __init__(
        self,
        computation_tree,
        parent=None,
        max_iter=None,
        description=None,
        tree_status_idx=0,
        idx=0,
    ):
        self.computation_tree = computation_tree
        self.parent = parent
        self.max_iter = max_iter
        self.description = description
        self.tree_status_idx = tree_status_idx
        self.idx = idx
        self.children = []
        self.depth = 0 if self.parent is None else self.parent.depth + 1

    def get_ancestors(self, include_ancestor_trees=True):
        """Get the list of all nodes in the path from the node to the root

        Parameters
        ----------
        include_ancestor_trees : bool, default=True
            If True, propagate to the tree of the `parent_node` of this tree if it
            exists and so on.

        Returns
        -------
        ancestors : list
            The list of ancestors of this node (included).
        """
        node = self
        ancestors = [node]

        while node.parent is not None:
            node = node.parent
            ancestors.append(node)

        if include_ancestor_trees:
            node_parent_tree = node.computation_tree.parent_node
            if node_parent_tree is not None:
                ancestors.extend(node_parent_tree.get_ancestors())

        return ancestors

    def __repr__(self):
        return (
            f"ComputationNode(description={self.description}, "
            f"depth={self.depth}, idx={self.idx})"
        )


class ComputationTree:
    """Data structure to store the computation tree of an estimator

    Parameters
    ----------
    estimator_name : str
        The name of the estimator.

    levels : list of dict
        A description of the nested levels of computation of the estimator to build the
        tree. It's a list of dict with "descr" and "max_iter" keys.

    parent_node : ComputationNode, default=None
        The node where the estimator is used in the computation tree of a
        meta-estimator. This node is not set to be the parent of the root of this tree.

    Attributes
    ----------
    depth : int
        The depth of the tree. It corresponds to the depth of its deepest leaf.

    root : ComputationNode instance
        The root of the computation tree.

    tree_dir : pathlib.Path instance
        The path of the directory where the computation tree is dumped during the fit of
        its estimator. If it has a parent tree, this is a sub-directory of the
        `tree_dir` of its parent.

    uid : uuid.UUID
        Unique indentifier for a ComputationTree instance.
    """

    def __init__(self, estimator_name, levels, *, parent_node=None):
        self.estimator_name = estimator_name
        self.parent_node = parent_node

        self.depth = len(levels) - 1
        self.root, self.n_nodes = self._build_tree(levels)

        self.uid = uuid4()

        parent_tree_dir = (
            None
            if self.parent_node is None
            else self.parent_node.computation_tree.tree_dir
        )
        if parent_tree_dir is None:
            self.tree_dir = Path(mkdtemp())
        else:
            # This tree has a parent tree. Place it in a subdir of its parent dir
            # and give it a name that allows from the parent tree to find the sub dir
            # of the sub tree of a given leaf.
            self.tree_dir = parent_tree_dir / str(parent_node.tree_status_idx)
            self.tree_dir.mkdir()
        self._filename = self.tree_dir / "tree_status.memmap"

        self._set_tree_status(mode="w+")
        self._tree_status[:] = False

    def _build_tree(self, levels):
        """Build the computation tree from the description of the levels"""
        root = ComputationNode(
            computation_tree=self,
            max_iter=levels[0]["max_iter"],
            description=levels[0]["descr"],
        )

        n_nodes = self._recursive_build_tree(root, levels)

        return root, n_nodes

    def _recursive_build_tree(self, parent, levels, n_nodes=1):
        """Recursively build the tree from the root the leaves"""
        if parent.depth == self.depth:
            return n_nodes

        for i in range(parent.max_iter):
            children_max_iter = levels[parent.depth + 1]["max_iter"]
            description = levels[parent.depth + 1]["descr"]

            node = ComputationNode(
                computation_tree=self,
                parent=parent,
                max_iter=children_max_iter,
                description=description,
                tree_status_idx=n_nodes,
                idx=i,
            )
            parent.children.append(node)

            n_nodes = self._recursive_build_tree(node, levels, n_nodes + 1)

        return n_nodes

    def _set_tree_status(self, mode):
        """Create a memory-map to the tree_status array stored on the disk"""
        # This has to be done each time we unpickle the tree
        self._tree_status = np.memmap(
            self._filename, dtype=bool, mode=mode, shape=(self.n_nodes,)
        )

    def get_progress(self, node):
        """Return the number of finished child nodes of this node"""
        if self._tree_status[node.tree_status_idx]:
            return node.max_iter

        # Since the children of a node are not ordered (to account for parallel
        # execution), we can't rely on the highest index for which the status is True.
        return sum(
            [self._tree_status[child.tree_status_idx] for child in node.children]
        )

    def get_child_computation_tree_dir(self, node):
        if node.children:
            raise ValueError("node is not a leaf")
        return self.tree_dir / str(node.tree_status_idx)

    def iterate(self, include_leaves=False):
        """Return an iterable over the nodes of the computation tree

        Nodes are discovered in a depth first search manner.

        Parameters
        ----------
        include_leaves : bool
            Whether or not to include the leaves of the tree in the iterable

        Returns
        -------
        nodes_list : list
            A list of the nodes of the computation tree.
        """
        return self._recursive_iterate(include_leaves=include_leaves)

    def _recursive_iterate(self, node=None, include_leaves=False, node_list=None):
        """Recursively constructs the iterable"""
        # TODO make it an iterator ?
        if node is None:
            node = self.root
            node_list = []

        if node.children or include_leaves:
            node_list.append(node)

        for child in node.children:
            self._recursive_iterate(child, include_leaves, node_list)

        return node_list

    def __repr__(self):
        res = (
            f"[{self.estimator_name}] {self.root.description} : progress "
            f"{self.get_progress(self.root)} / {self.root.max_iter}\n"
        )
        for node in self.iterate(include_leaves=False):
            if node is not self.root:
                res += (
                    f"{'  ' * node.depth}{node.description} {node.idx}: progress "
                    f"{self.get_progress(node)} / {node.max_iter}\n"
                )
        return res


def load_computation_tree(directory):
    """load the computation tree of a directory

    Parameters
    ----------
    directory : pathlib.Path instance
        The directory where the computation tree is dumped

    Returns
    -------
    computation_tree : ComputationTree instance
        The loaded computation tree
    """
    file_path = directory / "computation_tree.pkl"
    if not file_path.exists() or not os.path.getsize(file_path) > 0:
        # Do not try to load the tree when it's created but not yet written
        return

    with open(file_path, "rb") as f:
        computation_tree = pickle.load(f)

    computation_tree._set_tree_status(mode="r")

    return computation_tree
