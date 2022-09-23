# License: BSD 3 clause

import importlib
from threading import Thread, Event

from . import BaseCallback
from . import load_computation_tree


def _check_backend_support(backend, caller_name):
    """Raise ImportError with detailed error message if backend is not installed.

    Parameters
    ----------
    backend : {"rich", "tqdm"}
        The requested backend.

    caller_name : str
        The name of the caller that requires the backend.
    """
    try:
        importlib.import_module(backend)  # noqa
    except ImportError as e:
        raise ImportError(f"{caller_name} requires {backend} installed.") from e


class ProgressBar(BaseCallback):
    """Callback that displays progress bars for each iterative steps of the estimator

    Parameters
    ----------
    backend: {"rich"}, default="rich"
        The backend for the progress bars display.

    max_depth_show : int, default=None
        The maximum nested level of progress bars to display.

    max_depth_keep : int, default=None
        The maximum nested level of progress bars to keep displayed when they are
        finished.
    """

    auto_propagate = True

    def __init__(self, backend="rich", max_depth_show=None, max_depth_keep=None):
        if backend not in ("rich", "tqdm"):
            raise ValueError(f"backend should be 'rich' or 'tqdm', got {self.backend} instead.")
        _check_backend_support(backend, caller_name="Progressbar")
        self.backend = backend

        if max_depth_show is not None and max_depth_show < 0:
            raise ValueError(f"max_depth_show should be >= 0.")
        self.max_depth_show = max_depth_show

        if max_depth_keep is not None and max_depth_keep < 0:
            raise ValueError(f"max_depth_keep should be >= 0.")
        self.max_depth_keep = max_depth_keep

    def on_fit_begin(self, estimator, X=None, y=None):
        self._stop_event = Event()

        if self.backend == "rich":
            self.progress_monitor = _RichProgressMonitor(
                estimator=estimator,
                event=self._stop_event,
                max_depth_show=self.max_depth_show,
                max_depth_keep=self.max_depth_keep,
            )
        elif self.backend == "tqdm":
            self.progress_monitor = _TqdmProgressMonitor(
                estimator=estimator,
                event=self._stop_event,
            )

        self.progress_monitor.start()

    def on_fit_iter_end(self, *, estimator, node, **kwargs):
        pass

    def on_fit_end(self):
        self._stop_event.set()
        self.progress_monitor.join()

    def on_fit_exception(self):
        pass

    def __getstate__(self):
        state = self.__dict__.copy()
        if "_stop_event" in state:
            del state["_stop_event"]
        if "progress_monitor" in state:
            del state["progress_monitor"]
        return state


# Custom Progress class to allow showing the tasks in a given order (given by setting
# the _ordered_tasks attribute). In particular it allows to dynamically create and
# insert tasks between existing tasks.

try:
    from rich.progress import Progress
    class _Progress(Progress):
        def get_renderables(self):
            table = self.make_tasks_table(getattr(self, "_ordered_tasks", []))
            yield table
except:
    pass


class _RichProgressMonitor(Thread):
    """Thread monitoring the progress of an estimator with rich based display

    The display is a list of nested rich tasks using rich.Progress. There is one for
    each node in the computation tree of the estimator and in the computation trees of
    estimators used in the estimator.

    Parameters
    ----------
    estimator : estimator instance
        The estimator to monitor

    event : threading.Event instance
        This thread will run until event is set.

    max_depth_show : int, default=None
        The maximum nested level of progress bars to display.

    max_depth_keep : int, default=None
        The maximum nested level of progress bars to keep displayed when they are
        finished.
    """

    def __init__(self, estimator, event, max_depth_show=None, max_depth_keep=None):
        Thread.__init__(self)
        self.computation_tree = estimator._computation_tree
        self.event = event
        self.max_depth_show = max_depth_show
        self.max_depth_keep = max_depth_keep

        # _computation_trees is a dict `directory: tuple` where
        # - tuple[0] is the computation tree of the directory
        # - tuple[1] is a dict `node.tree_status_idx: task_id`
        self._computation_trees = {}

    def run(self):
        from rich.progress import BarColumn, TimeRemainingColumn, TextColumn
        from rich.style import Style

        with _Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(
                complete_style=Style(color="dark_orange"),
                finished_style=Style(color="cyan"),
            ),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            auto_refresh=False,
        ) as progress_ctx:
            self._progress_ctx = progress_ctx

            while not self.event.wait(0.05):
                self._recursive_update_tasks()
                self._progress_ctx.refresh()

            self._recursive_update_tasks()
            self._progress_ctx.refresh()

    def _recursive_update_tasks(self, this_dir=None, depth=0):
        """Recursively loop through directories and init or update tasks

        Parameters
        ----------
        this_dir : pathlib.Path instance
            The directory to

        depth : int
            The current depth
        """
        if self.max_depth_show is not None and depth > self.max_depth_show:
            # Fast exit if this dir is deeper than what we want to show anyway
            return

        if this_dir is None:
            this_dir = self.computation_tree.tree_dir
            # _ordered_tasks holds the list of the tasks in the order we want them to
            # be displayed.
            self._progress_ctx._ordered_tasks = []

        if this_dir not in self._computation_trees:
            # First time we discover this directory -> store the computation tree
            # If the computation tree is not readable yet, skip and try again next time
            computation_tree = load_computation_tree(this_dir)
            if computation_tree is None:
                return

            self._computation_trees[this_dir] = (computation_tree, {})

        computation_tree, task_ids = self._computation_trees[this_dir]

        for node in computation_tree.iterate(include_leaves=True):
            if node.children:
                # node is not a leaf, create or update its task
                if node.tree_status_idx not in task_ids:
                    visible = True
                    if (
                        self.max_depth_show is not None
                        and depth + node.depth > self.max_depth_show
                    ):
                        # If this node is deeper than what we want to show, we create
                        # the task anyway but make it not visible
                        visible = False

                    task_ids[node.tree_status_idx] = self._progress_ctx.add_task(
                        self._format_task_description(node, computation_tree, depth),
                        total=node.max_iter,
                        visible=visible,
                    )

                task_id = task_ids[node.tree_status_idx]
                task = self._progress_ctx.tasks[task_id]
                self._progress_ctx._ordered_tasks.append(task)

                parent_task = self._get_parent_task(node, computation_tree, task_ids)
                if parent_task is not None and parent_task.finished:
                    # If the task of the parent node is finished, make this task
                    # finished. It can happen if some computations are stopped
                    # before reaching max_iter.
                    visible = True
                    if (
                        self.max_depth_keep is not None
                        and depth + node.depth > self.max_depth_keep
                    ):
                        # If this node is deeper than what we want to keep in the output
                        # make it not visible
                        visible = False
                    self._progress_ctx.update(
                        task_id, completed=node.max_iter, visible=visible, refresh=False
                    )
                else:
                    node_progress = computation_tree.get_progress(node)
                    if node_progress != task.completed:
                        self._progress_ctx.update(
                            task_id, completed=node_progress, refresh=False
                        )
            else:
                # node is a leaf, look for tasks of its sub computation tree before
                # going to the next node
                child_dir = computation_tree.get_child_computation_tree_dir(node)
                # child_dir = this_dir / str(node.tree_status_idx)
                if child_dir.exists():
                    self._recursive_update_tasks(
                        child_dir, depth + computation_tree.depth
                    )

    def _format_task_description(self, node, computation_tree, depth):
        """Return a formatted description for the task of the node"""
        colors = ["red", "green", "blue", "yellow"]

        indent = f"{'  ' * (depth + node.depth)}"
        style = f"[{colors[(depth + node.depth)%len(colors)]}]"

        description = f"{computation_tree.estimator_name} - {node.description}"
        if node.parent is None and computation_tree.parent_node is not None:
            description = (
                f"{computation_tree.parent_node.description} {computation_tree.parent_node.idx} |"
                f" {description}"
            )
        if node.parent is not None:
            description = f"{description} {node.idx}"

        return f"{style}{indent}{description}"

    def _get_parent_task(self, node, computation_tree, task_ids):
        """Get the task of the parent node"""
        if node.parent is not None:
            # node is not the root, return the task of its parent
            task_id = task_ids[node.parent.tree_status_idx]
            return self._progress_ctx.tasks[task_id]
        if computation_tree.parent_node is not None:
            # node is the root, return the task of the parent of the parent_node of
            # its computation tree
            parent_dir = computation_tree.parent_node.computation_tree.tree_dir
            _, parent_tree_task_ids = self._computation_trees[parent_dir]
            task_id = parent_tree_task_ids[
                computation_tree.parent_node.parent.tree_status_idx
            ]
            return self._progress_ctx._tasks[task_id]
        return


class _TqdmProgressMonitor(Thread):
    def __init__(self, estimator, event):
        Thread.__init__(self)
        self.computation_tree = estimator._computation_tree
        self.event = event

    def run(self):
        from tqdm import tqdm

        root = self.computation_tree.root

        with tqdm(total=len(root.children)) as pbar:
            while not self.event.wait(0.05):
                node_progress = self.computation_tree.get_progress(root)
                if node_progress != pbar.total:
                    pbar.update(node_progress - pbar.n)

            pbar.update(pbar.total - pbar.n)
