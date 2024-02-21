# License: BSD 3 clause
# Authors: the scikit-learn developers

from multiprocessing import Manager
from threading import Thread

from ..utils._optional_dependencies import check_rich_support
from . import BaseCallback


class ProgressBar(BaseCallback):
    """Callback that displays progress bars for each iterative steps of an estimator."""

    auto_propagate = True

    def __init__(self):
        check_rich_support("Progressbar")

    def on_fit_begin(self, estimator, data):
        self._queue = Manager().Queue()
        self.progress_monitor = _RichProgressMonitor(queue=self._queue)
        self.progress_monitor.start()

    def on_fit_iter_end(self, *, estimator, node, **kwargs):
        self._queue.put(node)

    def on_fit_end(self):
        self._queue.put(None)
        self.progress_monitor.join()

    def __getstate__(self):
        state = self.__dict__.copy()
        if "progress_monitor" in state:
            del state["progress_monitor"]  # a thread is not picklable
        return state


try:
    from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn
    from rich.style import Style

    class _Progress(Progress):
        # Custom Progress class to allow showing the tasks in a given order (given by
        # setting the _ordered_tasks attribute). In particular it allows to dynamically
        # create and insert tasks between existing tasks.
        def get_renderables(self):
            table = self.make_tasks_table(getattr(self, "_ordered_tasks", []))
            yield table

except ImportError:
    pass


class _RichProgressMonitor(Thread):
    """Thread monitoring the progress of an estimator with rich based display.

    The display is a list of nested rich tasks using `rich.Progress`. There is one for
    each non-leaf node in the computation tree of the estimator.

    Parameters
    ----------
    queue : `multiprocessing.Manager.Queue` instance
        This thread will run until the queue is empty.
    """

    def __init__(self, *, queue):
        Thread.__init__(self)
        self.queue = queue

    def run(self):
        self.progress_ctx = _Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(
                complete_style=Style(color="dark_orange"),
                finished_style=Style(color="cyan"),
            ),
            TextColumn("[bright_magenta]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            auto_refresh=False,
        )

        # Holds the root of the tree of rich tasks (i.e. progress bars) that will be
        # created dynamically as the computation tree of the estimator is traversed.
        self.root_task = None

        with self.progress_ctx:
            while node := self.queue.get():
                self._update_task_tree(node)
                self._update_tasks()
                self.progress_ctx.refresh()

    def _update_task_tree(self, node):
        """Update the tree of tasks from a new node."""
        curr_task, parent_task = None, None

        for curr_node in node.path:
            if curr_node.parent is None:  # root node
                if self.root_task is None:
                    self.root_task = TaskNode(curr_node, progress_ctx=self.progress_ctx)
                curr_task = self.root_task
            elif curr_node.idx not in parent_task.children:
                curr_task = TaskNode(
                    curr_node, progress_ctx=self.progress_ctx, parent=parent_task
                )
                parent_task.children[curr_node.idx] = curr_task
            else:  # task already exists
                curr_task = parent_task.children[curr_node.idx]
            parent_task = curr_task

        # Mark the deepest task as finished (this is the one corresponding the
        # computation node that we just get from the queue).
        curr_task.finished = True

    def _update_tasks(self):
        """Loop through the tasks in their display order and update their progress."""
        self.progress_ctx._ordered_tasks = []

        for task_node in self.root_task:
            task = self.progress_ctx.tasks[task_node.task_id]

            if task_node.parent is not None and task_node.parent.finished:
                # If the parent task is finished, then mark the current task as
                # finished. It can happen if an estimator doesn't reach its max number
                # of iterations (e.g. early stopping).
                completed = task.total
            else:
                completed = sum(t.finished for t in task_node.children.values())

            if completed == task.total:
                task_node.finished = True

            self.progress_ctx.update(
                task_node.task_id, completed=completed, refresh=False
            )
            self.progress_ctx._ordered_tasks.append(task)


class TaskNode:
    """A node in the tree of rich tasks.

    Parameters
    ----------
    node : `ComputationNode` instance
        The computation node this task corresponds to.

    progress_ctx : `rich.Progress` instance
        The progress context to which this task belongs.

    parent : `TaskNode` instance
        The parent of this task.
    """

    def __init__(self, node, progress_ctx, parent=None):
        self.node_idx = node.idx
        self.parent = parent
        self.children = {}
        self.finished = False

        if node.n_children is not None:
            description = self._format_task_description(node)
            self.task_id = progress_ctx.add_task(description, total=node.n_children)

    def _format_task_description(self, node):
        """Return a formatted description for the task of the node."""
        colors = ["bright_magenta", "cyan", "dark_orange"]

        indent = f"{'  ' * (node.depth)}"
        style = f"[{colors[(node.depth)%len(colors)]}]"

        description = f"{node.estimator_name[0]} - {node.stage[0]}"
        if node.parent is not None:
            description += f" #{node.idx}"
        if len(node.estimator_name) == 2:
            description += f" | {node.estimator_name[1]} - {node.stage[1]}"

        return f"{style}{indent}{description}"

    def __iter__(self):
        """Pre-order depth-first traversal, excluding leaves."""
        if self.children:
            yield self
            for child in self.children.values():
                yield from child
