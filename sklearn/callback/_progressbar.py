# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import sys
import warnings
from threading import Thread

from sklearn.callback._callback_context import get_context_path
from sklearn.callback._callback_support import get_callback_manager
from sklearn.utils._optional_dependencies import check_rich_support


class ProgressBar:
    """Callback that displays progress bars for each iterative step of an estimator.

    Parameters
    ----------
    max_propagation_depth : int, default=1
        The maximum depth of nested levels of estimators to display progress bars for.
        0 means that the progress of only the outermost estimator is displayed.
        If set to None, all levels are displayed.

    Notes
    -----
    The use of this callback on python versions inferior to 3.12.8 might lead to
    unexpected crashes related to multiprocessing bugs on certain platforms.
    """

    def __init__(self, max_propagation_depth=1):
        if sys.version_info < (3, 12, 8):
            warnings.warn(
                "The use of the ProgressBar callback on python versions inferior "
                "to 3.12.8 might lead to unexpected crashes related to multiprocessing "
                "bugs on certain platforms."
            )

        check_rich_support("Progressbar")

        self.max_propagation_depth = max_propagation_depth
        self._manager = get_callback_manager()
        # Queue proxies need to be shared across callback copies in subprocesses,
        # while monitor threads must stay process-local (they are not picklable).
        self._run_queues = self._manager.dict()
        self._run_monitors = {}

    def setup(self, estimator, context):
        if not hasattr(self, "_manager"):
            # If the outer function call supports callback, it would typically have
            # initialized the manager and monitor in the same process and we can
            # directly reuse it. If the same callback is used to collect progress of
            # sub-estimators in subprocess parallel workers the setup/teardown is not
            # needed because it is performed only once, typically in the parent process.
            # However, if the outer function call does not support callbacks explicitly,
            # we need to reinitialize a working callback state in worker processes:
            # the callback will work in slightly degraded mode with redundant managers
            # and progress monitors but this should not crash.
            self._manager = get_callback_manager()
            self._run_queues = self._manager.dict()
            self._run_monitors = {}

        queue = self._manager.Queue()
        progress_monitor = RichProgressMonitor(queue=queue)
        progress_monitor.start()
        self._run_queues[context.root_uuid] = queue
        self._run_monitors[context.root_uuid] = progress_monitor

    def on_fit_task_begin(self, estimator, context):
        # A new progress bar is created at the beginning of each task that is not a
        # leaf, except if it's also the root task of an estimator.
        if (
            context.max_subtasks != 0
            or context.parent is None
            or context.source_estimator_name is not None
        ):
            # We pass the minimal information to the queue that is necessary to create a
            # progress bar and not the context to avoid pickling the whole context tree.
            path = [ctx.task_id for ctx in get_context_path(context)]
            self._run_queues[context.root_uuid].put(
                {
                    "event": "begin",
                    "path": path,
                    "task_name": context.task_name,
                    "task_id": context.task_id,
                    "max_subtasks": context.max_subtasks,
                    "estimator_name": context.estimator_name,
                    "source_estimator_name": context.source_estimator_name,
                    "source_task_name": context.source_task_name,
                }
            )

    def on_fit_task_end(self, estimator, context):
        # The path is enough to update the progress of the task and its ancestors.
        self._run_queues[context.root_uuid].put(
            {
                "event": "end",
                "path": [ctx.task_id for ctx in get_context_path(context)],
            }
        )

    def teardown(self, estimator, context):
        # Signal that the queue won't receive any more tasks.
        self._run_queues[context.root_uuid].put(None)
        self._run_monitors[context.root_uuid].join()

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop("_manager", None)
        state.pop("_run_monitors", None)
        # Note that run queues are pickleable and are expected to be shared between
        # the parent and worker processes.
        return state


try:
    from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn
    from rich.style import Style

    class _Progress(Progress):
        # Custom Progress class to allow showing the tasks in a given order (given by
        # setting the _ordered_tasks attribute). In particular it allows to dynamically
        # create and insert tasks between existing tasks.
        def get_renderables(self):
            yield self.make_tasks_table(getattr(self, "_ordered_tasks", []))

    class _StyledColumnMixin:
        """Apply finished/in-progress color style to rendered text."""

        def render(self, task):
            text = super().render(task)
            text.style = "#29ABE2" if task.finished else "#F7931E"
            return text

    class _StyledTimeRemainingColumn(_StyledColumnMixin, TimeRemainingColumn):
        """Time column with color styling."""

    class _StyledPercentageColumn(_StyledColumnMixin, TextColumn):
        """Percentage column with color styling."""

except ImportError:
    pass


class RichProgressMonitor(Thread):
    """Thread monitoring the progress of an estimator with rich based display.

    The display is a list of nested rich tasks using `rich.Progress`. There is one for
    each non-leaf node in the task tree of the estimator.

    Parameters
    ----------
    queue : `multiprocessing.Manager.Queue` instance
        This thread will run until the queue is empty.
    """

    def __init__(self, *, queue):
        super().__init__()
        self.queue = queue

    def run(self):
        self.progress_ctx = _Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(
                complete_style=Style(color="#F7931E"),
                finished_style=Style(color="#29ABE2"),
                pulse_style=Style(color="#F7931E"),
            ),
            _StyledPercentageColumn("{task.percentage:>3.0f}%"),
            _StyledTimeRemainingColumn(elapsed_when_finished=True),
        )

        # Holds the root of the tree of rich tasks (i.e. progress bars) that will be
        # created dynamically as the task tree of the estimator is traversed.
        self.root_rich_task = None

        with self.progress_ctx:
            while task_info := self.queue.get():
                if task_info.pop("event") == "begin":
                    self._on_task_begin(task_info)
                else:
                    self._on_task_end(task_info)

            self.progress_ctx.refresh()

    def _on_task_begin(self, task_info):
        """Create a progress bar for the task and update the list of ordered tasks."""
        path = task_info.pop("path")

        rich_task = RichTask(
            progress_ctx=self.progress_ctx, task_info=task_info, depth=len(path) - 1
        )
        if rich_task.depth == 0:
            self.root_rich_task = rich_task
        else:
            parent = self.root_rich_task.get_descendants(path)[-2]
            parent.children[path[-1]] = rich_task

        self.progress_ctx._ordered_tasks = [
            self.progress_ctx.tasks[task.id] for task in self.root_rich_task
        ]

    def _on_task_end(self, task_info):
        """Update the progress of the task and its ancestors recursively."""
        path = task_info.pop("path")
        *ancestors, task = self.root_rich_task.get_descendants(path)

        if task is None:
            # a leaf task of the estimator, no progress bar was created for it
            task = RichTask(self.progress_ctx, task_info, depth=len(ancestors))
            ancestors[-1].children[path[-1]] = task
        else:
            # Task is finished. Render the progress bar as 100% completed regardless of
            # its progress because a task may execute less tasks than its max_subtasks.
            task.completed = task.total
            self.progress_ctx.update(task.id, completed=1, total=1)

        for ancestor in reversed(ancestors):
            if ancestor.total is None:
                continue

            completed = ancestor.progress * ancestor.total
            self.progress_ctx.update(ancestor.id, completed=completed)


class RichTask:
    """A task, i.e. progressbar, in the tree of rich tasks.

    There is a rich task for each non-leaf task in the context tree of the estimator.

    Parameters
    ----------
    progress_ctx : `rich.Progress` instance
        The progress context to which this task belongs.

    task_info : dict
        Information about the task for which this rich task is created.

    depth : int
        The depth of the task in the tree of rich tasks.

    Attributes
    ----------
    completed : int
        The number of completed subtasks.

    total : int or None
        The total number of subtasks. None if the total number of subtasks is not known.

    progress : float
        The fraction, between 0 and 1, of the task that is completed.

    id : int
        The ID of the task in the Progress context.

    children : dict
        A mapping from the index of a child to the child node `{idx: RichTask}`.
        For a leaf, it's an empty dictionary.
    """

    def __init__(self, progress_ctx, task_info, *, depth):
        self.children = {}
        self.depth = depth
        self.completed = 0
        self.total = task_info.pop("max_subtasks", 0)

        if task_info:
            description = self._format_task_description(task_info)
            self.id = progress_ctx.add_task(description, total=self.total)

    @property
    def progress(self):
        """Return the fraction of the task that is completed."""
        if self.completed == self.total:
            return 1.0
        if self.total is None:
            return 0.0
        return sum(child.progress for child in self.children.values()) / self.total

    def _format_task_description(self, task_info):
        """Return a formatted description for the task."""
        indent = f"{'  ' * self.depth}"
        task_desc = f"{task_info['estimator_name']} - {task_info['task_name']}"
        id_mark = f" #{task_info['task_id']}" if self.depth > 0 else ""
        source_task_desc = (
            f"{task_info['source_estimator_name']} - {task_info['source_task_name']} | "
            if task_info["source_estimator_name"] is not None
            else ""
        )
        return f"{indent}{source_task_desc}{task_desc}{id_mark}"

    def get_descendants(self, path):
        """Return the descendants from this task along the given path."""
        if len(path) == 1:
            return [self]
        if path[1] not in self.children:
            return [self, None]
        return [self] + self.children[path[1]].get_descendants(path[1:])

    def __iter__(self):
        """Pre-order depth-first traversal, only including tasks with a progress bar."""
        if hasattr(self, "id"):
            yield self
            for child in self.children.values():
                yield from child
