# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from multiprocessing import Manager
from threading import Thread

from sklearn.callback._callback_context import get_task_info_path
from sklearn.utils._optional_dependencies import check_rich_support


class ProgressBar:
    """Callback that displays progress bars for each iterative steps of an estimator.

    Parameters
    ----------
    max_estimator_depth : int, default=1
        The maximum number of nested levels of estimators to display progress bars for.
        By default, only the progress bars of the outermost estimator are displayed.
        If set to None, all levels are displayed.
    """

    def __init__(self, max_estimator_depth=1):
        check_rich_support("Progressbar")

        self.max_estimator_depth = max_estimator_depth

    def _on_fit_begin(self, estimator):
        self._queue = Manager().Queue()
        self.progress_monitor = RichProgressMonitor(queue=self._queue)
        self.progress_monitor.start()

    def _on_fit_task_end(self, estimator, task_info, **kwargs):
        self._queue.put(task_info)

    def _on_fit_end(self, estimator, task_info):
        self._queue.put(task_info)
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
            TimeRemainingColumn(elapsed_when_finished=True),
            auto_refresh=False,
        )

        # Holds the root of the tree of rich tasks (i.e. progress bars) that will be
        # created dynamically as the computation tree of the estimator is traversed.
        self.root_rich_task = None

        with self.progress_ctx:
            while task_info := self.queue.get():
                task_info_path = get_task_info_path(task_info)
                self._update_task_tree(task_info_path)
                self._update_tasks()
                self.progress_ctx.refresh()

    def _update_task_tree(self, task_info_path):
        """Update the tree of rich tasks from the path of a new task.

        A new rich task is created for the task and all its ancestors if needed.
        """
        curr_rich_task, parent_rich_task = None, None

        for task_info in task_info_path:
            if task_info["parent_task_info"] is None:  # root node
                if self.root_rich_task is None:
                    self.root_rich_task = RichTask(
                        task_info, progress_ctx=self.progress_ctx
                    )
                curr_rich_task = self.root_rich_task
            elif task_info["task_id"] not in parent_rich_task.children:
                curr_rich_task = RichTask(
                    task_info, progress_ctx=self.progress_ctx, parent=parent_rich_task
                )
                parent_rich_task.children[task_info["task_id"]] = curr_rich_task
            else:  # task already exists
                curr_rich_task = parent_rich_task.children[task_info["task_id"]]
            parent_rich_task = curr_rich_task

        # Mark the deepest task as finished (this is the one corresponding to the
        # task that we just get from the queue).
        curr_rich_task.finished = True

    def _update_tasks(self):
        """Loop through the tasks in their display order and update their progress."""
        self.progress_ctx._ordered_tasks = []

        for rich_task_node in self.root_rich_task:
            task = self.progress_ctx.tasks[rich_task_node.task_id]

            total = task.total

            if rich_task_node.finished:
                # It's possible that a task finishes without reaching its total
                # (e.g. early stopping). We mark it as 100% completed.

                if task.total is None:
                    # Indeterminate task is finished. Set total to an arbitrary
                    # value to render its completion as 100%.
                    completed = total = 1
                else:
                    completed = total
            else:
                completed = sum(t.finished for t in rich_task_node.children.values())

            self.progress_ctx.update(
                rich_task_node.task_id, completed=completed, total=total, refresh=False
            )
            self.progress_ctx._ordered_tasks.append(task)


class RichTask:
    """A task, i.e. progressbar, in the tree of rich tasks.

    Parameters
    ----------
    task_info : dict
        Available information about the estimator task for which this rich task is
        created. See :meth:`~sklearn.callback.CallbackContext.task_info` for a detailed
        description of the keys of this dictionary.

    progress_ctx : `rich.Progress` instance
        The progress context to which this task belongs.

    parent : `RichTask` instance
        The parent of this task.

    Attributes
    ----------
    finished : bool
        Whether the task is finished.

    task_id : int
        The ID of the task in the Progress context.

    children : dict
        A mapping from the index of a child to the child node `{idx: RichTask}`.
        For a leaf, it's an empty dictionary.
    """

    def __init__(self, task_info, progress_ctx, parent=None):
        self.parent = parent
        self.children = {}
        self.finished = False
        self.depth = 0 if parent is None else parent.depth + 1

        if task_info["max_subtasks"] != 0:
            description = self._format_task_description(task_info)
            self.task_id = progress_ctx.add_task(
                description, total=task_info["max_subtasks"]
            )

    def _format_task_description(self, task_info):
        """Return a formatted description for the task."""
        colors = ["bright_magenta", "cyan", "dark_orange"]

        indent = f"{'  ' * self.depth}"
        style = f"[{colors[(self.depth) % len(colors)]}]"

        task_desc = f"{task_info['estimator_name']} - {task_info['task_name']}"
        id_mark = (
            f" #{task_info['task_id']}"
            if task_info["parent_task_info"] is not None
            else ""
        )
        prev_task_desc = (
            f"{task_info['prev_estimator_name']} - {task_info['prev_task_name']} | "
            if task_info["prev_estimator_name"] is not None
            else ""
        )

        return f"{style}{indent}{prev_task_desc}{task_desc}{id_mark}"

    def __iter__(self):
        """Pre-order depth-first traversal, excluding leaves."""
        if self.children:
            yield self
            for child in self.children.values():
                yield from child
