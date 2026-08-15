# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from datetime import timedelta
from numbers import Integral
from queue import Queue
from threading import Thread

from sklearn.callback._callback_context import get_context_path
from sklearn.callback._transport import close_listener, open_listener, send
from sklearn.utils._optional_dependencies import check_rich_support
from sklearn.utils._param_validation import Interval, validate_params

# Per-fit local queues and monitors, keyed by the run's `root_uuid`. Both are not
# picklable so they live here rather than on the callback instance. Entries are added in
# `setup` and removed in `teardown`.
_run_queues = {}
_run_monitors = {}


class ProgressBar:
    """Callback that displays progress bars for each iterative step of an estimator.

    Parameters
    ----------
    max_propagation_depth : int, default=1
        The maximum depth of nested levels of estimators to display progress bars for.
        0 means that the progress of only the outermost estimator is displayed.
        If set to None, all levels are displayed.
    """

    def __init__(self, max_propagation_depth=1):
        self.max_propagation_depth = max_propagation_depth

    def setup(self, estimator, context):
        check_rich_support("ProgressBar")

        if context.root_uuid not in _run_queues:
            _run_queues[context.root_uuid] = Queue()
            _run_monitors[context.root_uuid] = Thread(
                target=_monitor_progress,
                args=(_run_queues[context.root_uuid],),
                daemon=True,
            )
            _run_monitors[context.root_uuid].start()

        if not hasattr(self, "_listener_handles"):
            self._listener_handles = {}

        self._listener_handles[context.root_uuid] = open_listener(
            _run_queues[context.root_uuid]
        )

    def on_step(self, estimator, context):
        if (
            self.max_propagation_depth is not None
            and len(get_context_path(context)) > self.max_propagation_depth + 1
        ):
            return

        send(
            self._listener_handles[context.root_uuid],
            {
                "event": "step",
                "path": [ctx.task_id for ctx in get_context_path(context)],
            },
        )

    def teardown(self, estimator, context):
        # Fit is finished. Signal that the queue won't receive any more tasks, close
        # the monitor thread and the listener.
        queue = _run_queues.pop(context.root_uuid, None)
        if queue is not None:
            queue.put(None)
        monitor = _run_monitors.pop(context.root_uuid, None)
        if monitor is not None:
            monitor.join()
        if hasattr(self, "_listener_handles"):
            handle = self._listener_handles.pop(context.root_uuid, None)
            if handle is not None:
                close_listener(handle)
