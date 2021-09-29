"""Compatibility fixes for older version of threadpoolctl

This module is separated from utils.fixes because it need to not be imported too soon
when importing sklearn but utils.fixes is.
"""

import threadpoolctl


if hasattr(threadpoolctl, "ThreadpoolController"):
    _sklearn_threadpool_controller = threadpoolctl.ThreadpoolController()
    threadpool_limits = _sklearn_threadpool_controller.limit
else:
    _sklearn_threadpool_controller = None
    threadpool_limits = threadpoolctl.threadpool_limits
