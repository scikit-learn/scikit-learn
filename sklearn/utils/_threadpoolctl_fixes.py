"""Compatibility fixes for older version of threadpoolctl

This module is separated from utils.fixes because it need to not be imported too soon
when importing sklearn but utils.fixes is.
"""

"""Compatibility fixes for older version of threadpoolctl"""

import threadpoolctl


if hasattr(threadpoolctl, "ThreadpoolController"):
    _sklearn_threadpool_controller = threadpoolctl.ThreadpoolController()
    threadpool_limits = _sklearn_threadpool_controller.limit
    threadpool_info = _sklearn_threadpool_controller.info
else:
    _sklearn_threadpool_controller = None
    threadpool_limits = threadpoolctl.threadpool_limits
    threadpool_info = threadpoolctl.threadpool_info
