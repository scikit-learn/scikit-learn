import threading


def flag_current_thread_clean_exit():
    """Put a ``_clean_exit`` flag on the current thread"""
    thread = threading.current_thread()
    thread._clean_exit = True


def is_crashed(thread):
    """Check if a thread has terminated unexpectedly for any version of python.

    The thread should be flagged using ``flag_current_thread_clean_exit()`` for
    all the cases where there was no crash. This permits to avoid false
    positive. If this flag is not set, this function return True whenever a
    thread is stopped.
    """
    if thread is None:
        return False
    if hasattr(thread, "_started"):
        terminate = thread._started.is_set() and not thread.is_alive()
    else:
        terminate = thread._Thread__started.is_set() and not thread.is_alive()
    return terminate and not getattr(thread, "_clean_exit", False)
