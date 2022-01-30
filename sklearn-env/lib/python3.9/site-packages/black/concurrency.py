import asyncio
import logging
import sys
from typing import Any, Iterable

from black.output import err


def maybe_install_uvloop() -> None:
    """If our environment has uvloop installed we use it.

    This is called only from command-line entry points to avoid
    interfering with the parent process if Black is used as a library.

    """
    try:
        import uvloop

        uvloop.install()
    except ImportError:
        pass


def cancel(tasks: Iterable["asyncio.Task[Any]"]) -> None:
    """asyncio signal handler that cancels all `tasks` and reports to stderr."""
    err("Aborted!")
    for task in tasks:
        task.cancel()


def shutdown(loop: asyncio.AbstractEventLoop) -> None:
    """Cancel all pending tasks on `loop`, wait for them, and close the loop."""
    try:
        if sys.version_info[:2] >= (3, 7):
            all_tasks = asyncio.all_tasks
        else:
            all_tasks = asyncio.Task.all_tasks
        # This part is borrowed from asyncio/runners.py in Python 3.7b2.
        to_cancel = [task for task in all_tasks(loop) if not task.done()]
        if not to_cancel:
            return

        for task in to_cancel:
            task.cancel()
        loop.run_until_complete(
            asyncio.gather(*to_cancel, loop=loop, return_exceptions=True)
        )
    finally:
        # `concurrent.futures.Future` objects cannot be cancelled once they
        # are already running. There might be some when the `shutdown()` happened.
        # Silence their logger's spew about the event loop being closed.
        cf_logger = logging.getLogger("concurrent.futures")
        cf_logger.setLevel(logging.CRITICAL)
        loop.close()
