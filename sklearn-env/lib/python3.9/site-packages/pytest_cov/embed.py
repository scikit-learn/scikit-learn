"""Activate coverage at python startup if appropriate.

The python site initialisation will ensure that anything we import
will be removed and not visible at the end of python startup.  However
we minimise all work by putting these init actions in this separate
module and only importing what is needed when needed.

For normal python startup when coverage should not be activated the pth
file checks a single env var and does not import or call the init fn
here.

For python startup when an ancestor process has set the env indicating
that code coverage is being collected we activate coverage based on
info passed via env vars.
"""
import atexit
import os
import signal

_active_cov = None


def multiprocessing_start(_):
    global _active_cov
    cov = init()
    if cov:
        _active_cov = cov
        multiprocessing.util.Finalize(None, cleanup, exitpriority=1000)


try:
    import multiprocessing.util
except ImportError:
    pass
else:
    multiprocessing.util.register_after_fork(multiprocessing_start, multiprocessing_start)


def init():
    # Only continue if ancestor process has set everything needed in
    # the env.
    global _active_cov

    cov_source = os.environ.get('COV_CORE_SOURCE')
    cov_config = os.environ.get('COV_CORE_CONFIG')
    cov_datafile = os.environ.get('COV_CORE_DATAFILE')
    cov_branch = True if os.environ.get('COV_CORE_BRANCH') == 'enabled' else None
    cov_context = os.environ.get('COV_CORE_CONTEXT')

    if cov_datafile:
        if _active_cov:
            cleanup()
        # Import what we need to activate coverage.
        import coverage

        # Determine all source roots.
        if cov_source in os.pathsep:
            cov_source = None
        else:
            cov_source = cov_source.split(os.pathsep)
        if cov_config == os.pathsep:
            cov_config = True

        # Activate coverage for this process.
        cov = _active_cov = coverage.Coverage(
            source=cov_source,
            branch=cov_branch,
            data_suffix=True,
            config_file=cov_config,
            auto_data=True,
            data_file=cov_datafile
        )
        cov.load()
        cov.start()
        if cov_context:
            cov.switch_context(cov_context)
        cov._warn_no_data = False
        cov._warn_unimported_source = False
        return cov


def _cleanup(cov):
    if cov is not None:
        cov.stop()
        cov.save()
        cov._auto_save = False  # prevent autosaving from cov._atexit in case the interpreter lacks atexit.unregister
        try:
            atexit.unregister(cov._atexit)
        except Exception:
            pass


def cleanup():
    global _active_cov
    global _cleanup_in_progress
    global _pending_signal

    _cleanup_in_progress = True
    _cleanup(_active_cov)
    _active_cov = None
    _cleanup_in_progress = False
    if _pending_signal:
        pending_singal = _pending_signal
        _pending_signal = None
        _signal_cleanup_handler(*pending_singal)


multiprocessing_finish = cleanup  # in case someone dared to use this internal

_previous_handlers = {}
_pending_signal = None
_cleanup_in_progress = False


def _signal_cleanup_handler(signum, frame):
    global _pending_signal
    if _cleanup_in_progress:
        _pending_signal = signum, frame
        return
    cleanup()
    _previous_handler = _previous_handlers.get(signum)
    if _previous_handler == signal.SIG_IGN:
        return
    elif _previous_handler and _previous_handler is not _signal_cleanup_handler:
        _previous_handler(signum, frame)
    elif signum == signal.SIGTERM:
        os._exit(128 + signum)
    elif signum == signal.SIGINT:
        raise KeyboardInterrupt()


def cleanup_on_signal(signum):
    previous = signal.getsignal(signum)
    if previous is not _signal_cleanup_handler:
        _previous_handlers[signum] = previous
        signal.signal(signum, _signal_cleanup_handler)


def cleanup_on_sigterm():
    cleanup_on_signal(signal.SIGTERM)
