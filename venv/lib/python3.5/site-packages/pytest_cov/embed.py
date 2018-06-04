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
import os
import signal

active_cov = None


def multiprocessing_start(_):
    cov = init()
    if cov:
        multiprocessing.util.Finalize(None, cleanup, args=(cov,), exitpriority=1000)


try:
    import multiprocessing.util
except ImportError:
    pass
else:
    multiprocessing.util.register_after_fork(multiprocessing_start, multiprocessing_start)


def init():
    # Only continue if ancestor process has set everything needed in
    # the env.
    global active_cov

    cov_source = os.environ.get('COV_CORE_SOURCE')
    cov_config = os.environ.get('COV_CORE_CONFIG')
    cov_datafile = os.environ.get('COV_CORE_DATAFILE')
    cov_branch = True if os.environ.get('COV_CORE_BRANCH') == 'enabled' else None

    if cov_datafile:
        # Import what we need to activate coverage.
        import coverage

        # Determine all source roots.
        if not cov_source:
            cov_source = None
        else:
            cov_source = cov_source.split(os.pathsep)
        if not cov_config:
            cov_config = True

        # Activate coverage for this process.
        cov = active_cov = coverage.coverage(
            source=cov_source,
            branch=cov_branch,
            data_suffix=True,
            config_file=cov_config,
            auto_data=True,
            data_file=cov_datafile
        )
        cov.load()
        cov.start()
        cov._warn_no_data = False
        cov._warn_unimported_source = False
        return cov


def _cleanup(cov):
    if cov is not None:
        cov.stop()
        cov.save()


def cleanup(cov=None):
    global active_cov

    _cleanup(cov)
    if active_cov is not cov:
        _cleanup(active_cov)
    active_cov = None


multiprocessing_finish = cleanup  # in case someone dared to use this internal


def _sigterm_handler(*_):
    cleanup()


def cleanup_on_sigterm():
    signal.signal(signal.SIGTERM, _sigterm_handler)
