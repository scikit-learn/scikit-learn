"""Sphinx test fixtures for pytest"""

from __future__ import annotations

import shutil
import subprocess
import sys
from collections import namedtuple
from io import StringIO
from typing import TYPE_CHECKING

import pytest

from sphinx.testing.util import SphinxTestApp, SphinxTestAppWrapperForSkipBuilding

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from pathlib import Path
    from typing import Any

DEFAULT_ENABLED_MARKERS = [
    # The marker signature differs from the constructor signature
    # since the way it is processed assumes keyword arguments for
    # the 'testroot' and 'srcdir'.
    (
        'sphinx('
        'buildername="html", *, '
        'testroot="root", srcdir=None, '
        'confoverrides=None, freshenv=False, '
        'warningiserror=False, tags=None, verbosity=0, parallel=0, '
        'builddir=None, docutils_conf=None'
        '): arguments to initialize the sphinx test application.'
    ),
    'test_params(shared_result=...): test parameters.',
]


def pytest_configure(config: pytest.Config) -> None:
    """Register custom markers"""
    for marker in DEFAULT_ENABLED_MARKERS:
        config.addinivalue_line('markers', marker)


@pytest.fixture(scope='session')
def rootdir() -> Path | None:
    return None


class SharedResult:
    cache: dict[str, dict[str, str]] = {}

    def store(self, key: str, app_: SphinxTestApp) -> Any:
        if key in self.cache:
            return
        data = {
            'status': app_.status.getvalue(),
            'warning': app_.warning.getvalue(),
        }
        self.cache[key] = data

    def restore(self, key: str) -> dict[str, StringIO]:
        if key not in self.cache:
            return {}
        data = self.cache[key]
        return {
            'status': StringIO(data['status']),
            'warning': StringIO(data['warning']),
        }


@pytest.fixture
def app_params(
    request: Any,
    test_params: dict[str, Any],
    shared_result: SharedResult,
    sphinx_test_tempdir: str,
    rootdir: Path,
) -> _app_params:
    """
    Parameters that are specified by 'pytest.mark.sphinx' for
    sphinx.application.Sphinx initialization
    """
    # ##### process pytest.mark.sphinx

    pargs: dict[int, Any] = {}
    kwargs: dict[str, Any] = {}

    # to avoid stacking positional args
    for info in reversed(list(request.node.iter_markers('sphinx'))):
        pargs |= dict(enumerate(info.args))
        kwargs.update(info.kwargs)

    args = [pargs[i] for i in sorted(pargs.keys())]

    # ##### process pytest.mark.test_params
    if test_params['shared_result']:
        if 'srcdir' in kwargs:
            msg = 'You can not specify shared_result and srcdir in same time.'
            pytest.fail(msg)
        kwargs['srcdir'] = test_params['shared_result']
        restore = shared_result.restore(test_params['shared_result'])
        kwargs.update(restore)

    # ##### prepare Application params

    testroot = kwargs.pop('testroot', 'root')
    kwargs['srcdir'] = srcdir = sphinx_test_tempdir / kwargs.get('srcdir', testroot)

    # special support for sphinx/tests
    if rootdir and not srcdir.exists():
        testroot_path = rootdir / ('test-' + testroot)
        shutil.copytree(testroot_path, srcdir)

    return _app_params(args, kwargs)


_app_params = namedtuple('_app_params', 'args,kwargs')


@pytest.fixture
def test_params(request: Any) -> dict[str, Any]:
    """
    Test parameters that are specified by 'pytest.mark.test_params'

    :param Union[str] shared_result:
       If the value is provided, app._status and app._warning objects will be
       shared in the parametrized test functions and/or test functions that
       have same 'shared_result' value.
       **NOTE**: You can not specify both shared_result and srcdir.
    """
    env = request.node.get_closest_marker('test_params')
    kwargs = env.kwargs if env else {}
    result = {
        'shared_result': None,
    }
    result.update(kwargs)

    if result['shared_result'] and not isinstance(result['shared_result'], str):
        msg = 'You can only provide a string type of value for "shared_result"'
        raise pytest.Exception(msg)
    return result


@pytest.fixture
def app(
    test_params: dict[str, Any],
    app_params: _app_params,
    make_app: Callable[[], SphinxTestApp],
    shared_result: SharedResult,
) -> Iterator[SphinxTestApp]:
    """
    Provides the 'sphinx.application.Sphinx' object
    """
    args, kwargs = app_params
    app_ = make_app(*args, **kwargs)
    yield app_

    print('# testroot:', kwargs.get('testroot', 'root'))
    print('# builder:', app_.builder.name)
    print('# srcdir:', app_.srcdir)
    print('# outdir:', app_.outdir)
    print('# status:', '\n' + app_.status.getvalue())
    print('# warning:', '\n' + app_.warning.getvalue())

    if test_params['shared_result']:
        shared_result.store(test_params['shared_result'], app_)


@pytest.fixture
def status(app: SphinxTestApp) -> StringIO:
    """
    Back-compatibility for testing with previous @with_app decorator
    """
    return app.status


@pytest.fixture
def warning(app: SphinxTestApp) -> StringIO:
    """
    Back-compatibility for testing with previous @with_app decorator
    """
    return app.warning


@pytest.fixture
def make_app(test_params: dict[str, Any]) -> Iterator[Callable[[], SphinxTestApp]]:
    """
    Provides make_app function to initialize SphinxTestApp instance.
    if you want to initialize 'app' in your test function. please use this
    instead of using SphinxTestApp class directory.
    """
    apps = []
    syspath = sys.path.copy()

    def make(*args: Any, **kwargs: Any) -> SphinxTestApp:
        status, warning = StringIO(), StringIO()
        kwargs.setdefault('status', status)
        kwargs.setdefault('warning', warning)
        app_: SphinxTestApp
        if test_params['shared_result']:
            app_ = SphinxTestAppWrapperForSkipBuilding(*args, **kwargs)
        else:
            app_ = SphinxTestApp(*args, **kwargs)
        apps.append(app_)
        return app_

    yield make

    sys.path[:] = syspath
    for app_ in reversed(apps):  # clean up applications from the new ones
        app_.cleanup()


@pytest.fixture
def shared_result() -> SharedResult:
    return SharedResult()


@pytest.fixture(scope='module', autouse=True)
def _shared_result_cache() -> None:
    SharedResult.cache.clear()


@pytest.fixture
def if_graphviz_found(app: SphinxTestApp) -> None:  # NoQA: PT004
    """
    The test will be skipped when using 'if_graphviz_found' fixture and graphviz
    dot command is not found.
    """
    graphviz_dot = getattr(app.config, 'graphviz_dot', '')
    try:
        if graphviz_dot:
            # print the graphviz_dot version, to check that the binary is available
            subprocess.run([graphviz_dot, '-V'], capture_output=True, check=False)
            return
    except OSError:  # No such file or directory
        pass

    pytest.skip('graphviz "dot" is not available')


@pytest.fixture(scope='session')
def sphinx_test_tempdir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Temporary directory."""
    return tmp_path_factory.getbasetemp()


@pytest.fixture
def rollback_sysmodules() -> Iterator[None]:  # NoQA: PT004
    """
    Rollback sys.modules to its value before testing to unload modules
    during tests.

    For example, used in test_ext_autosummary.py to permit unloading the
    target module to clear its cache.
    """
    sysmodules = list(sys.modules)
    try:
        yield
    finally:
        for modname in list(sys.modules):
            if modname not in sysmodules:
                sys.modules.pop(modname)
