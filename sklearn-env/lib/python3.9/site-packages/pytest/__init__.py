# PYTHON_ARGCOMPLETE_OK
"""pytest: unit and functional testing with Python."""
from . import collect
from _pytest import __version__
from _pytest.assertion import register_assert_rewrite
from _pytest.cacheprovider import Cache
from _pytest.capture import CaptureFixture
from _pytest.config import cmdline
from _pytest.config import console_main
from _pytest.config import ExitCode
from _pytest.config import hookimpl
from _pytest.config import hookspec
from _pytest.config import main
from _pytest.config import UsageError
from _pytest.debugging import pytestPDB as __pytestPDB
from _pytest.fixtures import _fillfuncargs
from _pytest.fixtures import fixture
from _pytest.fixtures import FixtureLookupError
from _pytest.fixtures import FixtureRequest
from _pytest.fixtures import yield_fixture
from _pytest.freeze_support import freeze_includes
from _pytest.logging import LogCaptureFixture
from _pytest.main import Session
from _pytest.mark import MARK_GEN as mark
from _pytest.mark import param
from _pytest.monkeypatch import MonkeyPatch
from _pytest.nodes import Collector
from _pytest.nodes import File
from _pytest.nodes import Item
from _pytest.outcomes import exit
from _pytest.outcomes import fail
from _pytest.outcomes import importorskip
from _pytest.outcomes import skip
from _pytest.outcomes import xfail
from _pytest.pytester import Pytester
from _pytest.pytester import Testdir
from _pytest.python import Class
from _pytest.python import Function
from _pytest.python import Instance
from _pytest.python import Module
from _pytest.python import Package
from _pytest.python_api import approx
from _pytest.python_api import raises
from _pytest.recwarn import deprecated_call
from _pytest.recwarn import WarningsRecorder
from _pytest.recwarn import warns
from _pytest.tmpdir import TempdirFactory
from _pytest.tmpdir import TempPathFactory
from _pytest.warning_types import PytestAssertRewriteWarning
from _pytest.warning_types import PytestCacheWarning
from _pytest.warning_types import PytestCollectionWarning
from _pytest.warning_types import PytestConfigWarning
from _pytest.warning_types import PytestDeprecationWarning
from _pytest.warning_types import PytestExperimentalApiWarning
from _pytest.warning_types import PytestUnhandledCoroutineWarning
from _pytest.warning_types import PytestUnhandledThreadExceptionWarning
from _pytest.warning_types import PytestUnknownMarkWarning
from _pytest.warning_types import PytestUnraisableExceptionWarning
from _pytest.warning_types import PytestWarning

set_trace = __pytestPDB.set_trace

__all__ = [
    "__version__",
    "_fillfuncargs",
    "approx",
    "Cache",
    "CaptureFixture",
    "Class",
    "cmdline",
    "collect",
    "Collector",
    "console_main",
    "deprecated_call",
    "exit",
    "ExitCode",
    "fail",
    "File",
    "fixture",
    "FixtureLookupError",
    "FixtureRequest",
    "freeze_includes",
    "Function",
    "hookimpl",
    "hookspec",
    "importorskip",
    "Instance",
    "Item",
    "LogCaptureFixture",
    "main",
    "mark",
    "Module",
    "MonkeyPatch",
    "Package",
    "param",
    "PytestAssertRewriteWarning",
    "PytestCacheWarning",
    "PytestCollectionWarning",
    "PytestConfigWarning",
    "PytestDeprecationWarning",
    "PytestExperimentalApiWarning",
    "Pytester",
    "PytestUnhandledCoroutineWarning",
    "PytestUnhandledThreadExceptionWarning",
    "PytestUnknownMarkWarning",
    "PytestUnraisableExceptionWarning",
    "PytestWarning",
    "raises",
    "register_assert_rewrite",
    "Session",
    "set_trace",
    "skip",
    "TempPathFactory",
    "Testdir",
    "TempdirFactory",
    "UsageError",
    "WarningsRecorder",
    "warns",
    "xfail",
    "yield_fixture",
]
