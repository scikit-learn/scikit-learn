"""Deprecation messages and bits of code used elsewhere in the codebase that
is planned to be removed in the next pytest release.

Keeping it in a central location makes it easy to track what is deprecated and should
be removed when the time comes.

All constants defined in this module should be either instances of
:class:`PytestWarning`, or :class:`UnformattedWarning`
in case of warnings which need to format their messages.
"""
from warnings import warn

from _pytest.warning_types import PytestDeprecationWarning
from _pytest.warning_types import UnformattedWarning

# set of plugins which have been integrated into the core; we use this list to ignore
# them during registration to avoid conflicts
DEPRECATED_EXTERNAL_PLUGINS = {
    "pytest_catchlog",
    "pytest_capturelog",
    "pytest_faulthandler",
}


FILLFUNCARGS = UnformattedWarning(
    PytestDeprecationWarning,
    "{name} is deprecated, use "
    "function._request._fillfixtures() instead if you cannot avoid reaching into internals.",
)

PYTEST_COLLECT_MODULE = UnformattedWarning(
    PytestDeprecationWarning,
    "pytest.collect.{name} was moved to pytest.{name}\n"
    "Please update to the new name.",
)

YIELD_FIXTURE = PytestDeprecationWarning(
    "@pytest.yield_fixture is deprecated.\n"
    "Use @pytest.fixture instead; they are the same."
)

MINUS_K_DASH = PytestDeprecationWarning(
    "The `-k '-expr'` syntax to -k is deprecated.\nUse `-k 'not expr'` instead."
)

MINUS_K_COLON = PytestDeprecationWarning(
    "The `-k 'expr:'` syntax to -k is deprecated.\n"
    "Please open an issue if you use this and want a replacement."
)

WARNING_CAPTURED_HOOK = PytestDeprecationWarning(
    "The pytest_warning_captured is deprecated and will be removed in a future release.\n"
    "Please use pytest_warning_recorded instead."
)

FSCOLLECTOR_GETHOOKPROXY_ISINITPATH = PytestDeprecationWarning(
    "The gethookproxy() and isinitpath() methods of FSCollector and Package are deprecated; "
    "use self.session.gethookproxy() and self.session.isinitpath() instead. "
)

STRICT_OPTION = PytestDeprecationWarning(
    "The --strict option is deprecated, use --strict-markers instead."
)

PRIVATE = PytestDeprecationWarning("A private pytest class or function was used.")


# You want to make some `__init__` or function "private".
#
#   def my_private_function(some, args):
#       ...
#
# Do this:
#
#   def my_private_function(some, args, *, _ispytest: bool = False):
#       check_ispytest(_ispytest)
#       ...
#
# Change all internal/allowed calls to
#
#   my_private_function(some, args, _ispytest=True)
#
# All other calls will get the default _ispytest=False and trigger
# the warning (possibly error in the future).
def check_ispytest(ispytest: bool) -> None:
    if not ispytest:
        warn(PRIVATE, stacklevel=3)
