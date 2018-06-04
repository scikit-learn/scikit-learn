"""
This module contains deprecation messages and bits of code used elsewhere in the codebase
that is planned to be removed in the next pytest release.

Keeping it in a central location makes it easy to track what is deprecated and should
be removed when the time comes.
"""
from __future__ import absolute_import, division, print_function


class RemovedInPytest4Warning(DeprecationWarning):
    """warning class for features removed in pytest 4.0"""


MAIN_STR_ARGS = 'passing a string to pytest.main() is deprecated, ' \
    'pass a list of arguments instead.'

YIELD_TESTS = 'yield tests are deprecated, and scheduled to be removed in pytest 4.0'

FUNCARG_PREFIX = (
    '{name}: declaring fixtures using "pytest_funcarg__" prefix is deprecated '
    'and scheduled to be removed in pytest 4.0.  '
    'Please remove the prefix and use the @pytest.fixture decorator instead.')

CFG_PYTEST_SECTION = '[pytest] section in {filename} files is deprecated, use [tool:pytest] instead.'

GETFUNCARGVALUE = "use of getfuncargvalue is deprecated, use getfixturevalue"

RESULT_LOG = (
    '--result-log is deprecated and scheduled for removal in pytest 4.0.\n'
    'See https://docs.pytest.org/en/latest/usage.html#creating-resultlog-format-files for more information.'
)

MARK_INFO_ATTRIBUTE = RemovedInPytest4Warning(
    "MarkInfo objects are deprecated as they contain the merged marks"
)

MARK_PARAMETERSET_UNPACKING = RemovedInPytest4Warning(
    "Applying marks directly to parameters is deprecated,"
    " please use pytest.param(..., marks=...) instead.\n"
    "For more details, see: https://docs.pytest.org/en/latest/parametrize.html"
)

RECORD_XML_PROPERTY = (
    'Fixture renamed from "record_xml_property" to "record_property" as user '
    'properties are now available to all reporters.\n'
    '"record_xml_property" is now deprecated.'
)

COLLECTOR_MAKEITEM = RemovedInPytest4Warning(
    "pycollector makeitem was removed "
    "as it is an accidentially leaked internal api"
)

METAFUNC_ADD_CALL = (
    "Metafunc.addcall is deprecated and scheduled to be removed in pytest 4.0.\n"
    "Please use Metafunc.parametrize instead."
)

PYTEST_PLUGINS_FROM_NON_TOP_LEVEL_CONFTEST = RemovedInPytest4Warning(
    "Defining pytest_plugins in a non-top-level conftest is deprecated, "
    "because it affects the entire directory tree in a non-explicit way.\n"
    "Please move it to the top level conftest file instead."
)
