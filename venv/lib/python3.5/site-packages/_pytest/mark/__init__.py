""" generic mechanism for marking and selecting python functions. """
from __future__ import absolute_import, division, print_function
from _pytest.config import UsageError
from .structures import (
    ParameterSet, EMPTY_PARAMETERSET_OPTION, MARK_GEN,
    Mark, MarkInfo, MarkDecorator, MarkGenerator,
    transfer_markers, get_empty_parameterset_mark
)
from .legacy import matchkeyword, matchmark

__all__ = [
    'Mark', 'MarkInfo', 'MarkDecorator', 'MarkGenerator',
    'transfer_markers', 'get_empty_parameterset_mark'
]


class MarkerError(Exception):

    """Error in use of a pytest marker/attribute."""


def param(*values, **kw):
    """Specify a parameter in `pytest.mark.parametrize`_ calls or
    :ref:`parametrized fixtures <fixture-parametrize-marks>`.

    .. code-block:: python

        @pytest.mark.parametrize("test_input,expected", [
            ("3+5", 8),
            pytest.param("6*9", 42, marks=pytest.mark.xfail),
        ])
        def test_eval(test_input, expected):
            assert eval(test_input) == expected

    :param values: variable args of the values of the parameter set, in order.
    :keyword marks: a single mark or a list of marks to be applied to this parameter set.
    :keyword str id: the id to attribute to this parameter set.
    """
    return ParameterSet.param(*values, **kw)


def pytest_addoption(parser):
    group = parser.getgroup("general")
    group._addoption(
        '-k',
        action="store", dest="keyword", default='', metavar="EXPRESSION",
        help="only run tests which match the given substring expression. "
             "An expression is a python evaluatable expression "
             "where all names are substring-matched against test names "
             "and their parent classes. Example: -k 'test_method or test_"
             "other' matches all test functions and classes whose name "
             "contains 'test_method' or 'test_other', while -k 'not test_method' "
             "matches those that don't contain 'test_method' in their names. "
             "Additionally keywords are matched to classes and functions "
             "containing extra names in their 'extra_keyword_matches' set, "
             "as well as functions which have names assigned directly to them."
    )

    group._addoption(
        "-m",
        action="store", dest="markexpr", default="", metavar="MARKEXPR",
        help="only run tests matching given mark expression.  "
             "example: -m 'mark1 and not mark2'."
    )

    group.addoption(
        "--markers", action="store_true",
        help="show markers (builtin, plugin and per-project ones)."
    )

    parser.addini("markers", "markers for test functions", 'linelist')
    parser.addini(
        EMPTY_PARAMETERSET_OPTION,
        "default marker for empty parametersets")


def pytest_cmdline_main(config):
    import _pytest.config
    if config.option.markers:
        config._do_configure()
        tw = _pytest.config.create_terminal_writer(config)
        for line in config.getini("markers"):
            parts = line.split(":", 1)
            name = parts[0]
            rest = parts[1] if len(parts) == 2 else ''
            tw.write("@pytest.mark.%s:" % name, bold=True)
            tw.line(rest)
            tw.line()
        config._ensure_unconfigure()
        return 0


pytest_cmdline_main.tryfirst = True


def deselect_by_keyword(items, config):
    keywordexpr = config.option.keyword.lstrip()
    if keywordexpr.startswith("-"):
        keywordexpr = "not " + keywordexpr[1:]
    selectuntil = False
    if keywordexpr[-1:] == ":":
        selectuntil = True
        keywordexpr = keywordexpr[:-1]

    remaining = []
    deselected = []
    for colitem in items:
        if keywordexpr and not matchkeyword(colitem, keywordexpr):
            deselected.append(colitem)
        else:
            if selectuntil:
                keywordexpr = None
            remaining.append(colitem)

    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = remaining


def deselect_by_mark(items, config):
    matchexpr = config.option.markexpr
    if not matchexpr:
        return

    remaining = []
    deselected = []
    for item in items:
        if matchmark(item, matchexpr):
            remaining.append(item)
        else:
            deselected.append(item)

    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = remaining


def pytest_collection_modifyitems(items, config):
    deselect_by_keyword(items, config)
    deselect_by_mark(items, config)


def pytest_configure(config):
    config._old_mark_config = MARK_GEN._config
    if config.option.strict:
        MARK_GEN._config = config

    empty_parameterset = config.getini(EMPTY_PARAMETERSET_OPTION)

    if empty_parameterset not in ('skip', 'xfail', None, ''):
        raise UsageError(
            "{!s} must be one of skip and xfail,"
            " but it is {!r}".format(EMPTY_PARAMETERSET_OPTION, empty_parameterset))


def pytest_unconfigure(config):
    MARK_GEN._config = getattr(config, '_old_mark_config', None)
