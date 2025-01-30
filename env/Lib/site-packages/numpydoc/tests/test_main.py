import inspect
import io
import sys

import pytest

import numpydoc
import numpydoc.cli


def _capture_stdout(func_name, *args, **kwargs):
    """
    Return stdout of calling `func_name`.

    This docstring should be perfect, as it is used to test the
    validation with a docstring without errors.

    Parameters
    ----------
    func_name : callable
        Function to be called.
    *args, **kwargs
        Will be passed to `func_name`.

    Returns
    -------
    str
        The content that the function printed.

    See Also
    --------
    sys.stdout : Python's file handler for stdout.

    Examples
    --------
    >>> _capture_stdout(print, "hello world")
    'hello world'
    """
    f = io.StringIO()
    sys.stdout, old_stdout = f, sys.stdout
    try:
        func_name(*args, **kwargs)
        return f.getvalue().strip("\n\r")
    finally:
        sys.stdout = old_stdout


def _docstring_with_errors():
    """
    this docstring should report some errors

    Parameters
    ----------
    made_up_param : str
    """


def _invalid_docstring():
    """
    This docstring should break the parsing.

    See Also
    --------
    : this is invalid
    """


def test_renders_package_docstring():
    out = _capture_stdout(numpydoc.cli.render_object, "numpydoc")
    assert out.startswith("This package provides the numpydoc Sphinx")


def test_renders_module_docstring(capsys):
    numpydoc.cli.main(["render", "numpydoc.cli"])
    out = capsys.readouterr().out.strip("\n\r")
    assert out.startswith(numpydoc.cli.__doc__)


def test_renders_function_docstring():
    out = _capture_stdout(
        numpydoc.cli.render_object, "numpydoc.tests.test_main._capture_stdout"
    )
    assert out.startswith("Return stdout of calling")


def test_render_object_returns_correct_exit_status():
    exit_status = numpydoc.cli.render_object("numpydoc.tests.test_main._capture_stdout")
    assert exit_status == 0

    with pytest.raises(ValueError):
        numpydoc.cli.render_object("numpydoc.tests.test_main._invalid_docstring")


def test_validate_detects_errors():
    out = _capture_stdout(
        numpydoc.cli.validate_object,
        "numpydoc.tests.test_main._docstring_with_errors",
    )
    assert "SS02" in out
    assert "Summary does not start with a capital letter" in out

    exit_status = numpydoc.cli.validate_object(
        "numpydoc.tests.test_main._docstring_with_errors"
    )
    assert exit_status > 0


def test_validate_perfect_docstring():
    out = _capture_stdout(
        numpydoc.cli.validate_object, "numpydoc.tests.test_main._capture_stdout"
    )
    assert out == ""

    exit_status = numpydoc.cli.validate_object(
        "numpydoc.tests.test_main._capture_stdout"
    )
    assert exit_status == 0


@pytest.mark.parametrize("args", [[], ["--ignore", "ES01", "SA01", "EX01"]])
def test_lint(capsys, args):
    argv = ["lint", "numpydoc/__main__.py"] + args
    if args:
        expected = ""
        expected_status = 0
    else:
        expected = inspect.cleandoc(
            """
            +------------------------+----------+---------+----------------------------+
            | file                   | item     | check   | description                |
            +========================+==========+=========+============================+
            | numpydoc/__main__.py:1 | __main__ | ES01    | No extended summary found  |
            +------------------------+----------+---------+----------------------------+
            | numpydoc/__main__.py:1 | __main__ | SA01    | See Also section not found |
            +------------------------+----------+---------+----------------------------+
            | numpydoc/__main__.py:1 | __main__ | EX01    | No examples section found  |
            +------------------------+----------+---------+----------------------------+
        """
        )
        expected_status = 1

    return_status = numpydoc.cli.main(argv)
    err = capsys.readouterr().err.strip("\n\r")
    assert err == expected
    assert return_status == expected_status


def test_lint_help(capsys):
    """Test that lint help section is displaying."""

    with pytest.raises(SystemExit):
        return_code = numpydoc.cli.main(["lint", "--help"])
        assert return_code == 0

    out = capsys.readouterr().out
    assert "--ignore" in out
    assert "--config" in out
