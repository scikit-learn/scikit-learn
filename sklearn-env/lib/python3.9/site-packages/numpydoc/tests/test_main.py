import sys
import io
import pytest
import numpydoc
import numpydoc.__main__


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
    >>> _capture_stdout(print, 'hello world')
    'hello world'
    """
    f = io.StringIO()
    sys.stdout, old_stdout = f, sys.stdout
    try:
        func_name(*args, **kwargs)
        return f.getvalue().strip('\n\r')
    finally:
        sys.stdout = old_stdout


def _docstring_with_errors():
    """
    this docstring should report some errors

    Parameters
    ----------
    made_up_param : str
    """
    pass


def _invalid_docstring():
    """
    This docstring should break the parsing.

    See Also
    --------
    : this is invalid
    """
    pass


def test_renders_package_docstring():
    out = _capture_stdout(numpydoc.__main__.render_object,
                          'numpydoc')
    assert out.startswith('This package provides the numpydoc Sphinx')


def test_renders_module_docstring():
    out = _capture_stdout(numpydoc.__main__.render_object,
                          'numpydoc.__main__')
    assert out.startswith('Implementing `python -m numpydoc` functionality.')


def test_renders_function_docstring():
    out = _capture_stdout(numpydoc.__main__.render_object,
                          'numpydoc.tests.test_main._capture_stdout')
    assert out.startswith('Return stdout of calling')


def test_render_object_returns_correct_exit_status():
    exit_status = numpydoc.__main__.render_object(
        'numpydoc.tests.test_main._capture_stdout')
    assert exit_status == 0

    with pytest.raises(ValueError):
        numpydoc.__main__.render_object(
            'numpydoc.tests.test_main._invalid_docstring')


def test_validate_detects_errors():
    out = _capture_stdout(numpydoc.__main__.validate_object,
                          'numpydoc.tests.test_main._docstring_with_errors')
    assert 'SS02' in out
    assert 'Summary does not start with a capital letter' in out

    exit_status = numpydoc.__main__.validate_object(
        'numpydoc.tests.test_main._docstring_with_errors')
    assert exit_status > 0


def test_validate_perfect_docstring():
    out = _capture_stdout(numpydoc.__main__.validate_object,
                          'numpydoc.tests.test_main._capture_stdout')
    assert out == ''

    exit_status = numpydoc.__main__.validate_object(
        'numpydoc.tests.test_main._capture_stdout')
    assert exit_status == 0
