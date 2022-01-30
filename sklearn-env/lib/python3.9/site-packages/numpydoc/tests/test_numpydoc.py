# -*- encoding:utf-8 -*-
import pytest
from io import StringIO
from copy import deepcopy
from numpydoc.numpydoc import (
    mangle_docstrings, _clean_text_signature, update_config
)
from numpydoc.xref import DEFAULT_LINKS
from sphinx.ext.autodoc import ALL
from sphinx.util import logging


class MockConfig():
    numpydoc_use_plots = False
    numpydoc_use_blockquotes = True
    numpydoc_show_class_members = True
    numpydoc_show_inherited_class_members = True
    numpydoc_class_members_toctree = True
    numpydoc_xref_param_type = False
    numpydoc_xref_aliases = {}
    numpydoc_xref_aliases_complete = deepcopy(DEFAULT_LINKS)
    numpydoc_xref_ignore = set()
    templates_path = []
    numpydoc_edit_link = False
    numpydoc_citation_re = '[a-z0-9_.-]+'
    numpydoc_attributes_as_param_list = True
    numpydoc_validation_checks = set()
    numpydoc_validation_exclude = set()


class MockBuilder():
    config = MockConfig()


class MockApp():
    config = MockConfig()
    builder = MockBuilder()
    translator = None

    def __init__(self):
        self.builder.app = self
        # Attrs required for logging
        self.verbosity = 2
        self._warncount = 0
        self.warningiserror = False


def test_mangle_docstrings():
    s = '''
A top section before

.. autoclass:: str
    '''
    lines = s.split('\n')
    mangle_docstrings(MockApp(), 'class', 'str', str, {}, lines)
    assert 'rpartition' in [x.strip() for x in lines]

    lines = s.split('\n')
    mangle_docstrings(
        MockApp(), 'class', 'str', str, {'members': ['upper']}, lines)
    assert 'rpartition' not in [x.strip() for x in lines]
    assert 'upper' in [x.strip() for x in lines]

    lines = s.split('\n')
    mangle_docstrings(
        MockApp(), 'class', 'str', str, {'exclude-members': ALL}, lines)
    assert 'rpartition' not in [x.strip() for x in lines]
    assert 'upper' not in [x.strip() for x in lines]

    lines = s.split('\n')
    mangle_docstrings(
        MockApp(), 'class', 'str', str, {'exclude-members': ['upper']}, lines)
    assert 'rpartition' in [x.strip() for x in lines]
    assert 'upper' not in [x.strip() for x in lines]


def test_clean_text_signature():
    assert _clean_text_signature(None) is None
    assert _clean_text_signature('func($self)') == 'func()'
    assert (_clean_text_signature('func($self, *args, **kwargs)')
            == 'func(*args, **kwargs)')
    assert _clean_text_signature('($self)') == '()'
    assert _clean_text_signature('()') == '()'
    assert _clean_text_signature('func()') == 'func()'
    assert (_clean_text_signature('func($self, /, *args, **kwargs)')
            == 'func(*args, **kwargs)')
    assert (_clean_text_signature('func($self, other, /, *args, **kwargs)')
            == 'func(other, *args, **kwargs)')
    assert _clean_text_signature('($module)') == '()'
    assert _clean_text_signature('func($type)') == 'func()'
    assert (_clean_text_signature('func($self, foo="hello world")')
            == 'func(foo="hello world")')
    assert (_clean_text_signature("func($self, foo='hello world')")
            == "func(foo='hello world')")
    assert (_clean_text_signature('func(foo="hello world")')
            == 'func(foo="hello world")')
    assert (_clean_text_signature('func(foo="$self")')
            == 'func(foo="$self")')
    assert (_clean_text_signature('func($self, foo="$self")')
            == 'func(foo="$self")')
    assert _clean_text_signature('func(self, other)') == 'func(self, other)'
    assert _clean_text_signature('func($self, *args)') == 'func(*args)'


@pytest.fixture
def f():
    def _function_without_seealso_and_examples():
        """
        A function whose docstring has no examples or see also section.

        Expect SA01 and EX01 errors if validation enabled.
        """
        pass
    return _function_without_seealso_and_examples


@pytest.mark.parametrize(
    (
        'numpydoc_validation_checks',
        'expected_warn',
        'non_warnings',
    ),
    (
        # Validation configured off - expect no warnings
        (set(), [], []),
        # Validation on with expected warnings
        (set(['SA01', 'EX01']), ('SA01', 'EX01'), []),
        # Validation on with only one activated check
        (set(['SA01']), ('SA01',), ('EX01',)),
    ),
)
def test_mangle_docstring_validation_warnings(
    f,
    numpydoc_validation_checks,
    expected_warn,
    non_warnings,
):
    app = MockApp()
    # Set up config for test
    app.config.numpydoc_validation_checks = numpydoc_validation_checks
    # Update configuration
    update_config(app)
    # Set up logging
    status, warning = StringIO(), StringIO()
    logging.setup(app, status, warning)
    # Run mangle docstrings with the above configuration
    mangle_docstrings(app, 'function', 'f', f, None, f.__doc__.split('\n'))
    # Assert that all (and only) expected warnings are logged
    warnings = warning.getvalue()
    for w in expected_warn:
        assert w in warnings
    for w in non_warnings:
        assert w not in warnings


def test_mangle_docstring_validation_exclude():
    def function_with_bad_docstring():
        """
        This docstring will raise docstring validation warnings."""
    app = MockApp()
    app.config.numpydoc_validation_checks = {"all"}
    app.config.numpydoc_validation_exclude = [r"_bad_"]
    # Call update_config to construct regexp from config value
    update_config(app)
    # Setup for catching warnings
    status, warning = StringIO(), StringIO()
    logging.setup(app, status, warning)
    # Run mangle docstrings on function_with_bad_docstring
    mangle_docstrings(
        app,
        'function',
        function_with_bad_docstring.__name__,
        function_with_bad_docstring,
        None,
        function_with_bad_docstring.__doc__.split('\n'),
    )
    # Validation is skipped due to exclude pattern matching fn name, therefore
    # no warnings expected
    assert warning.getvalue() == ""


def test_update_config_invalid_validation_set():
    app = MockApp()
    # Results in {'a', 'l'} instead of {"all"}
    app.config.numpydoc_validation_checks = set("all")
    with pytest.raises(ValueError, match="Unrecognized validation code"):
        update_config(app)


def test_update_config_exclude_str():
    app = MockApp()
    app.config.numpydoc_validation_checks = set()
    app.config.numpydoc_validation_exclude = "shouldnt-be-a-str"
    with pytest.raises(ValueError, match="\['shouldnt-be-a-str'\]"):
        update_config(app)


if __name__ == "__main__":
    import pytest
    pytest.main()
