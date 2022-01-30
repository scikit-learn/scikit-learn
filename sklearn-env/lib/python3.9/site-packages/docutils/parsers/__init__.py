# $Id: __init__.py 8671 2021-04-07 12:09:51Z milde $
# Author: David Goodger <goodger@python.org>
# Copyright: This module has been placed in the public domain.

"""
This package contains Docutils parser modules.
"""

__docformat__ = 'reStructuredText'

import sys
from importlib import import_module

from docutils import Component, frontend


class Parser(Component):
    settings_spec = (
        'Generic Parser Options',
        None,
        (('Disable directives that insert the contents of an external file; '
          'replaced with a "warning" system message.',
          ['--no-file-insertion'],
          {'action': 'store_false', 'default': 1,
           'dest': 'file_insertion_enabled',
           'validator': frontend.validate_boolean}),
         ('Enable directives that insert the contents '
          'of an external file. (default)',
          ['--file-insertion-enabled'],
          {'action': 'store_true'}),
         ('Disable the "raw" directive; '
          'replaced with a "warning" system message.',
          ['--no-raw'],
          {'action': 'store_false', 'default': 1, 'dest': 'raw_enabled',
           'validator': frontend.validate_boolean}),
         ('Enable the "raw" directive. (default)',
          ['--raw-enabled'],
          {'action': 'store_true'}),
         ('Maximal number of characters in an input line. Default 10 000.',
          ['--line-length-limit'],
          {'metavar': '<length>', 'type': 'int', 'default': 10000,
           'validator': frontend.validate_nonnegative_int}),
        ))
    component_type = 'parser'
    config_section = 'parsers'

    def parse(self, inputstring, document):
        """Override to parse `inputstring` into document tree `document`."""
        raise NotImplementedError('subclass must override this method')

    def setup_parse(self, inputstring, document):
        """Initial parse setup.  Call at start of `self.parse()`."""
        self.inputstring = inputstring
        # provide fallbacks in case the document has only generic settings
        document.settings.setdefault('file_insertion_enabled', False)
        document.settings.setdefault('raw_enabled', False)
        document.settings.setdefault('line_length_limit', 10000)
        self.document = document
        document.reporter.attach_observer(document.note_parse_message)

    def finish_parse(self):
        """Finalize parse details.  Call at end of `self.parse()`."""
        self.document.reporter.detach_observer(
            self.document.note_parse_message)


_parser_aliases = {
      'restructuredtext': 'rst',
      'rest': 'rst',
      'restx': 'rst',
      'rtxt': 'rst',
      'recommonmark': 'recommonmark_wrapper',
      'commonmark': 'recommonmark_wrapper',
      'markdown': 'recommonmark_wrapper'}

def get_parser_class(parser_name):
    """Return the Parser class from the `parser_name` module."""
    parser_name = parser_name.lower()
    if parser_name in _parser_aliases:
        parser_name = _parser_aliases[parser_name]
    try:
        module = import_module('docutils.parsers.'+parser_name)
    except ImportError:
        module = import_module(parser_name)
    return module.Parser
