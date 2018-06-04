# Copyright 2018 Palantir Technologies, Inc.
import logging
from autopep8 import fix_code
from pyls import hookimpl

log = logging.getLogger(__name__)


@hookimpl(tryfirst=True)  # Prefer autopep8 over YAPF
def pyls_format_document(config, document):
    log.info("Formatting document %s with autopep8", document)
    return _format(config, document)


@hookimpl(tryfirst=True)  # Prefer autopep8 over YAPF
def pyls_format_range(config, document, range):  # pylint: disable=redefined-builtin
    log.info("Formatting document %s in range %s with autopep8", document, range)

    # First we 'round' the range up/down to full lines only
    range['start']['character'] = 0
    range['end']['line'] += 1
    range['end']['character'] = 0

    # Add 1 for 1-indexing vs LSP's 0-indexing
    line_range = (range['start']['line'] + 1, range['end']['line'] + 1)
    return _format(config, document, line_range=line_range)


def _format(config, document, line_range=None):
    options = _autopep8_config(config)
    if line_range:
        options['line_range'] = list(line_range)

    new_source = fix_code(document.source, options=options)

    if new_source == document.source:
        return []

    # I'm too lazy at the moment to parse diffs into TextEdit items
    # So let's just return the entire file...
    return [{
        'range': {
            'start': {'line': 0, 'character': 0},
            # End char 0 of the line after our document
            'end': {'line': len(document.lines), 'character': 0}
        },
        'newText': new_source
    }]


def _autopep8_config(config):
    # We user pycodestyle settings to avoid redefining things
    settings = config.plugin_settings('pycodestyle')
    options = {
        'exclude': settings.get('exclude'),
        'hang_closing': settings.get('hangClosing'),
        'ignore': settings.get('ignore'),
        'max_line_length': settings.get('maxLineLength'),
        'select': settings.get('select'),
    }

    # Filter out null options
    return {k: v for k, v in options.items() if v}
