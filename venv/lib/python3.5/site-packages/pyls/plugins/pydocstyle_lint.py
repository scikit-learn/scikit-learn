# Copyright 2017 Palantir Technologies, Inc.
import contextlib
import logging
import os
import re
import sys

import pydocstyle
from pyls import hookimpl, lsp

log = logging.getLogger(__name__)

# PyDocstyle is a little verbose in debug message
pydocstyle_logger = logging.getLogger(pydocstyle.utils.__name__)
pydocstyle_logger.setLevel(logging.INFO)

DEFAULT_MATCH_RE = pydocstyle.config.ConfigurationParser.DEFAULT_MATCH_RE
DEFAULT_MATCH_DIR_RE = pydocstyle.config.ConfigurationParser.DEFAULT_MATCH_DIR_RE


@hookimpl
def pyls_settings():
    # Default pydocstyle to disabled
    return {'plugins': {'pydocstyle': {'enabled': False}}}


@hookimpl
def pyls_lint(config, document):
    settings = config.plugin_settings('pydocstyle')
    log.debug("Got pydocstyle settings: %s", settings)

    # Explicitly passing a path to pydocstyle means it doesn't respect the --match flag, so do it ourselves
    filename_match_re = re.compile(settings.get('match', DEFAULT_MATCH_RE) + '$')
    if not filename_match_re.match(os.path.basename(document.path)):
        return []

    # Likewise with --match-dir
    dir_match_re = re.compile(settings.get('matchDir', DEFAULT_MATCH_DIR_RE) + '$')
    if not dir_match_re.match(os.path.basename(os.path.dirname(document.path))):
        return []

    args = [document.path]

    if settings.get('convention'):
        args.append('--convention=' + settings['convention'])

        if settings.get('addSelect'):
            args.append('--add-select=' + ','.join(settings['addSelect']))
        if settings.get('addIgnore'):
            args.append('--add-ignore=' + ','.join(settings['addIgnore']))

    elif settings.get('select'):
        args.append('--select=' + ','.join(settings['select']))
    elif settings.get('ignore'):
        args.append('--ignore=' + ','.join(settings['ignore']))

    log.info("Using pydocstyle args: %s", args)

    conf = pydocstyle.config.ConfigurationParser()
    with _patch_sys_argv(args):
        # TODO(gatesn): We can add more pydocstyle args here from our pyls config
        conf.parse()

    # Will only yield a single filename, the document path
    diags = []
    for filename, checked_codes, ignore_decorators in conf.get_files_to_check():
        errors = pydocstyle.checker.ConventionChecker().check_source(
            document.source, filename, ignore_decorators=ignore_decorators
        )

        try:
            for error in errors:
                if error.code not in checked_codes:
                    continue
                diags.append(_parse_diagnostic(document, error))
        except pydocstyle.parser.ParseError:
            # In the case we cannot parse the Python file, just continue
            pass

    log.debug("Got pydocstyle errors: %s", diags)
    return diags


def _parse_diagnostic(document, error):
    lineno = error.definition.start - 1
    line = document.lines[0] if document.lines else ""

    start_character = len(line) - len(line.lstrip())
    end_character = len(line)

    return {
        'source': 'pydocstyle',
        'code': error.code,
        'message': error.message,
        'severity': lsp.DiagnosticSeverity.Warning,
        'range': {
            'start': {
                'line': lineno,
                'character': start_character
            },
            'end': {
                'line': lineno,
                'character': end_character
            }
        }
    }


@contextlib.contextmanager
def _patch_sys_argv(arguments):
    old_args = sys.argv

    # Preserve argv[0] since it's the executable
    sys.argv = old_args[0:1] + arguments

    try:
        yield
    finally:
        sys.argv = old_args
