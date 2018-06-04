# Copyright 2017 Palantir Technologies, Inc.
import logging
from pyls import hookimpl, lsp, _utils

log = logging.getLogger(__name__)


@hookimpl
def pyls_completions(document, position):
    definitions = document.jedi_script(position).completions()
    return [{
        'label': _label(d),
        'kind': _kind(d),
        'detail': _detail(d),
        'documentation': _utils.format_docstring(d.docstring()),
        'sortText': _sort_text(d),
        'insertText': d.name
    } for d in definitions] or None


def _label(definition):
    if definition.type in ('function', 'method'):
        params = ', '.join(param.name for param in definition.params)
        return '{}({})'.format(definition.name, params)

    return definition.name


def _detail(definition):
    return definition.parent().full_name or ''


def _sort_text(definition):
    """ Ensure builtins appear at the bottom.
    Description is of format <type>: <module>.<item>
    """

    # If its 'hidden', put it next last
    prefix = 'z{}' if definition.name.startswith('_') else 'a{}'
    return prefix.format(definition.name)


def _kind(d):
    """ Return the VSCode type """
    MAP = {
        'none': lsp.CompletionItemKind.Value,
        'type': lsp.CompletionItemKind.Class,
        'tuple': lsp.CompletionItemKind.Class,
        'dict': lsp.CompletionItemKind.Class,
        'dictionary': lsp.CompletionItemKind.Class,
        'function': lsp.CompletionItemKind.Function,
        'lambda': lsp.CompletionItemKind.Function,
        'generator': lsp.CompletionItemKind.Function,
        'class': lsp.CompletionItemKind.Class,
        'instance': lsp.CompletionItemKind.Reference,
        'method': lsp.CompletionItemKind.Method,
        'builtin': lsp.CompletionItemKind.Class,
        'builtinfunction': lsp.CompletionItemKind.Function,
        'module': lsp.CompletionItemKind.Module,
        'file': lsp.CompletionItemKind.File,
        'xrange': lsp.CompletionItemKind.Class,
        'slice': lsp.CompletionItemKind.Class,
        'traceback': lsp.CompletionItemKind.Class,
        'frame': lsp.CompletionItemKind.Class,
        'buffer': lsp.CompletionItemKind.Class,
        'dictproxy': lsp.CompletionItemKind.Class,
        'funcdef': lsp.CompletionItemKind.Function,
        'property': lsp.CompletionItemKind.Property,
        'import': lsp.CompletionItemKind.Module,
        'keyword': lsp.CompletionItemKind.Keyword,
        'constant': lsp.CompletionItemKind.Variable,
        'variable': lsp.CompletionItemKind.Variable,
        'value': lsp.CompletionItemKind.Value,
        'param': lsp.CompletionItemKind.Variable,
        'statement': lsp.CompletionItemKind.Keyword,
    }

    return MAP.get(d.type)
