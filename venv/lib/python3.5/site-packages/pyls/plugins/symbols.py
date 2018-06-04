# Copyright 2017 Palantir Technologies, Inc.
import logging
from pyls import hookimpl
from pyls.lsp import SymbolKind

log = logging.getLogger(__name__)


@hookimpl
def pyls_document_symbols(config, document):
    all_scopes = config.plugin_settings('jedi_symbols').get('all_scopes', True)
    definitions = document.jedi_names(all_scopes=all_scopes)
    return [{
        'name': d.name,
        'containerName': _container(d),
        'location': {
            'uri': document.uri,
            'range': _range(d),
        },
        'kind': _kind(d),
    } for d in definitions if _include_def(d)]


def _include_def(definition):
    return (
        # Don't tend to include parameters as symbols
        definition.type != 'param' and
        # Unused vars should also be skipped
        definition.name != '_' and
        _kind(definition) is not None
    )


def _container(definition):
    try:
        # Jedi sometimes fails here.
        parent = definition.parent()
        # Here we check that a grand-parent exists to avoid declaring symbols
        # as children of the module.
        if parent.parent():
            return parent.name
    except:  # pylint: disable=bare-except
        return None

    return None


def _range(definition):
    # This gets us more accurate end position
    definition = definition._name.tree_name.get_definition()
    (start_line, start_column) = definition.start_pos
    (end_line, end_column) = definition.end_pos
    return {
        'start': {'line': start_line - 1, 'character': start_column},
        'end': {'line': end_line - 1, 'character': end_column}
    }


def _kind(d):
    """ Return the VSCode Symbol Type """
    MAP = {
        'none': SymbolKind.Variable,
        'type': SymbolKind.Class,
        'tuple': SymbolKind.Class,
        'dict': SymbolKind.Class,
        'dictionary': SymbolKind.Class,
        'function': SymbolKind.Function,
        'lambda': SymbolKind.Function,
        'generator': SymbolKind.Function,
        'class': SymbolKind.Class,
        'instance': SymbolKind.Class,
        'method': SymbolKind.Method,
        'builtin': SymbolKind.Class,
        'builtinfunction': SymbolKind.Function,
        'module': SymbolKind.Module,
        'file': SymbolKind.File,
        'xrange': SymbolKind.Array,
        'slice': SymbolKind.Class,
        'traceback': SymbolKind.Class,
        'frame': SymbolKind.Class,
        'buffer': SymbolKind.Array,
        'dictproxy': SymbolKind.Class,
        'funcdef': SymbolKind.Function,
        'property': SymbolKind.Property,
        'import': SymbolKind.Module,
        'keyword': SymbolKind.Variable,
        'constant': SymbolKind.Constant,
        'variable': SymbolKind.Variable,
        'value': SymbolKind.Variable,
        'param': SymbolKind.Variable,
        'statement': SymbolKind.Variable,
        'boolean': SymbolKind.Boolean,
        'int': SymbolKind.Number,
        'longlean': SymbolKind.Number,
        'float': SymbolKind.Number,
        'complex': SymbolKind.Number,
        'string': SymbolKind.String,
        'unicode': SymbolKind.String,
        'list': SymbolKind.Array,
    }

    return MAP.get(d.type)
