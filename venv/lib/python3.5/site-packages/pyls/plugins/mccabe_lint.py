# Copyright 2017 Palantir Technologies, Inc.
import ast
import logging
import mccabe
from pyls import hookimpl, lsp

log = logging.getLogger(__name__)

THRESHOLD = 'threshold'
DEFAULT_THRESHOLD = 15


@hookimpl
def pyls_lint(config, document):
    threshold = config.plugin_settings('mccabe').get(THRESHOLD, DEFAULT_THRESHOLD)
    log.debug("Running mccabe lint with threshold: %s", threshold)

    try:
        tree = compile(document.source, document.path, "exec", ast.PyCF_ONLY_AST)
    except SyntaxError:
        # We'll let the other linters point this one out
        return None

    visitor = mccabe.PathGraphingAstVisitor()
    visitor.preorder(tree, visitor)

    diags = []
    for graph in visitor.graphs.values():
        if graph.complexity() >= threshold:
            diags.append({
                'source': 'mccabe',
                'range': {
                    'start': {'line': graph.lineno, 'character': graph.column},
                    'end': {'line': graph.lineno, 'character': len(document.lines[graph.lineno])},
                },
                'message': 'Cyclomatic complexity too high: %s (threshold %s)' % (graph.complexity(), threshold),
                'severity': lsp.DiagnosticSeverity.Warning
            })

    return diags
