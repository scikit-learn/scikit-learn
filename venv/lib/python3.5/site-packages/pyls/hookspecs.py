# Copyright 2017 Palantir Technologies, Inc.
# pylint: disable=redefined-builtin, unused-argument
from pyls import hookspec


@hookspec
def pyls_code_actions(config, workspace, document, range, context):
    pass


@hookspec
def pyls_code_lens(config, workspace, document):
    pass


@hookspec
def pyls_commands(config, workspace):
    """The list of command strings supported by the server.

    Returns:
        List[str]: The supported commands.
    """


@hookspec
def pyls_completions(config, workspace, document, position):
    pass


@hookspec
def pyls_definitions(config, workspace, document, position):
    pass


@hookspec
def pyls_dispatchers(config, workspace):
    pass


@hookspec
def pyls_document_did_open(config, workspace, document):
    pass


@hookspec
def pyls_document_did_save(config, workspace, document):
    pass


@hookspec
def pyls_document_highlight(config, workspace, document, position):
    pass


@hookspec
def pyls_document_symbols(config, workspace, document):
    pass


@hookspec(firstresult=True)
def pyls_execute_command(config, workspace, command, arguments):
    pass


@hookspec
def pyls_experimental_capabilities(config, workspace):
    pass


@hookspec(firstresult=True)
def pyls_format_document(config, workspace, document):
    pass


@hookspec(firstresult=True)
def pyls_format_range(config, workspace, document, range):
    pass


@hookspec(firstresult=True)
def pyls_hover(config, workspace, document, position):
    pass


@hookspec
def pyls_initialize(config, workspace):
    pass


@hookspec
def pyls_lint(config, workspace, document):
    pass


@hookspec
def pyls_references(config, workspace, document, position, exclude_declaration):
    pass


@hookspec(firstresult=True)
def pyls_rename(config, workspace, document, position, new_name):
    pass


@hookspec
def pyls_settings(config):
    pass


@hookspec(firstresult=True)
def pyls_signature_help(config, workspace, document, position):
    pass
