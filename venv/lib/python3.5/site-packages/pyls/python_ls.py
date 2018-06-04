# Copyright 2017 Palantir Technologies, Inc.
import logging
import socketserver

from . import lsp, _utils, uris
from .config import config
from .jsonrpc.dispatchers import MethodDispatcher
from .jsonrpc.endpoint import Endpoint
from .jsonrpc.streams import JsonRpcStreamReader, JsonRpcStreamWriter
from .workspace import Workspace

log = logging.getLogger(__name__)


LINT_DEBOUNCE_S = 0.5  # 500 ms


class _StreamHandlerWrapper(socketserver.StreamRequestHandler, object):
    """A wrapper class that is used to construct a custom handler class."""

    delegate = None

    def setup(self):
        super(_StreamHandlerWrapper, self).setup()
        # pylint: disable=no-member
        self.delegate = self.DELEGATE_CLASS(self.rfile, self.wfile)

    def handle(self):
        self.delegate.start()


def start_tcp_lang_server(bind_addr, port, handler_class):
    if not issubclass(handler_class, PythonLanguageServer):
        raise ValueError('Handler class must be an instance of PythonLanguageServer')

    # Construct a custom wrapper class around the user's handler_class
    wrapper_class = type(
        handler_class.__name__ + 'Handler',
        (_StreamHandlerWrapper,),
        {'DELEGATE_CLASS': handler_class}
    )

    server = socketserver.TCPServer((bind_addr, port), wrapper_class)
    try:
        log.info('Serving %s on (%s, %s)', handler_class.__name__, bind_addr, port)
        server.serve_forever()
    finally:
        log.info('Shutting down')
        server.server_close()


def start_io_lang_server(rfile, wfile, handler_class):
    if not issubclass(handler_class, PythonLanguageServer):
        raise ValueError('Handler class must be an instance of PythonLanguageServer')
    log.info('Starting %s IO language server', handler_class.__name__)
    server = handler_class(rfile, wfile)
    server.start()


class PythonLanguageServer(MethodDispatcher):
    """ Implementation of the Microsoft VSCode Language Server Protocol
    https://github.com/Microsoft/language-server-protocol/blob/master/versions/protocol-1-x.md
    """

    # pylint: disable=too-many-public-methods,redefined-builtin

    def __init__(self, rx, tx):
        self.workspace = None
        self.config = None

        self._jsonrpc_stream_reader = JsonRpcStreamReader(rx)
        self._jsonrpc_stream_writer = JsonRpcStreamWriter(tx)
        self._endpoint = Endpoint(self, self._jsonrpc_stream_writer.write)
        self._dispatchers = []
        self._shutdown = False

    def start(self):
        """Entry point for the server."""
        self._jsonrpc_stream_reader.listen(self._endpoint.consume)

    def __getitem__(self, item):
        """Override getitem to fallback through multiple dispatchers."""
        if self._shutdown and item != 'exit':
            # exit is the only allowed method during shutdown
            log.debug("Ignoring non-exit method during shutdown: %s", item)
            raise KeyError

        try:
            return super(PythonLanguageServer, self).__getitem__(item)
        except KeyError:
            # Fallback through extra dispatchers
            for dispatcher in self._dispatchers:
                try:
                    return dispatcher[item]
                except KeyError:
                    continue

        raise KeyError()

    def m_shutdown(self, **_kwargs):
        self._shutdown = True
        return None

    def m_exit(self, **_kwargs):
        self._endpoint.shutdown()
        self._jsonrpc_stream_reader.close()
        self._jsonrpc_stream_writer.close()

    def _hook(self, hook_name, doc_uri=None, **kwargs):
        """Calls hook_name and returns a list of results from all registered handlers"""
        doc = self.workspace.get_document(doc_uri) if doc_uri else None
        hook_handlers = self.config.plugin_manager.subset_hook_caller(hook_name, self.config.disabled_plugins)
        return hook_handlers(config=self.config, workspace=self.workspace, document=doc, **kwargs)

    def capabilities(self):
        server_capabilities = {
            'codeActionProvider': True,
            'codeLensProvider': {
                'resolveProvider': False,  # We may need to make this configurable
            },
            'completionProvider': {
                'resolveProvider': False,  # We know everything ahead of time
                'triggerCharacters': ['.']
            },
            'documentFormattingProvider': True,
            'documentHighlightProvider': True,
            'documentRangeFormattingProvider': True,
            'documentSymbolProvider': True,
            'definitionProvider': True,
            'executeCommandProvider': {
                'commands': flatten(self._hook('pyls_commands'))
            },
            'hoverProvider': True,
            'referencesProvider': True,
            'renameProvider': True,
            'signatureHelpProvider': {
                'triggerCharacters': ['(', ',']
            },
            'textDocumentSync': lsp.TextDocumentSyncKind.INCREMENTAL,
            'experimental': merge(self._hook('pyls_experimental_capabilities'))
        }
        log.info('Server capabilities: %s', server_capabilities)
        return server_capabilities

    def m_initialize(self, processId=None, rootUri=None, rootPath=None, initializationOptions=None, **_kwargs):
        log.debug('Language server initialized with %s %s %s %s', processId, rootUri, rootPath, initializationOptions)
        if rootUri is None:
            rootUri = uris.from_fs_path(rootPath) if rootPath is not None else ''

        self.workspace = Workspace(rootUri, self._endpoint)
        self.config = config.Config(rootUri, initializationOptions or {})
        self._dispatchers = self._hook('pyls_dispatchers')
        self._hook('pyls_initialize')

        # Get our capabilities
        return {'capabilities': self.capabilities()}

    def m_initialized(self, **_kwargs):
        pass

    def code_actions(self, doc_uri, range, context):
        return flatten(self._hook('pyls_code_actions', doc_uri, range=range, context=context))

    def code_lens(self, doc_uri):
        return flatten(self._hook('pyls_code_lens', doc_uri))

    def completions(self, doc_uri, position):
        completions = self._hook('pyls_completions', doc_uri, position=position)
        return {
            'isIncomplete': False,
            'items': flatten(completions)
        }

    def definitions(self, doc_uri, position):
        return flatten(self._hook('pyls_definitions', doc_uri, position=position))

    def document_symbols(self, doc_uri):
        return flatten(self._hook('pyls_document_symbols', doc_uri))

    def execute_command(self, command, arguments):
        return self._hook('pyls_execute_command', command=command, arguments=arguments)

    def format_document(self, doc_uri):
        return self._hook('pyls_format_document', doc_uri)

    def format_range(self, doc_uri, range):
        return self._hook('pyls_format_range', doc_uri, range=range)

    def highlight(self, doc_uri, position):
        return flatten(self._hook('pyls_document_highlight', doc_uri, position=position)) or None

    def hover(self, doc_uri, position):
        return self._hook('pyls_hover', doc_uri, position=position) or {'contents': ''}

    @_utils.debounce(LINT_DEBOUNCE_S, keyed_by='doc_uri')
    def lint(self, doc_uri):
        # Since we're debounced, the document may no longer be open
        if doc_uri in self.workspace.documents:
            self.workspace.publish_diagnostics(doc_uri, flatten(self._hook('pyls_lint', doc_uri)))

    def references(self, doc_uri, position, exclude_declaration):
        return flatten(self._hook(
            'pyls_references', doc_uri, position=position,
            exclude_declaration=exclude_declaration
        ))

    def rename(self, doc_uri, position, new_name):
        return self._hook('pyls_rename', doc_uri, position=position, new_name=new_name)

    def signature_help(self, doc_uri, position):
        return self._hook('pyls_signature_help', doc_uri, position=position)

    def m_text_document__did_close(self, textDocument=None, **_kwargs):
        self.workspace.rm_document(textDocument['uri'])

    def m_text_document__did_open(self, textDocument=None, **_kwargs):
        self.workspace.put_document(textDocument['uri'], textDocument['text'], version=textDocument.get('version'))
        self._hook('pyls_document_did_open', textDocument['uri'])
        self.lint(textDocument['uri'])

    def m_text_document__did_change(self, contentChanges=None, textDocument=None, **_kwargs):
        for change in contentChanges:
            self.workspace.update_document(
                textDocument['uri'],
                change,
                version=textDocument.get('version')
            )
        self.lint(textDocument['uri'])

    def m_text_document__did_save(self, textDocument=None, **_kwargs):
        self.lint(textDocument['uri'])

    def m_text_document__code_action(self, textDocument=None, range=None, context=None, **_kwargs):
        return self.code_actions(textDocument['uri'], range, context)

    def m_text_document__code_lens(self, textDocument=None, **_kwargs):
        return self.code_lens(textDocument['uri'])

    def m_text_document__completion(self, textDocument=None, position=None, **_kwargs):
        return self.completions(textDocument['uri'], position)

    def m_text_document__definition(self, textDocument=None, position=None, **_kwargs):
        return self.definitions(textDocument['uri'], position)

    def m_text_document__document_highlight(self, textDocument=None, position=None, **_kwargs):
        return self.highlight(textDocument['uri'], position)

    def m_text_document__hover(self, textDocument=None, position=None, **_kwargs):
        return self.hover(textDocument['uri'], position)

    def m_text_document__document_symbol(self, textDocument=None, **_kwargs):
        return self.document_symbols(textDocument['uri'])

    def m_text_document__formatting(self, textDocument=None, _options=None, **_kwargs):
        # For now we're ignoring formatting options.
        return self.format_document(textDocument['uri'])

    def m_text_document__rename(self, textDocument=None, position=None, newName=None, **_kwargs):
        return self.rename(textDocument['uri'], position, newName)

    def m_text_document__range_formatting(self, textDocument=None, range=None, _options=None, **_kwargs):
        # Again, we'll ignore formatting options for now.
        return self.format_range(textDocument['uri'], range)

    def m_text_document__references(self, textDocument=None, position=None, context=None, **_kwargs):
        exclude_declaration = not context['includeDeclaration']
        return self.references(textDocument['uri'], position, exclude_declaration)

    def m_text_document__signature_help(self, textDocument=None, position=None, **_kwargs):
        return self.signature_help(textDocument['uri'], position)

    def m_workspace__did_change_configuration(self, settings=None):
        self.config.update((settings or {}).get('pyls', {}))
        for doc_uri in self.workspace.documents:
            self.lint(doc_uri)

    def m_workspace__did_change_watched_files(self, **_kwargs):
        # Externally changed files may result in changed diagnostics
        for doc_uri in self.workspace.documents:
            self.lint(doc_uri)

    def m_workspace__execute_command(self, command=None, arguments=None):
        return self.execute_command(command, arguments)


def flatten(list_of_lists):
    return [item for lst in list_of_lists for item in lst]


def merge(list_of_dicts):
    return {k: v for dictionary in list_of_dicts for k, v in dictionary.items()}
