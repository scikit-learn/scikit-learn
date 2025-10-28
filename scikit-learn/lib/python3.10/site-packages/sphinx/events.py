"""Sphinx core events.

Gracefully adapted from the TextPress system by Armin.
"""

from __future__ import annotations

from collections import defaultdict
from operator import attrgetter
from typing import TYPE_CHECKING, NamedTuple, overload

from sphinx.errors import ExtensionError, SphinxError
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.inspect import safe_getattr

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence, Set
    from pathlib import Path
    from typing import Any, Literal

    from docutils import nodes

    from sphinx import addnodes
    from sphinx.application import Sphinx
    from sphinx.builders import Builder
    from sphinx.config import Config
    from sphinx.domains import Domain
    from sphinx.environment import BuildEnvironment
    from sphinx.ext.autodoc import _AutodocProcessDocstringListener
    from sphinx.ext.todo import todo_node


logger = logging.getLogger(__name__)


class EventListener(NamedTuple):
    id: int
    handler: Callable
    priority: int


# List of all known core events. Maps name to arguments description.
core_events = {
    'config-inited': 'config',
    'builder-inited': '',
    'env-get-outdated': 'env, added, changed, removed',
    'env-before-read-docs': 'env, docnames',
    'env-purge-doc': 'env, docname',
    'source-read': 'docname, source text',
    'include-read': 'relative path, parent docname, source text',
    'doctree-read': 'the doctree before being pickled',
    'env-merge-info': 'env, read docnames, other env instance',
    'env-updated': 'env',
    'env-get-updated': 'env',
    'env-check-consistency': 'env',
    'write-started': 'builder',
    'doctree-resolved': 'doctree, docname',
    'missing-reference': 'env, node, contnode',
    'warn-missing-reference': 'domain, node',
    'build-finished': 'exception',
}


class EventManager:
    """Event manager for Sphinx."""

    def __init__(self, app: Sphinx) -> None:
        self.app = app
        self.events = core_events.copy()
        self.listeners: dict[str, list[EventListener]] = defaultdict(list)
        self.next_listener_id = 0

    def add(self, name: str) -> None:
        """Register a custom Sphinx event."""
        if name in self.events:
            raise ExtensionError(__('Event %r already present') % name)
        self.events[name] = ''

    # ---- Core events -------------------------------------------------------

    @overload
    def connect(
        self,
        name: Literal['config-inited'],
        callback: Callable[[Sphinx, Config], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['builder-inited'],
        callback: Callable[[Sphinx], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-get-outdated'],
        callback: Callable[
            [Sphinx, BuildEnvironment, Set[str], Set[str], Set[str]], Sequence[str]
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-before-read-docs'],
        callback: Callable[[Sphinx, BuildEnvironment, list[str]], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-purge-doc'],
        callback: Callable[[Sphinx, BuildEnvironment, str], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['source-read'],
        callback: Callable[[Sphinx, str, list[str]], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['include-read'],
        callback: Callable[[Sphinx, Path, str, list[str]], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['doctree-read'],
        callback: Callable[[Sphinx, nodes.document], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-merge-info'],
        callback: Callable[
            [Sphinx, BuildEnvironment, list[str], BuildEnvironment], None
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-updated'],
        callback: Callable[[Sphinx, BuildEnvironment], str],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-get-updated'],
        callback: Callable[[Sphinx, BuildEnvironment], Iterable[str]],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['env-check-consistency'],
        callback: Callable[[Sphinx, BuildEnvironment], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['write-started'],
        callback: Callable[[Sphinx, Builder], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['doctree-resolved'],
        callback: Callable[[Sphinx, nodes.document, str], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['missing-reference'],
        callback: Callable[
            [Sphinx, BuildEnvironment, addnodes.pending_xref, nodes.TextElement],
            nodes.reference | None,
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['warn-missing-reference'],
        callback: Callable[[Sphinx, Domain, addnodes.pending_xref], bool | None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['build-finished'],
        callback: Callable[[Sphinx, Exception | None], None],
        priority: int,
    ) -> int: ...

    # ---- Events from builtin builders --------------------------------------

    @overload
    def connect(
        self,
        name: Literal['html-collect-pages'],
        callback: Callable[[Sphinx], Iterable[tuple[str, dict[str, Any], str]]],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['html-page-context'],
        callback: Callable[
            [Sphinx, str, str, dict[str, Any], nodes.document], str | None
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['linkcheck-process-uri'],
        callback: Callable[[Sphinx, str], str | None],
        priority: int,
    ) -> int: ...

    # ---- Events from builtin extensions-- ----------------------------------

    @overload
    def connect(
        self,
        name: Literal['object-description-transform'],
        callback: Callable[[Sphinx, str, str, addnodes.desc_content], None],
        priority: int,
    ) -> int: ...

    # ---- Events from first-party extensions --------------------------------

    @overload
    def connect(
        self,
        name: Literal['autodoc-process-docstring'],
        callback: _AutodocProcessDocstringListener,
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['autodoc-before-process-signature'],
        callback: Callable[[Sphinx, Any, bool], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['autodoc-process-signature'],
        callback: Callable[
            [
                Sphinx,
                Literal[
                    'module', 'class', 'exception', 'function', 'method', 'attribute'
                ],
                str,
                Any,
                dict[str, bool],
                str | None,
                str | None,
            ],
            tuple[str | None, str | None] | None,
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['autodoc-process-bases'],
        callback: Callable[[Sphinx, str, Any, dict[str, bool], list[str]], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['autodoc-skip-member'],
        callback: Callable[
            [
                Sphinx,
                Literal[
                    'module', 'class', 'exception', 'function', 'method', 'attribute'
                ],
                str,
                Any,
                bool,
                dict[str, bool],
            ],
            bool,
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['todo-defined'],
        callback: Callable[[Sphinx, todo_node], None],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['viewcode-find-source'],
        callback: Callable[
            [Sphinx, str],
            tuple[str, dict[str, tuple[Literal['class', 'def', 'other'], int, int]]],
        ],
        priority: int,
    ) -> int: ...

    @overload
    def connect(
        self,
        name: Literal['viewcode-follow-imported'],
        callback: Callable[[Sphinx, str, str], str | None],
        priority: int,
    ) -> int: ...

    # ---- Catch-all ---------------------------------------------------------

    @overload
    def connect(
        self,
        name: str,
        callback: Callable[..., Any],
        priority: int,
    ) -> int: ...

    def connect(self, name: str, callback: Callable, priority: int) -> int:
        """Connect a handler to specific event."""
        if name not in self.events:
            raise ExtensionError(__('Unknown event name: %s') % name)

        listener_id = self.next_listener_id
        self.next_listener_id += 1
        self.listeners[name].append(EventListener(listener_id, callback, priority))
        return listener_id

    def disconnect(self, listener_id: int) -> None:
        """Disconnect a handler."""
        for listeners in self.listeners.values():
            for listener in listeners.copy():
                if listener.id == listener_id:
                    listeners.remove(listener)

    def emit(
        self,
        name: str,
        *args: Any,
        allowed_exceptions: tuple[type[Exception], ...] = (),
    ) -> list:
        """Emit a Sphinx event."""
        # not every object likes to be repr()'d (think
        # random stuff coming via autodoc)
        try:
            repr_args = repr(args)
        except Exception:
            pass
        else:
            logger.debug('[app] emitting event: %r%s', name, repr_args)

        results = []
        listeners = sorted(self.listeners[name], key=attrgetter('priority'))
        for listener in listeners:
            try:
                results.append(listener.handler(self.app, *args))
            except allowed_exceptions:
                # pass through the errors specified as *allowed_exceptions*
                raise
            except SphinxError:
                raise
            except Exception as exc:
                if self.app.pdb:
                    # Just pass through the error, so that it can be debugged.
                    raise
                modname = safe_getattr(listener.handler, '__module__', None)
                raise ExtensionError(
                    __('Handler %r for event %r threw an exception')
                    % (listener.handler, name),
                    exc,
                    modname=modname,
                ) from exc
        return results

    def emit_firstresult(
        self,
        name: str,
        *args: Any,
        allowed_exceptions: tuple[type[Exception], ...] = (),
    ) -> Any:
        """Emit a Sphinx event and returns first result.

        This returns the result of the first handler that doesn't return ``None``.
        """
        for result in self.emit(name, *args, allowed_exceptions=allowed_exceptions):
            if result is not None:
                return result
        return None
