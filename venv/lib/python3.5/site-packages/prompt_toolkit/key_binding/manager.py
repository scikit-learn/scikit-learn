"""
DEPRECATED:
Use `prompt_toolkit.key_binding.defaults.load_key_bindings` instead.

:class:`KeyBindingManager` is a utility (or shortcut) for loading all the key
bindings in a key binding registry, with a logic set of filters to quickly to
quickly change from Vi to Emacs key bindings at runtime.

You don't have to use this, but it's practical.

Usage::

    manager = KeyBindingManager()
    app = Application(key_bindings_registry=manager.registry)
"""
from __future__ import unicode_literals
from .defaults import load_key_bindings
from prompt_toolkit.filters import to_cli_filter
from prompt_toolkit.key_binding.registry import Registry, ConditionalRegistry, MergedRegistry

__all__ = (
    'KeyBindingManager',
)


class KeyBindingManager(object):
    """
    Utility for loading all key bindings into memory.

    :param registry: Optional `Registry` instance.
    :param enable_abort_and_exit_bindings: Filter to enable Ctrl-C and Ctrl-D.
    :param enable_system_bindings: Filter to enable the system bindings
            (meta-! prompt and Control-Z suspension.)
    :param enable_search: Filter to enable the search bindings.
    :param enable_open_in_editor: Filter to enable open-in-editor.
    :param enable_open_in_editor: Filter to enable open-in-editor.
    :param enable_extra_page_navigation: Filter for enabling extra page navigation.
        (Bindings for up/down scrolling through long pages, like in Emacs or Vi.)
    :param enable_auto_suggest_bindings: Filter to enable fish-style suggestions.

    :param enable_vi_mode: Deprecated!
    """
    def __init__(self,
                 registry=None,  # XXX: not used anymore.
                 enable_vi_mode=None,  # (`enable_vi_mode` is deprecated.)
                 enable_all=True,  #
                 get_search_state=None,
                 enable_abort_and_exit_bindings=False,
                 enable_system_bindings=False,
                 enable_search=False,
                 enable_open_in_editor=False,
                 enable_extra_page_navigation=False,
                 enable_auto_suggest_bindings=False):

        assert registry is None or isinstance(registry, Registry)
        assert get_search_state is None or callable(get_search_state)
        enable_all = to_cli_filter(enable_all)

        defaults = load_key_bindings(
             get_search_state=get_search_state,
             enable_abort_and_exit_bindings=enable_abort_and_exit_bindings,
             enable_system_bindings=enable_system_bindings,
             enable_search=enable_search,
             enable_open_in_editor=enable_open_in_editor,
             enable_extra_page_navigation=enable_extra_page_navigation,
             enable_auto_suggest_bindings=enable_auto_suggest_bindings)

        # Note, we wrap this whole thing again in a MergedRegistry, because we
        # don't want the `enable_all` settings to apply on items that were
        # added to the registry as a whole.
        self.registry = MergedRegistry([
            ConditionalRegistry(defaults, enable_all)
        ])

    @classmethod
    def for_prompt(cls, **kw):
        """
        Create a ``KeyBindingManager`` with the defaults for an input prompt.
        This activates the key bindings for abort/exit (Ctrl-C/Ctrl-D),
        incremental search and auto suggestions.

        (Not for full screen applications.)
        """
        kw.setdefault('enable_abort_and_exit_bindings', True)
        kw.setdefault('enable_search', True)
        kw.setdefault('enable_auto_suggest_bindings', True)

        return cls(**kw)

    def reset(self, cli):
        # For backwards compatibility.
        pass

    def get_vi_state(self, cli):
        # Deprecated!
        return cli.vi_state
