"""
Default key bindings.::

    registry = load_key_bindings()
    app = Application(key_bindings_registry=registry)
"""
from __future__ import unicode_literals
from prompt_toolkit.key_binding.registry import ConditionalRegistry, MergedRegistry
from prompt_toolkit.key_binding.bindings.basic import load_basic_bindings, load_abort_and_exit_bindings, load_basic_system_bindings, load_auto_suggestion_bindings, load_mouse_bindings
from prompt_toolkit.key_binding.bindings.emacs import load_emacs_bindings, load_emacs_system_bindings, load_emacs_search_bindings, load_emacs_open_in_editor_bindings, load_extra_emacs_page_navigation_bindings
from prompt_toolkit.key_binding.bindings.vi import load_vi_bindings, load_vi_system_bindings, load_vi_search_bindings, load_vi_open_in_editor_bindings, load_extra_vi_page_navigation_bindings
from prompt_toolkit.filters import to_cli_filter

__all__ = (
    'load_key_bindings',
    'load_key_bindings_for_prompt',
)


def load_key_bindings(
        get_search_state=None,
        enable_abort_and_exit_bindings=False,
        enable_system_bindings=False,
        enable_search=False,
        enable_open_in_editor=False,
        enable_extra_page_navigation=False,
        enable_auto_suggest_bindings=False):
    """
    Create a Registry object that contains the default key bindings.

    :param enable_abort_and_exit_bindings: Filter to enable Ctrl-C and Ctrl-D.
    :param enable_system_bindings: Filter to enable the system bindings (meta-!
            prompt and Control-Z suspension.)
    :param enable_search: Filter to enable the search bindings.
    :param enable_open_in_editor: Filter to enable open-in-editor.
    :param enable_open_in_editor: Filter to enable open-in-editor.
    :param enable_extra_page_navigation: Filter for enabling extra page
        navigation. (Bindings for up/down scrolling through long pages, like in
        Emacs or Vi.)
    :param enable_auto_suggest_bindings: Filter to enable fish-style suggestions.
    """

    assert get_search_state is None or callable(get_search_state)

    # Accept both Filters and booleans as input.
    enable_abort_and_exit_bindings = to_cli_filter(enable_abort_and_exit_bindings)
    enable_system_bindings = to_cli_filter(enable_system_bindings)
    enable_search = to_cli_filter(enable_search)
    enable_open_in_editor = to_cli_filter(enable_open_in_editor)
    enable_extra_page_navigation = to_cli_filter(enable_extra_page_navigation)
    enable_auto_suggest_bindings = to_cli_filter(enable_auto_suggest_bindings)

    registry = MergedRegistry([
        # Load basic bindings.
        load_basic_bindings(),
        load_mouse_bindings(),

        ConditionalRegistry(load_abort_and_exit_bindings(),
                            enable_abort_and_exit_bindings),

        ConditionalRegistry(load_basic_system_bindings(),
                            enable_system_bindings),

        # Load emacs bindings.
        load_emacs_bindings(),

        ConditionalRegistry(load_emacs_open_in_editor_bindings(),
                            enable_open_in_editor),

        ConditionalRegistry(load_emacs_search_bindings(get_search_state=get_search_state),
                            enable_search),

        ConditionalRegistry(load_emacs_system_bindings(),
                            enable_system_bindings),

        ConditionalRegistry(load_extra_emacs_page_navigation_bindings(),
                            enable_extra_page_navigation),

        # Load Vi bindings.
        load_vi_bindings(get_search_state=get_search_state),

        ConditionalRegistry(load_vi_open_in_editor_bindings(),
                            enable_open_in_editor),

        ConditionalRegistry(load_vi_search_bindings(get_search_state=get_search_state),
                            enable_search),

        ConditionalRegistry(load_vi_system_bindings(),
                            enable_system_bindings),

        ConditionalRegistry(load_extra_vi_page_navigation_bindings(),
                            enable_extra_page_navigation),

        # Suggestion bindings.
        # (This has to come at the end, because the Vi bindings also have an
        # implementation for the "right arrow", but we really want the
        # suggestion binding when a suggestion is available.)
        ConditionalRegistry(load_auto_suggestion_bindings(),
                            enable_auto_suggest_bindings),
    ])

    return registry


def load_key_bindings_for_prompt(**kw):
    """
    Create a ``Registry`` object with the defaults key bindings for an input
    prompt.

    This activates the key bindings for abort/exit (Ctrl-C/Ctrl-D),
    incremental search and auto suggestions.

    (Not for full screen applications.)
    """
    kw.setdefault('enable_abort_and_exit_bindings', True)
    kw.setdefault('enable_search', True)
    kw.setdefault('enable_auto_suggest_bindings', True)

    return load_key_bindings(**kw)
