"""
Completer for a regular grammar.
"""
from __future__ import unicode_literals

from prompt_toolkit.completion import Completer, Completion
from prompt_toolkit.document import Document

from .compiler import _CompiledGrammar

__all__ = (
    'GrammarCompleter',
)


class GrammarCompleter(Completer):
    """
    Completer which can be used for autocompletion according to variables in
    the grammar. Each variable can have a different autocompleter.

    :param compiled_grammar: `GrammarCompleter` instance.
    :param completers: `dict` mapping variable names of the grammar to the
                       `Completer` instances to be used for each variable.
    """
    def __init__(self, compiled_grammar, completers):
        assert isinstance(compiled_grammar, _CompiledGrammar)
        assert isinstance(completers, dict)

        self.compiled_grammar = compiled_grammar
        self.completers = completers

    def get_completions(self, document, complete_event):
        m = self.compiled_grammar.match_prefix(document.text_before_cursor)

        if m:
            completions = self._remove_duplicates(
                self._get_completions_for_match(m, complete_event))

            for c in completions:
                yield c

    def _get_completions_for_match(self, match, complete_event):
        """
        Yield all the possible completions for this input string.
        (The completer assumes that the cursor position was at the end of the
        input string.)
        """
        for match_variable in match.end_nodes():
            varname = match_variable.varname
            start = match_variable.start

            completer = self.completers.get(varname)

            if completer:
                text = match_variable.value

                # Unwrap text.
                unwrapped_text = self.compiled_grammar.unescape(varname, text)

                # Create a document, for the completions API (text/cursor_position)
                document = Document(unwrapped_text, len(unwrapped_text))

                # Call completer
                for completion in completer.get_completions(document, complete_event):
                    new_text = unwrapped_text[:len(text) + completion.start_position] + completion.text

                    # Wrap again.
                    yield Completion(
                        text=self.compiled_grammar.escape(varname, new_text),
                        start_position=start - len(match.string),
                        display=completion.display,
                        display_meta=completion.display_meta)

    def _remove_duplicates(self, items):
        """
        Remove duplicates, while keeping the order.
        (Sometimes we have duplicates, because the there several matches of the
        same grammar, each yielding similar completions.)
        """
        result = []
        for i in items:
            if i not in result:
                result.append(i)
        return result
