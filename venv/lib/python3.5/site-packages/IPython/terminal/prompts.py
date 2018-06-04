"""Terminal input and output prompts."""

from pygments.token import Token
import sys

from IPython.core.displayhook import DisplayHook

from prompt_toolkit.layout.utils import token_list_width

class Prompts(object):
    def __init__(self, shell):
        self.shell = shell

    def in_prompt_tokens(self, cli=None):
        return [
            (Token.Prompt, 'In ['),
            (Token.PromptNum, str(self.shell.execution_count)),
            (Token.Prompt, ']: '),
        ]

    def _width(self):
        return token_list_width(self.in_prompt_tokens())

    def continuation_prompt_tokens(self, cli=None, width=None):
        if width is None:
            width = self._width()
        return [
            (Token.Prompt, (' ' * (width - 5)) + '...: '),
        ]

    def rewrite_prompt_tokens(self):
        width = self._width()
        return [
            (Token.Prompt, ('-' * (width - 2)) + '> '),
        ]

    def out_prompt_tokens(self):
        return [
            (Token.OutPrompt, 'Out['),
            (Token.OutPromptNum, str(self.shell.execution_count)),
            (Token.OutPrompt, ']: '),
        ]

class ClassicPrompts(Prompts):
    def in_prompt_tokens(self, cli=None):
        return [
            (Token.Prompt, '>>> '),
        ]

    def continuation_prompt_tokens(self, cli=None, width=None):
        return [
            (Token.Prompt, '... ')
        ]

    def rewrite_prompt_tokens(self):
        return []

    def out_prompt_tokens(self):
        return []

class RichPromptDisplayHook(DisplayHook):
    """Subclass of base display hook using coloured prompt"""
    def write_output_prompt(self):
        sys.stdout.write(self.shell.separate_out)
        # If we're not displaying a prompt, it effectively ends with a newline,
        # because the output will be left-aligned.
        self.prompt_end_newline = True

        if self.do_full_cache:
            tokens = self.shell.prompts.out_prompt_tokens()
            prompt_txt = ''.join(s for t, s in tokens)
            if prompt_txt and not prompt_txt.endswith('\n'):
                # Ask for a newline before multiline output
                self.prompt_end_newline = False

            if self.shell.pt_cli:
                self.shell.pt_cli.print_tokens(tokens)
            else:
                sys.stdout.write(prompt_txt)
