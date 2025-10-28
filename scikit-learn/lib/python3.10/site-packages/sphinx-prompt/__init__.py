#!/usr/bin/python

from typing import Any

import sphinx.application
from docutils import nodes
from docutils.parsers import rst
from docutils.parsers.rst import directives
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import BashLexer, BatchLexer, PowerShellLexer, PythonLexer, ScalaLexer, TextLexer
from sphinx.util import logging

logger = logging.getLogger(__name__)


class PromptCache:
    """The cache of different prompt."""

    def __init__(self) -> None:
        """Initialize."""
        self.next_index = 1
        self.prompts: dict[str, int] = {}

    def clear(self, *args: Any) -> None:
        """Clear all cache."""
        del args
        self.next_index = 1
        self.prompts = {}

    def register_prompt(self, prompt: str) -> str:
        """Initialize the prompts."""
        if prompt in self.prompts:
            return ""
        else:
            index = self.next_index
            self.next_index = index + 1
            self.prompts[prompt] = index
            return f"""span.prompt{index}:before {{
  content: "{prompt} ";
}}
"""

    def get_prompt_class(self, prompt: str) -> str:
        """Get the CSS class name."""
        return f"prompt{self.prompts[prompt]}"


_cache = PromptCache()
PROMPTS = {
    "bash": "$",
    "batch": r"C:\\>",
    "powershell": r"PS C:\\>",
}
LEXERS = {
    "bash": BashLexer,
    "batch": BatchLexer,
    "powershell": PowerShellLexer,
    "python": PythonLexer,
    "scala": ScalaLexer,
}


class PromptDirective(rst.Directive):
    """The prompt directive."""

    optional_arguments = 3
    option_spec = {
        "language": directives.unchanged_required,
        "prompts": directives.unchanged_required,
        "modifiers": directives.unchanged_required,
    }
    has_content = True

    def run(self) -> list[nodes.raw]:
        """Run the directive."""
        self.assert_has_content()

        arg_count = len(self.arguments)

        for idx, option_name in enumerate(("language", "prompts", "modifiers")):
            if arg_count > idx:
                if self.options.get(option_name):
                    logger.warning(
                        "%s is already passed as an option, ignoring the value passed"
                        " as positional argument and all arguments that come after it.",
                        option_name,
                        location=(self.state.document.settings.env.docname, self.lineno),
                    )
                    break
                else:
                    self.options[option_name] = self.arguments[idx]

        language: str = self.options.get("language") or "text"
        prompt: str = self.options.get("prompts") or PROMPTS.get(language, "")
        modifiers: list[str] = self.options.get("modifiers", "").split(",")
        if "auto" in modifiers:
            prompts: list[str] = prompt.split(",")

        html = '<div class="highlight-default notranslate"><div class="highlight"><pre>'
        styles = ""
        if "auto" in modifiers:
            for prompt in prompts:
                styles += _cache.register_prompt(prompt)
        else:
            if prompt is not None:
                styles += _cache.register_prompt(prompt)
        if styles:
            html += '<style type="text/css">\n' + styles + "</style>"
        latex = "\\begin{Verbatim}[commandchars=\\\\\\{\\}]"

        Lexer = LEXERS.get(language, TextLexer)  # noqa: N806, pylint: disable=invalid-name

        statement: list[str] = []
        if "auto" in modifiers:
            prompt_class = ""
            for line in self.content:
                latex += "\n" + line

                for prompt in prompts:
                    if line.startswith(prompt):
                        if len(statement) > 0:
                            highlighted_line = highlight(
                                "\n".join(statement), Lexer(), HtmlFormatter(nowrap=True)
                            ).strip("\r\n")
                            html += f'<span class="{prompt_class}">{highlighted_line}</span>\n'
                            statement = []
                        line = line[len(prompt) + 1 :].rstrip()
                        prompt_class = _cache.get_prompt_class(prompt)
                        break

                statement.append(line)
            # Add last prompt
            if len(statement) > 0:
                highlighted_line = highlight("\n".join(statement), Lexer(), HtmlFormatter(nowrap=True)).strip(
                    "\r\n"
                )
                html += f'<span class="{prompt_class}">{highlighted_line}</span>\n'
        elif language in ["bash", "python"]:
            for line in self.content:
                statement.append(line)
                highlighted_line = highlight("\n".join(statement), Lexer(), HtmlFormatter(nowrap=True)).strip(
                    "\r\n"
                )
                if len(line) == 0 or not line[-1] == "\\":
                    html += f'<span class="{_cache.get_prompt_class(prompt)}">{highlighted_line}</span>\n'
                    if prompt is not None:
                        statements = "\n".join(statement)
                        latex += f"\n{prompt} {statements}"
                    else:
                        latex += "\n" + "\n".join(statement)
                    statement = []
        else:
            for line in self.content:
                highlighted_line = highlight(line, Lexer(), HtmlFormatter(nowrap=True)).strip("\r\n")
                html += f'<span class="{_cache.get_prompt_class(prompt)}">{highlighted_line}</span>\n'
                if prompt is not None:
                    latex += f"\n{prompt} {line}"
                else:
                    latex += "\n" + line

        html += "</pre></div></div>"
        latex += "\n\\end{Verbatim}"

        return [
            nodes.raw("\n".join(self.content), html, format="html"),
            nodes.raw("\n".join(self.content), latex, format="latex"),
        ]


def setup(app: sphinx.application.Sphinx) -> dict[str, bool]:
    """Register the plugin."""
    app.add_directive("prompt", PromptDirective)
    app.connect("env-purge-doc", _cache.clear)
    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
