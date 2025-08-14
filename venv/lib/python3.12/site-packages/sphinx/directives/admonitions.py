from __future__ import annotations

from typing import TYPE_CHECKING

from docutils import nodes
from docutils.parsers.rst.directives.admonitions import BaseAdmonition

from sphinx import addnodes
from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from typing import ClassVar

    from docutils.nodes import Node

    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata, OptionSpec


def _collapsible_arg(argument: str | None) -> str:
    if argument is None:
        return "open"
    if (value := argument.lower().strip()) in {"open", "closed"}:
        return value
    msg = f'"{argument}" unknown; choose from "open" or "closed".'
    raise ValueError(msg)


class SphinxAdmonition(BaseAdmonition, SphinxDirective):
    option_spec: ClassVar[OptionSpec] = BaseAdmonition.option_spec.copy()  # type: ignore[union-attr]
    option_spec |= {
        "collapsible": _collapsible_arg,
    }

    node_class: type[nodes.Admonition] = nodes.admonition
    """Subclasses must set this to the appropriate admonition node class."""

    def run(self) -> list[Node]:
        (admonition_node,) = super().run()
        return [admonition_node]


class Admonition(SphinxAdmonition):
    required_arguments = 1
    node_class = nodes.admonition


class Attention(SphinxAdmonition):
    node_class = nodes.attention


class Caution(SphinxAdmonition):
    node_class = nodes.caution


class Danger(SphinxAdmonition):
    node_class = nodes.danger


class Error(SphinxAdmonition):
    node_class = nodes.error


class Hint(SphinxAdmonition):
    node_class = nodes.hint


class Important(SphinxAdmonition):
    node_class = nodes.important


class Note(SphinxAdmonition):
    node_class = nodes.note


class Tip(SphinxAdmonition):
    node_class = nodes.tip


class Warning(SphinxAdmonition):
    node_class = nodes.warning


class SeeAlso(SphinxAdmonition):
    """An admonition mentioning things to look at as reference."""

    node_class = addnodes.seealso


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive("admonition", Admonition, override=True)
    app.add_directive("attention", Attention, override=True)
    app.add_directive("caution", Caution, override=True)
    app.add_directive("danger", Danger, override=True)
    app.add_directive("error", Error, override=True)
    app.add_directive("hint", Hint, override=True)
    app.add_directive("important", Important, override=True)
    app.add_directive("note", Note, override=True)
    app.add_directive("tip", Tip, override=True)
    app.add_directive("warning", Warning, override=True)
    app.add_directive("seealso", SeeAlso, override=True)

    return {
        "version": "builtin",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
