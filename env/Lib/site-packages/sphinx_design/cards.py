import re
from typing import NamedTuple, Optional

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import StringList
from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.logging import getLogger

from ._compat import findall
from .shared import (
    WARNING_TYPE,
    PassthroughTextElement,
    SdDirective,
    create_component,
    is_component,
    make_choice,
    margin_option,
    text_align,
)

LOGGER = getLogger(__name__)

DIRECTIVE_NAME_CARD = "card"
DIRECTIVE_NAME_CAROUSEL = "card-carousel"
REGEX_HEADER = re.compile(r"^\^{3,}\s*$")
REGEX_FOOTER = re.compile(r"^\+{3,}\s*$")


def setup_cards(app: Sphinx) -> None:
    """Setup the card components."""
    app.add_directive(DIRECTIVE_NAME_CARD, CardDirective)
    app.add_directive(DIRECTIVE_NAME_CAROUSEL, CardCarouselDirective)


class CardContent(NamedTuple):
    """Split card into header (optional), body, footer (optional).

    (offset, content)
    """

    body: tuple[int, StringList]
    header: Optional[tuple[int, StringList]] = None
    footer: Optional[tuple[int, StringList]] = None


class CardDirective(SdDirective):
    """A card component."""

    has_content = True
    required_arguments = 0
    optional_arguments = 1  # card title
    final_argument_whitespace = True
    option_spec = {
        "width": make_choice(["auto", "25%", "50%", "75%", "100%"]),
        "margin": margin_option,
        "text-align": text_align,
        "img-top": directives.uri,
        "img-bottom": directives.uri,
        "img-background": directives.uri,
        "img-alt": directives.unchanged,
        "link": directives.uri,
        "link-type": make_choice(["url", "any", "ref", "doc"]),
        "link-alt": directives.unchanged,
        "shadow": make_choice(["none", "sm", "md", "lg"]),
        "class-card": directives.class_option,
        "class-header": directives.class_option,
        "class-body": directives.class_option,
        "class-title": directives.class_option,
        "class-footer": directives.class_option,
        "class-img-top": directives.class_option,
        "class-img-bottom": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        return [self.create_card(self, self.arguments, self.options)]

    @classmethod
    def create_card(  # noqa: PLR0915
        cls, inst: SphinxDirective, arguments: Optional[list], options: dict
    ) -> nodes.Node:
        """Run the directive."""
        # TODO better degradation for latex
        card_classes = ["sd-card", "sd-sphinx-override"]
        if "width" in options:
            card_classes += [f'sd-w-{options["width"].rstrip("%")}']
        card_classes += options.get("margin", ["sd-mb-3"])
        card_classes += [f"sd-shadow-{options.get('shadow', 'sm')}"]
        if "link" in options:
            card_classes += ["sd-card-hover"]
        card = create_component(
            "card",
            card_classes
            + options.get("text-align", [])
            + options.get("class-card", []),
        )
        inst.set_source_info(card)

        img_alt = options.get("img-alt") or ""

        container = card
        if "img-background" in options:
            card.append(
                nodes.image(
                    uri=options["img-background"],
                    classes=["sd-card-img"],
                    alt=img_alt,
                )
            )
            overlay = create_component("card-overlay", ["sd-card-img-overlay"])
            inst.set_source_info(overlay)
            card += overlay
            container = overlay

        if "img-top" in options:
            image_top = nodes.image(
                "",
                uri=options["img-top"],
                alt=img_alt,
                classes=["sd-card-img-top", *options.get("class-img-top", [])],
            )
            container.append(image_top)

        components = cls.split_content(inst.content, inst.content_offset)

        if components.header:
            container.append(
                cls._create_component(
                    inst, "header", options, components.header[0], components.header[1]
                )
            )

        body = cls._create_component(
            inst, "body", options, components.body[0], components.body[1]
        )
        if arguments:
            title = create_component(
                "card-title",
                [
                    "sd-card-title",
                    "sd-font-weight-bold",
                    *options.get("class-title", []),
                ],
            )
            textnodes, _ = inst.state.inline_text(arguments[0], inst.lineno)
            title_container = PassthroughTextElement()
            title_container.extend(textnodes)
            inst.set_source_info(title_container)
            title.append(title_container)
            body.insert(0, title)
        container.append(body)

        if components.footer:
            container.append(
                cls._create_component(
                    inst, "footer", options, components.footer[0], components.footer[1]
                )
            )

        if "img-bottom" in options:
            image_bottom = nodes.image(
                "",
                uri=options["img-bottom"],
                alt=img_alt,
                classes=["sd-card-img-bottom", *options.get("class-img-bottom", [])],
            )
            container.append(image_bottom)

        if "link" in options:
            link_container = PassthroughTextElement()
            _classes = ["sd-stretched-link", "sd-hide-link-text"]
            _rawtext = options.get("link-alt") or options["link"]
            if options.get("link-type", "url") == "url":
                link = nodes.reference(
                    _rawtext,
                    "",
                    nodes.inline(_rawtext, _rawtext),
                    refuri=options["link"],
                    classes=_classes,
                )
            else:
                options = {
                    # TODO the presence of classes raises an error if the link cannot be found
                    "classes": _classes,
                    "reftarget": options["link"],
                    "refdoc": inst.env.docname,
                    "refdomain": "" if options["link-type"] == "any" else "std",
                    "reftype": options["link-type"],
                    "refexplicit": "link-alt" in options,
                    "refwarn": True,
                }
                link = addnodes.pending_xref(
                    _rawtext, nodes.inline(_rawtext, _rawtext), **options
                )
            inst.set_source_info(link)
            link_container += link
            container.append(link_container)

        return card

    @staticmethod
    def split_content(content: StringList, offset: int) -> CardContent:
        """Split the content into header, body and footer."""
        header_index, footer_index, header, footer = None, None, None, None
        body_offset = offset
        for index, line in enumerate(content):
            # match the first occurrence of a header regex
            if (header_index is None) and REGEX_HEADER.match(line):
                header_index = index
            # match the final occurrence of a footer regex
            if REGEX_FOOTER.match(line):
                footer_index = index
        if header_index is not None:
            header = (offset, content[:header_index])
            body_offset += header_index + 1
        if footer_index is not None:
            footer = (offset + footer_index + 1, content[footer_index + 1 :])
        body = (
            body_offset,
            content[
                (header_index + 1 if header_index is not None else None) : footer_index
            ],
        )
        return CardContent(body, header, footer)

    @classmethod
    def _create_component(
        cls,
        inst: SphinxDirective,
        name: str,
        options: dict,
        offset: int,
        content: StringList,
    ) -> nodes.container:
        """Create the header, body, or footer."""
        component = create_component(
            f"card-{name}", [f"sd-card-{name}", *options.get(f"class-{name}", [])]
        )
        inst.set_source_info(component)  # TODO set proper lines
        inst.state.nested_parse(content, offset, component)
        cls.add_card_child_classes(component)
        return component

    @staticmethod
    def add_card_child_classes(node):
        """Add classes to specific child nodes."""
        for para in findall(node)(nodes.paragraph):
            para["classes"] = [*para.get("classes", []), "sd-card-text"]
        # for title in findall(node)(nodes.title):
        #     title["classes"] = ([] if "classes" not in title else title["classes"]) + [
        #         "sd-card-title"
        #     ]


class CardCarouselDirective(SdDirective):
    """A component, which is a container for cards in a single scrollable row."""

    has_content = True
    required_arguments = 1  # columns
    optional_arguments = 0
    option_spec = {
        "class": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        self.assert_has_content()
        try:
            cols = make_choice([str(i) for i in range(1, 13)])(
                self.arguments[0].strip()
            )
        except ValueError as exc:
            raise self.error(f"Invalid directive argument: {exc}") from exc
        container = create_component(
            "card-carousel",
            [
                "sd-sphinx-override",
                "sd-cards-carousel",
                f"sd-card-cols-{cols}",
                *self.options.get("class", []),
            ],
        )
        self.set_source_info(container)
        self.state.nested_parse(self.content, self.content_offset, container)
        for item in container.children:
            if not is_component(item, "card"):
                LOGGER.warning(
                    "All children of a 'card-carousel' "
                    f"should be 'card' [{WARNING_TYPE}.card]",
                    location=item,
                    type=WARNING_TYPE,
                    subtype="card",
                )
                break
        return [container]
