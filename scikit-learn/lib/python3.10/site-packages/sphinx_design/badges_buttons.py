from typing import Optional

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.util.docutils import ReferenceRole, SphinxRole

from sphinx_design.shared import SEMANTIC_COLORS, SdDirective, make_choice, text_align

ROLE_NAME_BADGE_PREFIX = "bdg"
ROLE_NAME_LINK_PREFIX = "bdg-link"
ROLE_NAME_REF_PREFIX = "bdg-ref"
DIRECTIVE_NAME_BUTTON_LINK = "button-link"
DIRECTIVE_NAME_BUTTON_REF = "button-ref"

# TODO defining arbitrary classes for badges
# (maybe split text right of last `;`, then split that by comma)
# in particular for rounded-pill class etc


def setup_badges_and_buttons(app: Sphinx) -> None:
    """Setup the badge components."""
    app.add_role(ROLE_NAME_BADGE_PREFIX, BadgeRole())
    app.add_role(ROLE_NAME_LINK_PREFIX, LinkBadgeRole())
    app.add_role(ROLE_NAME_REF_PREFIX, XRefBadgeRole())
    for color in SEMANTIC_COLORS:
        app.add_role("-".join((ROLE_NAME_BADGE_PREFIX, color)), BadgeRole(color))
        app.add_role(
            "-".join((ROLE_NAME_BADGE_PREFIX, color, "line")),
            BadgeRole(color, outline=True),
        )
        app.add_role("-".join((ROLE_NAME_LINK_PREFIX, color)), LinkBadgeRole(color))
        app.add_role(
            "-".join((ROLE_NAME_LINK_PREFIX, color, "line")),
            LinkBadgeRole(color, outline=True),
        )
        app.add_role("-".join((ROLE_NAME_REF_PREFIX, color)), XRefBadgeRole(color))
        app.add_role(
            "-".join((ROLE_NAME_REF_PREFIX, color, "line")),
            XRefBadgeRole(color, outline=True),
        )

    app.add_directive(DIRECTIVE_NAME_BUTTON_LINK, ButtonLinkDirective)
    app.add_directive(DIRECTIVE_NAME_BUTTON_REF, ButtonRefDirective)


def create_bdg_classes(color: Optional[str], outline: bool) -> list[str]:
    """Create the badge classes."""
    classes = [
        "sd-sphinx-override",
        "sd-badge",
    ]
    if color is None:
        return classes
    if outline:
        classes.extend([f"sd-outline-{color}", f"sd-text-{color}"])
    else:
        classes.extend([f"sd-bg-{color}", f"sd-bg-text-{color}"])
    return classes


class BadgeRole(SphinxRole):
    """Role to display a badge."""

    def __init__(self, color: Optional[str] = None, *, outline: bool = False) -> None:
        super().__init__()
        self.color = color
        self.outline = outline

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        node = nodes.inline(
            self.rawtext,
            self.text,
            classes=create_bdg_classes(self.color, self.outline),
        )
        self.set_source_info(node)
        return [node], []


class LinkBadgeRole(ReferenceRole):
    """Role to display a badge with an external link."""

    def __init__(self, color: Optional[str] = None, *, outline: bool = False) -> None:
        super().__init__()
        self.color = color
        self.outline = outline

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        node = nodes.reference(
            self.rawtext,
            refuri=self.target,
            classes=create_bdg_classes(self.color, self.outline),
        )
        # TODO open in new tab
        self.set_source_info(node)
        # if self.target != self.title:
        #     node["reftitle"] = self.target
        node += nodes.inline(self.title, self.title)
        return [node], []


class XRefBadgeRole(ReferenceRole):
    """Role to display a badge with an internal link."""

    def __init__(self, color: Optional[str] = None, *, outline: bool = False) -> None:
        super().__init__()
        self.color = color
        self.outline = outline

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        options = {
            "classes": create_bdg_classes(self.color, self.outline),
            "reftarget": self.target,
            "refdoc": self.env.docname,
            "refdomain": "",
            "reftype": "any",
            "refexplicit": self.has_explicit_title,
            "refwarn": True,
        }
        node = addnodes.pending_xref(self.rawtext, **options)
        self.set_source_info(node)
        node += nodes.inline(self.title, self.title, classes=["xref", "any"])
        return [node], []


class _ButtonDirective(SdDirective):
    """A base button directive."""

    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    has_content = True
    option_spec = {
        "color": make_choice(SEMANTIC_COLORS),
        "outline": directives.flag,
        "align": text_align,
        # expand to fit parent width
        "expand": directives.flag,
        # make parent also clickable
        "click-parent": directives.flag,
        "tooltip": directives.unchanged_required,
        "shadow": directives.flag,
        # ref button only
        "ref-type": make_choice(["any", "ref", "doc", "myst"]),
        "class": directives.class_option,
    }

    def create_ref_node(
        self, rawtext: str, target: str, explicit_title: bool, classes: list[str]
    ) -> nodes.Node:
        """Create the reference node."""
        raise NotImplementedError

    def run_with_defaults(self) -> list[nodes.Node]:
        rawtext = self.arguments[0]
        target = directives.uri(rawtext)
        classes = ["sd-sphinx-override", "sd-btn", "sd-text-wrap"]
        if "color" in self.options:
            if "outline" in self.options:
                classes.append(f"sd-btn-outline-{self.options['color']}")
            else:
                classes.append(f"sd-btn-{self.options['color']}")
        if "click-parent" in self.options:
            classes.append("sd-stretched-link")
        if "shadow" in self.options:
            classes.append("sd-shadow-sm")
        if "class" in self.options:
            classes.extend(self.options["class"])
        node = self.create_ref_node(rawtext, target, bool(self.content), classes)
        # TODO open in new tab
        self.set_source_info(node)
        if "tooltip" in self.options:
            node["reftitle"] = self.options["tooltip"]  # TODO escape HTML

        if self.content:
            textnodes, _ = self.state.inline_text(
                "\n".join(self.content), self.lineno + self.content_offset
            )
            content = nodes.inline("", "")
            content.extend(textnodes)
        else:
            content = nodes.inline(target, target)
        node.append(content)

        if "expand" in self.options:
            grid_container = nodes.inline(classes=["sd-d-grid"])
            self.set_source_info(grid_container)
            grid_container += node
            node = grid_container

        # `visit_reference` requires that a reference be inside a `TextElement` parent
        container = nodes.paragraph(classes=self.options.get("align", []))
        self.set_source_info(container)
        container += node

        return [container]


class ButtonLinkDirective(_ButtonDirective):
    """A button directive with an external link."""

    def create_ref_node(
        self, rawtext: str, target: str, explicit_title: bool, classes: list[str]
    ) -> nodes.Node:
        """Create the reference node."""
        return nodes.reference(
            rawtext,
            refuri=target,
            classes=classes,
        )


class ButtonRefDirective(_ButtonDirective):
    """A button directive with an internal link."""

    def create_ref_node(
        self, rawtext: str, target: str, explicit_title: bool, classes: list[str]
    ) -> nodes.Node:
        """Create the reference node."""
        ref_type = self.options.get("ref-type", "any")
        options = {
            # TODO the presence of classes raises an error if the link cannot be found
            "classes": classes,
            "reftarget": target,
            "refdoc": self.env.docname,
            "refdomain": "std" if ref_type in {"ref", "doc"} else "",
            "reftype": ref_type,
            "refexplicit": explicit_title,
            "refwarn": True,
        }
        return addnodes.pending_xref(rawtext, **options)
