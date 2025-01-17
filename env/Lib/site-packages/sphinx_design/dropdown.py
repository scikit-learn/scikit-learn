"""Originally Adapted from sphinxcontrib.details.directive"""

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.application import Sphinx
from sphinx.transforms.post_transforms import SphinxPostTransform

from sphinx_design.shared import (
    SEMANTIC_COLORS,
    SdDirective,
    create_component,
    is_component,
    make_choice,
    margin_option,
)

from ._compat import findall
from .icons import get_octicon, list_octicons


def setup_dropdown(app: Sphinx) -> None:
    app.add_node(dropdown_main, html=(visit_dropdown_main, depart_dropdown_main))
    app.add_node(dropdown_title, html=(visit_dropdown_title, depart_dropdown_title))
    app.add_directive("dropdown", DropdownDirective)
    app.add_post_transform(DropdownHtmlTransform)


class dropdown_main(nodes.Element, nodes.General):  # noqa: N801
    pass


class dropdown_title(nodes.TextElement, nodes.General):  # noqa: N801
    pass


def visit_dropdown_main(self, node):
    if node.get("opened"):
        self.body.append(self.starttag(node, "details", open="open"))
    else:
        self.body.append(self.starttag(node, "details"))


def depart_dropdown_main(self, node):
    self.body.append("</details>")


def visit_dropdown_title(self, node):
    self.body.append(self.starttag(node, "summary"))


def depart_dropdown_title(self, node):
    self.body.append("</summary>")


class DropdownDirective(SdDirective):
    """A directive to generate a collapsible container.

    Note: This directive generates a single container,
    for the title (optional) and content::

        <container design_component="dropdown" has_title=True>
            <rubric>
                ...title nodes
        ...content nodes

    This allows for a default rendering in non-HTML outputs.

    The ``DropdownHtmlTransform`` then transforms this container
    into the HTML specific structure.
    """

    optional_arguments = 1  # title of dropdown
    final_argument_whitespace = True
    has_content = True
    option_spec = {
        "open": directives.flag,  # make open by default
        "color": make_choice(SEMANTIC_COLORS),
        "icon": make_choice(list_octicons()),
        "chevron": make_choice(
            ["right-down", "down-up"]
        ),  # chevron direction closed-open
        "animate": make_choice(("fade-in", "fade-in-slide-down")),
        "margin": margin_option,
        "name": directives.unchanged,
        "class-container": directives.class_option,
        "class-title": directives.class_option,
        "class-body": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        # default classes
        classes = {
            "container_classes": self.options.get("margin", ["sd-mb-3"])
            + self.options.get("class-container", []),
            "title_classes": self.options.get("class-title", []),
            "body_classes": self.options.get("class-body", []),
        }

        # add color classes
        title_color = self.options.get("color")
        if title_color:
            classes["title_classes"].extend(
                [f"sd-bg-{title_color}", f"sd-bg-text-{title_color}"]
            )

        # add animation classes
        if (
            "animate" in self.options
            and ("sd-" + self.options["animate"]) not in classes["container_classes"]
        ):
            classes["container_classes"].append("sd-" + self.options["animate"])

        container = create_component(
            "dropdown",
            opened="open" in self.options,
            type="dropdown",
            has_title=len(self.arguments) > 0,
            icon=self.options.get("icon"),
            chevron=self.options.get("chevron"),
            **classes,
        )
        self.set_source_info(container)
        if self.arguments:
            textnodes, messages = self.state.inline_text(self.arguments[0], self.lineno)
            title_node = nodes.rubric(self.arguments[0], "", *textnodes)
            container += title_node
            container += messages
            # where possible we add the target to the title node,
            # so that it can be used as the reference text
            self.add_name(title_node)
        else:
            self.add_name(container)
        self.state.nested_parse(self.content, self.content_offset, container)
        return [container]


class DropdownHtmlTransform(SphinxPostTransform):
    """Transform dropdown containers into the HTML specific AST structures::

    <details class="sd-sphinx-override sd-dropdown sd-card">
        <summary class="sd-summary-title sd-card-header">
            ...title nodes
        <div class="sd-summary-content sd-card-body">
            ...content nodes

    """

    default_priority = 199
    formats = ("html",)

    def run(self) -> None:
        """Run the transform"""
        document: nodes.document = self.document
        for node in findall(document)(lambda node: is_component(node, "dropdown")):
            # TODO option to not have card css (but requires more formatting)
            use_card = True

            marker_type = (
                "chevron-down" if node["chevron"] == "down-up" else "chevron-right"
            )
            state_marker = nodes.inline(
                "",
                "",
                nodes.raw(
                    "",
                    nodes.Text(get_octicon(marker_type, height="1.5em")),
                    format="html",
                ),
                classes=["sd-summary-state-marker", f"sd-summary-{marker_type}"],
            )

            newnode = dropdown_main(
                opened=node["opened"],
                classes=["sd-sphinx-override", "sd-dropdown"]
                + (["sd-card"] if use_card else ["sd-d-flex-column"])
                + node["container_classes"],
            )

            if node["has_title"]:
                title_text_children = node[0].children
                if node[0].get("ids"):
                    newnode["ids"] += node[0]["ids"]
                body_children = node[1:]
            else:
                title_text_children = [
                    nodes.raw(
                        "...",
                        nodes.Text(
                            get_octicon(
                                "kebab-horizontal", height="1.5em", classes=["no-title"]
                            )
                        ),
                        format="html",
                    )
                ]
                body_children = node.children
            title_text_node = nodes.inline(
                "",
                "",
                *title_text_children,
                classes=["sd-summary-text"],
            )
            title_children = [title_text_node, state_marker]
            if node["icon"]:
                title_children.insert(
                    0,
                    nodes.raw(
                        "",
                        get_octicon(node["icon"], height="1em"),
                        classes=["sd-summary-icon"],
                        format="html",
                    ),
                )

            newnode += dropdown_title(
                "",
                "",
                *title_children,
                classes=["sd-summary-title"]
                + (["sd-card-header"] if use_card else [])
                + node["title_classes"],
            )
            body_node = create_component(
                "dropdown-body",
                classes=["sd-summary-content"]
                + (["sd-card-body"] if use_card else [])
                + node["body_classes"],
                children=body_children,
            )
            if use_card:
                for para in findall(body_node)(nodes.paragraph):
                    para["classes"] = ([] if "classes" in para else para["classes"]) + [
                        "sd-card-text"
                    ]
            newnode += body_node
            # newnode += open_marker
            node.replace_self(newnode)
