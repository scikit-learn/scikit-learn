from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.application import Sphinx
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx.util.logging import getLogger

from ._compat import findall
from .shared import WARNING_TYPE, SdDirective, create_component, is_component

LOGGER = getLogger(__name__)


def setup_tabs(app: Sphinx) -> None:
    app.add_directive("tab-set", TabSetDirective)
    app.add_directive("tab-item", TabItemDirective)
    app.add_directive("tab-set-code", TabSetCodeDirective)
    app.add_post_transform(TabSetHtmlTransform)
    app.add_node(sd_tab_input, html=(visit_tab_input, depart_tab_input))
    app.add_node(sd_tab_label, html=(visit_tab_label, depart_tab_label))


class TabSetDirective(SdDirective):
    """A container for a set of tab items."""

    has_content = True
    option_spec = {
        "sync-group": directives.unchanged_required,
        "class": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        self.assert_has_content()
        tab_set = create_component(
            "tab-set", classes=["sd-tab-set", *self.options.get("class", [])]
        )
        self.set_source_info(tab_set)
        self.state.nested_parse(self.content, self.content_offset, tab_set)
        for item in tab_set.children:
            if not is_component(item, "tab-item"):
                LOGGER.warning(
                    f"All children of a 'tab-set' "
                    f"should be 'tab-item' [{WARNING_TYPE}.tab]",
                    location=item,
                    type=WARNING_TYPE,
                    subtype="tab",
                )
                break
            if "sync_id" in item.children[0]:
                item.children[0]["sync_group"] = self.options.get("sync-group", "tab")
        return [tab_set]


class TabItemDirective(SdDirective):
    """A single tab item in a tab set.

    Note: This directive generates a single container,
    for the label and content::

        <container design_component="tab-item" has_title=True>
            <rubric>
                ...title nodes
            <container design_component="tab-content">
                ...content nodes

    This allows for a default rendering in non-HTML outputs.

    The ``TabHtmlTransform`` then transforms this container
    into the HTML specific structure.
    """

    required_arguments = 1  # the tab label is the first argument
    final_argument_whitespace = True
    has_content = True
    option_spec = {
        "selected": directives.flag,
        "sync": directives.unchanged_required,
        "name": directives.unchanged,
        "class-container": directives.class_option,
        "class-label": directives.class_option,
        "class-content": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        self.assert_has_content()
        if not is_component(self.state_machine.node, "tab-set"):
            LOGGER.warning(
                f"The parent of a 'tab-item' should be a 'tab-set' [{WARNING_TYPE}.tab]",
                location=(self.env.docname, self.lineno),
                type=WARNING_TYPE,
                subtype="tab",
            )
        tab_item = create_component(
            "tab-item",
            classes=["sd-tab-item", *self.options.get("class-container", [])],
            selected=("selected" in self.options),
        )

        # add tab label
        textnodes, _ = self.state.inline_text(self.arguments[0], self.lineno)
        tab_label = nodes.rubric(
            self.arguments[0],
            "",
            *textnodes,
            classes=["sd-tab-label", *self.options.get("class-label", [])],
        )
        if "sync" in self.options:
            tab_label["sync_id"] = self.options["sync"]
        self.add_name(tab_label)
        tab_item += tab_label

        # add tab content
        tab_content = create_component(
            "tab-content",
            classes=["sd-tab-content", *self.options.get("class-content", [])],
        )
        self.state.nested_parse(self.content, self.content_offset, tab_content)
        tab_item += tab_content

        return [tab_item]


class TabSetCodeDirective(SdDirective):
    """A container for a set of tab items, generated from code blocks."""

    has_content = True
    option_spec = {
        "no-sync": directives.flag,
        "sync-group": directives.unchanged_required,
        "class-set": directives.class_option,
        "class-item": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        self.assert_has_content()
        tab_set = create_component(
            "tab-set", classes=["sd-tab-set", *self.options.get("class-set", [])]
        )
        self.set_source_info(tab_set)
        self.state.nested_parse(self.content, self.content_offset, tab_set)
        new_children = []
        for item in tab_set.children:
            if not isinstance(item, nodes.literal_block):
                LOGGER.warning(
                    f"All children of a 'tab-code-set' "
                    f"should be a 'literal_block' [{WARNING_TYPE}.tab_code]",
                    location=tab_set,
                    type=WARNING_TYPE,
                    subtype="tab_code",
                )
                continue
            language = item.get("language", "unknown")
            tab_label = nodes.rubric(
                language.upper(),
                nodes.Text(language.upper()),
                classes=["sd-tab-label", *self.options.get("class-label", [])],
            )
            if "no-sync" not in self.options:
                tab_label["sync_group"] = self.options.get("sync-group", "code")
                tab_label["sync_id"] = language
            tab_content = create_component(
                "tab-content",
                children=[item],
                classes=["sd-tab-content", *self.options.get("class-content", [])],
            )
            tab_item = create_component(
                "tab-item",
                children=[tab_label, tab_content],
                classes=["sd-tab-item", *self.options.get("class-item", [])],
            )
            new_children.append(tab_item)
        tab_set.children = new_children
        return [tab_set]


class sd_tab_input(nodes.Element, nodes.General):  # noqa: N801
    pass


class sd_tab_label(nodes.TextElement, nodes.General):  # noqa: N801
    pass


def visit_tab_input(self, node):
    attributes = {"ids": [node["id"]], "type": node["type"], "name": node["set_id"]}
    if node["checked"]:
        attributes["checked"] = "checked"
    self.body.append(self.starttag(node, "input", **attributes))


def depart_tab_input(self, node):
    # note do not add a closing tag, since this is a void element:
    # https://developer.mozilla.org/en-US/docs/Glossary/Void_element
    pass


def visit_tab_label(self, node):
    attributes = {"for": node["input_id"]}
    if "sync_id" in node and "sync_group" in node:
        attributes["data-sync-id"] = node["sync_id"]
        attributes["data-sync-group"] = node["sync_group"]
    self.body.append(self.starttag(node, "label", **attributes))


def depart_tab_label(self, node):
    self.body.append("</label>")


class TabSetHtmlTransform(SphinxPostTransform):
    """Transform tab-set to HTML specific AST structure."""

    default_priority = 200
    formats = ("html",)

    def run(self) -> None:
        """Run the transform."""

        # setup id generators
        tab_set_id_base = "sd-tab-set-"
        tab_item_id_base = "sd-tab-item-"
        tab_item_id_num = 0

        for tab_set_id_num, tab_set in enumerate(
            findall(self.document)(lambda node: is_component(node, "tab-set"))
        ):
            tab_set_identity = tab_set_id_base + str(tab_set_id_num)
            children = []
            # get the first selected node
            selected_idx = None
            for idx, tab_item in enumerate(tab_set.children):
                if tab_item.get("selected", False):
                    if selected_idx is None:
                        selected_idx = idx
                    else:
                        LOGGER.warning(
                            f"Multiple selected 'tab-item' directives [{WARNING_TYPE}.tab]",
                            location=tab_item,
                            type=WARNING_TYPE,
                            subtype="tab",
                        )
            selected_idx = 0 if selected_idx is None else selected_idx

            for idx, tab_item in enumerate(tab_set.children):
                try:
                    tab_label, tab_content = tab_item.children
                except ValueError:
                    raise
                tab_item_identity = tab_item_id_base + str(tab_item_id_num)
                tab_item_id_num += 1

                # create: <input checked="checked" id="id" type="radio">
                input_node = sd_tab_input(
                    "",
                    id=tab_item_identity,
                    set_id=tab_set_identity,
                    type="radio",
                    checked=(idx == selected_idx),
                )
                input_node.source, input_node.line = tab_item.source, tab_item.line
                children.append(input_node)

                # create: <label for="id">...</label>
                label_node = sd_tab_label(
                    "",
                    "",
                    *tab_label.children,
                    input_id=tab_item_identity,
                    classes=tab_label["classes"],
                )
                if tab_label.get("ids"):
                    label_node["ids"] += tab_label["ids"]
                if "sync_group" in tab_label and "sync_id" in tab_label:
                    label_node["sync_group"] = tab_label["sync_group"]
                    label_node["sync_id"] = tab_label["sync_id"]
                label_node.source, label_node.line = tab_item.source, tab_item.line
                children.append(label_node)

                # add content
                children.append(tab_content)

            tab_set.children = children
