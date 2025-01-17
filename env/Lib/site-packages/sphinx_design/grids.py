from typing import Optional

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.application import Sphinx
from sphinx.util.logging import getLogger

from .cards import CardDirective
from .shared import (
    WARNING_TYPE,
    SdDirective,
    create_component,
    is_component,
    make_choice,
    margin_option,
    padding_option,
    text_align,
)

LOGGER = getLogger(__name__)


DIRECTIVE_NAME_GRID = "grid"
DIRECTIVE_NAME_GRID_ITEM = "grid-item"
DIRECTIVE_NAME_GRID_ITEM_CARD = "grid-item-card"


def setup_grids(app: Sphinx):
    """Setup the grid components."""
    app.add_directive(DIRECTIVE_NAME_GRID, GridDirective)
    app.add_directive(DIRECTIVE_NAME_GRID_ITEM, GridItemDirective)
    app.add_directive(DIRECTIVE_NAME_GRID_ITEM_CARD, GridItemCardDirective)


def _media_option(
    argument: Optional[str],
    prefix: str,
    *,
    allow_auto: bool = False,
    min_num: int = 1,
    max_num: int = 12,
) -> list[str]:
    """Validate the number of columns (out of 12).

    One or four integers (for "xs sm md lg") between 1 and 12.
    """
    validate_error_msg = (
        "argument must be 1 or 4 (xs sm md lg) values, and each value should be "
        f"{'either auto or ' if allow_auto else ''}an integer from {min_num} to {max_num}"
    )
    if argument is None:
        raise ValueError(validate_error_msg)
    values = argument.strip().split()
    if len(values) == 1:
        values = [values[0], values[0], values[0], values[0]]
    if len(values) != 4:
        raise ValueError(validate_error_msg)
    for value in values:
        if allow_auto and value == "auto":
            continue
        try:
            int_value = int(value)
        except Exception as exc:
            raise ValueError(validate_error_msg) from exc
        if not (min_num <= int_value <= max_num):
            raise ValueError(validate_error_msg)
    return [f"{prefix}{values[0]}"] + [
        f"{prefix}{size}-{value}"
        for size, value in zip(["xs", "sm", "md", "lg"], values)
    ]


def row_columns_option(argument: Optional[str]) -> list[str]:
    """Validate the number of columns (out of 12) a grid row will have.

    One or four integers (for "xs sm md lg") between 1 and 12  (or 'auto').
    """
    return _media_option(argument, "sd-row-cols-", allow_auto=True)


def item_columns_option(argument: Optional[str]) -> list[str]:
    """Validate the number of columns (out of 12) a grid-item will take up.

    One or four integers (for "xs sm md lg") between 1 and 12 (or 'auto').
    """
    return _media_option(argument, "sd-col-", allow_auto=True)


def gutter_option(argument: Optional[str]) -> list[str]:
    """Validate the gutter size between grid items.

    One or four integers (for "xs sm md lg") between 0 and 5.
    """
    return _media_option(argument, "sd-g-", min_num=0, max_num=5)


class GridDirective(SdDirective):
    """A grid component, which is a container for grid items (i.e. columns)."""

    has_content = True
    required_arguments = 0
    optional_arguments = 1  # columns
    final_argument_whitespace = True
    option_spec = {
        "gutter": gutter_option,
        "margin": margin_option,
        "padding": padding_option,
        "outline": directives.flag,
        "reverse": directives.flag,
        "class-container": directives.class_option,
        "class-row": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        try:
            column_classes = (
                row_columns_option(self.arguments[0]) if self.arguments else []
            )
        except ValueError as exc:
            raise self.error(f"Invalid directive argument: {exc}") from exc
        self.assert_has_content()
        # container-fluid is 100% width for all breakpoints,
        # rather than the fixed width of the breakpoint (like container)
        grid_classes = ["sd-container-fluid", "sd-sphinx-override"]
        container = create_component(
            "grid-container",
            grid_classes
            + self.options.get("margin", ["sd-mb-4"])
            + self.options.get("padding", [])
            + (["sd-border-1"] if "outline" in self.options else [])
            + self.options.get("class-container", []),
        )
        self.set_source_info(container)
        row = create_component(
            "grid-row",
            ["sd-row"]
            + column_classes
            + self.options.get("gutter", [])
            + (["sd-flex-row-reverse"] if "reverse" in self.options else [])
            + self.options.get("class-row", []),
        )
        self.set_source_info(row)
        container += row
        self.state.nested_parse(self.content, self.content_offset, row)
        # each item in a row should be a column
        for item in row.children:
            if not is_component(item, "grid-item"):
                LOGGER.warning(
                    f"All children of a 'grid-row' "
                    f"should be 'grid-item' [{WARNING_TYPE}.grid]",
                    location=item,
                    type=WARNING_TYPE,
                    subtype="grid",
                )
                break
        return [container]


class GridItemDirective(SdDirective):
    """An item within a grid row.

    Can "occupy" 1 to 12 columns.
    """

    has_content = True
    option_spec = {
        "columns": item_columns_option,
        "margin": margin_option,
        "padding": padding_option,
        "child-direction": make_choice(["column", "row"]),
        "child-align": make_choice(["start", "end", "center", "justify", "spaced"]),
        "outline": directives.flag,
        "class": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        if not is_component(self.state_machine.node, "grid-row"):
            LOGGER.warning(
                f"The parent of a 'grid-item' should be a 'grid-row' [{WARNING_TYPE}.grid]",
                location=(self.env.docname, self.lineno),
                type=WARNING_TYPE,
                subtype="grid",
            )
        column = create_component(
            "grid-item",
            [
                "sd-col",
                f"sd-d-flex-{self.options.get('child-direction', 'column')}",
            ]
            + self.options.get("columns", [])
            + self.options.get("margin", [])
            + self.options.get("padding", [])
            + (
                [f'sd-align-major-{self.options["child-align"]}']
                if "child-align" in self.options
                else []
            )
            + (["sd-border-1"] if "outline" in self.options else [])
            + self.options.get("class", []),
        )
        self.set_source_info(column)
        self.state.nested_parse(self.content, self.content_offset, column)
        return [column]


class GridItemCardDirective(SdDirective):
    """An item within a grid row, with an internal card."""

    has_content = True
    required_arguments = 0
    optional_arguments = 1  # card title
    final_argument_whitespace = True
    option_spec = {
        "columns": item_columns_option,
        "margin": margin_option,
        "padding": padding_option,
        "class-item": directives.class_option,
        # The options below must be sync'ed with CardDirective.option_spec (minus margin).
        "width": make_choice(["auto", "25%", "50%", "75%", "100%"]),
        "text-align": text_align,
        "img-background": directives.uri,
        "img-top": directives.uri,
        "img-bottom": directives.uri,
        "img-alt": directives.unchanged,
        "link": directives.uri,
        "link-type": make_choice(["url", "any", "ref", "doc"]),
        "link-alt": directives.unchanged,
        "shadow": make_choice(["none", "sm", "md", "lg"]),
        "class-card": directives.class_option,
        "class-body": directives.class_option,
        "class-title": directives.class_option,
        "class-header": directives.class_option,
        "class-footer": directives.class_option,
        "class-img-top": directives.class_option,
        "class-img-bottom": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        if not is_component(self.state_machine.node, "grid-row"):
            LOGGER.warning(
                f"The parent of a 'grid-item' should be a 'grid-row' [{WARNING_TYPE}.grid]",
                location=(self.env.docname, self.lineno),
                type=WARNING_TYPE,
                subtype="grid",
            )
        column = create_component(
            "grid-item",
            [
                "sd-col",
                "sd-d-flex-row",
                *self.options.get("columns", []),
                *self.options.get("margin", []),
                *self.options.get("padding", []),
                *self.options.get("class-item", []),
            ],
        )
        card_options = {
            key: value
            for key, value in self.options.items()
            if key
            in [
                "width",
                "text-align",
                "img-background",
                "img-top",
                "img-bottom",
                "img-alt",
                "link",
                "link-type",
                "link-alt",
                "shadow",
                "class-card",
                "class-body",
                "class-title",
                "class-header",
                "class-footer",
                "class-img-top",
                "class-img-bottom",
            ]
        }
        if "width" not in card_options:
            card_options["width"] = "100%"
        card_options["margin"] = []
        card = CardDirective.create_card(self, self.arguments, card_options)
        column += card
        return [column]
