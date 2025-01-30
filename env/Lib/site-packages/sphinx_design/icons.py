from collections.abc import Sequence
from functools import lru_cache
import json
import re
from typing import Any, Optional

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.application import Sphinx
from sphinx.util import logging
from sphinx.util.docutils import SphinxRole

from . import compiled
from ._compat import read_text
from .shared import WARNING_TYPE, SdDirective

logger = logging.getLogger(__name__)

OCTICON_VERSION = "v19.8.0"

OCTICON_CSS = """\
.octicon {
  display: inline-block;
  vertical-align: text-top;
  fill: currentColor;
}"""


def setup_icons(app: Sphinx) -> None:
    app.add_role("octicon", OcticonRole())
    app.add_directive("_all-octicon", AllOcticons)
    for style in ["fa", "fas", "fab", "far"]:
        # note: fa is deprecated in v5, fas is the default and fab is the other free option
        app.add_role(style, FontawesomeRole(style))
    for style in ["regular", "outlined", "round", "sharp", "twotone"]:
        app.add_role("material-" + style, MaterialRole(style))
    app.add_config_value("sd_fontawesome_latex", False, "env")
    app.connect("config-inited", add_fontawesome_pkg)
    app.add_node(
        fontawesome,
        html=(visit_fontawesome_html, depart_fontawesome_html),
        latex=(visit_fontawesome_latex, None),
        man=(visit_fontawesome_warning, None),
        text=(visit_fontawesome_warning, None),
        texinfo=(visit_fontawesome_warning, None),
    )


@lru_cache(1)
def get_octicon_data() -> dict[str, Any]:
    """Load all octicon data."""
    content = read_text(compiled, "octicons.json")
    return json.loads(content)


def list_octicons() -> list[str]:
    """List available octicon names."""
    return list(get_octicon_data().keys())


HEIGHT_REGEX = re.compile(r"^(?P<value>\d+(\.\d+)?)(?P<unit>px|em|rem)$")


def get_octicon(
    name: str,
    height: str = "1em",
    classes: Sequence[str] = (),
    aria_label: Optional[str] = None,
) -> str:
    """Return the HTML for an GitHub octicon SVG icon.

    :height: the height of the octicon, with suffix unit 'px', 'em' or 'rem'.
    """
    try:
        data = get_octicon_data()[name]
    except KeyError as exc:
        raise KeyError(f"Unrecognised octicon: {name}") from exc

    match = HEIGHT_REGEX.match(height)
    if not match:
        raise ValueError(
            f"Invalid height: '{height}', must be format <integer><px|em|rem>"
        )
    height_value = round(float(match.group("value")), 3)
    height_unit = match.group("unit")

    original_height = 16
    if "16" not in data["heights"]:
        original_height = int(next(iter(data["heights"].keys())))
    elif "24" in data["heights"]:
        if height_unit == "px":
            if height_value >= 24:
                original_height = 24
        elif height_value >= 1.5:
            original_height = 24
    original_width = data["heights"][str(original_height)]["width"]
    width_value = round(original_width * height_value / original_height, 3)
    content = data["heights"][str(original_height)]["path"]
    options = {
        "version": "1.1",
        "width": f"{width_value}{height_unit}",
        "height": f"{height_value}{height_unit}",
        "class": " ".join(("sd-octicon", f"sd-octicon-{name}", *classes)),
    }

    options["viewBox"] = f"0 0 {original_width} {original_height}"

    if aria_label is not None:
        options["aria-label"] = aria_label
        options["role"] = "img"
    else:
        options["aria-hidden"] = "true"

    opt_string = " ".join(f'{k}="{v}"' for k, v in options.items())
    return f"<svg {opt_string}>{content}</svg>"


class OcticonRole(SphinxRole):
    """Role to display a GitHub octicon SVG.

    Additional classes can be added to the element after a semicolon.
    """

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        values = self.text.split(";") if ";" in self.text else [self.text]
        icon = values[0]
        height = "1em" if len(values) < 2 else values[1]
        classes = "" if len(values) < 3 else values[2]
        icon = icon.strip()
        try:
            svg = get_octicon(icon, height=height, classes=classes.split())
        except Exception as exc:
            msg = self.inliner.reporter.error(
                f"Invalid octicon content: {exc}",
                line=self.lineno,
            )
            prb = self.inliner.problematic(self.rawtext, self.rawtext, msg)
            return [prb], [msg]
        node = nodes.raw("", nodes.Text(svg), format="html")
        self.set_source_info(node)
        return [node], []


class AllOcticons(SdDirective):
    """Directive to generate all octicon icons.

    Primarily for self documentation.
    """

    option_spec = {
        "class": directives.class_option,
    }

    def run_with_defaults(self) -> list[nodes.Node]:
        classes = self.options.get("class", [])
        table = nodes.table()
        group = nodes.tgroup(cols=2)
        table += group
        group.extend(
            (
                nodes.colspec(colwidth=1),
                nodes.colspec(colwidth=1),
            )
        )
        body = nodes.tbody()
        group += body
        for icon in list_octicons():
            row = nodes.row()
            body += row
            cell = nodes.entry()
            row += cell
            cell += nodes.literal(icon, icon)
            cell = nodes.entry()
            row += cell
            cell += nodes.raw(
                "",
                get_octicon(icon, classes=classes),
                format="html",
            )
        return [table]


class fontawesome(nodes.Element, nodes.General):  # noqa: N801
    """Node for rendering fontawesome icon."""


class FontawesomeRole(SphinxRole):
    """Role to display a Fontawesome icon.

    Additional classes can be added to the element after a semicolon.
    """

    def __init__(self, style: str) -> None:
        super().__init__()
        self.style = style

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        icon, classes = self.text.split(";", 1) if ";" in self.text else [self.text, ""]
        icon = icon.strip()
        node = fontawesome(
            icon=icon, classes=[self.style, f"fa-{icon}", *classes.split()]
        )
        self.set_source_info(node)
        return [node], []


def visit_fontawesome_html(self, node):
    self.body.append(self.starttag(node, "span", ""))


def depart_fontawesome_html(self, node):
    self.body.append("</span>")


def add_fontawesome_pkg(app, config):
    if app.config.sd_fontawesome_latex:
        app.add_latex_package("fontawesome")


def visit_fontawesome_latex(self, node):
    """Add latex fonteawesome icon, if configured, else warn."""
    if self.config.sd_fontawesome_latex:
        self.body.append(f"\\faicon{{{node['icon']}}}")
    else:
        logger.warning(
            "Fontawesome icons not included in LaTeX output, "
            f"consider 'sd_fontawesome_latex=True' [{WARNING_TYPE}.fa-build]",
            location=node,
            type=WARNING_TYPE,
            subtype="fa-build",
        )
    raise nodes.SkipNode


def visit_fontawesome_warning(self, node: nodes.Element) -> None:
    """Warn that fontawesome is not supported for this builder."""
    logger.warning(
        "Fontawesome icons not supported for builder: "
        f"{self.builder.name} [{WARNING_TYPE}.fa-build]",
        location=node,
        type=WARNING_TYPE,
        subtype="fa-build",
    )
    raise nodes.SkipNode


@lru_cache(1)
def get_material_icon_data(style: str) -> dict[str, Any]:
    """Load all octicon data."""
    content = read_text(compiled, f"material_{style}.json")
    return json.loads(content)


def get_material_icon(
    style: str,
    name: str,
    height: str = "1em",
    classes: Sequence[str] = (),
    aria_label: Optional[str] = None,
) -> str:
    """Return the HTML for an Google material icon SVG icon.

    :height: the height of the material icon, with suffix unit 'px', 'em' or 'rem'.
    """
    try:
        data = get_material_icon_data(style)[name]
    except KeyError as exc:
        raise KeyError(f"Unrecognised material-{style} icon: {name}") from exc

    match = HEIGHT_REGEX.match(height)
    if not match:
        raise ValueError(
            f"Invalid height: '{height}', must be format <integer><px|em|rem>"
        )
    height_value = round(float(match.group("value")), 3)
    height_unit = match.group("unit")

    original_height = 20
    if "20" not in data["heights"]:
        original_height = int(next(iter(data["heights"].keys())))
    elif "24" in data["heights"]:
        if height_unit == "px":
            if height_value >= 24:
                original_height = 24
        elif height_value >= 1.5:
            original_height = 24
    original_width = data["heights"][str(original_height)]["width"]
    width_value = round(original_width * height_value / original_height, 3)
    content = data["heights"][str(original_height)]["path"]
    options = {
        "version": "4.0.0.63c5cb3",
        "width": f"{width_value}{height_unit}",
        "height": f"{height_value}{height_unit}",
        "class": " ".join(("sd-material-icon", f"sd-material-icon-{name}", *classes)),
    }

    options["viewBox"] = f"0 0 {original_width} {original_height}"

    if aria_label is not None:
        options["aria-label"] = aria_label
        options["role"] = "img"
    else:
        options["aria-hidden"] = "true"

    opt_string = " ".join(f'{k}="{v}"' for k, v in options.items())
    return f"<svg {opt_string}>{content}</svg>"


class MaterialRole(SphinxRole):
    """Role to display a Material-* icon.

    Additional classes can be added to the element after a semicolon.
    """

    def __init__(self, style: str) -> None:
        super().__init__()
        self.style = style

    def run(self) -> tuple[list[nodes.Node], list[nodes.system_message]]:
        """Run the role."""
        values = self.text.split(";") if ";" in self.text else [self.text]
        icon = values[0]
        height = "1em" if len(values) < 2 else values[1]
        classes = "" if len(values) < 3 else values[2]
        icon = icon.strip()
        try:
            svg = get_material_icon(
                self.style, icon, height=height, classes=classes.split()
            )
        except Exception as exc:
            msg = self.inliner.reporter.error(
                f"Invalid material-{self.style} icon content: {type(exc)} {exc}",
                line=self.lineno,
            )
            prb = self.inliner.problematic(self.rawtext, self.rawtext, msg)
            return [prb], [msg]
        node = nodes.raw("", nodes.Text(svg), format="html")
        self.set_source_info(node)
        return [node], []
