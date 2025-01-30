"""General helpers for the management of config parameters."""

import copy
import os
import re

from typing import Any, Callable, Dict, Iterable, List, Optional, Union

from docutils.nodes import Node
from sphinx.application import Sphinx
from sphinx.util import logging, matching


def get_theme_options_dict(app: Sphinx) -> Dict[str, Any]:
    """Return theme options for the application w/ a fallback if they don't exist.

    The "top-level" mapping (the one we should usually check first, and modify
    if desired) is ``app.builder.theme_options``. It is created by Sphinx as a
    copy of ``app.config.html_theme_options`` (containing user-configs from
    their ``conf.py``); sometimes that copy never occurs though which is why we
    check both.
    """
    if hasattr(app.builder, "theme_options"):
        return app.builder.theme_options
    elif hasattr(app.config, "html_theme_options"):
        return app.config.html_theme_options
    else:
        return {}


def config_provided_by_user(app: Sphinx, key: str) -> bool:
    """Check if the user has manually provided the config."""
    return any(key in ii for ii in [app.config.overrides, app.config._raw_config])


def traverse_or_findall(
    node: Node, condition: Union[Callable, type], **kwargs
) -> Iterable[Node]:
    """Triage node.traverse (docutils <0.18.1) vs node.findall.

    TODO: This check can be removed when the minimum supported docutils version
    for numpydoc is docutils>=0.18.1.
    """
    return (
        node.findall(condition, **kwargs)
        if hasattr(node, "findall")
        else node.traverse(condition, **kwargs)
    )


def escape_ansi(string: str) -> str:
    """Helper function to remove ansi coloring from sphinx warnings."""
    ansi_escape = re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]")
    return ansi_escape.sub("", string)


SPHINX_LOGGER = logging.getLogger(__name__)


def maybe_warn(app: Sphinx, msg, *args, **kwargs):
    """Wraps the Sphinx logger to allow warning suppression."""
    theme_options = get_theme_options_dict(app)
    should_warn = theme_options.get("surface_warnings", False)
    if should_warn:
        SPHINX_LOGGER.warning(msg, *args, **kwargs)
    else:
        SPHINX_LOGGER.info(msg, *args, **kwargs)


def set_secondary_sidebar_items(
    app: Sphinx, pagename: str, templatename: str, context, doctree
) -> None:
    """Set the secondary sidebar items to render for the given pagename."""
    if "theme_secondary_sidebar_items" in context:
        templates = context["theme_secondary_sidebar_items"]
        if isinstance(templates, dict):
            templates = _get_matching_sidebar_items(pagename, templates)

        context["secondary_sidebar_items"] = _update_and_remove_templates(
            app,
            context,
            templates,
            "theme_secondary_sidebar_items",
        )


def _update_and_remove_templates(
    app: Sphinx,
    context: Dict[str, Any],
    templates: Union[List, str],
    section: str,
    templates_skip_empty_check: Optional[List[str]] = None,
) -> List[str]:
    """
    Update templates to include html suffix if needed; remove templates which render
    empty.

    Args:
        app: Sphinx application passed to the html page context
        context: The html page context; dictionary of values passed to the templating
            engine
        templates: A list of template names, or a string of comma separated template
            names
        section: Name of the template section where the templates are to be rendered.
            Valid section names include any of the ``sphinx`` or ``html_theme_options``
            that take templates or lists of templates as arguments, for example:
            ``theme_navbar_start``, ``theme_primary_sidebar_end``,
            ``theme_secondary_sidebar_items``, ``sidebars``, etc.
            For a complete list of valid section names, see the source for
            :py:func:`pydata_sphinx_theme.update_and_remove_templates` and
            :py:func:`pydata_sphinx_theme.utils.set_secondary_sidebar_items`,
            both of which call this function.
        templates_skip_empty_check: Names of any templates which should never be removed
            from the list of filtered templates returned by this function.
            These templates aren't checked if they render empty, which can save time if
            the template is slow to render.

    Returns:
        A list of template names (including '.html' suffix) to render into the section
    """
    if templates_skip_empty_check is None:
        templates_skip_empty_check = []

    # Break apart `,` separated strings so we can use , in the defaults
    if isinstance(templates, str):
        templates = [template.strip() for template in templates.split(",")]

    # Add `.html` to templates with no suffix
    suffixed_templates = []
    for template in templates:
        if os.path.splitext(template)[1]:
            suffixed_templates.append(template)
        else:
            suffixed_templates.append(f"{template}.html")

    ctx = copy.copy(context)
    ctx.update({section: suffixed_templates})

    # Check whether the template renders to an empty string; remove if this is the case
    # Skip templates that are slow to render with templates_skip_empty_check
    filtered_templates = []
    for template in suffixed_templates:
        if any(template.endswith(item) for item in templates_skip_empty_check):
            filtered_templates.append(template)
        else:
            rendered = app.builder.templates.render(template, ctx)
            if len(rendered.strip()) != 0:
                filtered_templates.append(template)

    return filtered_templates


def _get_matching_sidebar_items(
    pagename: str, sidebars: Dict[str, List[str]]
) -> List[str]:
    """Get the matching sidebar templates to render for the given pagename.

    If a page matches more than one pattern, a warning is emitted, and the templates
    for the last matching pattern are used.

    This function was adapted from
    sphinx.builders.html.StandaloneHTMLBuilder.add_sidebars.
    """
    matched = None
    secondary_sidebar_items = []
    for pattern, sidebar_items in sidebars.items():
        if matching.patmatch(pagename, pattern):
            if matched and _has_wildcard(pattern) and _has_wildcard(matched):
                (
                    SPHINX_LOGGER.warning(
                        "Page %s matches two wildcard patterns in secondary_sidebar_items: %s and %s",  # noqa: E501
                        pagename,
                        matched,
                        pattern,
                    ),
                )
            matched = pattern
            secondary_sidebar_items = sidebar_items
    return secondary_sidebar_items


def _has_wildcard(pattern: str) -> bool:
    """Check whether the pattern contains a wildcard.

    Taken from sphinx.builders.StandaloneHTMLBuilder.add_sidebars.
    """
    return any(char in pattern for char in "*?[")
