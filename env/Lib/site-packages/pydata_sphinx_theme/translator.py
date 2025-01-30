"""A custom Sphinx HTML Translator for Bootstrap layout."""

import types

from sphinx.application import Sphinx
from sphinx.ext.autosummary import autosummary_table
from sphinx.util import logging


logger = logging.getLogger(__name__)


class BootstrapHTML5TranslatorMixin:
    """Mixin HTML Translator for a Bootstrap-ified Sphinx layout.

    Only a couple of functions have been overridden to produce valid HTML to be
    directly styled with Bootstrap, and fulfill acessibility best practices.
    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.settings.table_style = "table"

    def starttag(self, *args, **kwargs):
        """Perform small modifications to tags.

        - ensure aria-level is set for any tag with heading role
        """
        if kwargs.get("ROLE") == "heading" and "ARIA-LEVEL" not in kwargs:
            kwargs["ARIA-LEVEL"] = "2"

        return super().starttag(*args, **kwargs)

    def visit_table(self, node):
        """
        Custom visit table method.

        Copy of sphinx source to *not* add 'docutils' and 'align-default' classes but
        add 'table' class.
        """
        # init the attributes
        atts = {}

        self._table_row_indices.append(0)

        # get the classes
        classes = [cls.strip(" \t\n") for cls in self.settings.table_style.split(",")]

        # we're looking at the 'real_table', which is wrapped by an autosummary
        if isinstance(node.parent, autosummary_table):
            classes += ["autosummary"]

        # add the width if set in a style attribute
        if "width" in node:
            atts["style"] = f'width: {node["width"]}'

        # add specific class if align is set
        if "align" in node:
            classes.append(f'table-{node["align"]}')

        # put table within a scrollable container (for tables that are too wide)
        self.body.append('<div class="pst-scrollable-table-container">')

        tag = self.starttag(node, "table", CLASS=" ".join(classes), **atts)
        self.body.append(tag)

    def depart_table(self, node):
        """
        Custom depart_table method to close the scrollable div we add in
        visit_table.
        """
        super().depart_table(node)
        self.body.append("</div>\n")


def setup_translators(app: Sphinx):
    """
    Add bootstrap HTML functionality if we are using an HTML translator.

    This re-uses the pre-existing Sphinx translator and adds extra functionality
    defined in ``BootstrapHTML5TranslatorMixin``.
    This way we can retain the original translator's
    behavior and configuration, and _only_ add the extra bootstrap rules.
    If we don't detect an HTML-based translator, then we do nothing.
    """
    if not app.registry.translators.items():
        try:
            default_translator_class = app.builder.default_translator_class
        except AttributeError:
            # some builders, e.g. linkcheck, do not define 'default_translator_class'
            return

        translator = types.new_class(
            "BootstrapHTML5Translator",
            (
                BootstrapHTML5TranslatorMixin,
                default_translator_class,
            ),
            {},
        )
        app.set_translator(app.builder.name, translator, override=True)
    else:
        for name, klass in app.registry.translators.items():
            if app.builder.format != "html":
                # Skip translators that are not HTML
                continue

            translator = types.new_class(
                "BootstrapHTML5Translator",
                (
                    BootstrapHTML5TranslatorMixin,
                    klass,
                ),
                {},
            )
            app.set_translator(name, translator, override=True)
