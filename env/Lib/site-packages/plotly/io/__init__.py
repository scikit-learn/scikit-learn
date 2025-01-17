from _plotly_utils.importers import relative_import
import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._kaleido import to_image, write_image, full_figure_for_development
    from . import orca, kaleido
    from . import json
    from ._json import to_json, from_json, read_json, write_json
    from ._templates import templates, to_templated
    from ._html import to_html, write_html
    from ._renderers import renderers, show
    from . import base_renderers

    __all__ = [
        "to_image",
        "write_image",
        "orca",
        "json",
        "to_json",
        "from_json",
        "read_json",
        "write_json",
        "templates",
        "to_templated",
        "to_html",
        "write_html",
        "renderers",
        "show",
        "base_renderers",
        "full_figure_for_development",
    ]
else:
    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".orca", ".kaleido", ".json", ".base_renderers"],
        [
            "._kaleido.to_image",
            "._kaleido.write_image",
            "._kaleido.full_figure_for_development",
            "._json.to_json",
            "._json.from_json",
            "._json.read_json",
            "._json.write_json",
            "._templates.templates",
            "._templates.to_templated",
            "._html.to_html",
            "._html.write_html",
            "._renderers.renderers",
            "._renderers.show",
        ],
    )

    # Set default template (for < 3.7 this is done in ploty/__init__.py)
    from plotly.io import templates

    templates._default = "plotly"
