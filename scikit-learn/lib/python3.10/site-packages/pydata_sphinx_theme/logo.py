"""customize events for logo management.

we use one event to copy over custom logo images to _static
and another even to link them in the html context
"""

from functools import partial
from pathlib import Path

from docutils.nodes import Node
from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx.util import isurl
from sphinx.util.fileutil import copy_asset_file

from .utils import get_theme_options_dict, maybe_warn


def setup_logo_path(
    app: Sphinx, pagename: str, templatename: str, context: dict, doctree: Node
) -> None:
    """Set up relative paths to logos in our HTML templates.

    In Sphinx, the context["logo"] is a path to the `html_logo` image now in the output
    `_static` folder.

    If logo["image_light"] and logo["image_dark"] are given, we must modify them to
    follow the same pattern. They have already been copied to the output folder
    in the `update_config` event.
    """
    # get information from the context "logo_url" for sphinx>=6, "logo" sphinx<6
    pathto = context.get("pathto")
    logo = context.get("logo_url") or context.get("logo")
    theme_logo = context.get("theme_logo", {})

    # Define the final path to logo images in the HTML context
    theme_logo["image_relative"] = {}
    for kind in ["light", "dark"]:
        image_kind_logo = theme_logo.get(f"image_{kind}")

        # If it's a URL the "relative" path is just the URL
        # else we need to calculate the relative path to a local file
        if image_kind_logo:
            if not isurl(image_kind_logo):
                image_kind_name = Path(image_kind_logo).name
                image_kind_logo = pathto(f"_static/{image_kind_name}", resource=True)
            theme_logo["image_relative"][kind] = image_kind_logo

        # If there's no custom logo for this kind, just use `html_logo`
        # If `logo` is also None, then do not add this key to context.
        elif isinstance(logo, str) and len(logo) > 0:
            theme_logo["image_relative"][kind] = logo

    # Update our context logo variables with the new image paths
    context["theme_logo"] = theme_logo


def copy_logo_images(app: Sphinx, exception=None) -> None:
    """
    Copy logo image to the _static directory.

    If logo image paths are given, copy them to the `_static` folder.
    Then we can link to them directly in an html_page_context event.
    """
    warning = partial(maybe_warn, app)
    logo = get_theme_options_dict(app).get("logo", {})
    staticdir = Path(app.builder.outdir) / "_static"
    assert staticdir.is_absolute()
    for kind in ["light", "dark"]:
        path_image = logo.get(f"image_{kind}")
        if not path_image or isurl(path_image):
            continue
        if (staticdir / Path(path_image).name).exists():
            # file already exists in static dir e.g. because a theme has
            # bundled the logo and installed it there
            continue
        full_logo_path = Path(app.srcdir) / path_image
        assert full_logo_path.is_absolute()
        if not full_logo_path.exists():
            warning(f"Path to {kind} image logo does not exist: {path_image}")
        # Ensure templates cannot be passed for logo path to avoid security
        # vulnerability
        if path_image.lower().endswith("_t"):
            raise ExtensionError(
                f"The {kind} logo path '{path_image}' looks like a Sphinx template; "
                "please provide a static logo image."
            )
        copy_asset_file(str(full_logo_path), staticdir)
