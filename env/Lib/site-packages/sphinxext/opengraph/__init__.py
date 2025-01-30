from typing import Any, Dict
from urllib.parse import urljoin, urlparse, urlunparse
from pathlib import Path

import docutils.nodes as nodes
from sphinx.application import Sphinx

from .descriptionparser import get_description
from .metaparser import get_meta_description
from .titleparser import get_title

try:
    import matplotlib
except ImportError:
    print("matplotlib is not installed, social cards will not be generated")
    create_social_card = None
    DEFAULT_SOCIAL_CONFIG = {}
else:
    from .socialcards import create_social_card, DEFAULT_SOCIAL_CONFIG

import os


DEFAULT_DESCRIPTION_LENGTH = 200
DEFAULT_DESCRIPTION_LENGTH_SOCIAL_CARDS = 160
DEFAULT_PAGE_LENGTH_SOCIAL_CARDS = 80

# A selection from https://www.iana.org/assignments/media-types/media-types.xhtml#image
IMAGE_MIME_TYPES = {
    "gif": "image/gif",
    "apng": "image/apng",
    "webp": "image/webp",
    "jpeg": "image/jpeg",
    "jpg": "image/jpeg",
    "png": "image/png",
    "bmp": "image/bmp",
    "heic": "image/heic",
    "heif": "image/heif",
    "tiff": "image/tiff",
}


def make_tag(property: str, content: str, type_: str = "property") -> str:
    # Parse quotation, so they won't break html tags if smart quotes are disabled
    content = content.replace('"', "&quot;")
    return f'<meta {type_}="{property}" content="{content}" />'


def get_tags(
    app: Sphinx,
    context: Dict[str, Any],
    doctree: nodes.document,
    config: Dict[str, Any],
) -> str:
    # Get field lists for per-page overrides
    fields = context["meta"]
    if fields is None:
        fields = {}

    if "ogp_disable" in fields:
        return ""

    tags = {}
    meta_tags = {}  # For non-og meta tags

    # Set length of description
    try:
        desc_len = int(
            fields.get("ogp_description_length", config["ogp_description_length"])
        )
    except ValueError:
        desc_len = DEFAULT_DESCRIPTION_LENGTH

    # Get the title and parse any html in it
    title = get_title(context["title"], skip_html_tags=False)
    title_excluding_html = get_title(context["title"], skip_html_tags=True)

    # Parse/walk doctree for metadata (tag/description)
    description = get_description(doctree, desc_len, [title, title_excluding_html])

    # title tag
    tags["og:title"] = title

    # type tag
    tags["og:type"] = config["ogp_type"]

    if os.getenv("READTHEDOCS") and not config["ogp_site_url"]:
        # readthedocs uses html_baseurl for sphinx > 1.8
        parse_result = urlparse(config["html_baseurl"])

        if config["html_baseurl"] is None:
            raise OSError("ReadTheDocs did not provide a valid canonical URL!")

        # Grab root url from canonical url
        config["ogp_site_url"] = urlunparse(
            (
                parse_result.scheme,
                parse_result.netloc,
                parse_result.path,
                "",
                "",
                "",
            )
        )

    # url tag
    # Get the URL of the specific page
    page_url = urljoin(
        config["ogp_site_url"], app.builder.get_target_uri(context["pagename"])
    )
    tags["og:url"] = page_url

    # site name tag, False disables, default to project if ogp_site_name not
    # set.
    if config["ogp_site_name"] is False:
        site_name = None
    elif config["ogp_site_name"] is None:
        site_name = config["project"]
    else:
        site_name = config["ogp_site_name"]
    if site_name:
        tags["og:site_name"] = site_name

    # description tag
    if description:
        tags["og:description"] = description

        if config["ogp_enable_meta_description"] and not get_meta_description(
            context["metatags"]
        ):
            meta_tags["description"] = description

    # image tag
    # Get basic values from config
    if "og:image" in fields:
        image_url = fields["og:image"]
        ogp_use_first_image = False
        ogp_image_alt = fields.get("og:image:alt")
        fields.pop("og:image", None)
    else:
        image_url = config["ogp_image"]
        ogp_use_first_image = config["ogp_use_first_image"]
        ogp_image_alt = fields.get("og:image:alt", config["ogp_image_alt"])

    # Decide whether to add social media card images for each page.
    # Only do this as a fallback if the user hasn't given any configuration
    # to add other images.
    config_social = DEFAULT_SOCIAL_CONFIG.copy()
    social_card_user_options = app.config.ogp_social_cards or {}
    config_social.update(social_card_user_options)
    if (
        not (image_url or ogp_use_first_image)
        and config_social.get("enable") is not False
        and create_social_card is not None
    ):
        # Description
        description_max_length = config_social.get(
            "description_max_length", DEFAULT_DESCRIPTION_LENGTH_SOCIAL_CARDS - 3
        )
        if len(description) > description_max_length:
            description = description[:description_max_length].strip() + "..."

        # Page title
        pagetitle = title
        if len(pagetitle) > DEFAULT_PAGE_LENGTH_SOCIAL_CARDS:
            pagetitle = pagetitle[:DEFAULT_PAGE_LENGTH_SOCIAL_CARDS] + "..."

        # Site URL
        site_url = config_social.get("site_url", True)
        if site_url is True:
            url_text = app.config.ogp_site_url.split("://")[-1]
        elif isinstance(site_url, str):
            url_text = site_url

        # Plot an image with the given metadata to the output path
        image_path = create_social_card(
            app,
            config_social,
            site_name,
            pagetitle,
            description,
            url_text,
            context["pagename"],
        )
        ogp_use_first_image = False

        # Alt text is taken from description unless given
        if "og:image:alt" in fields:
            ogp_image_alt = fields.get("og:image:alt")
        else:
            ogp_image_alt = description

        # Link the image in our page metadata
        # We use os.path.sep to standardize behavior acros *nix and Windows
        url = app.config.ogp_site_url.strip("/")
        image_path = str(image_path).replace(os.path.sep, "/").strip("/")
        image_url = f"{url}/{image_path}"

        # If the social card objects have been added we add special metadata for them
        # These are the dimensions *in pixels* of the card
        # They were chosen by looking at the image pixel dimensions on disk
        tags["og:image:width"] = "1146"
        tags["og:image:height"] = "600"
        meta_tags["twitter:card"] = "summary_large_image"

    fields.pop("og:image:alt", None)

    first_image = None
    if ogp_use_first_image:
        # Use the first image that is defined in the current page
        first_image = doctree.next_node(nodes.image)
        if (
            first_image
            and Path(first_image.get("uri", "")).suffix[1:].lower() in IMAGE_MIME_TYPES
        ):
            image_url = first_image["uri"]
            ogp_image_alt = first_image.get("alt", None)
        else:
            first_image = None

    if image_url:
        # temporarily disable relative image paths with field lists
        if "og:image" not in fields:
            image_url_parsed = urlparse(image_url)
            if not image_url_parsed.scheme:
                # Relative image path detected, relative to the source. Make absolute.
                if first_image:
                    root = page_url
                else:  # ogp_image is set
                    # ogp_image is defined as being relative to the site root.
                    # This workaround is to keep that functionality from breaking.
                    root = config["ogp_site_url"]

                image_url = urljoin(root, image_url_parsed.path)
            tags["og:image"] = image_url

        # Add image alt text (either provided by config or from site_name)
        if isinstance(ogp_image_alt, str):
            tags["og:image:alt"] = ogp_image_alt
        elif ogp_image_alt is None and site_name:
            tags["og:image:alt"] = site_name
        elif ogp_image_alt is None and title:
            tags["og:image:alt"] = title

    # arbitrary tags and overrides
    tags.update({k: v for k, v in fields.items() if k.startswith("og:")})

    return (
        "\n".join(
            [make_tag(p, c) for p, c in tags.items()]
            + [make_tag(p, c, "name") for p, c in meta_tags.items()]
            + config["ogp_custom_meta_tags"]
        )
        + "\n"
    )


def html_page_context(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: Dict[str, Any],
    doctree: nodes.document,
) -> None:
    if doctree:
        context["metatags"] += get_tags(app, context, doctree, app.config)


def setup(app: Sphinx) -> Dict[str, Any]:
    # ogp_site_url="" allows relative by default, even though it's not
    # officially supported by OGP.
    app.add_config_value("ogp_site_url", "", "html")
    app.add_config_value("ogp_description_length", DEFAULT_DESCRIPTION_LENGTH, "html")
    app.add_config_value("ogp_image", None, "html")
    app.add_config_value("ogp_image_alt", None, "html")
    app.add_config_value("ogp_use_first_image", False, "html")
    app.add_config_value("ogp_type", "website", "html")
    app.add_config_value("ogp_site_name", None, "html")
    app.add_config_value("ogp_social_cards", None, "html")
    app.add_config_value("ogp_custom_meta_tags", [], "html")
    app.add_config_value("ogp_enable_meta_description", True, "html")

    # Main Sphinx OpenGraph linking
    app.connect("html-page-context", html_page_context)

    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
