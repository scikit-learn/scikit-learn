from __future__ import annotations

import os
import posixpath
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import urljoin, urlparse, urlsplit, urlunsplit

from docutils import nodes

from sphinxext.opengraph._description_parser import get_description
from sphinxext.opengraph._meta_parser import get_meta_description
from sphinxext.opengraph._title_parser import get_title

try:
    from types import NoneType
except ImportError:
    NoneType = type(None)

if TYPE_CHECKING:
    from typing import Any

    from sphinx.application import Sphinx
    from sphinx.builders import Builder
    from sphinx.config import Config
    from sphinx.environment import BuildEnvironment
    from sphinx.util.typing import ExtensionMetadata

try:
    from sphinxext.opengraph._social_cards import (
        DEFAULT_SOCIAL_CONFIG,
        create_social_card,
    )
except ImportError:
    print('matplotlib is not installed, social cards will not be generated')
    create_social_card = None
    DEFAULT_SOCIAL_CONFIG = {}

__version__ = '0.13.0'
version_info = (0, 13, 0)

DEFAULT_DESCRIPTION_LENGTH = 200
DEFAULT_DESCRIPTION_LENGTH_SOCIAL_CARDS = 160
DEFAULT_PAGE_LENGTH_SOCIAL_CARDS = 80

# A selection from https://www.iana.org/assignments/media-types/media-types.xhtml#image
IMAGE_MIME_TYPES = {
    'gif': 'image/gif',
    'apng': 'image/apng',
    'webp': 'image/webp',
    'jpeg': 'image/jpeg',
    'jpg': 'image/jpeg',
    'png': 'image/png',
    'bmp': 'image/bmp',
    'heic': 'image/heic',
    'heif': 'image/heif',
    'tiff': 'image/tiff',
}


def html_page_context(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: dict[str, Any],
    doctree: nodes.document,
) -> None:
    if app.builder.name == 'epub':
        return

    if doctree:
        context['metatags'] += get_tags(
            context,
            doctree,
            srcdir=app.srcdir,
            outdir=app.outdir,
            config=app.config,
            builder=app.builder,
            env=app.env,
        )


def get_tags(
    context: dict[str, Any],
    doctree: nodes.document,
    *,
    srcdir: str | Path,
    outdir: str | Path,
    config: Config,
    builder: Builder,
    env: BuildEnvironment,
) -> str:
    # Get field lists for per-page overrides
    fields = context['meta']
    if fields is None:
        fields = {}

    if 'ogp_disable' in fields:
        return ''

    tags = {}
    meta_tags = {}  # For non-og meta tags

    # Set length of description
    try:
        desc_len = int(
            fields.get('ogp_description_length', config.ogp_description_length)
        )
    except ValueError:
        desc_len = DEFAULT_DESCRIPTION_LENGTH

    # Get the title and parse any html in it
    title, title_excluding_html = get_title(context['title'])

    # Parse/walk doctree for metadata (tag/description)
    description = get_description(doctree, desc_len, {title, title_excluding_html})

    # title tag
    tags['og:title'] = title

    # type tag
    tags['og:type'] = config.ogp_type

    if not config.ogp_site_url and os.getenv('READTHEDOCS'):
        ogp_site_url = ambient_site_url()
    else:
        ogp_site_url = config.ogp_site_url

    # If ogp_canonical_url is not set, default to the value of ogp_site_url
    ogp_canonical_url = config.ogp_canonical_url or ogp_site_url

    # url tag
    # Get the URL of the specific page
    page_url = urljoin(ogp_canonical_url, builder.get_target_uri(context['pagename']))
    tags['og:url'] = page_url

    # site name tag, False disables, default to project if ogp_site_name not
    # set.
    if config.ogp_site_name is False:
        site_name = None
    elif config.ogp_site_name is None:
        site_name = config.project
    else:
        site_name = config.ogp_site_name
    if site_name:
        tags['og:site_name'] = site_name

    # description tag
    if description:
        tags['og:description'] = description

        if config.ogp_enable_meta_description and not get_meta_description(
            context['metatags']
        ):
            meta_tags['description'] = description

    # image tag
    # Get basic values from config
    if 'og:image' in fields:
        image_url = fields['og:image']
        ogp_use_first_image = False
        ogp_image_alt = fields.get('og:image:alt')
        fields.pop('og:image', None)
    else:
        image_url = config.ogp_image
        ogp_use_first_image = config.ogp_use_first_image
        ogp_image_alt = fields.get('og:image:alt', config.ogp_image_alt)

    # Decide whether to add social media card images for each page.
    # Only do this as a fallback if the user hasn't given any configuration
    # to add other images.
    config_social = DEFAULT_SOCIAL_CONFIG.copy()
    social_card_user_options = config.ogp_social_cards or {}
    config_social.update(social_card_user_options)
    if (
        not (image_url or ogp_use_first_image)
        and config_social.get('enable') is not False
        and create_social_card is not None
    ):
        image_url = social_card_for_page(
            config_social=config_social,
            site_name=site_name,
            title=title,
            description=description,
            pagename=context['pagename'],
            ogp_site_url=ogp_site_url,
            ogp_canonical_url=ogp_canonical_url,
            srcdir=srcdir,
            outdir=outdir,
            config=config,
            env=env,
        )
        ogp_use_first_image = False

        # Alt text is taken from description unless given
        if 'og:image:alt' in fields:
            ogp_image_alt = fields.get('og:image:alt')
        else:
            ogp_image_alt = description

        # If the social card objects have been added we add special metadata for them
        # These are the dimensions *in pixels* of the card
        # They were chosen by looking at the image pixel dimensions on disk
        tags['og:image:width'] = '1146'
        tags['og:image:height'] = '600'
        meta_tags['twitter:card'] = 'summary_large_image'

    fields.pop('og:image:alt', None)

    first_image = None
    if ogp_use_first_image:
        # Use the first image that is defined in the current page
        first_image = doctree.next_node(nodes.image)
        if (
            first_image
            and Path(first_image.get('uri', '')).suffix[1:].lower() in IMAGE_MIME_TYPES
        ):
            image_url = first_image['uri']
            ogp_image_alt = first_image.get('alt', None)
        else:
            first_image = None

    if image_url:
        # temporarily disable relative image paths with field lists
        if 'og:image' not in fields:
            image_url_parsed = urlparse(image_url)
            if not image_url_parsed.scheme:
                # Relative image path detected, relative to the source. Make absolute.
                if first_image:  # NoQA: SIM108
                    root = page_url
                else:  # ogp_image is set
                    # ogp_image is defined as being relative to the site root.
                    # This workaround is to keep that functionality from breaking.
                    root = ogp_site_url

                image_url = urljoin(root, image_url_parsed.path)
            tags['og:image'] = image_url

        # Add image alt text (either provided by config or from site_name)
        if isinstance(ogp_image_alt, str):
            tags['og:image:alt'] = ogp_image_alt
        elif ogp_image_alt is None and site_name:
            tags['og:image:alt'] = site_name
        elif ogp_image_alt is None and title:
            tags['og:image:alt'] = title

    # arbitrary tags and overrides
    tags.update({k: v for k, v in fields.items() if k.startswith('og:')})

    return (
        '\n'.join(
            [make_tag(p, c) for p, c in tags.items()]
            + [make_tag(p, c, 'name') for p, c in meta_tags.items()]
            + list(config.ogp_custom_meta_tags)
        )
        + '\n'
    )


def ambient_site_url() -> str:
    # readthedocs addons sets the READTHEDOCS_CANONICAL_URL variable
    if rtd_canonical_url := os.getenv('READTHEDOCS_CANONICAL_URL'):
        parse_result = urlsplit(rtd_canonical_url)
    else:
        msg = 'ReadTheDocs did not provide a valid canonical URL!'
        raise RuntimeError(msg)

    # Grab root url from canonical url
    return urlunsplit(
        (parse_result.scheme, parse_result.netloc, parse_result.path, '', '')
    )


def social_card_for_page(
    config_social: dict[str, bool | str],
    site_name: str,
    title: str,
    description: str,
    pagename: str,
    ogp_site_url: str,
    ogp_canonical_url: str,
    *,
    srcdir: str | Path,
    outdir: str | Path,
    config: Config,
    env: BuildEnvironment,
) -> str:
    # Description
    description_max_length = config_social.get(
        'description_max_length', DEFAULT_DESCRIPTION_LENGTH_SOCIAL_CARDS - 3
    )
    if len(description) > description_max_length:
        description = description[:description_max_length].strip() + '...'

    # Page title
    pagetitle = title
    if len(pagetitle) > DEFAULT_PAGE_LENGTH_SOCIAL_CARDS:
        pagetitle = pagetitle[:DEFAULT_PAGE_LENGTH_SOCIAL_CARDS] + '...'

    # Site URL
    site_url = config_social.get('site_url', True)
    if site_url is True:
        url_text = ogp_canonical_url.split('://')[-1]
    elif isinstance(site_url, str):
        url_text = site_url

    # Plot an image with the given metadata to the output path
    image_path = create_social_card(
        config_social,
        site_name,
        pagetitle,
        description,
        url_text,
        pagename,
        srcdir=srcdir,
        outdir=outdir,
        env=env,
        html_logo=config.html_logo,
    )

    # Link the image in our page metadata
    return posixpath.join(ogp_site_url, image_path.as_posix())


def make_tag(property: str, content: str, type_: str = 'property') -> str:
    # Parse quotation, so they won't break html tags if smart quotes are disabled
    content = content.replace('"', '&quot;')
    return f'<meta {type_}="{property}" content="{content}" />'


def setup(app: Sphinx) -> ExtensionMetadata:
    # ogp_site_url="" allows relative by default, even though it's not
    # officially supported by OGP.
    app.add_config_value('ogp_site_url', '', 'html', types=frozenset({str}))
    app.add_config_value('ogp_canonical_url', '', 'html', types=frozenset({str}))
    app.add_config_value(
        'ogp_description_length',
        DEFAULT_DESCRIPTION_LENGTH,
        'html',
        types=frozenset({int}),
    )
    app.add_config_value('ogp_image', None, 'html', types=frozenset({str, NoneType}))
    app.add_config_value(
        'ogp_image_alt', None, 'html', types=frozenset({str, bool, NoneType})
    )
    app.add_config_value('ogp_use_first_image', False, 'html', types=frozenset({bool}))
    app.add_config_value('ogp_type', 'website', 'html', types=frozenset({str}))
    app.add_config_value(
        'ogp_site_name', None, 'html', types=frozenset({str, bool, NoneType})
    )
    app.add_config_value(
        'ogp_social_cards', None, 'html', types=frozenset({dict, NoneType})
    )
    app.add_config_value(
        'ogp_custom_meta_tags', (), 'html', types=frozenset({list, tuple})
    )
    app.add_config_value(
        'ogp_enable_meta_description', True, 'html', types=frozenset({bool})
    )

    # Main Sphinx OpenGraph linking
    app.connect('html-page-context', html_page_context)

    return {
        'version': __version__,
        'env_version': 1,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
