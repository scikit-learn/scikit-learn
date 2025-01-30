"""Bootstrap-based sphinx theme from the PyData community."""

import json

from functools import partial
from pathlib import Path
from typing import Dict
from urllib.parse import urlparse

import requests

from requests.exceptions import ConnectionError, HTTPError, RetryError
from sphinx.application import Sphinx
from sphinx.builders.dirhtml import DirectoryHTMLBuilder
from sphinx.errors import ExtensionError

from . import edit_this_page, logo, pygments, short_link, toctree, translator, utils


__version__ = "0.16.1"


def update_config(app):
    """Update config with new default values and handle deprecated keys."""
    # By the time `builder-inited` happens, `app.builder.theme_options` already exists.
    # At this point, modifying app.config.html_theme_options will NOT update the
    # page's HTML context (e.g. in jinja, `theme_keyword`).
    # To do this, you must manually modify `app.builder.theme_options`.
    theme_options = utils.get_theme_options_dict(app)
    warning = partial(utils.maybe_warn, app)

    # TODO: DEPRECATE after v1.0
    themes = ["light", "dark"]
    for theme in themes:
        if style := theme_options.get(f"pygment_{theme}_style"):
            theme_options[f"pygments_{theme}_style"] = style
            warning(
                f'The parameter "pygment_{theme}_style" was renamed to '
                f'"pygments_{theme}_style" (note the "s" on "pygments").'
            )

    # Validate icon links
    if not isinstance(theme_options.get("icon_links", []), list):
        raise ExtensionError(
            "`icon_links` must be a list of dictionaries, you provided "
            f"type {type(theme_options.get('icon_links'))}."
        )

    # Set the anchor link default to be # if the user hasn't provided their own
    if not utils.config_provided_by_user(app, "html_permalinks_icon"):
        app.config.html_permalinks_icon = "#"

    # check the validity of the theme switcher file
    is_dict = isinstance(theme_options.get("switcher"), dict)
    should_test = theme_options.get("check_switcher", True)
    if is_dict and should_test:
        theme_switcher = theme_options.get("switcher")

        # raise an error if one of these compulsory keys is missing
        json_url = theme_switcher["json_url"]
        theme_switcher["version_match"]

        # try to read the json file. If it's a url we use request,
        # else we simply read the local file from the source directory
        # display a log warning if the file cannot be reached
        reading_error = None
        if urlparse(json_url).scheme in ["http", "https"]:
            try:
                request = requests.get(json_url)
                request.raise_for_status()
                content = request.text
            except (ConnectionError, HTTPError, RetryError) as e:
                reading_error = repr(e)
        else:
            try:
                content = Path(app.srcdir, json_url).read_text()
            except FileNotFoundError as e:
                reading_error = repr(e)

        if reading_error is not None:
            warning(
                f'The version switcher "{json_url}" file cannot be read due to '
                f"the following error:\n{reading_error}"
            )
        else:
            # check that the json file is not illformed,
            # throw a warning if the file is ill formed and an error if it's not json
            switcher_content = json.loads(content)
            missing_url = any(["url" not in e for e in switcher_content])
            missing_version = any(["version" not in e for e in switcher_content])
            if missing_url or missing_version:
                warning(
                    f'The version switcher "{json_url}" file is malformed; '
                    'at least one of the items is missing the "url" or "version" key'
                )

    # Add an analytics ID to the site if provided
    analytics = theme_options.get("analytics", {})
    if analytics:
        # Plausible analytics
        plausible_domain = analytics.get("plausible_analytics_domain")
        plausible_url = analytics.get("plausible_analytics_url")

        # Ref: https://plausible.io/docs/plausible-script
        if plausible_domain and plausible_url:
            kwargs = {
                "loading_method": "defer",
                "data-domain": plausible_domain,
                "filename": plausible_url,
            }
            app.add_js_file(**kwargs)

        # Google Analytics
        gid = analytics.get("google_analytics_id")
        if gid:
            gid_js_path = f"https://www.googletagmanager.com/gtag/js?id={gid}"
            gid_script = f"""
                window.dataLayer = window.dataLayer || [];
                function gtag(){{ dataLayer.push(arguments); }}
                gtag('js', new Date());
                gtag('config', '{gid}');
            """

            # Link the JS files
            app.add_js_file(gid_js_path, loading_method="async")
            app.add_js_file(None, body=gid_script)

    # Update ABlog configuration default if present
    fa_provided = utils.config_provided_by_user(app, "fontawesome_included")
    if "ablog" in app.config.extensions and not fa_provided:
        app.config.fontawesome_included = True

    # Handle icon link shortcuts
    shortcuts = [
        ("twitter_url", "fa-brands fa-square-twitter", "Twitter"),
        ("bitbucket_url", "fa-brands fa-bitbucket", "Bitbucket"),
        ("gitlab_url", "fa-brands fa-square-gitlab", "GitLab"),
        ("github_url", "fa-brands fa-square-github", "GitHub"),
    ]
    # Add extra icon links entries if there were shortcuts present
    # TODO: Deprecate this at some point in the future?
    icon_links = theme_options.get("icon_links", [])
    for url, icon, name in shortcuts:
        if theme_options.get(url):
            # This defaults to an empty list so we can always insert
            icon_links.insert(
                0,
                {
                    "url": theme_options.get(url),
                    "icon": icon,
                    "name": name,
                    "type": "fontawesome",
                },
            )
    theme_options["icon_links"] = icon_links

    # Prepare the logo config dictionary
    theme_logo = theme_options.get("logo")
    if not theme_logo:
        # In case theme_logo is an empty string
        theme_logo = {}
    if not isinstance(theme_logo, dict):
        raise ValueError(f"Incorrect logo config type: {type(theme_logo)}")
    theme_logo_link = theme_options.get("theme_logo_link")
    if theme_logo_link:
        theme_logo["link"] = theme_logo_link
    theme_options["logo"] = theme_logo


def update_and_remove_templates(
    app: Sphinx, pagename: str, templatename: str, context, doctree
) -> None:
    """Update template names and assets for page build."""
    # Allow for more flexibility in template names
    template_sections = [
        "theme_navbar_start",
        "theme_navbar_center",
        "theme_navbar_persistent",
        "theme_navbar_end",
        "theme_article_header_start",
        "theme_article_header_end",
        "theme_article_footer_items",
        "theme_content_footer_items",
        "theme_footer_start",
        "theme_footer_center",
        "theme_footer_end",
        "theme_primary_sidebar_end",
        "sidebars",
    ]
    for section in template_sections:
        if context.get(section):
            context[section] = utils._update_and_remove_templates(
                app=app,
                context=context,
                templates=context.get(section, []),
                section=section,
                templates_skip_empty_check=["sidebar-nav-bs.html", "navbar-nav.html"],
            )

    # Remove a duplicate entry of the theme CSS. This is because it is in both:
    # - theme.conf
    # - manually linked in `webpack-macros.html`
    if "css_files" in context:
        theme_css_name = "_static/styles/pydata-sphinx-theme.css"
        for i in range(len(context["css_files"])):
            asset = context["css_files"][i]
            # TODO: eventually the contents of context['css_files'] etc should probably
            # only be _CascadingStyleSheet etc. For now, assume mixed with strings.
            asset_path = getattr(asset, "filename", str(asset))
            if asset_path == theme_css_name:
                del context["css_files"][i]
                break
    # Add links for favicons in the topbar
    for favicon in context.get("theme_favicons", []):
        icon_type = Path(favicon["href"]).suffix.strip(".")
        opts = {
            "rel": favicon.get("rel", "icon"),
            "sizes": favicon.get("sizes", "16x16"),
            "type": f"image/{icon_type}",
        }
        if "color" in favicon:
            opts["color"] = favicon["color"]
        # Sphinx will auto-resolve href if it's a local file
        app.add_css_file(favicon["href"], **opts)

    # Add metadata to DOCUMENTATION_OPTIONS so that we can re-use later
    # Pagename to current page
    app.add_js_file(None, body=f"DOCUMENTATION_OPTIONS.pagename = '{pagename}';")
    if isinstance(context.get("theme_switcher"), dict):
        theme_switcher = context["theme_switcher"]
        json_url = theme_switcher["json_url"]
        version_match = theme_switcher["version_match"]

        # Add variables to our JavaScript for re-use in our main JS script
        js = f"""
        DOCUMENTATION_OPTIONS.theme_version = '{__version__}';
        DOCUMENTATION_OPTIONS.theme_switcher_json_url = '{json_url}';
        DOCUMENTATION_OPTIONS.theme_switcher_version_match = '{version_match}';
        DOCUMENTATION_OPTIONS.show_version_warning_banner =
            {str(context["theme_show_version_warning_banner"]).lower()};
        """
        app.add_js_file(None, body=js)

    # Update version number for the "made with version..." component
    context["theme_version"] = __version__


def _fix_canonical_url(
    app: Sphinx, pagename: str, templatename: str, context: dict, doctree
) -> None:
    """Fix the canonical URL when using the dirhtml builder.

    Sphinx builds a canonical URL if ``html_baseurl`` config is set. However,
    it builds a URL ending with ".html" when using the dirhtml builder, which is
    incorrect. Detect this and generate the correct URL for each page.

    Workaround for https://github.com/sphinx-doc/sphinx/issues/9730; can be removed
    when that is fixed, released, and available in our minimum supported Sphinx version.
    """
    if (
        not app.config.html_baseurl
        or not isinstance(app.builder, DirectoryHTMLBuilder)
        or not context["pageurl"]
        or not context["pageurl"].endswith(".html")
    ):
        return

    target = app.builder.get_target_uri(pagename)
    context["pageurl"] = app.config.html_baseurl + target


def setup(app: Sphinx) -> Dict[str, str]:
    """Setup the Sphinx application."""
    here = Path(__file__).parent.resolve()
    theme_path = here / "theme" / "pydata_sphinx_theme"

    app.add_html_theme("pydata_sphinx_theme", str(theme_path))

    app.add_post_transform(short_link.ShortenLinkTransform)

    app.connect("builder-inited", translator.setup_translators)
    app.connect("builder-inited", update_config)
    app.connect("html-page-context", _fix_canonical_url)
    app.connect("html-page-context", edit_this_page.setup_edit_url)
    app.connect("html-page-context", toctree.add_toctree_functions)
    app.connect("html-page-context", update_and_remove_templates)
    app.connect("html-page-context", logo.setup_logo_path)
    app.connect("html-page-context", utils.set_secondary_sidebar_items)
    app.connect("build-finished", pygments.overwrite_pygments_css)
    app.connect("build-finished", logo.copy_logo_images)

    # https://www.sphinx-doc.org/en/master/extdev/i18n.html#extension-internationalization-i18n-and-localization-l10n-using-i18n-api
    app.add_message_catalog("sphinx", here / "locale")

    # Include component templates
    app.config.templates_path.append(str(theme_path / "components"))

    return {"parallel_read_safe": True, "parallel_write_safe": True}
