"""Handle Pygments css.

inspired by the Furo theme
https://github.com/pradyunsg/furo/blob/main/src/furo/__init__.py
"""

from functools import partial
from pathlib import Path

from pygments.formatters import HtmlFormatter
from pygments.styles import get_all_styles
from sphinx.application import Sphinx

from .utils import get_theme_options_dict, maybe_warn


def _get_styles(formatter: HtmlFormatter, prefix: str) -> None:
    """Get styles out of a formatter, where everything has the correct prefix."""
    for line in formatter.get_linenos_style_defs():
        yield f"{prefix} {line}"
    yield from formatter.get_background_style_defs(prefix)
    yield from formatter.get_token_style_defs(prefix)


def get_pygments_stylesheet(light_style: str, dark_style: str) -> str:
    """Generate the theme-specific pygments.css.

    There is no way to tell Sphinx how the theme handles modes.
    """
    light_formatter = HtmlFormatter(style=light_style)
    dark_formatter = HtmlFormatter(style=dark_style)

    lines = []

    light_prefix = 'html[data-theme="light"] .highlight'
    lines.extend(_get_styles(light_formatter, prefix=light_prefix))

    dark_prefix = 'html[data-theme="dark"] .highlight'
    lines.extend(_get_styles(dark_formatter, prefix=dark_prefix))

    return "\n".join(lines)


def overwrite_pygments_css(app: Sphinx, exception=None):
    """Overwrite pygments.css to allow dynamic light/dark switching.

    Sphinx natively supports config variables `pygments_style` and
    `pygments_dark_style`. However, quoting from
    www.sphinx-doc.org/en/master/development/theming.html#creating-themes

        The pygments_dark_style setting [...is used] when the CSS media query
        (prefers-color-scheme: dark) evaluates to true.

    This does not allow for dynamic switching by the user, so at build time we
    overwrite the pygment.css file so that it embeds 2 versions:

    - the light theme prefixed with "[data-theme="light"]"
    - the dark theme prefixed with "[data-theme="dark"]"

    Fallbacks are defined in this function in case the user-requested (or our
    theme-specified) pygments theme is not available.
    """
    if exception is not None:
        return

    assert app.builder
    theme_options = get_theme_options_dict(app)
    warning = partial(maybe_warn, app)
    pygments_styles = list(get_all_styles())
    fallbacks = dict(light="tango", dark="monokai")

    for light_or_dark, fallback in fallbacks.items():
        # make sure our fallbacks work; if not fall(further)back to "default"
        if fallback not in pygments_styles:
            fallback = pygments_styles[0]  # should resolve to "default"

        # see if user specified a light/dark pygments theme:
        style_key = f"pygments_{light_or_dark}_style"
        style_name = theme_options.get(style_key, None)
        # if not, use the one we set in `theme.conf`:
        if style_name is None and hasattr(app.builder, "theme"):
            style_name = app.builder.theme.get_options()[style_key]

        # make sure we can load the style
        if style_name not in pygments_styles:
            # only warn if user asked for a highlight theme that we can't find
            if style_name is not None:
                warning(
                    f"Highlighting style {style_name} not found by pygments, "
                    f"falling back to {fallback}."
                )
            style_name = fallback

        # assign to the appropriate variable
        if light_or_dark == "light":
            light_theme = style_name
        else:
            dark_theme = style_name

    # re-write pygments.css
    pygments_css = Path(app.builder.outdir) / "_static" / "pygments.css"
    with pygments_css.open("w") as f:
        f.write(get_pygments_stylesheet(light_theme, dark_theme))
