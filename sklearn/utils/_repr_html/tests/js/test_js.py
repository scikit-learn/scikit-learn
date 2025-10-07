from pathlib import Path

import pytest


def _make_page(body):
    """Helper to create a HTML page that includes `estimator.js` and the given body."""

    # use importskip to avoir importing BufferedIOBase
    # useful in CI step: debian_32bit -> Test Library
    js_path = Path(__file__).parent.parent.parent / "estimator.js"
    with open(js_path, "r", encoding="utf-8") as f:
        script = f.read()

    return f"""
    <html>
      <head>
      <script>{script}</script>
      </head>
      <body>
        {body}
      </body>
    </html>
    """


def test_copy_paste(page, local_server):
    """Test that copyToClipboard copies the right text to the clipboard.

    Test requires clipboard permissions, which are granted through page's context.
    Assertion is done by reading back the clipboard content from the browser.
    This is easier than writing a cross platform clipboard reader.
    """
    url, set_html_response = local_server

    copy_paste_html = (
        _make_page('<div class="sk-toggleable__content" data-param-prefix="prefix"/>'),
    )
    set_html_response(copy_paste_html)
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.goto(url)
    page.evaluate(
        "copyToClipboard('test', document.querySelector('.sk-toggleable__content'))"
    )
    clipboard_content = page.evaluate("navigator.clipboard.readText()")
    assert clipboard_content == "prefixtest"


@pytest.mark.parametrize(
    "color,expected_theme",
    [
        (
            "black",
            "light",
        ),
        (
            "white",
            "dark",
        ),
        (
            "#828282",
            "light",
        ),
    ],
)
def test_force_theme(page, local_server, color, expected_theme):
    """Test that forceTheme applies the right theme class to the element.

    A light color must lead to a dark theme and vice-versa.
    """
    url, set_html_response = local_server

    html = _make_page('<div id="test" style="color: ${color};"></div>')
    set_html_response(html.replace("${color}", color))
    page.goto(url)
    page.evaluate("forceTheme('test')")
    assert page.locator("#test").evaluate(
        f"el => el.classList.contains('{expected_theme}')"
    )
