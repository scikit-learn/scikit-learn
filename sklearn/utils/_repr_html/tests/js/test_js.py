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


@pytest.mark.parametrize(
    "local_server",
    [
        _make_page(
            '<div class="sk-toggleable__content fitted" data-param-prefix="prefix"/>'
        )
    ],
    indirect=True,
)
def test_copy_paste(page, local_server):
    """Test that copyToClipboard copies the right text to the clipboard.

    Test requires clipboard permissions, which are granted through page's context.
    Assertion is done by reading back the clipboard content from the browser.
    This is easier than writing a cross platform clipboard reader.
    """
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.goto(local_server)
    page.evaluate(
        "copyToClipboard('test', document.querySelector('.sk-toggleable__content'))"
    )
    clipboard_content = page.evaluate("navigator.clipboard.readText()")
    assert clipboard_content == "prefixtest"


@pytest.mark.parametrize(
    "local_server,theme",
    [
        (
            _make_page('<div id="test" style="color: black;"></div>'),
            "light",
        ),
        (
            _make_page('<div id="test" style="color: white;"></div>'),
            "dark",
        ),
        (
            _make_page('<div id="test" style="color: #828282;"></div>'),
            "light",
        ),
    ],
    indirect=["local_server"],
)
def test_force_theme(page, local_server, theme):
    """Test that forceTheme applies the right theme class to the element.

    A light color must lead to a dark theme and vice-versa.
    """
    page.goto(local_server)
    page.evaluate("forceTheme('test')")
    assert page.locator("#test").evaluate(f"el => el.classList.contains('{theme}')")
