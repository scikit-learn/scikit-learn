from pathlib import Path

import pytest


def _make_page(body: str) -> str:
    """Helper to create a HTML page that includes `estimator.js` and the given body."""
    current_dir = Path(__file__).parent
    js_path = current_dir / ".." / "estimator.js"
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
def test_copy_paste(page, local_server) -> None:
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.goto(local_server)
    page.evaluate(
        "copyToClipboard('test', document.querySelector('.sk-toggleable__content'))"
    )
    clipboard_content = page.evaluate("navigator.clipboard.readText()")
    assert clipboard_content == "prefixtest"
