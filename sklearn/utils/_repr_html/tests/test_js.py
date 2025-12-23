import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest


@pytest.fixture(scope="session", autouse=True)
def check_playwright():
    """Skip tests if playwright is not installed.

    This fixture is used by the next fixture (which is autouse) to skip all tests
    if playwright is not installed."""
    return pytest.importorskip("playwright")


@pytest.fixture
def local_server(request):
    """Start a simple HTTP server that serves custom HTML per test.

    Usage :

    ```python
    def test_something(page, local_server):
        url, set_html_response = local_server
        set_html_response("<html>...</html>")
        page.goto(url)
        ...
    ```
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        PORT = s.getsockname()[1]

    html_content = "<html><body>Default</body></html>"

    def set_html_response(content):
        nonlocal html_content
        html_content = content

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            self.wfile.write(html_content.encode("utf-8"))

        # suppress logging
        def log_message(self, format, *args):
            return

    httpd = HTTPServer(("127.0.0.1", PORT), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    yield f"http://127.0.0.1:{PORT}", set_html_response

    httpd.shutdown()


def _make_page(body):
    """Helper to create a HTML page that includes `estimator.js` and the given body."""

    js_path = Path(__file__).parent.parent / "estimator.js"
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

    copy_paste_html = _make_page(
        '<div class="sk-toggleable__content" data-param-prefix="prefix"/>'
    )

    set_html_response(copy_paste_html)
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.goto(url)
    page.evaluate(
        "copyToClipboard('test', document.querySelector('.sk-toggleable__content'))"
    )
    clipboard_content = page.evaluate("navigator.clipboard.readText()")

    # `copyToClipboard` function concatenates the `data-param-prefix` attribute
    #  with the first argument. Hence we expect "prefixtest" and not just test.
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

    html = _make_page('<div style="color: ${color};"><div id="test"></div></div>')
    set_html_response(html.replace("${color}", color))
    page.goto(url)
    page.evaluate("forceTheme('test')")
    assert page.locator("#test").evaluate(
        f"el => el.classList.contains('{expected_theme}')"
    )
