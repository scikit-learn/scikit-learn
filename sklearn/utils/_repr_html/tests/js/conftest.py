import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest


@pytest.fixture(scope="session")
def check_playwright():
    """Skip tests if playwright is not installed.

    This fixture is used by the next fixture (which is autouse) to skip all tests
    if playwright is not installed."""
    return pytest.importorskip("playwright")


@pytest.fixture(scope="session", autouse=True)
def browser_type_launch_args(check_playwright, browser_type_launch_args):
    """Ensure that the browser is launched in headless mode."""
    return {**browser_type_launch_args, "headless": True}


@pytest.fixture
def local_server(request):
    """Start a simple HTTP server that serves custom HTML per test.

    Usage :

    ```python
    def test_something(page, local_server):
        url, set_html_response = local_server
        set_html_response("<html>...</html>")
        page.goto(local_server)
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
