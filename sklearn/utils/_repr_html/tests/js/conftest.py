import socket
import threading

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
    @pytest.mark.parametrize(
        "local_server",
        ["<html><body>test body</body></html>",],
        indirect=True,
    )
    def test_something(page, local_server):
        page.goto(local_server)
        ...
    ```
    """
    # use importskip to avoir importing BufferedIOBase
    # useful in CI step: debian_32bit -> Test Library
    http_server = pytest.importorskip("http.server")
    html_content = getattr(request, "param", "<html><body>Default</body></html>")
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        PORT = s.getsockname()[1]

    class Handler(http_server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            self.wfile.write(html_content.encode("utf-8"))

        # suppress logging
        def log_message(self, format, *args):
            return

    httpd = http_server.HTTPServer(("127.0.0.1", PORT), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    yield f"http://127.0.0.1:{PORT}"

    httpd.shutdown()
