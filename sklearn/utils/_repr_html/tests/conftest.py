import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest


@pytest.fixture(scope="session")
def browser_type_launch_args(browser_type_launch_args):
    """Ensure that the browser is launched in headless mode."""
    return {**browser_type_launch_args, "headless": True}


@pytest.fixture
def local_server(request):
    """
    Start a simple HTTP server that serves custom HTML per test.
    Usage: pass `html_content` as a test parameter.
    """
    html_content = getattr(request, "param", "<html><body>Default</body></html>")
    PORT = 8000

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            self.wfile.write(html_content.encode("utf-8"))

        # suppress logging
        def log_message(self, format, *args):
            return

    httpd = HTTPServer(("localhost", PORT), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    yield f"http://localhost:{PORT}"

    httpd.shutdown()
