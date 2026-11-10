import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import numpy as np
import pytest

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils._repr_html.estimator import estimator_html_repr


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
    """Helper to create an HTML page that includes `estimator.js` and the given body."""

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


FEATURE_NAMES_HTML = """
                        <div class="features">
                            <details>
                                <summary>
                                    <span class="image-container">
                                        <button type="button" class="copy-paste-icon"
                                            title="Copy output features (max 100)"
                                            aria-label="Copy output features (max 100)">
                                        </button>
                                    </span>
                                </summary>
                                <div class="features-container">
                                    <table class="features-table">
                                        <tbody>
                                            <tr><td>feature1</td></tr>
                                            <tr><td>feature2</td></tr>
                                        </tbody>
                                    </table>
                                </div>
                            </details>
                        </div>
                     """


def test_copy_paste_feature_names(page, local_server):
    """Test that copyFeatureNamesToClipboard copies the right text to the clipboard.

    Test requires clipboard permissions, which are granted through page's context.
    Assertion is done by reading back the clipboard content from the browser.
    This is easier than writing a cross platform clipboard reader.

    Test adapted from test_copy_paste
    """
    url, set_html_response = local_server

    copy_paste_html = _make_page(FEATURE_NAMES_HTML)

    set_html_response(copy_paste_html)
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.goto(url)
    page.evaluate(
        "copyFeatureNamesToClipboard(document.querySelector('.copy-paste-icon'))"
    )
    clipboard_content = page.evaluate("navigator.clipboard.readText()")

    assert clipboard_content == '[\n    "feature1",\n    "feature2",\n]'


def _fitted_pipeline():
    X = np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
    y = np.array([0, 1, 0, 1])
    return make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)


def _diagram_document(*estimators, prelude=""):
    """Wrap estimator diagrams in one document.

    `estimator_html_repr` already embeds `estimator.js`, so this does not go
    through `_make_page`. `prelude` is placed before the diagrams and can hold
    a control that is outside every `.sk-top-container`.
    """
    fragments = []
    for estimator in estimators:
        fragment = estimator_html_repr(estimator)
        fragment = fragment.replace("<body>", "", 1).replace("</body>", "", 1)
        fragments.append(fragment)
    return (
        "<!doctype html><html><head><meta charset='utf-8'></head>"
        f"<body>{prelude}{''.join(fragments)}</body></html>"
    )


def _open_diagram(page, local_server, html):
    """Serve `html` and record how keydown events propagate.

    The diagram's guard listens on `window` during capture. It runs before a
    notebook listener on `document`, and before the listener installed here on
    `window`, so `__keydowns` can see whether the guard called `preventDefault`.
    """
    url, set_html_response = local_server
    set_html_response(html)
    page.goto(url)
    page.evaluate(
        """() => {
            window.__seenKeys = [];
            window.__keydowns = [];
            window.addEventListener('keydown', (event) => {
                window.__keydowns.push({
                    key: event.key,
                    prevented: event.defaultPrevented,
                });
            }, true);
            document.addEventListener('keydown', (event) => {
                window.__seenKeys.push(event.key);
            }, true);
        }"""
    )


def test_escape_returns_focus_to_the_containing_diagram(page, local_server):
    """Escape inside a diagram moves focus back to that diagram only."""
    _open_diagram(
        page,
        local_server,
        _diagram_document(
            LogisticRegression(),
            _fitted_pipeline(),
            # Focusable element placed before the diagrams, so it is not inside
            # a `.sk-top-container`. Escape here must not be captured.
            prelude="<div id='outside' tabindex='0'>outside</div>",
        ),
    )
    page.locator(".sk-top-container").nth(1).locator("summary").filter(
        has_text="StandardScaler"
    ).focus()
    page.keyboard.press("Escape")

    assert page.evaluate(
        "document.activeElement === document.querySelectorAll('.sk-top-container')[1]"
    )
    # Stopped on window, so a document listener (as in a notebook) does not see it.
    assert page.evaluate("window.__seenKeys") == []

    page.locator("#outside").focus()
    page.keyboard.press("Escape")
    assert page.evaluate("document.activeElement.id") == "outside"
    assert page.evaluate("window.__seenKeys") == ["Escape"]


def test_summary_toggles_from_the_keyboard(page, local_server):
    """Enter and Space still toggle a collapsible element.

    Neither key reaches the document.
    """
    _open_diagram(page, local_server, _diagram_document(_fitted_pipeline()))
    summary = page.locator("summary").filter(has_text="StandardScaler")
    details = summary.locator("xpath=..")

    assert details.evaluate("el => el.open") is False
    summary.focus()
    page.keyboard.press("Enter")
    assert details.evaluate("el => el.open") is True
    page.keyboard.press("Space")
    assert details.evaluate("el => el.open") is False
    assert page.evaluate("window.__seenKeys") == []


def test_space_inside_diagram_is_cancelled(page, local_server):
    """Space on the diagram is cancelled and does not reach the document.

    Space on a summary or button still has its native action. Anywhere else in
    the diagram it would scroll the surrounding page, so the guard cancels it.
    The same key outside a diagram is left alone.
    """
    _open_diagram(
        page,
        local_server,
        _diagram_document(
            LogisticRegression(),
            prelude="<div id='outside' tabindex='0'>outside</div>",
        ),
    )
    page.locator(".sk-top-container").focus()
    page.keyboard.press("Space")
    assert page.evaluate("window.__keydowns") == [{"key": " ", "prevented": True}]
    assert page.evaluate("window.__seenKeys") == []

    page.locator("#outside").focus()
    page.keyboard.press("Space")
    assert page.evaluate("window.__keydowns") == [
        {"key": " ", "prevented": True},
        {"key": " ", "prevented": False},
    ]
    assert page.evaluate("window.__seenKeys") == [" "]


def test_doc_link_activation_does_not_toggle_summary(page, local_server):
    """Activating a documentation link must not toggle its collapsible element."""
    _open_diagram(page, local_server, _diagram_document(_fitted_pipeline()))
    summary = page.locator("summary").filter(has_text="LogisticRegression")
    details = summary.locator("xpath=..")
    link = summary.locator("a.sk-estimator-doc-link")
    page.evaluate(
        """() => {
            window.__linkClicks = 0;
            document.querySelectorAll('a.sk-estimator-doc-link').forEach((link) => {
                link.addEventListener('click', (event) => {
                    // Keep the test on this page while preserving click propagation.
                    event.preventDefault();
                    window.__linkClicks += 1;
                });
            });
        }"""
    )

    # The nested estimator starts collapsed. Neither click nor Enter may open it.
    assert details.evaluate("el => el.open") is False
    link.click()
    assert details.evaluate("el => el.open") is False
    link.focus()
    page.keyboard.press("Enter")
    assert details.evaluate("el => el.open") is False
    assert page.evaluate("window.__linkClicks") == 2


def test_feature_copy_button_does_not_toggle_details(page, local_server):
    """The feature-name copy button does not open its collapsible element."""
    _open_diagram(page, local_server, _diagram_document(_fitted_pipeline()))
    button = page.locator(".features button.copy-paste-icon").first
    details = button.locator("xpath=ancestor::details[1]")
    page.context.grant_permissions(["clipboard-read", "clipboard-write"])
    page.evaluate(
        """() => {
            window.__copyClicks = 0;
            document.querySelectorAll('.features button.copy-paste-icon')
                .forEach((button) => {
                    button.addEventListener('click', () => {
                        window.__copyClicks += 1;
                    });
                });
        }"""
    )

    # The feature list starts collapsed. Neither Enter nor Space may open it.
    assert details.evaluate("el => el.open") is False
    button.focus()
    page.keyboard.press("Enter")
    page.keyboard.press("Space")

    assert details.evaluate("el => el.open") is False
    assert page.evaluate("window.__copyClicks") == 2
    # StandardScaler names the two numpy columns x0 and x1.
    expected_clipboard = '[\n    "x0",\n    "x1",\n]'
    page.wait_for_function(
        """(expected) => navigator.clipboard.readText().then(
            (text) => text === expected
        )""",
        arg=expected_clipboard,
    )
    assert page.evaluate("navigator.clipboard.readText()") == expected_clipboard


def _focus_marker(page):
    """Describe the focused control, its ring, and whether it is expanded."""
    return page.evaluate(
        """() => {
            const el = document.activeElement;
            const style = getComputedStyle(el);
            const details = el.closest('details');
            let kind = 'other';
            if (el.classList.contains('sk-top-container')) {
                kind = 'diagram';
            } else if (el.matches('a.sk-estimator-doc-link')) {
                kind = 'doc-link';
            } else if (el.matches('summary.sk-toggleable__label')) {
                kind = 'estimator';
            } else if (el.matches('summary')) {
                kind = 'parameters';
            }
            let tipDisplay = '';
            if (kind === 'doc-link') {
                tipDisplay = getComputedStyle(el.querySelector('span')).display;
            }
            return {
                kind,
                focusVisible: el.matches(':focus-visible'),
                width: style.outlineWidth,
                lineStyle: style.outlineStyle,
                color: style.outlineColor,
                open: details ? details.open : null,
                tip_display: tipDisplay,
            };
        }"""
    )


def test_keyboard_focus_ring_is_visible(page, local_server):
    """Tab, Shift+Tab, and Enter move through one diagram and toggle it.

    A single estimator starts expanded. The purple ring stays on the focused
    control, and the documentation link shows its tooltip while focused.
    """
    _open_diagram(page, local_server, _diagram_document(LogisticRegression()))

    # step(key, focused control, details.open, tooltip CSS display).
    # The purple ring itself is checked inside `step` for every stop.
    # `expected_tip_display` is the CSS `display` of the "?" link's tooltip
    # span ("Documentation for LogisticRegression"). That span is `display:
    # none` until the link is focused, which sets it to `display: block`.
    # Other stops pass "" because they have no such tooltip.
    def step(key, expected_kind, expected_open, expected_tip_display=""):
        page.keyboard.press(key)
        assert _focus_marker(page) == {
            "kind": expected_kind,
            "focusVisible": True,
            "width": "2px",
            "lineStyle": "solid",
            "color": "rgb(143, 60, 224)",
            "open": expected_open,
            "tip_display": expected_tip_display,
        }

    # The diagram is the first stop. It is not collapsible, so `open` is None.
    step("Tab", "diagram", None)
    # A single estimator starts expanded. Enter collapses it, then expands it.
    step("Tab", "estimator", True)
    step("Enter", "estimator", False)
    step("Enter", "estimator", True)
    # Focusing the "?" link reveals its tooltip (`display: block`). The
    # collapsible estimator stays open.
    step("Tab", "doc-link", True, expected_tip_display="block")
    # Parameters starts collapsed. Enter opens it, then closes it.
    step("Tab", "parameters", False)
    step("Enter", "parameters", True)
    step("Enter", "parameters", False)
    # Shift+Tab walks back. Collapsing Parameters did not collapse the estimator.
    # Back on the "?" link, the tooltip is visible again.
    step("Shift+Tab", "doc-link", True, expected_tip_display="block")
    step("Shift+Tab", "estimator", True)
    step("Shift+Tab", "diagram", None)
