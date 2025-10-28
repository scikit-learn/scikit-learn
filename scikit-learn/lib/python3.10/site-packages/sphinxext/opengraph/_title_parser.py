from __future__ import annotations

from html.parser import HTMLParser


def get_title(title: str) -> tuple[str, str]:
    htp = HTMLTextParser()
    htp.feed(title)
    htp.close()

    return htp.text, htp.text_outside_tags


class HTMLTextParser(HTMLParser):
    """Parse HTML into text."""

    def __init__(self) -> None:
        super().__init__()
        # All text found
        self.text = ''
        # Only text outside of html tags
        self.text_outside_tags = ''
        self.level = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.level += 1

    def handle_endtag(self, tag: str) -> None:
        self.level -= 1

    def handle_data(self, data: str) -> None:
        self.text += data
        if self.level == 0:
            self.text_outside_tags += data
