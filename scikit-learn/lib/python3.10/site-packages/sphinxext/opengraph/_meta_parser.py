from __future__ import annotations

from html.parser import HTMLParser


def get_meta_description(meta_tags: str) -> bool:
    htp = HTMLTextParser()
    htp.feed(meta_tags)
    htp.close()

    return htp.meta_description


class HTMLTextParser(HTMLParser):
    """Parse HTML into text."""

    def __init__(self) -> None:
        super().__init__()
        self.meta_description = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        # For example:
        # attrs = [("content", "My manual description"), ("name", "description")]
        if ('name', 'description') in attrs:
            self.meta_description = True
            for name, value in attrs:
                if name == 'content':
                    self.meta_description = value
                    break
