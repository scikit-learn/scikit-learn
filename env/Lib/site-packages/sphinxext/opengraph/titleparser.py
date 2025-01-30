from html.parser import HTMLParser


class HTMLTextParser(HTMLParser):
    """
    Parse HTML into text
    """

    def __init__(self):
        super().__init__()
        # All text found
        self.text = ""
        # Only text outside of html tags
        self.text_outside_tags = ""
        self.level = 0

    def handle_starttag(self, tag, attrs) -> None:
        self.level += 1

    def handle_endtag(self, tag) -> None:
        self.level -= 1

    def handle_data(self, data) -> None:
        self.text += data
        if self.level == 0:
            self.text_outside_tags += data


def get_title(title: str, skip_html_tags: bool = False):
    htp = HTMLTextParser()
    htp.feed(title)
    htp.close()

    if skip_html_tags:
        return htp.text_outside_tags
    else:
        return htp.text
