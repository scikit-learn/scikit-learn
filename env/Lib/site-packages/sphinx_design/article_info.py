from typing import Optional

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.application import Sphinx

from .icons import get_octicon
from .shared import SEMANTIC_COLORS, SdDirective, create_component, make_choice


def setup_article_info(app: Sphinx):
    """Setup the article information components."""
    app.add_directive("article-info", ArticleInfoDirective)


class ArticleInfoDirective(SdDirective):
    """ """

    has_content = False
    required_arguments = 0
    optional_arguments = 0
    option_spec = {
        "avatar": directives.uri,
        "avatar-alt": directives.unchanged,
        "avatar-link": directives.uri,
        "avatar-outline": make_choice(SEMANTIC_COLORS),
        "author": directives.unchanged_required,
        "date": directives.unchanged_required,
        "read-time": directives.unchanged_required,
        "class-container": directives.class_option,
        "class-avatar": directives.class_option,
    }

    def _parse_text(
        self, text: str, icon: Optional[nodes.Node] = None, parse: bool = False
    ) -> nodes.Node:
        """Parse the text."""
        if not parse:
            output = ([icon] if icon else []) + [nodes.Text(text)]
        else:
            text_nodes, _ = self.state.inline_text(text, self.lineno)
            text_nodes = ([icon] if icon else []) + text_nodes
            # note certain nodes (like references) need to be nested in a TextElement node
            # (e.g. a pargraph)
            para = nodes.paragraph("", "", *text_nodes, classes=["sd-p-0", "sd-m-0"])
            self.set_source_info(para)
            output = [para]
        return output

    def run_with_defaults(self) -> list[nodes.Node]:  # noqa: PLR0915
        parse_fields = True  # parse field text

        top_grid = create_component(
            "grid-container",
            [
                "sd-container-fluid",
                "sd-sphinx-override",
                "sd-p-0",
                "sd-mt-2",
                "sd-mb-4",
                *self.options.get("class-container", []),
            ],
        )
        self.set_source_info(top_grid)

        top_row = create_component(
            "grid-row",
            ["sd-row", "sd-row-cols-2", "sd-gx-2", "sd-gy-1"],
        )
        self.set_source_info(top_row)
        top_grid += top_row

        avatar_uri = self.options.get("avatar")
        if avatar_uri:
            # TODO only in html (hide in latex)
            avatar_column = create_component(
                "grid-item",
                ["sd-col", "sd-col-auto", "sd-d-flex-row", "sd-align-minor-center"],
            )
            self.set_source_info(avatar_column)
            avatar_classes = ["sd-avatar-sm"]
            if "avatar-outline" in self.options:
                avatar_classes.append(f"sd-outline-{self.options['avatar-outline']}")
            if "class-avatar" in self.options:
                avatar_classes += self.options["class-avatar"]
            avatar_image = nodes.image(
                "",
                uri=avatar_uri,
                alt=self.options.get("avatar-alt", ""),
                classes=avatar_classes,
            )
            self.set_source_info(avatar_image)
            if self.options.get("avatar-link"):
                avatar_link = nodes.reference(
                    "", "", refuri=self.options.get("avatar-link")
                )
                avatar_link += avatar_image
                avatar_image = avatar_link
            avatar_column += avatar_image
            top_row += avatar_column

        info_column = create_component(
            "grid-item",
            ["sd-col", "sd-d-flex-row", "sd-align-minor-center"],
        )
        self.set_source_info(info_column)
        top_row += info_column

        info_grid = create_component(
            "grid-container",
            [
                "sd-container-fluid",
                "sd-sphinx-override",
            ],
        )
        self.set_source_info(info_grid)
        info_column += info_grid

        info_row = create_component(
            "grid-row",
            [
                "sd-row",
                "sd-row-cols-2",
                "sd-row-cols-xs-2",
                "sd-row-cols-sm-3",
                "sd-row-cols-md-3",
                "sd-row-cols-lg-3",
                "sd-gx-3",
                "sd-gy-1",
            ],
        )
        self.set_source_info(info_row)
        info_grid += info_row

        author_text = self.options.get("author")
        if author_text:
            author_column = create_component(
                "grid-item",
                ["sd-col", "sd-col-auto", "sd-d-flex-row", "sd-align-minor-center"],
            )
            self.set_source_info(author_column)
            author_nodes = self._parse_text(author_text, parse=parse_fields)
            author_column.extend(author_nodes)
            info_row += author_column

        date_text = self.options.get("date")
        if date_text:
            date_column = create_component(
                "grid-item",
                ["sd-col", "sd-col-auto", "sd-d-flex-row", "sd-align-minor-center"],
            )
            self.set_source_info(date_column)
            date_icon = nodes.raw(
                "",
                nodes.Text(get_octicon("calendar", height="16px")),
                classes=["sd-pr-2"],
                format="html",
            )
            date_nodes = self._parse_text(date_text, icon=date_icon, parse=parse_fields)
            date_column.extend(date_nodes)
            info_row += date_column

        read_time_text = self.options.get("read-time")
        if read_time_text:
            read_time_column = create_component(
                "grid-item",
                ["sd-col", "sd-col-auto", "sd-d-flex-row", "sd-align-minor-center"],
            )
            self.set_source_info(read_time_column)
            read_time_icon = nodes.raw(
                "",
                nodes.Text(get_octicon("clock", height="16px")),
                classes=["sd-pr-2"],
                format="html",
            )
            read_time_nodes = self._parse_text(
                read_time_text, icon=read_time_icon, parse=parse_fields
            )
            read_time_column.extend(read_time_nodes)
            info_row += read_time_column

        return [top_grid]
