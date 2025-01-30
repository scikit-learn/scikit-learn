"""Methods to build the toctree used in the html pages."""

from dataclasses import dataclass
from functools import cache
from itertools import count
from textwrap import dedent
from typing import Iterator, List, Tuple, Union
from urllib.parse import urlparse

import sphinx

from bs4 import BeautifulSoup
from docutils import nodes
from docutils.nodes import Node
from sphinx.addnodes import toctree as TocTreeNodeClass
from sphinx.application import Sphinx
from sphinx.environment.adapters.toctree import TocTree
from sphinx.locale import _

from .utils import traverse_or_findall


def add_inline_math(node: Node) -> str:
    """Render a node with HTML tags that activate MathJax processing.

    This is meant for use with rendering section titles with math in them, because
    math outputs are ignored by pydata-sphinx-theme's header.

    related to the behaviour of a normal math node from:
    https://github.com/sphinx-doc/sphinx/blob/master/sphinx/ext/mathjax.py#L28
    """
    return (
        '<span class="math notranslate nohighlight">' rf"\({node.astext()}\)" "</span>"
    )


def _get_ancestor_pagename(app: Sphinx, pagename: str, startdepth: int) -> str:
    """
    Get the name of `pagename`'s ancestor that is rooted `startdepth` levels below the
    global root.
    """
    toctree = TocTree(app.env)
    if sphinx.version_info[:2] >= (7, 2):
        from sphinx.environment.adapters.toctree import _get_toctree_ancestors

        ancestors = [*_get_toctree_ancestors(app.env.toctree_includes, pagename)]
    else:
        ancestors = toctree.get_toctree_ancestors(pagename)
    try:
        out = ancestors[-startdepth]
    except IndexError:
        # eg for index.rst, but also special pages such as genindex, py-modindex, search
        # those pages don't have a "current" element in the toctree, so we can
        # directly return None instead of using the default sphinx
        # toctree.get_toctree_for(pagename, app.builder, collapse, **kwargs)
        out = None
    return out, toctree


@dataclass
class LinkInfo:
    """Dataclass to generate toctree data."""

    is_current: bool
    href: str
    title: str
    is_external: bool


def add_toctree_functions(
    app: Sphinx, pagename: str, templatename: str, context, doctree
) -> None:
    """Add functions so Jinja templates can add toctree objects."""

    def suppress_sidebar_toctree(startdepth: int = 1, **kwargs):
        """Check if there's a sidebar TocTree that needs to be rendered.

        Parameters:
            startdepth : The level of the TocTree at which to start. 0 includes the
            entire TocTree for the site; 1 (default) gets the TocTree for the current
            top-level section.

            kwargs : passed to the Sphinx `toctree` template function.
        """
        ancestorname, toctree_obj = _get_ancestor_pagename(
            app=app, pagename=pagename, startdepth=startdepth
        )
        if ancestorname is None:
            return True  # suppress
        if kwargs.get("includehidden", False):
            # if ancestor is found and `includehidden=True` we're guaranteed there's a
            # TocTree to be shown, so don't suppress
            return False

        # we've found an ancestor page, but `includehidden=False` so we can't be sure if
        # there's a TocTree fragment that should be shown on this page; unfortunately we
        # must resolve the whole TOC subtree to find out
        toctree = get_nonroot_toctree(
            app, pagename, ancestorname, toctree_obj, **kwargs
        )
        return toctree is None

    @cache
    def get_or_create_id_generator(base_id: str) -> Iterator[str]:
        for n in count(start=1):
            if n == 1:
                yield base_id
            else:
                yield f"{base_id}-{n}"

    def unique_html_id(base_id: str):
        """
        Create an id that is unique from other ids created by this function at build
        time.

        The function works by sequentially returning "<base_id>", "<base_id>-2",
        "<base_id>-3", etc. each time it is called.
        """
        return next(get_or_create_id_generator(base_id))

    @cache
    def _generate_nav_info() -> List[LinkInfo]:
        """Generate informations necessary to generate nav.

        Instead of messing with html later, having this as a util function
        should make it slightly easier to generate different html snippet for
        sidebar or navbar.
        """
        toctree = TocTree(app.env)

        # Find the active header navigation item so we decide whether to highlight
        # Will be empty if there is no active page (root_doc, or genindex etc)
        if sphinx.version_info[:2] >= (7, 2):
            from sphinx.environment.adapters.toctree import _get_toctree_ancestors

            # NOTE: `env.toctree_includes` is a dict mapping pagenames to any (possibly
            # hidden) TocTree directives on that page (i.e., the "child" pages nested
            # under `pagename`).
            header_pages = [*_get_toctree_ancestors(app.env.toctree_includes, pagename)]
        else:
            header_pages = toctree.get_toctree_ancestors(pagename)
        if header_pages:
            # The final list item will be the top-most ancestor
            active_header_page = header_pages[-1]
        else:
            active_header_page = None

        # NOTE: `env.tocs` is a dict mapping pagenames to hierarchical bullet-lists
        # ("nodetrees" in Sphinx parlance) of in-page headings (including `toctree::`
        # directives). Thus the `tocs` of `root_doc` yields the top-level pages that sit
        # just below the root of our site
        root_toc = app.env.tocs[app.config.root_doc]

        links_data = []

        # Iterate through each node in the root document toc.
        # Grab the toctree pages and find the relative link + title.
        for toc in traverse_or_findall(root_toc, TocTreeNodeClass):
            # TODO: ↑↑↑ use `root_toc.findall(TocTreeNodeClass)` ↑↑↑
            #              once docutils min version >=0.18.1
            for title, page in toc.attributes["entries"]:
                # if the page is using "self" use the correct link
                page = toc.attributes["parent"] if page == "self" else page

                # If this is the active ancestor page, add a class so we highlight it

                # sanitize page title for use in the html output if needed
                if title is None:
                    title = ""
                    for node in app.env.titles[page].children:
                        if isinstance(node, nodes.math):
                            title += add_inline_math(node)
                        else:
                            title += node.astext()

                # set up the status of the link and the path
                # if the path is relative then we use the context for the path
                # resolution and the internal class.
                # If it's an absolute one then we use the external class and
                # the complete url.
                is_absolute = bool(urlparse(page).netloc)
                link_href = page if is_absolute else context["pathto"](page)

                links_data.append(
                    LinkInfo(
                        is_current=(page == active_header_page),
                        href=link_href,
                        title=title,
                        is_external=is_absolute,
                    )
                )

        # Add external links defined in configuration as sibling list items
        for external_link in context["theme_external_links"]:
            links_data.append(
                LinkInfo(
                    is_current=False,
                    href=external_link["url"],
                    title=external_link["name"],
                    is_external=True,
                )
            )

        return links_data

    @cache
    def _generate_header_nav_before_dropdown(
        n_links_before_dropdown,
    ) -> Tuple[str, List[str]]:
        """Return html for navbar and dropdown.

        Given the number of links before the dropdown, return the html for the navbar,
        as well as the list of links to put in a dropdown.

        Returns:
            - HTML str for the navbar
            - list of HTML str for the dropdown
        """
        try:
            n_links_before_dropdown = int(n_links_before_dropdown)
        except Exception:
            raise ValueError(
                f"n_links_before_dropdown is not an int: {n_links_before_dropdown}"
            )
        links_data = _generate_nav_info()

        links_html = []
        links_dropdown = []
        boilerplate = dedent(
            """
            <li class="{nav_item} {active}">
              <a class="{nav_link} nav-{ext_int}" href="{href}">
                {title}
              </a>
            </li>
            """
        )
        nav_item = "nav-item"
        nav_link = "nav-link"
        dropdown_item = "dropdown-item"
        for link in links_data[:n_links_before_dropdown]:
            links_html.append(
                boilerplate.format(
                    active="current active" if link.is_current else "",
                    nav_link=nav_link,
                    nav_item=nav_item,
                    ext_int="external" if link.is_external else "internal",
                    href=link.href,
                    title=link.title,
                )
            )
        for link in links_data[n_links_before_dropdown:]:
            links_dropdown.append(
                boilerplate.format(
                    active="current active" if link.is_current else "",
                    nav_link=nav_link + " " + dropdown_item,
                    nav_item="",
                    ext_int="external" if link.is_external else "internal",
                    href=link.href,
                    title=link.title,
                )
            )

        # The first links will always be visible
        return "\n".join(links_html), links_dropdown

    def generate_header_nav_html(
        n_links_before_dropdown: int = 5, dropdown_text: str = "More"
    ) -> str:
        """Generate top-level links that are meant for the header navigation.

        We use this function instead of the TocTree-based one used for the
        sidebar because this one is much faster for generating the links and
        we don't need the complexity of the full Sphinx TocTree.

        This includes two kinds of links:

        - Links to pages described listed in the root_doc TocTrees
        - External links defined in theme configuration

        Additionally it will create a dropdown list for several links after
        a cutoff.

        Parameters:
            n_links_before_dropdown:The number of links to show before nesting the
                remaining links in a Dropdown element.
            dropdown_text:Text of the dropdown element button.
        """
        out, links_dropdown = _generate_header_nav_before_dropdown(
            n_links_before_dropdown
        )

        if links_dropdown:
            dropdown_id = unique_html_id("pst-nav-more-links")
            links_dropdown_html = "\n".join(links_dropdown)
            out += f"""
            <li class="nav-item dropdown">
                <button class="btn dropdown-toggle nav-item" type="button"
                data-bs-toggle="dropdown" aria-expanded="false"
                aria-controls="{dropdown_id}">
                    {_(dropdown_text)}
                </button>
                <ul id="{dropdown_id}" class="dropdown-menu">
                    {links_dropdown_html}
                </ul>
            </li>
            """

        return out

    # Cache this function because it is expensive to run, and because Sphinx
    # somehow runs this twice in some circumstances in unpredictable ways.
    @cache
    def generate_toctree_html(
        kind: str, startdepth: int = 1, show_nav_level: int = 1, **kwargs
    ) -> Union[BeautifulSoup, str]:
        """Return the navigation link structure in HTML.

        This is similar to Sphinx's own default TocTree generation, but it is modified
        to generate TocTrees for *second*-level pages and below (not supported
        by default in Sphinx).
        This is used for our sidebar, which starts at the second-level page.

        It also modifies the generated TocTree slightly for Bootstrap classes
        and structure (via BeautifulSoup).

        Arguments are passed to Sphinx "toctree" function (context["toctree"] below).

        ref: https://www.sphinx-doc.org/en/master/templating.html#toctree

        Parameters:
            kind : "sidebar" or "raw". Whether to generate HTML meant for sidebar
                navigation ("sidebar") or to return the raw BeautifulSoup object
                ("raw").
            startdepth : The level of the toctree at which to start. By default,
                for the navbar uses the normal toctree (`startdepth=0`), and for the
                sidebar starts from the second level (`startdepth=1`).
            show_nav_level : The level of the navigation bar to toggle as visible on
                page load. By default, this level is 1, and only top-level pages are
                shown, with drop-boxes to reveal children. Increasing `show_nav_level`
                will show child levels as well.
            kwargs : passed to the Sphinx `toctree` template function.

        Returns:
            HTML string (if kind == "sidebar") OR BeautifulSoup object
                (if kind == "raw")
        """
        if startdepth == 0:
            html_toctree = context["toctree"](**kwargs)
        else:
            # find relevant ancestor page; some pages (search, genindex) won't have one
            ancestorname, toctree_obj = _get_ancestor_pagename(
                app=app, pagename=pagename, startdepth=startdepth
            )
            if ancestorname is None:
                raise RuntimeError(
                    "Template requested to generate a TocTree fragment but no suitable "
                    "ancestor found to act as root node. Please report this to theme "
                    "developers."
                )
            # select the "active" subset of the navigation tree for the sidebar
            toctree_element = get_nonroot_toctree(
                app, pagename, ancestorname, toctree_obj, **kwargs
            )
            html_toctree = app.builder.render_partial(toctree_element)["fragment"]

        soup = BeautifulSoup(html_toctree, "html.parser")

        # pair "current" with "active" since that's what we use w/ bootstrap
        for li in soup("li", {"class": "current"}):
            li["class"].append("active")

        # Remove sidebar links to sub-headers on the page
        for li in soup.select("li"):
            # Remove
            if li.find("a"):
                href = li.find("a")["href"]
                if "#" in href and href != "#":
                    li.decompose()

        if kind == "sidebar":
            # Add bootstrap classes for first `ul` items
            for ul in soup("ul", recursive=False):
                ul.attrs["class"] = [*ul.attrs.get("class", []), "nav", "bd-sidenav"]

            # Add collapse boxes for parts/captions.
            # Wraps the TOC part in an extra <ul> to behave like chapters with toggles
            # show_nav_level: 0 means make parts collapsible.
            if show_nav_level == 0:
                partcaptions = soup.find_all("p", attrs={"class": "caption"})
                if len(partcaptions):
                    new_soup = BeautifulSoup(
                        "<ul class='list-caption'></ul>", "html.parser"
                    )
                    for caption in partcaptions:
                        # Assume that the next <ul> element is the TOC list
                        # for this part
                        for sibling in caption.next_siblings:
                            if sibling.name == "ul":
                                toclist = sibling
                                break
                        li = soup.new_tag("li", attrs={"class": "toctree-l0"})
                        li.extend([caption, toclist])
                        new_soup.ul.append(li)
                    soup = new_soup

            # Add icons and labels for collapsible nested sections
            add_collapse_checkboxes(soup)

            # Open the sidebar navigation to the proper depth
            for ii in range(int(show_nav_level)):
                for details in soup.select(f"li.toctree-l{ii} > details"):
                    details["open"] = "open"

        return soup

    @cache
    def generate_toc_html(kind: str = "html") -> BeautifulSoup:
        """Return the within-page TOC links in HTML."""
        if "toc" not in context:
            return ""

        soup = BeautifulSoup(context["toc"], "html.parser")

        # Add toc-hN + visible classes
        def add_header_level_recursive(ul, level):
            if ul is None:
                return
            if level <= (context["theme_show_toc_level"] + 1):
                ul["class"] = [*ul.get("class", []), "visible"]
            for li in ul("li", recursive=False):
                li["class"] = [*li.get("class", []), f"toc-h{level}"]
                add_header_level_recursive(li.find("ul", recursive=False), level + 1)

        add_header_level_recursive(soup.find("ul"), 1)

        # Add in CSS classes for bootstrap
        for ul in soup("ul"):
            ul["class"] = [*ul.get("class", []), "nav", "section-nav", "flex-column"]

        for li in soup("li"):
            li["class"] = [*li.get("class", []), "nav-item", "toc-entry"]
            if li.find("a"):
                a = li.find("a")
                a["class"] = [*a.get("class", []), "nav-link"]

        # If we only have one h1 header, assume it's a title
        h1_headers = soup.select(".toc-h1")
        if len(h1_headers) == 1:
            title = h1_headers[0]
            # If we have no sub-headers of a title then we won't have a TOC
            if not title.select(".toc-h2"):
                out = ""
            else:
                out = title.find("ul")
        # Else treat the h1 headers as sections
        else:
            out = soup

        # Return the toctree object
        if kind == "html":
            return out
        else:
            return soup

    def navbar_align_class() -> List[str]:
        """Return the class that aligns the navbar based on config."""
        align = context.get("theme_navbar_align", "content")
        align_options = {
            "content": ("col-lg-3", "col-lg-9", "me-auto"),
            "left": ("", "", "me-auto"),
            "right": ("", "", "ms-auto"),
        }
        if align not in align_options:
            raise ValueError(
                "Theme option navbar_align must be one of"
                f"{align_options.keys()}, got: {align}"
            )
        return align_options[align]

    context["unique_html_id"] = unique_html_id
    context["generate_header_nav_html"] = generate_header_nav_html
    context["suppress_sidebar_toctree"] = suppress_sidebar_toctree
    context["generate_toctree_html"] = generate_toctree_html
    context["generate_toc_html"] = generate_toc_html
    context["navbar_align_class"] = navbar_align_class


def add_collapse_checkboxes(soup: BeautifulSoup) -> None:
    """Add checkboxes to collapse children in a toctree."""
    # based on https://github.com/pradyunsg/furo

    for element in soup.find_all("li", recursive=True):
        # We check all "li" elements, to add a "current-page" to the correct li.
        classes = element.get("class", [])

        # expanding the parent part explicitly, if present
        if "current" in classes:
            parentli = element.find_parent("li", class_="toctree-l0")
            if parentli:
                parentli.find("details")["open"] = None

        # Nothing more to do, unless this has "children"
        if not element.find("ul"):
            continue

        # Add a class to indicate that this has children.
        element["class"] = [*classes, "has-children"]

        if soup.new_tag is None:
            continue

        # For table of contents nodes that have subtrees, we modify the HTML so
        # that the subtree can be expanded or collapsed in the browser.
        #
        # The HTML markup tree at the parent node starts with this structure:
        #
        # - li.has-children
        #   - a.reference or p.caption
        #   - ul
        #
        # Note the first child of li.has-children is p.caption only if this node
        # is a section heading. (This only happens when show_nav_level is set to
        # 0.)
        #
        # Now we modify the tree structure in one of two ways.
        #
        # (1) If the node holds a section heading, the HTML tree will be
        # modified like so:
        #
        # - li.has-children
        #   - details
        #     - summary
        #       - p.caption
        #       - .toctree-toggle
        #     - ul
        #
        # (2) Otherwise, if the node holds a link to a page in the docs:
        #
        # - li.has-children
        #   - a.reference
        #   - details
        #     - summary
        #       - .toctree-toggle
        #   - ul
        #
        # Why the difference? In the first case, the TOC section heading is not
        # a link, but in the second case it is. So in the first case it makes
        # sense to put the (non-link) text inside the summary tag so that the
        # user can click either the text or the .toctree-toggle chevron icon to
        # expand/collapse the TOC subtree. But in the second case, putting the
        # link in the summary tag would make it unclear whether clicking on the
        # link should expand the subtree or take you to the link.

        # Create <details> and put the entire subtree into it
        details = soup.new_tag("details")
        details.extend(element.contents)
        element.append(details)

        # Hoist the link to the top if there is one
        toc_link = element.select_one("details > a.reference")
        if toc_link:
            element.insert(0, toc_link)

        # Create <summary> with chevron icon
        summary = soup.new_tag("summary")
        span = soup.new_tag(
            "span",
            attrs={
                "class": "toctree-toggle",
                # This element and the chevron it contains are purely decorative;
                # the actual expand/collapse functionality is delegated to the
                # <summary> tag
                "role": "presentation",
            },
        )
        span.append(soup.new_tag("i", attrs={"class": "fa-solid fa-chevron-down"}))
        summary.append(span)

        # Prepend section heading (if there is one) to <summary>
        collapsible_section_heading = element.select_one("details > p.caption")
        if collapsible_section_heading:
            # Put heading inside summary so that the heading text (and chevron) are both
            # clickable
            summary.insert(0, collapsible_section_heading)

        # Prepend <summary> to <details>
        details.insert(0, summary)

        # If this TOC node has a "current" class, be expanded by default
        # (by opening the details/summary disclosure widget)
        if "current" in classes:
            details["open"] = "open"


def get_nonroot_toctree(
    app: Sphinx, pagename: str, ancestorname: str, toctree, **kwargs
):
    """Get the partial TocTree (rooted at `ancestorname`) that dominates `pagename`.

    Parameters:
    app : Sphinx app.
    pagename : Name of the current page (as Sphinx knows it; i.e., its relative path
    from the documentation root).
    ancestorname : Name of a page that dominates `pagename` and that will serve as the
    root of the TocTree fragment.
    toctree : A Sphinx TocTree object. Since this is always needed when finding the
    ancestorname (see _get_ancestor_pagename), it's more efficient to pass it here to
    re-use it.
    kwargs : passed to the Sphinx `toctree` template function.

    This is similar to `context["toctree"](**kwargs)` (AKA `toctree(**kwargs)` within a
    Jinja template), or `TocTree.get_toctree_for()`, which always uses the "root"
    doctree (i.e., `doctree = self.env.get_doctree(self.env.config.root_doc)`).
    """
    kwargs.setdefault("collapse", True)
    if "maxdepth" not in kwargs or not kwargs["maxdepth"]:
        kwargs["maxdepth"] = 0
    kwargs["maxdepth"] = int(kwargs["maxdepth"])
    # starting from ancestor page, recursively parse `toctree::` elements
    ancestor_doctree = toctree.env.tocs[ancestorname].deepcopy()
    toctrees = []

    # for each `toctree::` directive in the ancestor page...
    for toctree_node in traverse_or_findall(ancestor_doctree, TocTreeNodeClass):
        # TODO: ↑↑↑↑↑↑ use `ancestor_doctree.findall(TocTreeNodeClass)` ↑↑↑↑↑↑
        #              once docutils min version >=0.18.1

        # ... resolve that `toctree::` (recursively get children, prune, collapse, etc)
        resolved_toctree = toctree.resolve(
            docname=pagename,
            builder=app.builder,
            toctree=toctree_node,
            **kwargs,
        )
        # ... keep the non-empty ones
        if resolved_toctree:
            toctrees.append(resolved_toctree)
    if not toctrees:
        return None
    # ... and merge them into a single entity
    result = toctrees[0]
    for resolved_toctree in toctrees[1:]:
        result.extend(resolved_toctree.children)
    return result
