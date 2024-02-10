from bs4 import BeautifulSoup


def tweak_toctree_functions(app, pagename, templatename, context, doctree):
    def generate_header_nav_html_with_substitution(
        n_links_before_dropdown, dropdown_text
    ):
        """Generate the HTML for the header navigation bar with substititions.

        This is a wrapper around the `generate_header_nav_html` function that is added
        to the context by pydata-sphinx-theme. It makes certain substitutions to the
        generated HTML that cannot be achieved by naive sphinx configurations.

        Substitute the "Development" link to point to the development documentation if
        the current build is not a development release.
        https://github.com/scikit-learn/scikit-learn/pull/22550
        """
        generate_header_nav_html = context["generate_header_nav_html"]
        out = generate_header_nav_html(n_links_before_dropdown, dropdown_text)

        if not context["is_devrelease"]:
            soup = BeautifulSoup(out, "html.parser")
            target = soup.find(
                "a",
                class_="nav-link",
                string=lambda text: text.strip() == "Development",
            )

            # Link to the development documentation and change link type to external
            target["href"] = "https://scikit-learn.org/dev/developers/index.html"
            if "nav-internal" in target["class"]:
                target["class"].remove("nav-internal")
            target["class"].append("nav-external")

            return str(soup)

        return out

    context["generate_header_nav_html_with_substitution"] = (
        generate_header_nav_html_with_substitution
    )


def setup(app):
    app.connect("html-page-context", tweak_toctree_functions)
