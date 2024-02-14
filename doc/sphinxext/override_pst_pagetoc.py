from functools import cache


def override_pst_pagetoc(app, pagename, templatename, context, doctree):
    """Overrides the `generate_toc_html` function of pydata-sphinx-theme for API."""

    @cache
    def generate_api_toc_html(kind="html"):
        """Generate the in-page toc for an API page.

        This relies on the `generate_toc_html` function added by pydata-sphinx-theme
        into the context. We save the original function into `pst_generate_toc_html`
        and override `generate_toc_html` with this function for generated API pages.
        """
        soup = context["pst_generate_toc_html"](kind="soup")

        # On an API page, the pagetoc would look like the following:
        #
        # <ul class="visible nav section-nav flex-column">
        #  <li class="toc-h1 nav-item toc-entry">
        #   <a class="reference internal nav-link" href="#">{{obj}}</a>
        #   <ul class="visible nav section-nav flex-column">
        #    <li class="toc-h2 nav-item toc-entry">
        #     <a class="reference internal nav-link" href="#xxx">
        #      <code class="docutils literal notranslate">
        #       <span class="pre">{{obj}}</span>
        #      </code>
        #     </a>
        #     <ul class="nav section-nav flex-column">  <-- May or may not exist
        #      <li class="toc-h3 nav-item toc-entry">  <-- May have one or more
        #       <a class="reference internal nav-link" href="#xxx">
        #        <code class="docutils literal notranslate">
        #         <span class="pre">{{obj}}.{{method}}</span>
        #        </code>
        #       </a>
        #      </li>
        #     </ul>
        #    </li>
        #    <li class="toc-h2 nav-item toc-entry">    <-- May or may not exist
        #     <a class="reference internal nav-link" href="#xxx">
        #      Examples using
        #      <code class="docutils literal notranslate">
        #       <span class="pre">{{mod}}.{{obj}}</span>
        #      </code>
        #     </a>
        #    </li>
        #   </ul>
        #  </li>
        # </ul>

        try:
            # Unwrap the outermost level level
            soup.ul.unwrap()
            soup.li.unwrap()
            soup.a.decompose()

            # Get all toc-h2 level entries, where the first one should be the function
            # or class, and the second one, if exists, should be the examples; there
            # should theoretically be no more than two entries at this level
            lis = soup.ul.select("li.toc-h2")
            main_li = lis[0]
            meth_list = main_li.ul

            if meth_list is not None:
                # This is a class API page, we remove the class name from the method
                # names to make them shorter; also we make the toc-h3 level entries
                # always visible because of the bootstrap issue that scrollspy cannot
                # correct target elements with IDs containing dots; see
                # https://github.com/twbs/bootstrap/issues/39248
                # TODO: no need to add "visible" once the bootstrap issue is resolved
                # and when pydata-sphinx-theme updates to use that new version; see
                # https://github.com/pydata/pydata-sphinx-theme/issues/1435
                meth_list["class"].append("visible")
                main_li.a.code.span.unwrap()
                main_li.a.code.unwrap()
                to_strip = len(main_li.a.string)
                for meth in meth_list.find_all("li", {"class": "toc-h3"}):
                    meth.a.code.span.unwrap()
                    meth.a.code.unwrap()
                    meth.a.string = meth.a.string[(to_strip + 1) :]

            if len(lis) >= 2:
                # The examples entry exists; shorten the long "Examples using ..."
                # expression to better fit into the sidebar
                example_li = lis[1]
                example_li.a.clear()
                example_li.a.string = "Examples"

            # This corresponds to the behavior of `generate_toc_html`
            return str(soup) if kind == "html" else soup

        except:  # noqa: E722
            # Upon any failure we return an empty pagetoc
            return ""

    # Override the pydata-sphinx-theme implementation for generate API pages
    if pagename.startswith("modules/generated/"):
        context["pst_generate_toc_html"] = context["generate_toc_html"]
        context["generate_toc_html"] = generate_api_toc_html


def setup(app):
    # Need to be triggered after `pydata_sphinx_theme.toctree.add_toctree_functions`,
    # and since default priority is 500 we set 900 for safety
    app.connect("html-page-context", override_pst_pagetoc, priority=900)
