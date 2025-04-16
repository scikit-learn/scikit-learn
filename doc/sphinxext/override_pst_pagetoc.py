from functools import cache

from sphinx.util.logging import getLogger

logger = getLogger(__name__)


def override_pst_pagetoc(app, pagename, templatename, context, doctree):
    """Overrides the `generate_toc_html` function of pydata-sphinx-theme for API."""

    @cache
    def generate_api_toc_html(kind="html"):
        """Generate the in-page toc for an API page.

        This relies on the `generate_toc_html` function added by pydata-sphinx-theme
        into the context. We save the original function into `pst_generate_toc_html`
        and override `generate_toc_html` with this function for generated API pages.

        The pagetoc of an API page would look like the following:

        <ul class="visible ...">               <-- Unwrap
         <li class="toc-h1 ...">               <-- Unwrap
          <a class="..." href="#">{{obj}}</a>  <-- Decompose

          <ul class="visible ...">
           <li class="toc-h2 ...">
            ...object
            <ul class="...">                          <-- Set visible if exists
             <li class="toc-h3 ...">...method 1</li>  <-- Shorten
             <li class="toc-h3 ...">...method 2</li>  <-- Shorten
             ...more methods                          <-- Shorten
            </ul>
           </li>
           <li class="toc-h2 ...">...gallery examples</li>
          </ul>

         </li>                                 <-- Unwrapped
        </ul>                                  <-- Unwrapped
        """
        soup = context["pst_generate_toc_html"](kind="soup")

        try:
            # Unwrap the outermost level
            soup.ul.unwrap()
            soup.li.unwrap()
            soup.a.decompose()

            # Get all toc-h2 level entries, where the first one should be the function
            # or class, and the second one, if exists, should be the examples; there
            # should be no more than two entries at this level for generated API pages
            lis = soup.ul.select("li.toc-h2")
            main_li = lis[0]
            meth_list = main_li.ul

            if meth_list is not None:
                # This is a class API page, we remove the class name from the method
                # names to make them better fit into the secondary sidebar; also we
                # make the toc-h3 level entries always visible to more easily navigate
                # through the methods
                meth_list["class"].append("visible")
                for meth in meth_list.find_all("li", {"class": "toc-h3"}):
                    target = meth.a.code.span
                    target.string = target.string.split(".", 1)[1]

            # This corresponds to the behavior of `generate_toc_html`
            return str(soup) if kind == "html" else soup

        except Exception as e:
            # Upon any failure we return the original pagetoc
            logger.warning(
                f"Failed to generate API pagetoc for {pagename}: {e}; falling back"
            )
            return context["pst_generate_toc_html"](kind=kind)

    # Override the pydata-sphinx-theme implementation for generate API pages
    if pagename.startswith("modules/generated/"):
        context["pst_generate_toc_html"] = context["generate_toc_html"]
        context["generate_toc_html"] = generate_api_toc_html


def setup(app):
    # Need to be triggered after `pydata_sphinx_theme.toctree.add_toctree_functions`,
    # and since default priority is 500 we set 900 for safety
    app.connect("html-page-context", override_pst_pagetoc, priority=900)
