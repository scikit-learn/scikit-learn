/**
 * This script is used initialize Algolia DocSearch on each page. It will hydrate the
 * container with ID `docsearch` (see `doc/templates/algolia-searchbox.html`) with the
 * Algolia search widget.
 */

document.addEventListener("DOMContentLoaded", () => {
  // Figure out how to route to the search page from the current page where we will show
  // all search results
  const pagename = DOCUMENTATION_OPTIONS.pagename;
  let searchPageHref = "./";
  for (let i = 0; i < pagename.split("/").length - 1; i++) {
    searchPageHref += "../";
  }
  searchPageHref += "algolia-search.html";

  // Initialize the Algolia DocSearch widget
  docsearch({
    container: "#docsearch",
    appId: SKLEARN_ALGOLIA_APP_ID,
    apiKey: SKLEARN_ALGOLIA_API_KEY,
    indexName: SKLEARN_ALGOLIA_INDEX_NAME,
    placeholder: "Search the docs ...",
    // Redirect to the search page with the corresponding query
    resultsFooterComponent: ({ state }) => ({
      type: "a",
      ref: undefined,
      constructor: undefined,
      key: state.query,
      props: {
        id: "sk-search-all-results-link",
        href: `${searchPageHref}?q=${state.query}`,
        children: `Check all results...`,
      },
      __v: null,
    }),
    navigator: {
      // Hack implementation to navigate to the search page instead of navigating to the
      // corresponding search result page; `navigateNewTab` and `navigateNewWindow` are
      // still left as the default behavior
      navigate: () => {
        const link = document.getElementById("sk-search-all-results-link");
        if (link) {
          link.click();
        } else {
          window.location.assign(searchPageHref);
        }
      },
    },
  });
});
