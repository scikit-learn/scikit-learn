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
    placeholder: "Search the docs ... (Alt+Enter to go to search page)",
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
  });

  // Ctrl-Alt to navigate to the all results page
  document.addEventListener("keydown", (e) => {
    if (e.altKey && e.key === "Enter") {
      e.preventDefault();

      // Click the link if it exists so the query is preserved; otherwise navigate to
      // the search page without query
      const link = document.getElementById("sk-search-all-results-link");
      if (link) {
        link.click();
      } else {
        window.location.href = searchPageHref;
      }
    }
  });
});
