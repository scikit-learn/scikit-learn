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

  // Function to navigate to the all results page
  const navigateToResultsPage = () => {
    const link = document.getElementById("sk-search-all-results-link");
    if (link !== null) {
      // If there is the "see all results" link, just click it
      link.click();
      return;
    }

    const inputBox = document.getElementById("docsearch-input");
    if (inputBox === null || inputBox.value === "") {
      // If we cannot get the input box or the input box is empty, navigate to the
      // all results page with no query
      window.location.assign(searchPageHref);
      return;
    }
    // Navigate to the all results page with query constructed from the input
    const query = new URLSearchParams({ q: inputBox.value });
    window.location.assign(`${searchPageHref}?${query}`);
  };

  // Initialize the Algolia DocSearch widget
  docsearch({
    container: "#docsearch",
    appId: SKLEARN_ALGOLIA_APP_ID,
    apiKey: SKLEARN_ALGOLIA_API_KEY,
    indexName: SKLEARN_ALGOLIA_INDEX_NAME,
    placeholder: "Search the docs ...",
    searchParameters: { attributesToHighlight: ["hierarchy.lvl0"] },
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
      navigate: navigateToResultsPage,
    },
  });

  // The navigator API only works when there are search results; there are cases where
  // there are no hits, e.g. empty query + no history, in which case we need to manually
  // listen to the Enter key
  document.addEventListener("keydown", (e) => {
    if (
      e.key === "Enter" &&
      !e.shiftKey &&
      !e.ctrlKey &&
      !e.metaKey &&
      !e.altKey
    ) {
      const container = document.querySelector(
        ".DocSearch.DocSearch-Container"
      );
      if (container === null) {
        return;
      }
      e.preventDefault();
      e.stopPropagation();
      navigateToResultsPage();
    }
  });
});
