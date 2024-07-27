/**
 * This script is used initialize Algolia DocSearch on each page. It will hydrate the
 * container with ID `docsearch` (see `doc/templates/algolia_docsearch.html`) with the
 * Algolia search widget.
 */

// Figure out how to route to the search page from the current page where we will show
// all search results
let searchPageHref = "./";
for (let i = 0; i < DOCUMENTATION_OPTIONS.pagename.split("/").length - 1; i++) {
  searchPageHref += "../";
}
searchPageHref += "search.html";

docsearch({
  container: "#docsearch",
  appId: "R2IYF7ETH7",
  apiKey: "599cec31baffa4868cae4e79f180729b",
  indexName: "docsearch",
  placeholder: "Search the docs...",
  // Redirect to the search page with the corresponding query
  resultsFooterComponent: ({ state }) => ({
    type: "a",
    ref: undefined,
    constructor: undefined,
    key: state.query,
    props: {
      href: `${searchPageHref}?q=${state.query}`,
      children: `Check all ${state.context.nbHits} results...`,
    },
    __v: null,
  }),
});

// Ctrl-Alt to navigate to the all results page
document.addEventListener("keydown", (e) => {
  if (e.altKey && e.key === "Enter") {
    e.preventDefault();

    const link = document.getElementById("sk-search-all-results-link");
    if (link) {
      link.click();
    }
  }
});
