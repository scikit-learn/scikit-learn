/**
 * This script is used initialize Algolia DocSearch on the Algolia search page. It will
 * hydrate the search page (see `doc/templates/search.html`) and activate the search
 * functionalities.
 */

document.addEventListener("DOMContentLoaded", () => {
  let timer;
  const timeout = 500; // Debounce search-as-you-type

  const searchClient = algoliasearch(
    SKLEARN_ALGOLIA_APP_ID,
    SKLEARN_ALGOLIA_API_KEY
  );

  const search = instantsearch({
    indexName: SKLEARN_ALGOLIA_INDEX_NAME,
    initialUiState: {
      [SKLEARN_ALGOLIA_INDEX_NAME]: {
        query: new URLSearchParams(window.location.search).get("q") || "",
      },
    },
    searchClient,
  });

  search.addWidgets([
    // The powered-by widget as the heading
    instantsearch.widgets.poweredBy({
      container: "#docsearch-powered-by-light",
      theme: "light",
    }),
    instantsearch.widgets.poweredBy({
      container: "#docsearch-powered-by-dark",
      theme: "dark",
    }),
    // The search input box
    instantsearch.widgets.searchBox({
      container: "#docsearch-container",
      placeholder: "Search the docs ...",
      autofocus: true,
      // Debounce the search input to avoid making too many requests
      queryHook(query, refine) {
        clearTimeout(timer);
        timer = setTimeout(() => refine(query), timeout);
      },
    }),
    // The search statistics before the list of results
    instantsearch.widgets.stats({
      container: "#docsearch-stats",
      templates: {
        text: (data, { html }) => {
          if (data.query === "") {
            return "";
          }

          let count;
          if (data.hasManyResults) {
            count = `${data.nbHits} results`;
          } else if (data.hasOneResult) {
            count = "1 result";
          } else {
            count = "no results";
          }

          const stats = `Search finished, found ${count} matching the search query in ${data.processingTimeMS}ms.`;
          return html`
            <div class="sk-search-stats-heading">Search Results</div>
            <p class="sk-search-stats">${stats}</p>
          `;
        },
      },
    }),
    // The list of search results
    instantsearch.widgets.infiniteHits({
      container: "#docsearch-hits",
      transformItems: (items, { results }) => {
        if (results.query === "") {
          return [];
        }
        return items;
      },
      templates: {
        item: (hit, { html, components }) => {
          const hierarchy = Object.entries(hit._highlightResult.hierarchy);
          const lastKey = hierarchy[hierarchy.length - 1][0];

          const sharedHTML = html`
            <a class="sk-search-item-header" href="${hit.url}">
              ${components.Highlight({
                hit,
                attribute: `hierarchy.${lastKey}`,
              })}
            </a>
            <div class="sk-search-item-path">
              ${components.Highlight({ hit, attribute: "hierarchy.lvl0" })}
              ${hierarchy.slice(1, -1).map(([key, _]) => {
                return html`
                  <span class="sk-search-item-path-divider">»</span>
                  ${components.Highlight({
                    hit,
                    attribute: `hierarchy.${key}`,
                  })}
                `;
              })}
            </div>
          `;

          if (hit.type === "content") {
            return html`
              ${sharedHTML}
              <p class="sk-search-item-context">
                ${components.Snippet({ hit, attribute: "content" })}
              </p>
            `;
          } else {
            return sharedHTML;
          }
        },
        // We have stats widget that can imply "no results"
        empty: () => {
          return "";
        },
      },
    }),
    // Additional configuration of the widgets
    instantsearch.widgets.configure({
      hitsPerPage: 50,
      attributesToSnippet: ["content:60"], // Lengthen snippets to show more context
    }),
  ]);

  search.start();

  // Apart from the loading indicator in the search form, also show loading information
  // at the bottom so when clicking on "load more" we also have some feedback
  search.on("render", () => {
    const container = document.getElementById("docsearch-loading-indicator");
    if (search.status === "stalled") {
      container.innerText = "Loading search results...";
      container.style.marginTop = "0.4rem";
    } else {
      container.innerText = "";
      container.style.marginTop = "0";
    }
  });
});
