// Function to create permalink into <details> elements to be able to link them
// The assumption is that such a block will be defined as follows:
//     <details id="summary-anchor">
//     <summary class="btn btn-light">
//     Some title
//     <span class="tooltiptext">Click for more details</span>
//     <a class="headerlink" href="#summary-anchor" title="Permalink to this heading">¶</a>
//     </summary>
//     <div class="card">
//     Some details
//     </div>
//     </details>
// We seek to replace `#summary-anchor` with a unique identifier based on the
// summary text.
// This syntax is defined in `doc/conf.py` in the `rst_prolog` variable.
function updateIdAndHrefBasedOnSummaryText() {
    var allDetailsElements = document.querySelectorAll('details');
    // Counter to store the duplicated summary text to add it as a suffix in the
    // anchor ID
    var anchorIDCounters = {};

    allDetailsElements.forEach(function (detailsElement) {
        // Get the <summary> element within the current <details>
        var summaryElement = detailsElement.querySelector('summary');

        // The ID uses the first line, lowercased, and spaces replaced with dashes
        var anchorID = summaryElement.textContent.trim().split("\n")[0].replace(/\s+/g, '-').toLowerCase();

        // Suffix the anchor ID with a counter if it already exists
        if (anchorIDCounters[anchorID]) {
            anchorIDCounters[anchorID] += 1;
            anchorID = anchorID + '-' + anchorIDCounters[anchorID];
        } else {
            anchorIDCounters[anchorID] = 1;
        }

        detailsElement.setAttribute('id', anchorID);

        var anchorElement = summaryElement.querySelector('a.headerlink');
        anchorElement.setAttribute('href', '#' + anchorID);
    });
}

// Add an event listener to execute the function when the page is loaded
document.addEventListener('DOMContentLoaded', function () {
    updateIdAndHrefBasedOnSummaryText();
});
