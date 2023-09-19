// Function to create permalink into <details> elements to be able to link them
// The assumption is that such a block will be defined as follows:
//     <details>
//     <summary class="btn btn-light" id="summary-anchor">
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
function updateHrefBasedOnSummaryText() {
    var allSummaryElements = document.querySelectorAll('details summary');

    allSummaryElements.forEach(function (summaryElement) {
        // The ID uses the first line, lower the case and replace spaces with
        // dashes
        var anchorID = summaryElement.textContent.trim().split("\n")[0].replace(/\s+/g, '-').toLowerCase();
        summaryElement.setAttribute('id', anchorID);

        var hrefID = '#' + anchorID;
        var anchorElement = summaryElement.querySelector('a.headerlink');
        if (anchorElement) {
            anchorElement.setAttribute('href', hrefID);
        }
    });
}

// Add an event listener to execute the function when the page is loaded
document.addEventListener('DOMContentLoaded', function () {
    updateHrefBasedOnSummaryText();
});
