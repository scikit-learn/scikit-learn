/**
 * This script is for initializing the search table on the API index page. See
 * DataTables documentation for more information: https://datatables.net/
 */

document.addEventListener("DOMContentLoaded", function () {
  new DataTable("table.apisearch-table", {
    order: [], // Keep original order
    lengthMenu: [10, 25, 50, 100, { label: "All", value: -1 }],
    pageLength: -1, // Show all entries by default
  });
});
