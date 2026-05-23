/**
 * This script is used to add the functionality of collapsing/expanding all dropdowns
 * on the page to the sphinx-design dropdowns. This is because some browsers cannot
 * search into collapsed <details> (such as Firefox).
 *
 * The reason why the buttons are added to the page with JS (dynamic) instead of with
 * sphinx (static) is that the button will not work without JS activated, so we do not
 * want them to show up in that case.
 */

document.addEventListener("DOMContentLoaded", () => {
  // Get all sphinx-design dropdowns
  const allDropdowns = document.querySelectorAll("details.sd-dropdown");

  allDropdowns.forEach((dropdown) => {
    // Get the summary element of the dropdown, where we will place the buttons
    const summaryTitle = dropdown.querySelector("summary.sd-summary-title");

    // The state marker with the toggle all icon inside
    const newStateMarker = document.createElement("span");
    const newIcon = document.createElement("i");
    newIcon.classList.add("fa-solid", "fa-angles-right");
    newStateMarker.appendChild(newIcon);

    // Classes for styling; `sd-summary-state-marker` and `sd-summary-chevron-right` are
    // implemented by sphinx-design; `sk-toggle-all` is implemented by us
    newStateMarker.classList.add(
      "sd-summary-state-marker",
      "sd-summary-chevron-right",
      "sk-toggle-all"
    );

    // Bootstrap tooltip configurations
    newStateMarker.setAttribute("data-bs-toggle", "tooltip");
    newStateMarker.setAttribute("data-bs-placement", "top");
    newStateMarker.setAttribute("data-bs-offset", "0,10");
    newStateMarker.setAttribute("data-bs-title", "Toggle all dropdowns");
    // Enable the tooltip
    new bootstrap.Tooltip(newStateMarker);

    // Assign the collapse/expand action to the state marker
    newStateMarker.addEventListener("click", () => {
      if (dropdown.open) {
        console.log("[SK] Collapsing all dropdowns...");
        allDropdowns.forEach((node) => {
          if (node !== dropdown) {
            node.removeAttribute("open");
          }
        });
      } else {
        console.log("[SK] Expanding all dropdowns...");
        allDropdowns.forEach((node) => {
          if (node !== dropdown) {
            node.setAttribute("open", "");
          }
        });
      }
    });

    // Append the state marker to the summary element
    summaryTitle.insertBefore(newStateMarker, summaryTitle.lastElementChild);
  });
});
