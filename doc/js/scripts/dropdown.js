/**
 * This script is used to add the functionality of collapsing/expanding all dropdowns
 * on the page to the sphinx-design dropdowns. This is because some browsers cannot
 * search into collapsed <details> (such as Firefox).
 *
 * The reason why the buttons are added to the page with JS (dynamic) instead of with
 * sphinx (static) is that the button will not work without JS activated, so we do not
 * want them to show up in that case.
 */

function addToggleAllButtons() {
  // Get all sphinx-design dropdowns
  const allDropdowns = document.querySelectorAll("details.sd-dropdown");

  function collapseAll() {
    // Function to collapse all dropdowns on the page
    console.log("[SK] Collapsing all dropdowns...");
    allDropdowns.forEach((dropdown) => {
      dropdown.removeAttribute("open");
    });
  }

  function expandAll() {
    // Function to expand all dropdowns on the page
    console.log("[SK] Expanding all dropdowns...");
    allDropdowns.forEach((dropdown) => {
      dropdown.setAttribute("open", "");
    });
  }

  const buttonConfigs = new Map([
    ["up", { desc: "Collapse", action: collapseAll }],
    ["down", { desc: "Expand", action: expandAll }],
  ]);

  allDropdowns.forEach((dropdown) => {
    // Get the summary element of the dropdown, where we will place the buttons
    const summaryTitle = dropdown.querySelector("summary.sd-summary-title");
    for (const [direction, config] of buttonConfigs) {
      // Button with icon inside
      var newButton = document.createElement("button");
      var newIcon = document.createElement("i");
      newIcon.classList.add("fa-solid", `fa-angles-${direction}`);
      newButton.appendChild(newIcon);
      // Class for styling; `sd-summary-up/down` is implemented by sphinx-design;
      // `sk-toggle-all` is implemented by us
      newButton.classList.add(`sd-summary-${direction}`, `sk-toggle-all`);
      // Bootstrap tooltip configurations
      newButton.setAttribute("data-bs-toggle", "tooltip");
      newButton.setAttribute("data-bs-placement", "top");
      newButton.setAttribute("data-bs-offset", "0,10");
      newButton.setAttribute("data-bs-title", `${config.desc} all dropdowns`);
      // Assign the collapse/expand action to the button
      newButton.onclick = config.action;
      // Append the button to the summary element
      summaryTitle.appendChild(newButton);
    }
  });
}

document.addEventListener("DOMContentLoaded", addToggleAllButtons);
