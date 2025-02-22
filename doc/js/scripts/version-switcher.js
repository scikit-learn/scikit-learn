/**
 * Adds the link to available documentation page as the last entry in the version
 * switcher dropdown. Since other entries in the dropdown are also added dynamically,
 * we only add the link when the user clicks on some version switcher button to make
 * sure that this entry is the last one.
 */

function addVersionSwitcherAvailDocsLink() {
  var availDocsLinkAdded = false;

  // There can be multiple version switcher buttons because there is at least one for
  // laptop size and one for mobile size (in the sidebar)
  document
    .querySelectorAll(".version-switcher__button")
    .forEach(function (btn) {
      btn.addEventListener("click", function () {
        if (!availDocsLinkAdded) {
          // All version switcher dropdowns are updated once any button is clicked
          document
            .querySelectorAll(".version-switcher__menu")
            .forEach(function (menu) {
              var availDocsLink = document.createElement("a");
              availDocsLink.setAttribute(
                "href",
                "https://scikit-learn.org/dev/versions.html"
              );
              availDocsLink.innerHTML = "More";
              // We use the same class as the last entry to be safe
              availDocsLink.className = menu.lastChild.className;
              availDocsLink.classList.add("sk-avail-docs-link");
              menu.appendChild(availDocsLink);
            });
          // Set the flag so we do not add again
          availDocsLinkAdded = true;
        }
      });
    });
}

document.addEventListener("DOMContentLoaded", addVersionSwitcherAvailDocsLink);
