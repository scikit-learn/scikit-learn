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

function addVersionWarningNightlyLink() {
  var version = DOCUMENTATION_OPTIONS.VERSION || "";
  var isDevelopmentVersion = version.includes("dev");

  if (!isDevelopmentVersion) {
    return;
  }

  var banner = document.querySelector("#bd-header-version-warning");

  if (!banner) {
    return;
  }

  function addNightlyLink() {
    if (
      banner.classList.contains("d-none") ||
      banner.querySelector(".sk-nightly-install-link")
    ) {
      return;
    }

    var stableLink = banner.querySelector(".pst-button-link-to-stable-version");

    if (!stableLink) {
      return;
    }

    var nightlyLink = document.createElement("a");
    nightlyLink.href =
      "https://scikit-learn.org/dev/install.html#installing-nightly-builds";
    nightlyLink.innerText = "Install nightly build";
    nightlyLink.className = stableLink.className;
    nightlyLink.classList.add("sk-nightly-install-link");

    stableLink.insertAdjacentElement("afterend", nightlyLink);
  }

  addNightlyLink();

  var observer = new MutationObserver(addNightlyLink);
  observer.observe(banner, {
    attributes: true,
    attributeFilter: ["class"],
    childList: true,
    subtree: true,
  });
}

document.addEventListener("DOMContentLoaded", addVersionWarningNightlyLink);
