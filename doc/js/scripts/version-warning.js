/**
 * Adds the version warning banner to the top of the page. This relies on the version
 * switcher JSON file used by pydata-sphinx-theme, being available at the URL specified
 * by `DOCUMENTATION_OPTIONS.theme_switcher_json_url`.
 */

async function addVersionWarningBanner() {
  try {
    const version = DOCUMENTATION_OPTIONS.VERSION;

    // Fetch the version switcher JSON file and get the current entry; see
    // https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/version-dropdown.html
    const response = await fetch(DOCUMENTATION_OPTIONS.theme_switcher_json_url);
    const versionList = await response.json();
    const currentEntry = versionList.find((item) => item.version === version);
    if (currentEntry.name.includes("stable")) {
      return;
    }

    // Get the latest stable and dev entries to be used in messages
    const stableEntry = versionList.find((item) => item.preferred);
    const devEntry = versionList.find((item) => item.version.includes("dev"));

    // The banner div, using the same class as the pydata-sphinx-theme version warning
    // banner, so that it can be styled in the same way
    var versionWarningBanner = document.createElement("div");
    versionWarningBanner.className =
      "bd-header-version-warning container-fluid init";

    // The banner has flex display, so we need to wrap the message in a div
    var bannerContainer = document.createElement("div");
    if (currentEntry.name.includes("dev")) {
      bannerContainer.innerHTML = `This is documentation for the <strong>unstable development version</strong> of scikit-learn. To use it, <a href="https://scikit-learn.org/dev/developers/advanced_installation.html#installing-nightly-builds">install the nightly build</a>. The latest stable release is <a href="${stableEntry.url}">version ${stableEntry.version}</a>.`;
    } else {
      bannerContainer.innerHTML = `This is documentation for an <strong>old release</strong> of scikit-learn (version ${version}). Try the <a href="${stableEntry.url}">latest stable release</a> (version ${stableEntry.version}) or the <a href="${devEntry.url}">unstable development version</a>.`;
    }

    // Insert the banner into the top of the body
    versionWarningBanner.appendChild(bannerContainer);
    document.body.prepend(versionWarningBanner);

    // This is the animation for the banner, same as the following:
    // https://github.com/pydata/pydata-sphinx-theme/blob/4c7572b7c0e805a7585d5916ffbd2de5e6f90ac7/src/pydata_sphinx_theme/assets/scripts/pydata-sphinx-theme.js#L520
    const autoHeight = Math.max(
      versionWarningBanner.offsetHeight,
      3 * parseFloat(getComputedStyle(document.documentElement).fontSize)
    );
    versionWarningBanner.style.setProperty("height", 0);
    versionWarningBanner.style.setProperty("padding-top", 0);
    versionWarningBanner.style.setProperty("padding-bottom", 0);
    versionWarningBanner.classList.remove("init");
    setTimeout(() => {
      versionWarningBanner.style.setProperty("height", `${autoHeight}px`);
      setTimeout(() => {
        versionWarningBanner.style.removeProperty("height");
        versionWarningBanner.style.removeProperty("padding-top");
        versionWarningBanner.style.removeProperty("padding-bottom");
        versionWarningBanner.style.setProperty("min-height", "3rem");
      }, 320);
    }, 10);
  } catch (error) {
    console.error(error);
  }
}

document.addEventListener("DOMContentLoaded", addVersionWarningBanner);
