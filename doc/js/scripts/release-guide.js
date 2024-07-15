/**
 * This script is for substituting the versions in multiple places in the release guide
 * according to input, on the `developers/maintainer` page.
 */

// Parsing version strings with regular expressions according to PyPA
// https://packaging.python.org/en/latest/specifications/version-specifiers/#version-specifiers-regex
const VERSION_PATTERN = new RegExp(
  `^\\s*v?(?:(?<epoch>\\d+)!)?(?<release>\\d+(?:\\.\\d+)*)(?<pre>[-_\\.]?(?<pre_l>alpha|a|beta|b|preview|pre|c|rc)[-_\\.]?(?<pre_n>\\d+)?)?(?<post>(?:-(?<post_n1>\\d+))|(?:[-_\\.]?(?<post_l>post|rev|r)[-_\\.]?(?<post_n2>\\d+)?))?(?<dev>[-_\\.]?(?<dev_l>dev)[-_\\.]?(?<dev_n>\\d+)?)?(?:\\+(?<local>[a-z0-9]+(?:[-_\\.]?[a-z0-9]+)*))?\\s*$`,
  "i"
);

function validateVersion(version, key) {
  const match = version.match(VERSION_PATTERN);
  if (!match) {
    return "Invalid version format.";
  }

  if (match.groups.epoch || match.groups.dev || match.groups.local) {
    return "Version should not contain epoch, dev, or local parts.";
  }

  if (key == "rc") {
    if (
      match.groups.pre_l != "rc" ||
      !match.groups.pre_n ||
      Number(match.groups.pre_n) == 0
    ) {
      return "Release candidate should have a non-zero release candidate number.";
    }
  } else {
    if (match.groups.pre) {
      return "Non-RC releases should not have the release candidate part.";
    }
  }

  try {
    // Release number should follow major.minor.micro format
    const rel = match.groups.release.split(".").map(Number);
    if (rel.length == 3) {
      if (key == "bf") {
        if (rel[2] >= 1) {
          return rel;
        }
        return "Bug-fix releases should have non-zero micro version.";
      } else {
        if (rel[2] == 0) {
          return rel;
        }
        return "Non-bug-fix releases should have zero micro version.";
      }
    }
    return "Release should be in major.minor.micro format.";
  } catch {
    return "Release should be in major.minor.micro format.";
  }
}

document.addEventListener("DOMContentLoaded", () => {
  fetch(DOCUMENTATION_OPTIONS.theme_switcher_json_url).then((response) => {
    response.json().then((data) => {
      for (const key of ["rc", "finale", "bf"]) {
        const versionFullTargets = document.querySelectorAll(
          `.${key}-version-full`
        );
        const versionMainTargets = document.querySelectorAll(
          `.${key}-version-main`
        );
        const previousTagTargets = document.querySelectorAll(
          `.${key}-previous-tag`
        );
        const warningContainer = document.getElementById(
          `${key}-invalid-release-version-warning`
        );
        const versionFullInput = document.getElementById(
          `${key}-version-full-input`
        );
        const previousTagInput = document.getElementById(
          `${key}-previous-tag-input`
        );

        // Update target texts, styles, and the warning message
        const updateContents = (
          version,
          major,
          minor,
          warningMessage = "",
          color = "",
          cursor = "",
          title = ""
        ) => {
          versionFullTargets.forEach((target) => {
            target.textContent = version;
            target.style.color = color;
            target.style.cursor = cursor;
            target.title = title;
          });
          versionMainTargets.forEach((target) => {
            target.textContent = `${major}.${minor}`;
            target.style.color = color;
            target.style.cursor = cursor;
            target.title = title;
          });
          warningContainer.textContent = warningMessage;
        };

        // Update the previous tag text
        const updatePreviousTag = (previousTag) => {
          previousTagTargets.forEach((target) => {
            target.textContent = previousTag;
          });
        };

        // Static fallback as given in the rst
        let defaultVersionFull = versionFullTargets[0].textContent;
        const defaultVersionFullParts = defaultVersionFull.split(".");
        let defaultMajor = Number(defaultVersionFullParts[0]);
        let defaultMinor = Number(defaultVersionFullParts[1]);
        let defaultPreviousTag =
          key == "rc" ? null : previousTagTargets[0].textContent;

        // Try to get a more reasonable default
        if (key == "bf") {
          // Bug-fix release; should be micro+1 from the latest stable
          const stableVersionIndex = data.findIndex((v) => v.preferred);
          const stableVersionParts =
            data[stableVersionIndex].version.split(".");
          defaultMajor = Number(stableVersionParts[0]);
          defaultMinor = Number(stableVersionParts[1]);
          defaultVersionFull = `${defaultMajor}.${defaultMinor}.${
            Number(stableVersionParts[2]) + 1
          }`;
          defaultPreviousTag = data[stableVersionIndex - 1].version;
        } else if (key == "rc") {
          const rcVersionIndex = data.findIndex((v) =>
            v.version.includes(".0rc")
          );
          if (rcVersionIndex != -1) {
            // Not the first release candidate; should be rc+1 from the current rc
            const rcVersionParts = data[rcVersionIndex].version.split(".");
            defaultMajor = Number(rcVersionParts[0]);
            defaultMinor = Number(rcVersionParts[1]);
            const rcVersionPart = rcVersionParts[2];
            defaultVersionFull = `${defaultMajor}.${defaultMinor}.0rc${
              Number(rcVersionPart[rcVersionPart.length - 1]) + 1
            }`;
          } else {
            // First release candidate; should be minor+1 from the latest stable
            const stableVersionIndex = data.findIndex((v) => v.preferred);
            const stableVersionParts =
              data[stableVersionIndex].version.split(".");
            defaultMajor = Number(stableVersionParts[0]);
            defaultMinor = Number(stableVersionParts[1]) + 1;
            defaultVersionFull = `${defaultMajor}.${defaultMinor}.0rc1`;
          }
        } else {
          // Finale major/minor release; should be minor+1 from the latest stable
          const stableVersionIndex = data.findIndex((v) => v.preferred);
          const stableVersionParts =
            data[stableVersionIndex].version.split(".");
          defaultMajor = Number(stableVersionParts[0]);
          defaultMinor = Number(stableVersionParts[1]) + 1;
          defaultVersionFull = `${defaultMajor}.${defaultMinor}.0`;
          defaultPreviousTag = data[stableVersionIndex].version;
        }

        // Update once with the updated default
        versionFullInput.placeholder = defaultVersionFull;
        updateContents(defaultVersionFull, defaultMajor, defaultMinor);
        if (key != "rc") {
          previousTagInput.placeholder = defaultPreviousTag;
          updatePreviousTag(defaultPreviousTag);
        }

        // Event listener for the full version input
        versionFullInput.addEventListener("input", (e) => {
          const version = e.target.value || defaultVersionFull;

          // Validate the version string
          const validationResult = validateVersion(version, key);
          if (typeof validationResult !== "string") {
            // Valid version; major/minor/micro numbers returned
            updateContents(version, validationResult[0], validationResult[1]);
          } else {
            // Invalid version; invalid reason returned
            updateContents(
              defaultVersionFull,
              defaultMajor,
              defaultMinor,
              validationResult,
              "var(--pst-color-danger)",
              "help",
              "Invalid version; fallback to default"
            );
          }
        });

        // Event listener for the previous tag input
        if (key != "rc") {
          previousTagInput.addEventListener("input", (e) => {
            const previousTag = e.target.value || defaultPreviousTag;
            updatePreviousTag(previousTag);
          });
        }
      }
    });
  });
});
