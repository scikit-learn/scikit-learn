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

for (const key of ["rc", "finale", "bf"]) {
  document.addEventListener("DOMContentLoaded", () => {
    const versionFullTargets = document.querySelectorAll(
      `.${key}-version-full`
    );
    const versionMainTargets = document.querySelectorAll(
      `.${key}-version-main`
    );
    const warningContainer = document.getElementById(
      `${key}-invalid-release-version-warning`
    );

    const defaultVersionFull = versionFullTargets[0].textContent;
    const [defaultMajor, defaultMinor, _] = validateVersion(
      defaultVersionFull,
      key
    );

    document
      .getElementById(`${key}-version-to-release`)
      .addEventListener("input", (e) => {
        const version = e.target.value || defaultVersionFull;

        let currentVersionFull = defaultVersionFull;
        let currentMajor = defaultMajor;
        let currentMinor = defaultMinor;
        let warningMessage = "";
        let color = "";
        let cursor = "";
        let title = "";

        // Validate the version string
        const validationResult = validateVersion(version, key);
        if (typeof validationResult !== "string") {
          // Valid version; major/minor/micro numbers returned
          currentVersionFull = version;
          currentMajor = validationResult[0];
          currentMinor = validationResult[1];
        } else {
          // Invalid version; invalid reason returned
          warningMessage = validationResult;
          color = "var(--pst-color-danger)";
          cursor = "help";
          title = "Invalid version; fallback to default";
        }

        // Set warning message
        warningContainer.textContent = warningMessage;

        // Replace target texts
        versionFullTargets.forEach((target) => {
          target.textContent = currentVersionFull;
          target.style.color = color;
          target.style.cursor = cursor;
          target.title = title;
        });
        versionMainTargets.forEach((target) => {
          target.textContent = `${currentMajor}.${currentMinor}`;
          target.style.color = color;
          target.style.cursor = cursor;
          target.title = title;
        });
      });
  });
}
