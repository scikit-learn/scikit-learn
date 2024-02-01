/**
 * Script for the installation page, see `install.rst`. It refers to the following
 * pytorch page: https://pytorch.org/get-started/locally/, subject to modifications.
 * This script should be included only for the install page.
 */

function setupInstallInstructions() {
  const supportedOperatingSystems = new Map([
    ["linux", "linux"],
    ["mac", "macos"],
    ["win", "windows"],
  ]);

  /* Determine the default selection of OS */

  function getPlatformOS() {
    const platform = navigator.platform.toLowerCase();
    for (var [
      navPlatformSubstring,
      os,
    ] of supportedOperatingSystems.entries()) {
      if (platform.indexOf(navPlatformSubstring) !== -1) {
        return os;
      }
    }
    return supportedOperatingSystems.values().next().value; // just return something
  }

  function getDefaultSelectedOS() {
    const anchor = location.hash;
    const ANCHOR_REGEX = /^#[^ ]+$/;

    // Look for anchor in the href
    if (!ANCHOR_REGEX.test(anchor)) {
      return getPlatformOS();
    }

    // Look for anchor with OS in the first portion
    const testOS = anchor.slice(1).split("-")[0];
    for (var [
      navPlatformSubstring,
      os,
    ] of supportedOperatingSystems.entries()) {
      if (testOS.indexOf(navPlatformSubstring) !== -1) {
        return os;
      }
    }
    return getPlatformOS();
  }

  const selectedOptions = {
    os: getDefaultSelectedOS(),
    packager: "pip",
    virtualenv: "venv",
  };

  const osOptions = document.querySelectorAll("#osRow .option");
  const packagerOptions = document.querySelectorAll("#packagerRow .option");
  const virtualenvOptions = document.querySelectorAll("#virtualenvRow .option");

  /* Update the install instructions based on the selected options */

  const instructionBlock = document.getElementById("skInstallInstructions");

  const pipMapping = {
    windows: "pip",
    macos: "pip",
    linux: "pip3",
  };
  const pythonMapping = {
    windows: "python",
    macos: "python",
    linux: "python3",
  };

  function validateOptions(os, packager, virtualenv) {
    if (packager === "conda") {
      // strike out the venv option
      virtualenvOptions.forEach(function (option) {
        option.style.textDecoration =
          option.id === "venv" ? "line-through" : "none";
      });
      if (virtualenv === "venv") {
        return `Conda is not compatible with venv. Please select a different option.`;
      } else {
        return false;
      }
    } else {
      // unstrike all options in virtualenv
      virtualenvOptions.forEach(function (option) {
        option.style.textDecoration = "none";
      });
      return false;
    }
  }

  function getUpdatedInstruction() {
    const curOs = selectedOptions["os"];
    const curPackager = selectedOptions["packager"];
    const curVirtualenv = selectedOptions["virtualenv"];

    const validationError = validateOptions(curOs, curPackager, curVirtualenv);
    if (validationError) {
      return `<div class="admonition error"><p class="admonition-title">Error</p><p>${validationError}</p></div>`;
    }
    var instruction = "";

    // The opening paragraph
    instruction += `<p>`;
    switch (curPackager) {
      case "pip":
        switch (curOs) {
          case "windows":
            instruction += `Install the 64-bit version of Python 3, for instance from <a href="https://www.python.org/">https://www.python.org</a>.`;
            break;
          case "macos":
            instruction += `Install Python 3 using <a href="https://brew.sh/">homebrew</a> (<code>brew install python</code>) or by manually installing the package from <a href="https://www.python.org">https://www.python.org</a>.`;
            break;
          case "linux":
            instruction += `Install python3 and python3-pip using the package manager of the Linux Distribution.`;
            break;
          default:
            break;
        }
        break;
      case "conda":
        instruction += `Install conda using the <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/">Anaconda or miniconda</a> installers or the <a href="https://github.com/conda-forge/miniforge#miniforge">miniforge</a> installers (no administrator permission required for any of those).`;
        break;
      default:
        break;
    }
    instruction += ` Then run:</p>\n`;

    // The code block for installation
    instruction += `<div class="highlight-default"><div class="highlight"><pre class="sk-install-prompt">`;
    switch (curPackager) {
      case "pip":
        switch (curVirtualenv) {
          case "venv":
            instruction += `<span>${pythonMapping[curOs]} -m venv sklearn-env</span>\n`;
            switch (curOs) {
              case "windows":
                instruction += `<span>sklearn-env\\Scripts\\activate</span>\n`;
                break;
              case "macos":
                instruction += `<span>source sklearn-env/bin/activate</span>\n`;
                break;
              case "linux":
                instruction += `<span>source sklearn-env/bin/activate</span>\n`;
                break;
              default:
                break;
            }
          // no break here because we always need the line in the default case
          default:
            instruction += `<span>${pipMapping[curOs]} install -U scikit-learn</span>`;
            break;
        }
        break;
      case "conda":
        instruction += `<span>conda create -n sklearn-env -c conda-forge scikit-learn</span>\n`;
        instruction += `<span>conda activate sklearn-env</span>`;
        break;
      default:
        break;
    }
    instruction += `</pre></div></div>`;

    // The code block for checking installation
    instruction += `<p>In order to check your installation you can use</p>`;
    instruction += `<div class="highlight-default"><div class="highlight"><pre class="sk-install-prompt">`;
    switch (curPackager) {
      case "pip":
        instruction += `<span>${pythonMapping[curOs]} -m pip show scikit-learn  # show which version and where scikit-learn is installed</span>\n`;
        instruction += `<span>${pythonMapping[curOs]} -m pip freeze             # show all packages installed in the current environment</span>\n`;
        instruction += `<span>${pythonMapping[curOs]} -c "import sklearn; sklearn.show_versions()"</span>`;
        break;
      case "conda":
        instruction += `<span>conda list scikit-learn  # show which version and where scikit-learn is installed</span>\n`;
        instruction += `<span>conda list               # show all packages installed in the current environment</span>\n`;
        instruction += `<span>python -c "import sklearn; sklearn.show_versions()"</span>`;
        break;
      default:
        break;
    }
    instruction += `</pre></div></div>`;
    return instruction;
  }

  function updateInstruction() {
    instructionBlock.innerHTML = getUpdatedInstruction();
    instructionBlock.setAttribute("data-os", selectedOptions["os"]);
    addCopyButtonToCodeCells(); // See copybutton.js from sphinx-copybutton
  }

  /* Initialize the options and their click handlers */

  function initOptions(options, category) {
    options.forEach(function (option) {
      option.addEventListener("click", function () {
        selectedOption(options, this, category);
      });
      if (option.id === selectedOptions[category]) {
        option.classList.add("selected");
      }
    });
  }

  function selectedOption(options, selection, category) {
    options.forEach(function (option) {
      if (option === selection) {
        option.classList.add("selected");
      } else {
        option.classList.remove("selected");
      }
    });
    selectedOptions[category] = selection.id;
    updateInstruction();
  }

  initOptions(osOptions, "os");
  initOptions(packagerOptions, "packager");
  initOptions(virtualenvOptions, "virtualenv");
  updateInstruction();
}

document.addEventListener("DOMContentLoaded", setupInstallInstructions);
