// @ts-check

// Extra JS capability for selected tabs to be synced
// The selection is stored in local storage so that it persists across page loads.

/**
 * @type {Record<string, HTMLElement[]>}
 */
let sd_id_to_elements = {};
const storageKeyPrefix = "sphinx-design-tab-id-";

/**
 * Create a key for a tab element.
 * @param {HTMLElement} el - The tab element.
 * @returns {[string, string, string] | null} - The key.
 *
 */
function create_key(el) {
  let syncId = el.getAttribute("data-sync-id");
  let syncGroup = el.getAttribute("data-sync-group");
  if (!syncId || !syncGroup) return null;
  return [syncGroup, syncId, syncGroup + "--" + syncId];
}

/**
 * Initialize the tab selection.
 *
 */
function ready() {
  // Find all tabs with sync data

  /** @type {string[]} */
  let groups = [];

  document.querySelectorAll(".sd-tab-label").forEach((label) => {
    if (label instanceof HTMLElement) {
      let data = create_key(label);
      if (data) {
        let [group, id, key] = data;

        // add click event listener
        // @ts-ignore
        label.onclick = onSDLabelClick;

        // store map of key to elements
        if (!sd_id_to_elements[key]) {
          sd_id_to_elements[key] = [];
        }
        sd_id_to_elements[key].push(label);

        if (groups.indexOf(group) === -1) {
          groups.push(group);
          // Check if a specific tab has been selected via URL parameter
          const tabParam = new URLSearchParams(window.location.search).get(
            group
          );
          if (tabParam) {
            console.log(
              "sphinx-design: Selecting tab id for group '" +
                group +
                "' from URL parameter: " +
                tabParam
            );
            window.sessionStorage.setItem(storageKeyPrefix + group, tabParam);
          }
        }

        // Check is a specific tab has been selected previously
        let previousId = window.sessionStorage.getItem(
          storageKeyPrefix + group
        );
        if (previousId === id) {
          // console.log(
          //   "sphinx-design: Selecting tab from session storage: " + id
          // );
          // @ts-ignore
          label.previousElementSibling.checked = true;
        }
      }
    }
  });
}

/**
 *  Activate other tabs with the same sync id.
 *
 * @this {HTMLElement} - The element that was clicked.
 */
function onSDLabelClick() {
  let data = create_key(this);
  if (!data) return;
  let [group, id, key] = data;
  for (const label of sd_id_to_elements[key]) {
    if (label === this) continue;
    // @ts-ignore
    label.previousElementSibling.checked = true;
  }
  window.sessionStorage.setItem(storageKeyPrefix + group, id);
}

document.addEventListener("DOMContentLoaded", ready, false);
