/**
 * Aligns the card texts in the landing page. This should be loaded only for the
 * landing page.
 */

function alignLandingCards() {
  var skCardCollections = { title: {}, desc: {} };

  // Get the groups to align
  document.querySelectorAll(".card").forEach(function (card) {
    var groupId = card.getAttribute("sk-align-group");
    card.querySelectorAll(".sk-vert-align").forEach(function (item) {
      var alignName = item.getAttribute("sk-align-name");
      if (!skCardCollections[alignName][groupId]) {
        skCardCollections[alignName][groupId] = [];
      }
      skCardCollections[alignName][groupId].push(item);
    });
  });

  function align(group) {
    // Align items in each group by setting them to same height
    var maxHeight = 0;

    if (window.innerWidth >= 720) {
      // Get the maximum auto height as the final height to set
      group.forEach(function (item) {
        item.style.setProperty("height", "auto");
        maxHeight = Math.max(maxHeight, item.offsetHeight);
      });
      group.forEach(function (item) {
        item.style.setProperty("height", `${maxHeight}px`);
      });
    } else {
      // We do not need to align the cards if the screen is sufficiently small
      // such that each row contains only one card; set height back to auto
      group.forEach(function (item) {
        item.style.setProperty("height", "auto");
      });
    }
  }

  function alignAll() {
    // Align all groups from the collection
    for (var alignName in skCardCollections) {
      for (var groupId in skCardCollections[alignName]) {
        align(skCardCollections[alignName][groupId]);
      }
    }
  }

  // Align on load and realign on resize
  window.onload = alignAll;
  window.addEventListener("resize", alignAll);
}

document.addEventListener("DOMContentLoaded", alignLandingCards);
