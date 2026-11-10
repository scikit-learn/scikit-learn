(function () {
  const observer = new MutationObserver((mutationsList) => {
    for (const mutation of mutationsList) {
      if (
        mutation.type === "attributes" &&
        mutation.attributeName === "data-theme"
      ) {
        document
          .querySelectorAll(".sk-top-container")
          .forEach((estimatorElement) => {
            const newTheme = detectTheme(estimatorElement);
            estimatorElement.classList.remove("light", "dark");
            estimatorElement.classList.add(newTheme);
          });
      }
    }
  });

  observer.observe(document.documentElement, {
    attributes: true,
    attributeFilter: ["data-theme"],
  });
})();
