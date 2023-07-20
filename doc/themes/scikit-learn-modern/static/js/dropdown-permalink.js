// JavaScript for dynamically setting the id attribute
document.addEventListener("DOMContentLoaded", function () {
    const detailsElements = document.querySelectorAll("details");
    detailsElements.forEach((detailsElement) => {
        const summaryElement = detailsElement.querySelector("summary");
        if (summaryElement) {
            const strongElement = summaryElement.querySelector("strong");
            if (strongElement) {
                const summaryText = strongElement.textContent.trim();
                const id = summaryText.replace(/\s+/g, "-").toLowerCase();
                detailsElement.setAttribute("id", id);
            }
        }
    });

    const url = new URL(window.location.href);
    const fragment = url.hash.slice(1); // Get the URL fragment (without the '#' symbol)

    if (fragment) {
        const targetDetailsElement = document.getElementById(fragment);
        if (targetDetailsElement) {
            targetDetailsElement.open = true;
        }
    }
});
  