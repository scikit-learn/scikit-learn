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

                // Create the <a> element and set its attributes
                const anchorElement = document.createElement("a");
                anchorElement.setAttribute("class", "headerlink");
                anchorElement.setAttribute("href", `#${id}`);
                anchorElement.setAttribute("title", summaryText);
                anchorElement.textContent = "¶"

                // Insert the <a> element after the text inside the <strong> tag
                strongElement.appendChild(anchorElement);

                // Add event listener to the anchor element to toggle the 'open' attribute
                anchorElement.addEventListener("click", function (event) {
                    // event.preventDefault();
                    detailsElement.open = true;
                });
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
  