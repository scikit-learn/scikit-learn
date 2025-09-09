:html_theme.sidebar_secondary.remove:

.. _ml_map:

Choosing the right estimator
============================

Often the hardest part of solving a machine learning problem can be finding the right
estimator for the job. Different estimators are better suited for different types of
data and different problems.

The flowchart below is designed to give users a bit of a rough guide on how to approach
problems with regard to which estimators to try on your data. Click on any estimator in
the chart below to see its documentation. The **Try next** orange arrows are to be read as
"if this estimator does not achieve the desired outcome, then follow the arrow and try
the next one". Use scroll wheel to zoom in and out, and click and drag to pan around.
You can also download the chart: :download:`ml_map.svg <images/ml_map.svg>`.

.. raw:: html

  <style>
    #sk-ml-map {
      height: 80vh;
      margin: 1.5rem 0;
    }

    #sk-ml-map svg {
      height: 100%;
      width: 100%;
      border: 2px solid var(--pst-color-border);
      border-radius: 0.5rem;
    }

    html[data-theme="dark"] #sk-ml-map svg {
      filter: invert(90%) hue-rotate(180deg);
    }
  </style>

  <script src="_static/scripts/vendor/svg-pan-zoom.min.js"></script>
  <script>
    document.addEventListener("DOMContentLoaded", function () {
      const beforePan = function (oldPan, newPan) {
        const gutterWidth = 100, gutterHeight = 100;
        const sizes = this.getSizes();

        // Compute pan limits
        const leftLimit = -((sizes.viewBox.x + sizes.viewBox.width) * sizes.realZoom) + gutterWidth;
        const rightLimit = sizes.width - gutterWidth - (sizes.viewBox.x * sizes.realZoom);
        const topLimit = -((sizes.viewBox.y + sizes.viewBox.height) * sizes.realZoom) + gutterHeight;
        const bottomLimit = sizes.height - gutterHeight - (sizes.viewBox.y * sizes.realZoom);

        return {
          x: Math.max(leftLimit, Math.min(rightLimit, newPan.x)),
          y: Math.max(topLimit, Math.min(bottomLimit, newPan.y))
        };
      };

      // Limit the pan
      svgPanZoom("#sk-ml-map svg", {
        zoomEnabled: true,
        controlIconsEnabled: true,
        fit: 1,
        center: 1,
        beforePan: beforePan,
      });
    });
  </script>

  <div id="sk-ml-map">

.. raw:: html
  :file: images/ml_map.svg

.. raw:: html

  </div>
