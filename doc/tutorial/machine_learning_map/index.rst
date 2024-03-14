:html_theme.sidebar_secondary.remove:

.. _ml_map:

Choosing the right estimator
============================

Often the hardest part of solving a machine learning problem can be finding the right
estimator for the job. Different estimators are better suited for different types of
data and different problems.

The flowchart below is designed to give users a bit of a rough guide on how to approach
problems with regard to which estimators to try on your data. Click on any estimator in
the chart below to see its documentation. Use scroll wheel to zoom in and out, and click
and drag to pan around. You can also download the chart:
:download:`ml_map.svg <../../images/ml_map.svg>`.

.. raw:: html

  <style>
    .ml-map {
      height: 80vh;
      margin: 1.5rem 0;
    }

    .ml-map svg {
      height: 100%;
      width: 100%;
    }

    html[data-theme="dark"] .ml-map svg {
      filter: invert(90%) hue-rotate(180deg);
    }
  </style>

  <script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom-container@0.6.1"></script>

  <div class="ml-map" data-zoom-on-wheel data-pan-on-drag>

.. raw:: html
  :file: ../../images/ml_map.svg

.. raw:: html

  </div>
