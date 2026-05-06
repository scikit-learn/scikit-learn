The scikit-learn machine learning cheat sheet was originally created by Andreas Mueller:
https://peekaboo-vision.blogspot.de/2013/01/machine-learning-cheat-sheet-for-scikit.html

The current version of the chart is located at `doc/images/ml_map.svg` in SVG+XML
format, created using [draw.io](https://draw.io/). To edit the chart, open the file in
draw.io, make changes, and save. This should update the chart in-place. Another option
would be to re-export the chart as SVG and replace the existing file. The options used
for exporting the chart are:

- Zoom: 100%
- Border width: 15
- Size: Diagram
- Transparent Background: False
- Appearance: Light

Note that estimators nodes are clickable and should go to the estimator
documentation. After updating or re-exporting the SVG with draw.io, the links
may be prefixed with e.g. `https://app.diagrams.net/`. Remember to check and
remove them, for instance by replacing all occurrences of
`https://app.diagrams.net/./` with `./` with the following command:

.. prompt:: bash

  perl -pi -e 's@https://app.diagrams.net/\./@./@g' doc/images/ml_map.svg
