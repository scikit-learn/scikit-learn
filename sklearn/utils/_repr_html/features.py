# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


def _features_html(features):
    FEATURES_TABLE_TEMPLATE = """
        <div class="features">
          <details>
            <summary class="features-title">{total_features} features</summary>
              <ul>
               {rows}
              </ul>
          </details>
        </div>
        <br>
    """

    FEATURES_ROW_TEMPLATE = """
        <li>{feature}</li>
    """
    total_features = len(features)
    rows = []
    for feature in features:
        rows.append(FEATURES_ROW_TEMPLATE.format(feature=feature))
    return FEATURES_TABLE_TEMPLATE.format(
        total_features=total_features, rows="".join(rows)
    )
