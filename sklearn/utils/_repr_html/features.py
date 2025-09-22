# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


def _features_html(features):
    FEATURES_TABLE_TEMPLATE = """
        <div class="features">
          <details>
            <summary>{total_features} features
            </summary>
            <div class="features-container">
              <table>
                <tbody>
                  {rows}
                </tbody>
              </table>
            </div>
          </details>
        </div>
        <br>
    """

    FEATURES_ROW_TEMPLATE = """
        <tr><td>{feature}</td></tr>
    """
    total_features = len(features)
    rows = []
    for feature in features:
        rows.append(FEATURES_ROW_TEMPLATE.format(feature=feature))
    return FEATURES_TABLE_TEMPLATE.format(
        total_features=total_features, rows="".join(rows)
    )
