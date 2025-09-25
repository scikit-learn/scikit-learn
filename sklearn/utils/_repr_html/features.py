# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html


def _features_html(features):
    FEATURES_TABLE_TEMPLATE = """
        <div class="features">
          <details>
            <summary>{total_features} features </summary>
            <div class="features-container">
              <i class="copy-paste-icon">
              </i>
              <table = "features-table">
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
        escaped_feature = html.escape(feature)
        rows.append(FEATURES_ROW_TEMPLATE.format(feature=escaped_feature))

    return FEATURES_TABLE_TEMPLATE.format(
        total_features=total_features, rows="".join(rows)
    )
