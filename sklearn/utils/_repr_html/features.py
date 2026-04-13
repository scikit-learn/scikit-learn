# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html


def _features_html(features, is_fitted_css_class=""):
    """Generate HTML representation of feature names.

    Creates a collapsible HTML details element containing a table of feature
    names with a summary line showing the total count. Includes a copy-to-clipboard
    button for all feature names.
    """
    FEATURES_TABLE_TEMPLATE = """
        <div class="features {is_fitted_css_class}">
          <details>
            <summary>
              <div class="arrow"></div>
              <div>{total_features_line}</div>
              <div class="image-container" title="Copy all output features">
                <i class="copy-paste-icon"
                  onclick="
                  event.stopPropagation();
                  event.preventDefault();
                  copyFeatureNamesToClipboard(this);
                  "
                >
                </i>
              </div>
            </summary>
            <div class="features-container">
                <table class="features-table">
                  <tbody>
                    {rows}
                  </tbody>
                </table>
            </div>
          </details>
        </div>
    """

    FEATURES_ROW_TEMPLATE = """
        <tr>
          <td>{feature}</td>
        </tr>

    """
    total_features = len(features)
    total_features_line = (
        f"{total_features} {'feature' if total_features == 1 else 'features'}"
    )

    rows = [
        FEATURES_ROW_TEMPLATE.format(feature=html.escape(feature))
        for feature in features
    ]
    return FEATURES_TABLE_TEMPLATE.format(
        total_features_line=total_features_line,
        is_fitted_css_class=html.escape(is_fitted_css_class),
        rows="".join(rows),
    )
