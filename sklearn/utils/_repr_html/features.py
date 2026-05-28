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


def _name_details_html(name_details, is_fitted_css_class=""):
    """Generate collapsible table HTML for list-like name details."""
    if len(name_details) == 0:
        return ""

    NAME_DETAILS_TABLE_TEMPLATE = """
        <div class="name-details {is_fitted_css_class}">
          <div class="name-details-container">
            <table class="name-details-table">
              <tbody>
                {rows}
              </tbody>
            </table>
          </div>
        </div>
    """

    NAME_DETAILS_ROW_TEMPLATE = """
        <tr>
          <td>{item}</td>
        </tr>

    """
    rows = [
        NAME_DETAILS_ROW_TEMPLATE.format(item=html.escape(str(item)))
        for item in name_details
    ]
    return NAME_DETAILS_TABLE_TEMPLATE.format(
        is_fitted_css_class=html.escape(is_fitted_css_class),
        rows="".join(rows),
    )
