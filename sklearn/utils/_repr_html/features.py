# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html

_MAX_DISPLAY_FEATURES = 100


def _features_html(features, is_fitted_css_class=""):
    """Generate HTML representation of feature names.

    Creates a collapsible HTML details element containing a table of feature
    names with a summary line showing the total count. Includes a copy-to-clipboard
    button for all feature names.

    Only the first ``_MAX_DISPLAY_FEATURES`` features are rendered as table
    rows to keep the HTML lightweight.
    """
    FEATURES_TABLE_TEMPLATE = """
        <div class="features {is_fitted_css_class}">
          <details>
            <summary>
              <div class="arrow"></div>
              <div>{total_features_line}</div>
              <div class="image-container">
                <button type="button" class="copy-paste-icon"
                  title="{copy_features_label}"
                  aria-label="{copy_features_label}"
                  onclick="
                  event.stopPropagation();
                  event.preventDefault();
                  copyFeatureNamesToClipboard(this);
                  "
                >
                </button>
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
    copy_features_label = f"Copy output features (max {_MAX_DISPLAY_FEATURES})"
    total_features = len(features)
    display_features = features[:_MAX_DISPLAY_FEATURES]

    if total_features > _MAX_DISPLAY_FEATURES:
        total_features_line = f"{_MAX_DISPLAY_FEATURES} of {total_features:,} features"
    else:
        total_features_line = (
            f"{total_features} {'feature' if total_features == 1 else 'features'}"
        )

    rows = [
        FEATURES_ROW_TEMPLATE.format(feature=html.escape(feature))
        for feature in display_features
    ]

    return FEATURES_TABLE_TEMPLATE.format(
        total_features_line=total_features_line,
        is_fitted_css_class=html.escape(is_fitted_css_class),
        copy_features_label=copy_features_label,
        rows="".join(rows),
    )
