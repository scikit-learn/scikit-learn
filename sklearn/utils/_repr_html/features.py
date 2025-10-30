# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html


def _features_html(features, is_fitted_css_class=""):
    FEATURES_TABLE_TEMPLATE = """
        <div class="sk-estimator {is_fitted_css_class}">
          <details>
            <summary>
              <div>{total_features} features</div>
              <div class="image-container"
                   title="Copy all output features">
                <i class="copy-paste-icon"
                  onclick="
                  event.stopPropagation();
                  event.preventDefault();

                  var detailsElem = this.closest('.features')
                  .querySelector('details');
                  var wasOpen = detailsElem.open;
                  detailsElem.open = true;
                  var content = this.closest('.features').querySelector('tbody')
                                .innerText.trim();
                  if (!wasOpen) detailsElem.open = false;
                  copyRowsToClipboard(content, this);
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
          </div>
    """

    FEATURES_ROW_TEMPLATE = """
        <tr>
            <td>{feature}</td>
        </tr>

    """
    total_features = len(features)
    rows = []
    for feature in features:
        escaped_feature = html.escape(feature)
        rows.append(FEATURES_ROW_TEMPLATE.format(feature=escaped_feature))

    return FEATURES_TABLE_TEMPLATE.format(
        total_features=total_features,
        is_fitted_css_class=is_fitted_css_class,
        rows="".join(rows),
    )
