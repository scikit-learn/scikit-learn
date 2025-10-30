# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html

from sklearn.utils._repr_html.base import _IDCounter

_FEATURES_ID_COUNTER = _IDCounter("sk-features-id")


def _features_html(features, is_fitted_css_class=""):
    FEATURES_TABLE_TEMPLATE = """
        <div class="sk-item">
          <div class="sk-estimator {is_fitted_css_class} sk-toggleable">
            <input class="sk-toggleable__control sk-hidden--visually" id="{feature_id}"
            type="checkbox">
            <label for="{feature_id}" class="sk-toggleable__label {is_fitted_css_class}
            sk-toggleable__label-arrow">
              <div>{total_features} features</div>
              <div class="image-container" title="Copy all output features" onclick="
                event.preventDefault();
                event.stopPropagation();
                var checkboxElem = this.closest('.sk-toggleable')
                                  .querySelector('input');
                var wasChecked = checkboxElem.checked;
                checkboxElem.checked = true;
                var content = this.closest('.sk-item').querySelector('tbody')
                              .innerText.trim();
                if (!wasChecked) checkboxElem.checked = false;
                copyRowsToClipboard(content, this);
                ">
                  <i class="copy-paste-icon"></i>
              </div>
            </label>
            <div class="sk-toggleable__content {is_fitted_css_class}">
              <div class="estimator-table">
                <table class="parameters-table">
                  <tbody>
                    {rows}
                  </tbody>
                </table>
              </div>
            </div>
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
        feature_id=_FEATURES_ID_COUNTER.get_id(),
        rows="".join(rows),
    )
