# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn._config import get_config


class ReprHTMLMixin:
    @property
    def _repr_html_(self):
        # Taken from sklearn.base.BaseEstimator
        """HTML representation of estimator.
        This is redundant with the logic of `_repr_mimebundle_`. The latter
        should be favored in the long term, `_repr_html_` is only
        implemented for consumers who do not interpret `_repr_mimbundle_`.
        """
        if get_config()["display"] != "diagram":
            raise AttributeError(
                "_repr_html_ is only defined when the "
                "'display' configuration option is set to "
                "'diagram'"
            )
        return self._repr_html_inner

    def _repr_html_inner(self):
        from sklearn.utils._repr_html.params import _html_template

        return _html_template(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        from sklearn.utils._repr_html.params import _html_template

        output = {"text/plain": repr(self)}
        if get_config()["display"] == "diagram":
            output["text/html"] = _html_template(self)
        return output
