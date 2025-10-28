"""
The `trendline_functions` module contains functions which are called by Plotly Express
when the `trendline` argument is used. Valid values for `trendline` are the names of the
functions in this module, and the value of the `trendline_options` argument to PX
functions is passed in as the first argument to these functions when called.

Note that the functions in this module are not meant to be called directly, and are
exposed as part of the public API for documentation purposes.
"""

__all__ = ["ols", "lowess", "rolling", "ewm", "expanding"]


def ols(trendline_options, x_raw, x, y, x_label, y_label, non_missing):
    """Ordinary Least Squares (OLS) trendline function

    Requires `statsmodels` to be installed.

    This trendline function causes fit results to be stored within the figure,
    accessible via the `plotly.express.get_trendline_results` function. The fit results
    are the output of the `statsmodels.api.OLS` function.

    Valid keys for the `trendline_options` dict are:

    - `add_constant` (`bool`, default `True`): if `False`, the trendline passes through
    the origin but if `True` a y-intercept is fitted.

    - `log_x` and `log_y` (`bool`, default `False`): if `True` the OLS is computed with
    respect to the base 10 logarithm of the input. Note that this means no zeros can
    be present in the input.
    """
    import numpy as np

    valid_options = ["add_constant", "log_x", "log_y"]
    for k in trendline_options.keys():
        if k not in valid_options:
            raise ValueError(
                "OLS trendline_options keys must be one of [%s] but got '%s'"
                % (", ".join(valid_options), k)
            )

    import statsmodels.api as sm

    add_constant = trendline_options.get("add_constant", True)
    log_x = trendline_options.get("log_x", False)
    log_y = trendline_options.get("log_y", False)

    if log_y:
        if np.any(y <= 0):
            raise ValueError(
                "Can't do OLS trendline with `log_y=True` when `y` contains non-positive values."
            )
        y = np.log10(y)
        y_label = "log10(%s)" % y_label
    if log_x:
        if np.any(x <= 0):
            raise ValueError(
                "Can't do OLS trendline with `log_x=True` when `x`  contains non-positive values."
            )
        x = np.log10(x)
        x_label = "log10(%s)" % x_label
    if add_constant:
        x = sm.add_constant(x)
    fit_results = sm.OLS(y, x, missing="drop").fit()
    y_out = fit_results.predict()
    if log_y:
        y_out = np.power(10, y_out)
    hover_header = "<b>OLS trendline</b><br>"
    if len(fit_results.params) == 2:
        hover_header += "%s = %g * %s + %g<br>" % (
            y_label,
            fit_results.params[1],
            x_label,
            fit_results.params[0],
        )
    elif not add_constant:
        hover_header += "%s = %g * %s<br>" % (y_label, fit_results.params[0], x_label)
    else:
        hover_header += "%s = %g<br>" % (y_label, fit_results.params[0])
    hover_header += "R<sup>2</sup>=%f<br><br>" % fit_results.rsquared
    return y_out, hover_header, fit_results


def lowess(trendline_options, x_raw, x, y, x_label, y_label, non_missing):
    """LOcally WEighted Scatterplot Smoothing (LOWESS) trendline function

    Requires `statsmodels` to be installed.

    Valid keys for the `trendline_options` dict are:

    - `frac` (`float`, default `0.6666666`): the `frac` parameter from the
    `statsmodels.api.nonparametric.lowess` function
    """

    valid_options = ["frac"]
    for k in trendline_options.keys():
        if k not in valid_options:
            raise ValueError(
                "LOWESS trendline_options keys must be one of [%s] but got '%s'"
                % (", ".join(valid_options), k)
            )

    import statsmodels.api as sm

    frac = trendline_options.get("frac", 0.6666666)
    y_out = sm.nonparametric.lowess(y, x, missing="drop", frac=frac)[:, 1]
    hover_header = "<b>LOWESS trendline</b><br><br>"
    return y_out, hover_header, None


def _pandas(mode, trendline_options, x_raw, y, non_missing):
    import numpy as np

    try:
        import pandas as pd
    except ImportError:
        msg = "Trendline requires pandas to be installed"
        raise ImportError(msg)

    modes = dict(rolling="Rolling", ewm="Exponentially Weighted", expanding="Expanding")
    trendline_options = trendline_options.copy()
    function_name = trendline_options.pop("function", "mean")
    function_args = trendline_options.pop("function_args", dict())

    series = pd.Series(np.copy(y), index=x_raw.to_pandas())

    # TODO: Narwhals Series/DataFrame do not support rolling, ewm nor expanding, therefore
    # it fallbacks to pandas Series independently of the original type.
    # Plotly issue: https://github.com/plotly/plotly.py/issues/4834
    # Narwhals issue: https://github.com/narwhals-dev/narwhals/issues/1254
    agg = getattr(series, mode)  # e.g. series.rolling
    agg_obj = agg(**trendline_options)  # e.g. series.rolling(**opts)
    function = getattr(agg_obj, function_name)  # e.g. series.rolling(**opts).mean
    y_out = function(**function_args)  # e.g. series.rolling(**opts).mean(**opts)
    y_out = y_out[non_missing]
    hover_header = "<b>%s %s trendline</b><br><br>" % (modes[mode], function_name)
    return y_out, hover_header, None


def rolling(trendline_options, x_raw, x, y, x_label, y_label, non_missing):
    """Rolling trendline function

    The value of the `function` key of the `trendline_options` dict is the function to
    use (defaults to `mean`) and the value of the `function_args` key are taken to be
    its arguments as a dict. The remainder of  the `trendline_options` dict is passed as
    keyword arguments into the `pandas.Series.rolling` function.
    """
    return _pandas("rolling", trendline_options, x_raw, y, non_missing)


def expanding(trendline_options, x_raw, x, y, x_label, y_label, non_missing):
    """Expanding trendline function

    The value of the `function` key of the `trendline_options` dict is the function to
    use (defaults to `mean`) and the value of the `function_args` key are taken to be
    its arguments as a dict. The remainder of  the `trendline_options` dict is passed as
    keyword arguments into the `pandas.Series.expanding` function.
    """
    return _pandas("expanding", trendline_options, x_raw, y, non_missing)


def ewm(trendline_options, x_raw, x, y, x_label, y_label, non_missing):
    """Exponentially Weighted Moment (EWM) trendline function

    The value of the `function` key of the `trendline_options` dict is the function to
    use (defaults to `mean`) and the value of the `function_args` key are taken to be
    its arguments as a dict. The remainder of  the `trendline_options` dict is passed as
    keyword arguments into the `pandas.Series.ewm` function.
    """
    return _pandas("ewm", trendline_options, x_raw, y, non_missing)
