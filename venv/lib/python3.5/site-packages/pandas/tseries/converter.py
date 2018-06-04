# flake8: noqa
import warnings

from pandas.plotting._converter import (time2num,
                                        TimeConverter, TimeFormatter,
                                        PeriodConverter, get_datevalue,
                                        DatetimeConverter,
                                        PandasAutoDateFormatter,
                                        PandasAutoDateLocator,
                                        MilliSecondLocator, get_finder,
                                        TimeSeries_DateLocator,
                                        TimeSeries_DateFormatter)


def register():
    from pandas.plotting._converter import register as register_
    msg = ("'pandas.tseries.converter.register' has been moved and renamed to "
           "'pandas.plotting.register_matplotlib_converters'. ")
    warnings.warn(msg, FutureWarning, stacklevel=2)
    register_()
