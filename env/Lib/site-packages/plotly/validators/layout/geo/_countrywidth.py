import _plotly_utils.basevalidators


class CountrywidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="countrywidth", parent_name="layout.geo", **kwargs):
        super(CountrywidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
