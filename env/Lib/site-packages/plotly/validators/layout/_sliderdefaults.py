import _plotly_utils.basevalidators


class SliderdefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="sliderdefaults", parent_name="layout", **kwargs):
        super(SliderdefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Slider"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
