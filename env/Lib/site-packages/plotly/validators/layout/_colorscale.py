import _plotly_utils.basevalidators


class ColorscaleValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="colorscale", parent_name="layout", **kwargs):
        super(ColorscaleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Colorscale"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            diverging
                Sets the default diverging colorscale. Note
                that `autocolorscale` must be true for this
                attribute to work.
            sequential
                Sets the default sequential colorscale for
                positive values. Note that `autocolorscale`
                must be true for this attribute to work.
            sequentialminus
                Sets the default sequential colorscale for
                negative values. Note that `autocolorscale`
                must be true for this attribute to work.
""",
            ),
            **kwargs,
        )
