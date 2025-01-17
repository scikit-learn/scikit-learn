import _plotly_utils.basevalidators


class AutorangeoptionsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self,
        plotly_name="autorangeoptions",
        parent_name="layout.polar.radialaxis",
        **kwargs,
    ):
        super(AutorangeoptionsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Autorangeoptions"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            clipmax
                Clip autorange maximum if it goes beyond this
                value. Has no effect when
                `autorangeoptions.maxallowed` is provided.
            clipmin
                Clip autorange minimum if it goes beyond this
                value. Has no effect when
                `autorangeoptions.minallowed` is provided.
            include
                Ensure this value is included in autorange.
            includesrc
                Sets the source reference on Chart Studio Cloud
                for `include`.
            maxallowed
                Use this value exactly as autorange maximum.
            minallowed
                Use this value exactly as autorange minimum.
""",
            ),
            **kwargs,
        )
