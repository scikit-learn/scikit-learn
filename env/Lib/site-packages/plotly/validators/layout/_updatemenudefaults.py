import _plotly_utils.basevalidators


class UpdatemenudefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="updatemenudefaults", parent_name="layout", **kwargs
    ):
        super(UpdatemenudefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Updatemenu"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
