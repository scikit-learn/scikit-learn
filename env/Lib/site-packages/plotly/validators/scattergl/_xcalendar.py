import _plotly_utils.basevalidators


class XcalendarValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="xcalendar", parent_name="scattergl", **kwargs):
        super(XcalendarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values",
                [
                    "chinese",
                    "coptic",
                    "discworld",
                    "ethiopian",
                    "gregorian",
                    "hebrew",
                    "islamic",
                    "jalali",
                    "julian",
                    "mayan",
                    "nanakshahi",
                    "nepali",
                    "persian",
                    "taiwan",
                    "thai",
                    "ummalqura",
                ],
            ),
            **kwargs,
        )
