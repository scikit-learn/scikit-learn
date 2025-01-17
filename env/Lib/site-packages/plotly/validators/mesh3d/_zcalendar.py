import _plotly_utils.basevalidators


class ZcalendarValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="zcalendar", parent_name="mesh3d", **kwargs):
        super(ZcalendarValidator, self).__init__(
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
