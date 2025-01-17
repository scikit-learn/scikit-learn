from _plotly_utils.basevalidators import EnumeratedValidator, NumberValidator


class EasingValidator(EnumeratedValidator):
    def __init__(self, plotly_name="easing", parent_name="batch_animate", **_):
        super(EasingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            values=[
                "linear",
                "quad",
                "cubic",
                "sin",
                "exp",
                "circle",
                "elastic",
                "back",
                "bounce",
                "linear-in",
                "quad-in",
                "cubic-in",
                "sin-in",
                "exp-in",
                "circle-in",
                "elastic-in",
                "back-in",
                "bounce-in",
                "linear-out",
                "quad-out",
                "cubic-out",
                "sin-out",
                "exp-out",
                "circle-out",
                "elastic-out",
                "back-out",
                "bounce-out",
                "linear-in-out",
                "quad-in-out",
                "cubic-in-out",
                "sin-in-out",
                "exp-in-out",
                "circle-in-out",
                "elastic-in-out",
                "back-in-out",
                "bounce-in-out",
            ],
        )


class DurationValidator(NumberValidator):
    def __init__(self, plotly_name="duration"):
        super(DurationValidator, self).__init__(
            plotly_name=plotly_name, parent_name="batch_animate", min=0
        )
