import _plotly_utils.basevalidators


class BoundsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="bounds", parent_name="layout.mapbox", **kwargs):
        super(BoundsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Bounds"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            east
                Sets the maximum longitude of the map (in
                degrees East) if `west`, `south` and `north`
                are declared.
            north
                Sets the maximum latitude of the map (in
                degrees North) if `east`, `west` and `south`
                are declared.
            south
                Sets the minimum latitude of the map (in
                degrees North) if `east`, `west` and `north`
                are declared.
            west
                Sets the minimum longitude of the map (in
                degrees East) if `east`, `south` and `north`
                are declared.
""",
            ),
            **kwargs,
        )
