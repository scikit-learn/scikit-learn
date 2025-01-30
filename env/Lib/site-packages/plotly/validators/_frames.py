import _plotly_utils.basevalidators


class FramesValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="frames", parent_name="", **kwargs):
        super(FramesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Frame"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            baseframe
                The name of the frame into which this frame's
                properties are merged before applying. This is
                used to unify properties and avoid needing to
                specify the same values for the same properties
                in multiple frames.
            data
                A list of traces this frame modifies. The
                format is identical to the normal trace
                definition.
            group
                An identifier that specifies the group to which
                the frame belongs, used by animate to select a
                subset of frames.
            layout
                Layout properties which this frame modifies.
                The format is identical to the normal layout
                definition.
            name
                A label by which to identify the frame
            traces
                A list of trace indices that identify the
                respective traces in the data attribute
""",
            ),
            **kwargs,
        )
