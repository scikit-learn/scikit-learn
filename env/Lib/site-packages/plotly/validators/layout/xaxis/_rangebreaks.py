import _plotly_utils.basevalidators


class RangebreaksValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="rangebreaks", parent_name="layout.xaxis", **kwargs):
        super(RangebreaksValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Rangebreak"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bounds
                Sets the lower and upper bounds of this axis
                rangebreak. Can be used with `pattern`.
            dvalue
                Sets the size of each `values` item. The
                default is one day in milliseconds.
            enabled
                Determines whether this axis rangebreak is
                enabled or disabled. Please note that
                `rangebreaks` only work for "date" axis type.
            name
                When used in a template, named items are
                created in the output figure in addition to any
                items the figure already has in this array. You
                can modify these items in the output figure by
                making your own item with `templateitemname`
                matching this `name` alongside your
                modifications (including `visible: false` or
                `enabled: false` to hide it). Has no effect
                outside of a template.
            pattern
                Determines a pattern on the time line that
                generates breaks. If *day of week* - days of
                the week in English e.g. 'Sunday' or `sun`
                (matching is case-insensitive and considers
                only the first three characters), as well as
                Sunday-based integers between 0 and 6. If
                "hour" - hour (24-hour clock) as decimal
                numbers between 0 and 24. for more info.
                Examples: - { pattern: 'day of week', bounds:
                [6, 1] }  or simply { bounds: ['sat', 'mon'] }
                breaks from Saturday to Monday (i.e. skips the
                weekends). - { pattern: 'hour', bounds: [17, 8]
                }   breaks from 5pm to 8am (i.e. skips non-work
                hours).
            templateitemname
                Used to refer to a named item in this array in
                the template. Named items from the template
                will be created even without a matching item in
                the input figure, but you can modify one by
                making an item with `templateitemname` matching
                its `name`, alongside your modifications
                (including `visible: false` or `enabled: false`
                to hide it). If there is no template or no
                matching item, this item will be hidden unless
                you explicitly show it with `visible: true`.
            values
                Sets the coordinate values corresponding to the
                rangebreaks. An alternative to `bounds`. Use
                `dvalue` to set the size of the values along
                the axis.
""",
            ),
            **kwargs,
        )
