import _plotly_utils.basevalidators


class GridValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="grid", parent_name="layout", **kwargs):
        super(GridValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Grid"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            columns
                The number of columns in the grid. If you
                provide a 2D `subplots` array, the length of
                its longest row is used as the default. If you
                give an `xaxes` array, its length is used as
                the default. But it's also possible to have a
                different length, if you want to leave a row at
                the end for non-cartesian subplots.
            domain
                :class:`plotly.graph_objects.layout.grid.Domain
                ` instance or dict with compatible properties
            pattern
                If no `subplots`, `xaxes`, or `yaxes` are given
                but we do have `rows` and `columns`, we can
                generate defaults using consecutive axis IDs,
                in two ways: "coupled" gives one x axis per
                column and one y axis per row. "independent"
                uses a new xy pair for each cell, left-to-right
                across each row then iterating rows according
                to `roworder`.
            roworder
                Is the first row the top or the bottom? Note
                that columns are always enumerated from left to
                right.
            rows
                The number of rows in the grid. If you provide
                a 2D `subplots` array or a `yaxes` array, its
                length is used as the default. But it's also
                possible to have a different length, if you
                want to leave a row at the end for non-
                cartesian subplots.
            subplots
                Used for freeform grids, where some axes may be
                shared across subplots but others are not. Each
                entry should be a cartesian subplot id, like
                "xy" or "x3y2", or "" to leave that cell empty.
                You may reuse x axes within the same column,
                and y axes within the same row. Non-cartesian
                subplots and traces that support `domain` can
                place themselves in this grid separately using
                the `gridcell` attribute.
            xaxes
                Used with `yaxes` when the x and y axes are
                shared across columns and rows. Each entry
                should be an x axis id like "x", "x2", etc., or
                "" to not put an x axis in that column. Entries
                other than "" must be unique. Ignored if
                `subplots` is present. If missing but `yaxes`
                is present, will generate consecutive IDs.
            xgap
                Horizontal space between grid cells, expressed
                as a fraction of the total width available to
                one cell. Defaults to 0.1 for coupled-axes
                grids and 0.2 for independent grids.
            xside
                Sets where the x axis labels and titles go.
                "bottom" means the very bottom of the grid.
                "bottom plot" is the lowest plot that each x
                axis is used in. "top" and "top plot" are
                similar.
            yaxes
                Used with `yaxes` when the x and y axes are
                shared across columns and rows. Each entry
                should be an y axis id like "y", "y2", etc., or
                "" to not put a y axis in that row. Entries
                other than "" must be unique. Ignored if
                `subplots` is present. If missing but `xaxes`
                is present, will generate consecutive IDs.
            ygap
                Vertical space between grid cells, expressed as
                a fraction of the total height available to one
                cell. Defaults to 0.1 for coupled-axes grids
                and 0.3 for independent grids.
            yside
                Sets where the y axis labels and titles go.
                "left" means the very left edge of the grid.
                *left plot* is the leftmost plot that each y
                axis is used in. "right" and *right plot* are
                similar.
""",
            ),
            **kwargs,
        )
