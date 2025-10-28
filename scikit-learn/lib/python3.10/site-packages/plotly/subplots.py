import plotly.graph_objects as go
from . import _subplots as _sub
from ._subplots import SubplotXY, SubplotDomain, SubplotRef  # noqa: F401


def make_subplots(
    rows=1,
    cols=1,
    shared_xaxes=False,
    shared_yaxes=False,
    start_cell="top-left",
    print_grid=False,
    horizontal_spacing=None,
    vertical_spacing=None,
    subplot_titles=None,
    column_widths=None,
    row_heights=None,
    specs=None,
    insets=None,
    column_titles=None,
    row_titles=None,
    x_title=None,
    y_title=None,
    figure=None,
    **kwargs,
) -> go.Figure:
    """
    Return an instance of plotly.graph_objs.Figure with predefined subplots
    configured in 'layout'.

    Parameters
    ----------
    rows: int (default 1)
        Number of rows in the subplot grid. Must be greater than zero.

    cols: int (default 1)
        Number of columns in the subplot grid. Must be greater than zero.

    shared_xaxes: boolean or str (default False)
        Assign shared (linked) x-axes for 2D cartesian subplots

          - True or 'columns': Share axes among subplots in the same column
          - 'rows': Share axes among subplots in the same row
          - 'all': Share axes across all subplots in the grid.

    shared_yaxes: boolean or str (default False)
        Assign shared (linked) y-axes for 2D cartesian subplots

          - 'columns': Share axes among subplots in the same column
          - True or 'rows': Share axes among subplots in the same row
          - 'all': Share axes across all subplots in the grid.

    start_cell: 'bottom-left' or 'top-left' (default 'top-left')
        Choose the starting cell in the subplot grid used to set the
        domains_grid of the subplots.

          - 'top-left': Subplots are numbered with (1, 1) in the top
                        left corner
          - 'bottom-left': Subplots are numbererd with (1, 1) in the bottom
                           left corner

    print_grid: boolean (default False):
        If True, prints a string representation of the plot grid. Grid may
        also be printed using the `Figure.print_grid()` method on the
        resulting figure.

    horizontal_spacing: float (default 0.2 / cols)
        Space between subplot columns in normalized plot coordinates. Must be
        a float between 0 and 1.

        Applies to all columns (use 'specs' subplot-dependents spacing)

    vertical_spacing: float (default 0.3 / rows)
        Space between subplot rows in normalized plot coordinates. Must be
        a float between 0 and 1.

        Applies to all rows (use 'specs' subplot-dependents spacing)

    subplot_titles: list of str or None (default None)
        Title of each subplot as a list in row-major ordering.

        Empty strings ("") can be included in the list if no subplot title
        is desired in that space so that the titles are properly indexed.

    specs: list of lists of dict or None (default None)
        Per subplot specifications of subplot type, row/column spanning, and
        spacing.

        ex1: specs=[[{}, {}], [{'colspan': 2}, None]]

        ex2: specs=[[{'rowspan': 2}, {}], [None, {}]]

        - Indices of the outer list correspond to subplot grid rows
          starting from the top, if start_cell='top-left',
          or bottom, if start_cell='bottom-left'.
          The number of rows in 'specs' must be equal to 'rows'.

        - Indices of the inner lists correspond to subplot grid columns
          starting from the left. The number of columns in 'specs'
          must be equal to 'cols'.

        - Each item in the 'specs' list corresponds to one subplot
          in a subplot grid. (N.B. The subplot grid has exactly 'rows'
          times 'cols' cells.)

        - Use None for a blank a subplot cell (or to move past a col/row span).

        - Note that specs[0][0] has the specs of the 'start_cell' subplot.

        - Each item in 'specs' is a dictionary.
            The available keys are:
            * type (string, default 'xy'): Subplot type. One of
                - 'xy': 2D Cartesian subplot type for scatter, bar, etc.
                - 'scene': 3D Cartesian subplot for scatter3d, cone, etc.
                - 'polar': Polar subplot for scatterpolar, barpolar, etc.
                - 'ternary': Ternary subplot for scatterternary
                - 'map': Map subplot for scattermap
                - 'mapbox': Mapbox subplot for scattermapbox
                - 'domain': Subplot type for traces that are individually
                            positioned. pie, parcoords, parcats, etc.
                - trace type: A trace type which will be used to determine
                              the appropriate subplot type for that trace

            * secondary_y (bool, default False): If True, create a secondary
                y-axis positioned on the right side of the subplot. Only valid
                if type='xy'.
            * colspan (int, default 1): number of subplot columns
                for this subplot to span.
            * rowspan (int, default 1): number of subplot rows
                for this subplot to span.
            * l (float, default 0.0): padding left of cell
            * r (float, default 0.0): padding right of cell
            * t (float, default 0.0): padding right of cell
            * b (float, default 0.0): padding bottom of cell

        - Note: Use 'horizontal_spacing' and 'vertical_spacing' to adjust
          the spacing in between the subplots.

    insets: list of dict or None (default None):
        Inset specifications.  Insets are subplots that overlay grid subplots

        - Each item in 'insets' is a dictionary.
            The available keys are:

            * cell (tuple, default=(1,1)): (row, col) index of the
                subplot cell to overlay inset axes onto.
            * type (string, default 'xy'): Subplot type
            * l (float, default=0.0): padding left of inset
                  in fraction of cell width
            * w (float or 'to_end', default='to_end') inset width
                  in fraction of cell width ('to_end': to cell right edge)
            * b (float, default=0.0): padding bottom of inset
                  in fraction of cell height
            * h (float or 'to_end', default='to_end') inset height
                  in fraction of cell height ('to_end': to cell top edge)

    column_widths: list of numbers or None (default None)
        list of length `cols` of the relative widths of each column of subplots.
        Values are normalized internally and used to distribute overall width
        of the figure (excluding padding) among the columns.

        For backward compatibility, may also be specified using the
        `column_width` keyword argument.

    row_heights: list of numbers or None (default None)
        list of length `rows` of the relative heights of each row of subplots.
        If start_cell='top-left' then row heights are applied top to bottom.
        Otherwise, if start_cell='bottom-left' then row heights are applied
        bottom to top.

        For backward compatibility, may also be specified using the
        `row_width` kwarg. If specified as `row_width`, then the width values
        are applied from bottom to top regardless of the value of start_cell.
        This matches the legacy behavior of the `row_width` argument.

    column_titles: list of str or None (default None)
        list of length `cols` of titles to place above the top subplot in
        each column.

    row_titles: list of str or None (default None)
        list of length `rows` of titles to place on the right side of each
        row of subplots. If start_cell='top-left' then row titles are
        applied top to bottom. Otherwise, if start_cell='bottom-left' then
        row titles are applied bottom to top.

    x_title: str or None (default None)
        Title to place below the bottom row of subplots,
        centered horizontally

    y_title: str or None (default None)
        Title to place to the left of the left column of subplots,
        centered vertically

    figure: go.Figure or None (default None)
        If None, a new go.Figure instance will be created and its axes will be
        populated with those corresponding to the requested subplot geometry and
        this new figure will be returned.
        If a go.Figure instance, the axes will be added to the
        layout of this figure and this figure will be returned. If the figure
        already contains axes, they will be overwritten.

    Examples
    --------

    Example 1:

    >>> # Stack two subplots vertically, and add a scatter trace to each
    >>> from plotly.subplots import make_subplots
    >>> import plotly.graph_objects as go
    >>> fig = make_subplots(rows=2)

    This is the format of your plot grid:
    [ (1,1) xaxis1,yaxis1 ]
    [ (2,1) xaxis2,yaxis2 ]

    >>> fig.add_scatter(y=[2, 1, 3], row=1, col=1) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_scatter(y=[1, 3, 2], row=2, col=1) # doctest: +ELLIPSIS
    Figure(...)

    or see Figure.add_trace

    Example 2:

    >>> # Stack a scatter plot
    >>> fig = make_subplots(rows=2, shared_xaxes=True)

    This is the format of your plot grid:
    [ (1,1) xaxis1,yaxis1 ]
    [ (2,1) xaxis2,yaxis2 ]

    >>> fig.add_scatter(y=[2, 1, 3], row=1, col=1) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_scatter(y=[1, 3, 2], row=2, col=1) # doctest: +ELLIPSIS
    Figure(...)

    Example 3:

    >>> # irregular subplot layout (more examples below under 'specs')
    >>> fig = make_subplots(rows=2, cols=2,
    ...                     specs=[[{}, {}],
    ...                     [{'colspan': 2}, None]])

    This is the format of your plot grid:
    [ (1,1) xaxis1,yaxis1 ]  [ (1,2) xaxis2,yaxis2 ]
    [ (2,1) xaxis3,yaxis3           -              ]

    >>> fig.add_trace(go.Scatter(x=[1,2,3], y=[2,1,2]), row=1, col=1) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_trace(go.Scatter(x=[1,2,3], y=[2,1,2]), row=1, col=2) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_trace(go.Scatter(x=[1,2,3], y=[2,1,2]), row=2, col=1) # doctest: +ELLIPSIS
    Figure(...)

    Example 4:

    >>> # insets
    >>> fig = make_subplots(insets=[{'cell': (1,1), 'l': 0.7, 'b': 0.3}])

    This is the format of your plot grid:
    [ (1,1) xaxis1,yaxis1 ]

    With insets:
    [ xaxis2,yaxis2 ] over [ (1,1) xaxis1,yaxis1 ]

    >>> fig.add_scatter(x=[1,2,3], y=[2,1,1]) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_scatter(x=[1,2,3], y=[2,1,2], xaxis='x2', yaxis='y2') # doctest: +ELLIPSIS
    Figure(...)

    Example 5:

    >>> # include subplot titles
    >>> fig = make_subplots(rows=2, subplot_titles=('Plot 1','Plot 2'))

    This is the format of your plot grid:
    [ (1,1) x1,y1 ]
    [ (2,1) x2,y2 ]

    >>> fig.add_scatter(x=[1,2,3], y=[2,1,2], row=1, col=1) # doctest: +ELLIPSIS
    Figure(...)
    >>> fig.add_bar(x=[1,2,3], y=[2,1,2], row=2, col=1) # doctest: +ELLIPSIS
    Figure(...)

    Example 6:

    Subplot with mixed subplot types

    >>> fig = make_subplots(rows=2, cols=2,
    ...                     specs=[[{'type': 'xy'},    {'type': 'polar'}],
    ...                            [{'type': 'scene'}, {'type': 'ternary'}]])

    >>> fig.add_traces(
    ...     [go.Scatter(y=[2, 3, 1]),
    ...      go.Scatterpolar(r=[1, 3, 2], theta=[0, 45, 90]),
    ...      go.Scatter3d(x=[1, 2, 1], y=[2, 3, 1], z=[0, 3, 5]),
    ...      go.Scatterternary(a=[0.1, 0.2, 0.1],
    ...                        b=[0.2, 0.3, 0.1],
    ...                        c=[0.7, 0.5, 0.8])],
    ...     rows=[1, 1, 2, 2],
    ...     cols=[1, 2, 1, 2]) # doctest: +ELLIPSIS
    Figure(...)
    """

    return _sub.make_subplots(
        rows,
        cols,
        shared_xaxes,
        shared_yaxes,
        start_cell,
        print_grid,
        horizontal_spacing,
        vertical_spacing,
        subplot_titles,
        column_widths,
        row_heights,
        specs,
        insets,
        column_titles,
        row_titles,
        x_title,
        y_title,
        figure,
        **kwargs,
    )
