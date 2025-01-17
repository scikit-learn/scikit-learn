def _swatches(module_names, module_contents, template=None):
    """
    Parameters
    ----------
    template : str or dict or plotly.graph_objects.layout.Template instance
        The figure template name or definition.

    Returns
    -------
    fig : graph_objects.Figure containing the displayed image
        A `Figure` object. This figure demonstrates the color scales and
        sequences in this module, as stacked bar charts.
    """
    import plotly.graph_objs as go
    from plotly.express._core import apply_default_cascade

    args = dict(template=template)
    apply_default_cascade(args)

    sequences = [
        (k, v)
        for k, v in module_contents.items()
        if not (k.startswith("_") or k.startswith("swatches") or k.endswith("_r"))
    ]

    return go.Figure(
        data=[
            go.Bar(
                orientation="h",
                y=[name] * len(colors),
                x=[1] * len(colors),
                customdata=list(range(len(colors))),
                marker=dict(color=colors),
                hovertemplate="%{y}[%{customdata}] = %{marker.color}<extra></extra>",
            )
            for name, colors in reversed(sequences)
        ],
        layout=dict(
            title="plotly.colors." + module_names.split(".")[-1],
            barmode="stack",
            barnorm="fraction",
            bargap=0.5,
            showlegend=False,
            xaxis=dict(range=[-0.02, 1.02], showticklabels=False, showgrid=False),
            height=max(600, 40 * len(sequences)),
            template=args["template"],
            margin=dict(b=10),
        ),
    )


def _swatches_continuous(module_names, module_contents, template=None):
    """
    Parameters
    ----------
    template : str or dict or plotly.graph_objects.layout.Template instance
        The figure template name or definition.

    Returns
    -------
    fig : graph_objects.Figure containing the displayed image
        A `Figure` object. This figure demonstrates the color scales and
        sequences in this module, as stacked bar charts.
    """
    import plotly.graph_objs as go
    from plotly.express._core import apply_default_cascade

    args = dict(template=template)
    apply_default_cascade(args)

    sequences = [
        (k, v)
        for k, v in module_contents.items()
        if not (k.startswith("_") or k.startswith("swatches") or k.endswith("_r"))
    ]

    n = 100

    return go.Figure(
        data=[
            go.Bar(
                orientation="h",
                y=[name] * n,
                x=[1] * n,
                customdata=[(x + 1) / n for x in range(n)],
                marker=dict(color=list(range(n)), colorscale=name, line_width=0),
                hovertemplate="%{customdata}",
                name=name,
            )
            for name, colors in reversed(sequences)
        ],
        layout=dict(
            title="plotly.colors." + module_names.split(".")[-1],
            barmode="stack",
            barnorm="fraction",
            bargap=0.3,
            showlegend=False,
            xaxis=dict(range=[-0.02, 1.02], showticklabels=False, showgrid=False),
            height=max(600, 40 * len(sequences)),
            width=500,
            template=args["template"],
            margin=dict(b=10),
        ),
    )


def _swatches_cyclical(module_names, module_contents, template=None):
    """
    Parameters
    ----------
    template : str or dict or plotly.graph_objects.layout.Template instance
        The figure template name or definition.

    Returns
    -------
    fig : graph_objects.Figure containing the displayed image
        A `Figure` object. This figure demonstrates the color scales and
        sequences in this module, as polar bar charts.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plotly.express._core import apply_default_cascade

    args = dict(template=template)
    apply_default_cascade(args)

    rows = 2
    cols = 4
    scales = [
        (k, v)
        for k, v in module_contents.items()
        if not (k.startswith("_") or k.startswith("swatches") or k.endswith("_r"))
    ]
    names = [name for name, colors in scales]
    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=names,
        specs=[[{"type": "polar"}] * cols] * rows,
    )

    for i, (name, scale) in enumerate(scales):
        fig.add_trace(
            go.Barpolar(
                r=[1] * int(360 / 5),
                theta=list(range(0, 360, 5)),
                marker_color=list(range(0, 360, 5)),
                marker_cmin=0,
                marker_cmax=360,
                marker_colorscale=name,
                name=name,
            ),
            row=int(i / cols) + 1,
            col=i % cols + 1,
        )
    fig.update_traces(width=5.2, marker_line_width=0, base=0.5, showlegend=False)
    fig.update_polars(angularaxis_visible=False, radialaxis_visible=False)
    fig.update_layout(
        title="plotly.colors." + module_names.split(".")[-1], template=args["template"]
    )
    return fig
