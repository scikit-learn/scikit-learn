import _plotly_utils.basevalidators


class UpdatemenusValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="updatemenus", parent_name="layout", **kwargs):
        super(UpdatemenusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Updatemenu"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            active
                Determines which button (by index starting from
                0) is considered active.
            bgcolor
                Sets the background color of the update menu
                buttons.
            bordercolor
                Sets the color of the border enclosing the
                update menu.
            borderwidth
                Sets the width (in px) of the border enclosing
                the update menu.
            buttons
                A tuple of :class:`plotly.graph_objects.layout.
                updatemenu.Button` instances or dicts with
                compatible properties
            buttondefaults
                When used in a template (as layout.template.lay
                out.updatemenu.buttondefaults), sets the
                default property values to use for elements of
                layout.updatemenu.buttons
            direction
                Determines the direction in which the buttons
                are laid out, whether in a dropdown menu or a
                row/column of buttons. For `left` and `up`, the
                buttons will still appear in left-to-right or
                top-to-bottom order respectively.
            font
                Sets the font of the update menu button text.
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
            pad
                Sets the padding around the buttons or dropdown
                menu.
            showactive
                Highlights active dropdown item or active
                button if true.
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
            type
                Determines whether the buttons are accessible
                via a dropdown menu or whether the buttons are
                stacked horizontally or vertically
            visible
                Determines whether or not the update menu is
                visible.
            x
                Sets the x position (in normalized coordinates)
                of the update menu.
            xanchor
                Sets the update menu's horizontal position
                anchor. This anchor binds the `x` position to
                the "left", "center" or "right" of the range
                selector.
            y
                Sets the y position (in normalized coordinates)
                of the update menu.
            yanchor
                Sets the update menu's vertical position anchor
                This anchor binds the `y` position to the
                "top", "middle" or "bottom" of the range
                selector.
""",
            ),
            **kwargs,
        )
