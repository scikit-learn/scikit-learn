import weakref


def _check_axes_has_been_used(ax):
    """Return true if the axes has been used"""
    used_attrs = ['lines', 'patches', 'texts', 'tables', 'artists',
                  'tables', 'images']
    msg = "The ax was already used in a matplotlib plot function"
    if any(getattr(ax, attr, None) for attr in used_attrs):
        raise ValueError(msg)


_SKLEARN_AX_DISP_OBJ_REF_KEY = "_sklearn_display_object_ref"


def _check_axes_has_display_object(display_obj, ax):
    """Check if axes has a weak ref to a display_obj

    If the weak ref does not exist or points to None, the axes will be assigned
    a reference to the passed in display_obj.

    Used when the display object needs to use the axes to define a space to
    create multiple plots on the axes.
    """
    if not hasattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY):
        setattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY, weakref.ref(display_obj))
        ax.set_axis_off()
        return display_obj

    ax_display_obj = getattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY)()

    if ax_display_obj is None:  # display obj was deleted
        setattr(ax, _SKLEARN_AX_DISP_OBJ_REF_KEY, weakref.ref(display_obj))
        ax.set_axis_off()
        return display_obj
    elif not isinstance(ax_display_obj, display_obj.__class__):
        raise ValueError("The ax was already used by another "
                         "display object which is not an "
                         "instance of {}".format(
                             display_obj.__class__.__name__))

    # ax._sklearn_display_object is an instance of display_obj.__class__
    return ax_display_obj
