
def _check_axes_has_been_used(ax):
    """Return true if the axes has been used"""
    used_attrs = ['lines', 'patches', 'texts', 'tables', 'artists',
                  'tables', 'images']
    msg = "The ax was already used in another plot function"
    for used_attr in used_attrs:
        if hasattr(ax, used_attr) and getattr(ax, used_attr):
            raise ValueError(msg)
