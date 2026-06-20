"""Utility for issuing warnings that respect scikit-learn global config."""

import warnings


def sklearn_warn(message, category=UserWarning, stacklevel=2):
    """Issue a warning, respecting the `disable_warnings` global config.
    
    Parameters
    ----------
    message : str
        The warning message.
    category : Warning class or tuple of Warning classes
        The category of the warning.
    stacklevel : int
        The stack level for the warning.
    """
    # Import here to avoid circular imports
    from sklearn import get_config
    
    config = get_config()
    disable_warnings = config.get("disable_warnings", False)
    
    # If disable_warnings is True, suppress all warnings
    if disable_warnings is True:
        return
    
    # If it's a class or tuple of classes, check if our category matches
    if disable_warnings:
        if isinstance(disable_warnings, type) and issubclass(
            category, disable_warnings
        ):
            return
        if isinstance(disable_warnings, tuple) and any(
            issubclass(category, cat) for cat in disable_warnings
        ):
            return
    
    warnings.warn(message, category=category, stacklevel=stacklevel)