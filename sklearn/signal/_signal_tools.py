def envelope_detector(signal):
    """
    Compute the envelope of a 1D signal using a naive method.

    Parameters
    ----------
    signal : array-like
          Input signal.
    Returns
    -------
    envelope : array-like
          Envelope of the input signal.      
    """
    import numpy as np
    return np.abs(signal)