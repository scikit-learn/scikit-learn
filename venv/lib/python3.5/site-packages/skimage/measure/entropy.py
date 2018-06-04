from numpy import unique
from scipy.stats import entropy as scipy_entropy


def shannon_entropy(image, base=2):
    """Calculate the Shannon entropy of an image.

    The Shannon entropy is defined as S = -sum(pk * log(pk)),
    where pk are frequency/probability of pixels of value k.

    Parameters
    ----------
    image : (N, M) ndarray
        Grayscale input image.
    base : float, optional
        The logarithmic base to use.

    Returns
    -------
    entropy : float

    Notes
    -----
    The returned value is measured in bits or shannon (Sh) for base=2, natural
    unit (nat) for base=np.e and hartley (Hart) for base=10.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Entropy_(information_theory)
    .. [2] https://en.wiktionary.org/wiki/Shannon_entropy

    Examples
    --------
    >>> from skimage import data
    >>> shannon_entropy(data.camera())
    7.0479552324230861
    """

    _, counts = unique(image, return_counts=True)
    return scipy_entropy(counts, base=base)
