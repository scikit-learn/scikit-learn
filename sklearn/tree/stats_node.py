from __future__ import division, print_function


class StatsNode(object):
    """Container of the statistic for a specific node.

    It aims at improving readability in the other classes while
    manipulating statistics.

    Parameters
    ----------
    sum_residuals: float,
        Sum of the residuals.

    sum_sq_residuals: float,
        Sum of the squared residuals.

    n_samples: int,
        Number of samples.

    sum_weighted_samples: float,
        Sum of the weights associated to samples.

    Attributes
    ----------
    sum_residuals: float,
        Sum of the residuals.

    sum_sq_residuals: float,
        Sum of the squared residuals.

    n_samples: int,
        Number of samples.

    sum_weighted_samples: float,
        Sum of the weights associated to samples.
    """

    def __init__(self, sum_residuals, sum_sq_residuals, n_samples,
                 sum_weighted_samples):
        self._sum_residuals = float(sum_residuals)
        self._sum_sq_residuals = float(sum_sq_residuals)
        self._n_samples = int(n_samples)
        self._sum_weighted_samples = float(sum_weighted_samples)

    def reset(self, sum_residuals, sum_sq_residuals, n_samples,
              sum_weighted_samples):
        """Reset the different stats"""
        self._sum_residuals = float(sum_residuals)
        self._sum_sq_residuals = float(sum_sq_residuals)
        self._n_samples = int(n_samples)
        self._sum_weighted_samples = float(sum_weighted_samples)

    def clear(self):
        """Clear the stats"""
        self.__init__()

    @property
    def sum_residuals(self):
        return self._sum_residuals

    @sum_residuals.setter
    def sum_residuals(self, val):
        self._sum_residuals = float(val)

    @property
    def sum_sq_residuals(self):
        return self._sum_sq_residuals

    @sum_sq_residuals.setter
    def sum_sq_residuals(self, val):
        self._sum_sq_residuals = float(val)

    @property
    def n_samples(self):
        return self._n_samples

    @n_samples.setter
    def n_samples(self, val):
        self._n_samples = int(val)

    @property
    def sum_weighted_samples(self):
        return self._sum_weighted_samples

    @sum_weighted_samples.setter
    def sum_weighted_samples(self, val):
        return self._sum_weighted_samples

    def __iadd__(self, val):
        """In place add to StatsNodes together."""
        if isinstance(val, StatsNode):
            self._sum_residuals += val.sum_residuals
            self._sum_sq_residuals += val.sum_sq_residuals
            self._n_samples += val.n_samples
            self._sum_weighted_samples += val.sum_weighted_samples
        else:
            return NotImplemented

        return self

    def __add__(self, val):
        """Add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self._sum_residuals + val.sum_residuals,
                self._sum_sq_residuals + val.sum_sq_residuals,
                self._n_samples + val.n_samples,
                self._sum_weighted_samples + val.sum_weighted_samples)
        else:
            return NotImplemented

    def __isub__(self, val):
        """In place subtract to StatsNodes together."""
        if isinstance(val, StatsNode):
            self._sum_residuals -= val.sum_residuals
            self._sum_sq_residuals -= val.sum_sq_residuals
            self._n_samples -= val.n_samples
            self._sum_weighted_samples -= val.sum_weighted_samples
        else:
            return NotImplemented

        return self

    def __sub__(self, val):
        """Add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self._sum_residuals - val.sum_residuals,
                self._sum_sq_residuals - val.sum_sq_residuals,
                self._n_samples - val.n_samples,
                self._sum_weighted_samples - val.sum_weighted_samples)
        else:
            return NotImplemented

    def __str__(self):
        info = ("Sum of residuals: {}\n"
                "Sum of the squared residuals: {}\n"
                "Number of samples: {}\n"
                "Sum of weights associated to samples: {}\n".format(
                    self._sum_residuals,
                    self._sum_sq_residuals,
                    self._n_samples,
                    self._sum_weighted_samples))

        return info
