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
        self.sum_residuals = float(sum_residuals)
        self.sum_sq_residuals = float(sum_sq_residuals)
        self.n_samples = int(n_samples)
        self.sum_weighted_samples = float(sum_weighted_samples)

    def reset(self, sum_residuals, sum_sq_residuals, n_samples,
              sum_weighted_samples):
        """Reset the different stats"""
        self.sum_residuals = float(sum_residuals)
        self.sum_sq_residuals = float(sum_sq_residuals)
        self.n_samples = int(n_samples)
        self.sum_weighted_samples = float(sum_weighted_samples)

    def clear(self):
        """Clear the stats"""
        self.sum_residuals = 0
        self.sum_sq_residuals = 0
        self.n_samples = 0
        self.sum_weighted_samples = 0

    def __iadd__(self, val):
        """In place add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_residuals + val.sum_residuals,
                self.sum_sq_residuals + val.sum_sq_residuals,
                self.n_samples + val.n_samples,
                self.sum_weighted_samples + val.sum_weighted_samples)
        else:
            return NotImplemented

    def __add__(self, val):
        """Add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_residuals + val.sum_residuals,
                self.sum_sq_residuals + val.sum_sq_residuals,
                self.n_samples + val.n_samples,
                self.sum_weighted_samples + val.sum_weighted_samples)
        else:
            return NotImplemented

    def __isub__(self, val):
        """In place subtract to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_residuals - val.sum_residuals,
                self.sum_sq_residuals - val.sum_sq_residuals,
                self.n_samples - val.n_samples,
                self.sum_weighted_samples - val.sum_weighted_samples)
        else:
            return NotImplemented

    def __sub__(self, val):
        """Subtract to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_residuals - val.sum_residuals,
                self.sum_sq_residuals - val.sum_sq_residuals,
                self.n_samples - val.n_samples,
                self.sum_weighted_samples - val.sum_weighted_samples)
        else:
            return NotImplemented

    def __str__(self):
        info = ("Sum of residuals: {}\n"
                "Sum of the squared residuals: {}\n"
                "Number of samples: {}\n"
                "Sum of weights associated to samples: {}\n".format(
                    self.sum_residuals,
                    self.sum_sq_residuals,
                    self.n_samples,
                    self.sum_weighted_samples))

        return info
