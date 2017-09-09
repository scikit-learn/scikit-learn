#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

from __future__ import division, print_function

from numpy import all


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
    sum_y: float,
        Sum of the y.

    sum_sq_y: float,
        Sum of the squared y.

    n_samples: int,
        Number of samples.

    sum_weighted_samples: float,
        Sum of the weights associated to samples.
    """

    def __init__(self, sum_y, sum_sq_y, n_samples,
                 sum_weighted_samples):
        self.sum_y = float(sum_y)
        self.sum_sq_y = float(sum_sq_y)
        self.n_samples = int(n_samples)
        self.sum_weighted_samples = float(sum_weighted_samples)

    def reset(self, sum_y, sum_sq_y, n_samples,
              sum_weighted_samples):
        """Reset the different stats"""
        self.sum_y = float(sum_y)
        self.sum_sq_y = float(sum_sq_y)
        self.n_samples = int(n_samples)
        self.sum_weighted_samples = float(sum_weighted_samples)

    def clear(self):
        """Clear the stats"""
        self.sum_y = float(0.)
        self.sum_sq_y = float(0.)
        self.n_samples = int(0)
        self.sum_weighted_samples = float(0.)

    def __iadd__(self, val):
        """In place add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_y + val.sum_y,
                self.sum_sq_y + val.sum_sq_y,
                self.n_samples + val.n_samples,
                self.sum_weighted_samples + val.sum_weighted_samples)
        else:
            return NotImplemented

    def __add__(self, val):
        """Add to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_y + val.sum_y,
                self.sum_sq_y + val.sum_sq_y,
                self.n_samples + val.n_samples,
                self.sum_weighted_samples + val.sum_weighted_samples)
        else:
            return NotImplemented

    def __isub__(self, val):
        """In place subtract to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_y - val.sum_y,
                self.sum_sq_y - val.sum_sq_y,
                self.n_samples - val.n_samples,
                self.sum_weighted_samples - val.sum_weighted_samples)
        else:
            return NotImplemented

    def __sub__(self, val):
        """Subtract to StatsNodes together."""
        if isinstance(val, StatsNode):
            return StatsNode(
                self.sum_y - val.sum_y,
                self.sum_sq_y - val.sum_sq_y,
                self.n_samples - val.n_samples,
                self.sum_weighted_samples - val.sum_weighted_samples)
        else:
            return NotImplemented

    def __eq__(self, val):
        """Compare to StatsNode."""
        if isinstance(val, StatsNode):
            return all([
                self.sum_y == val.sum_y,
                self.sum_sq_y == val.sum_sq_y,
                self.n_samples == val.n_samples,
                self.sum_weighted_samples == val.sum_weighted_samples
            ])
        else:
            return NotImplemented

    # FIXME for debugging purpose
    def __str__(self):
        info = ("Sum of y: {}\n"
                "Sum of the squared y: {}\n"
                "Number of samples: {}\n"
                "Sum of weights associated to samples: {}\n".format(
                    self.sum_y,
                    self.sum_sq_y,
                    self.n_samples,
                    self.sum_weighted_samples))

        return info

    def copy_to(self, dest_stats_node):
        dest_stats_node.sum_y = self.sum_y
        dest_stats_node.sum_sq_y = self.sum_sq_y
        dest_stats_node.n_samples = self.n_samples
        dest_stats_node.sum_weighted_samples = self.sum_weighted_samples
