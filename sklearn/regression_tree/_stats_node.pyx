cdef class StatsNode:
    """Container of the statistic for a specific node.

    It aims at improving readability in the other classes while
    manipulating statistics.

    Parameters
    ----------
    sum_residuals: double,
        Sum of the residuals.

    sum_sq_residuals: double,
        Sum of the squared residuals.

    n_samples: int,
        Number of samples.

    sum_weighted_samples: double,
        Sum of the weights associated to samples.

    Attributes
    ----------
    sum_y: double,
        Sum of the y.

    sum_sq_y: double,
        Sum of the squared y.

    n_samples: int,
        Number of samples.

    sum_weighted_samples: double,
        Sum of the weights associated to samples.
    """

    def __init__(self, double sum_y, double sum_sq_y, int n_samples,
                 double sum_weighted_samples):
        self.sum_y = sum_y
        self.sum_sq_y = sum_sq_y
        self.n_samples = n_samples
        self.sum_weighted_samples = sum_weighted_samples

    cdef void reset(self, double sum_y, double sum_sq_y, int n_samples,
                     double sum_weighted_samples):
        """Reset the different stats"""
        self.sum_y = sum_y
        self.sum_sq_y = sum_sq_y
        self.n_samples = n_samples
        self.sum_weighted_samples = sum_weighted_samples

    cdef void clear(self):
        """Clear the stats"""
        self.sum_y = 0.
        self.sum_sq_y = 0.
        self.n_samples = 0
        self.sum_weighted_samples = 0.

    def __iadd__(self, StatsNode val):
        """In place add to StatsNodes together."""
        self.sum_y += val.sum_y
        self.sum_sq_y += val.sum_sq_y
        self.n_samples += val.n_samples
        self.sum_weighted_samples += val.sum_weighted_samples

        return self

    def __add__(self, StatsNode val):
        """Add to StatsNodes together."""
        return StatsNode(
            self.sum_y + val.sum_y,
            self.sum_sq_y + val.sum_sq_y,
            self.n_samples + val.n_samples,
            self.sum_weighted_samples + val.sum_weighted_samples)

    def __isub__(self, StatsNode val):
        """In place subtract to StatsNodes together."""
        self.sum_y -= val.sum_y
        self.sum_sq_y -= val.sum_sq_y
        self.n_samples -= val.n_samples
        self.sum_weighted_samples -= val.sum_weighted_samples

        return self

    def __sub__(self, StatsNode val):
        """Subtract to StatsNodes together."""
        return StatsNode(
            self.sum_y - val.sum_y,
            self.sum_sq_y - val.sum_sq_y,
            self.n_samples - val.n_samples,
            self.sum_weighted_samples - val.sum_weighted_samples)

    # def __eq__(self, val):
    #     """Compare to StatsNode."""
    #     if isinstance(val, StatsNode):
    #         return all([
    #             self.sum_y == val.sum_y,
    #             self.sum_sq_y == val.sum_sq_y,
    #             self.n_samples == val.n_samples,
    #             self.sum_weighted_samples == val.sum_weighted_samples
    #         ])
    #     else:
    #         return NotImplemented

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

    cdef void copy_to(self, StatsNode dest_stats_node):
        dest_stats_node.sum_y = self.sum_y
        dest_stats_node.sum_sq_y = self.sum_sq_y
        dest_stats_node.n_samples = self.n_samples
        dest_stats_node.sum_weighted_samples = self.sum_weighted_samples
