"""Adaptive strategies for Polynomial Chaos Expansions.

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

In Polynomial Chaos expansions, the choice of the basis is crucial, as it
determines the accuracy of the resulting estimator. However, for problems with
a large number of features, the required number of basis terms to reach a
certain accuracy may be prohibitively expensive. When the number of basis terms
is larger than the number of samples (i.e., there are more unknowns than
equations) the system is underdetermined. In that case, it is beneficial to use
linear solvers that promote sparsity in the solution (such as `Lasso`,
`Lars`...) as a way to regularize the problem. Since the required number of
basis terms is unknown a priori, it can be beneficial to use adaptive basis
incremental strategies. In these approaches, additional basis terms are added
to the orthogonal polynomial basis in an incremental fashion, according to a
certain refinement criterium.

This file contains classes and files that can be used to implement adaptive
basis growth strategies for Polynomial Chaos expansions. The main (abstract)
class is `BasisIncrementStrategy`. Concrete basis increment strategies can be
added by inheriting from this base class. One adaptive basis growth strategy is
provided as `GerstnerGriebel`, which implements the basis growth strategy by
Gerstner & Griebel, originally targeting sparse grids, outlined in their 2020
paper.

The main method is `propose`. Given a `PolynomialChaosRegression`, this method
proposes a new set of multiindices, to be used as the basis terms in the next
Polynomial Chaos fitting procedure.

New incremental basis strategies can easily be added by inheriting from the
abstract `BasisIncrementStrategy` class and overriding the `propose` method.
"""

# import statements
from abc import ABC, abstractmethod  # for abstract classes

import numpy as np


class BasisIncrementStrategy(ABC):
    """An abstract base class for incremental Polynomial Chaos basis growth
    strategies."""

    def __init__(self):
        pass

    @abstractmethod
    def propose(self, pce):
        """Propose the next multiindex set given the current Polynomial Chaos
        expansion.

        Parameters
        ----------
        pce : PolynomialChaosRegressor
            The current Polynomial Chaos regressor

        Returns
        -------
        multiindices : array-like of shape (`n_terms`, `dimension`)
            The proposed multiindices.
        """

    @staticmethod
    def from_string(name):
        """Returns the adaptive basis growth strategy with the given name, if it
        exists.

        Parameters
        ----------
        name : str
            The name of the adaptive basis growth strategy (lower case and
            using underscore separators).

        Returns
        -------
        strategy : BasisIncrementStrategy
            The basis growth strategy with the given name. If no multiindex set
            with the given name can be found, a `ValueError` is raised.
        """
        name_ = "".join(word.capitalize() for word in name.split("_"))
        for strategy in BasisIncrementStrategy.__subclasses__():
            if strategy.__name__ == name_:
                return strategy()
        raise ValueError(f"unknown strategy type '{name}'")


class GerstnerGriebel(BasisIncrementStrategy):
    """Adaptive basis growth strategy from Gerstner & Griebel, 2010.

    In this algorithm, the set of multiindices are split into an 'old' and
    'active' set. The 'old' multiindex set contains all indices for which all
    forward neighbors (i.e., all multiindices that have one of the coordinates
    increased by 1) are included in the multiindex set. The remaining indices
    in the active set is composed of those indices for which at least one
    forward neighbor is not included in the multiindex set. At each iteration,
    (that is, each call to 'propose'), we select the multiindex from the active
    set for which the variance contribution (coefficient squared times
    polynomial norm squared) is the largest. Next, that index is moved from the
    'active' to the 'old' set, and all forward neighbors of that index are
    considered. If the forward neighbor is admissible in the old set (i.e., all
    of its backward neighbors are in the 'old' set), that index is added to the
    active set. This procedure is then repeated.
    """

    def propose(self, pce):
        # interpret multiindices as list of lists for convenience
        multiindices = pce.multiindices_.tolist()

        # split in active and old set if this is the first time this
        # method is called
        if not hasattr(self, "active"):
            self.active = list()
            self.old = list()
            for index in multiindices:
                # check if all forward neighbors are in the set
                is_old = True
                for d in range(pce.n_features_in_):
                    index_to_test = index.copy()
                    index_to_test[d] += 1
                    if index_to_test not in multiindices:
                        is_old = False
                        break
                # if so, add multiindex to the old set
                if is_old:
                    self.old.append(index)
                else:  # otherwise, add to active set
                    self.active.append(index)

        # find the multiindex with the maximum variance contribution
        max_contribution = 0
        best_index = self.active[0]
        for index in self.active:
            j = multiindices.index(index)
            contribution = (
                np.amax(np.atleast_2d(pce.coef_)[:, j]) ** 2 * pce.norms_[j] ** 2
            )
            if contribution > max_contribution:
                max_contribution = contribution
                best_index = index

        # move that index from active to old
        self.old.append(best_index)
        self.active.remove(best_index)

        # add all forward neighbors that are admissible in the old set to the
        # active set
        for d in range(pce.n_features_in_):
            new_index = best_index.copy()
            new_index[d] += 1
            # if not new_index in self.active:
            valid = True
            for j in range(pce.n_features_in_):
                backward_index = new_index.copy()
                backward_index[j] -= 1
                if backward_index[j] < 0:
                    continue
                if backward_index not in self.old:
                    valid = False
                    break
            if valid:
                self.active.append(new_index)

        # return new multiindex set
        return (*self.old, *self.active)
