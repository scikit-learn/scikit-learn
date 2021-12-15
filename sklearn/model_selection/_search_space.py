"""
The :mod:`sklearn.model_selection._search_space` module includes tools to build
parameter spaces for model selection.
"""

# Author: Joel Nothman <joel.nothman@gmail.com>
# License: BSD 3 clause

from itertools import product
from typing import Mapping
from pprint import pformat

__all__ = ["GridFactory"]


class GridFactory:
    """Utility to construct powerful grid searches

    This tools allows you to set the parameters for searching on each
    individual estimator. This way, if you change the structure of a composite
    estimator (e.g. wrapping it in a Pipeline or TransformedTargetRegressor)
    you don't need to rework the parameter grid. With this tool, it is also
    easy to build searches over multiple grids, for example with an optional
    step or alternative estimators for a step in a Pipeline.

    Parameters
    ----------
    estimator_grids : Mapping[Estimator, Mapping[str, List[object]]]
        Initial settings. For each estimator, this stores a mapping from
        parameter name to a list of its candidate values.

    Examples
    --------
    >>> from sklearn.model_selection import GridFactory
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.decomposition import PCA
    >>> fac = GridFactory()
    >>> kbest = fac.set(SelectKBest(), k=[5, 10, 20])
    >>> pca = fac.set(PCA(), n_components=[5, 10])
    >>> lr = fac.set(LogisticRegression(), C=[.1, 1, 10])
    >>> pipe = fac.set(Pipeline([("reduce", None),
    ...                          ("clf", lr)]),
    ...                reduce=[kbest, pca])
    >>> fac.get_grid(pipe)
    [{'reduce': [SelectKBest()],
      'clf__C': [0.1, 1, 10],
      'reduce__k': [5, 10, 20]},
     {'reduce': [PCA()],
      'clf__C': [0.1, 1, 10],
      'reduce__n_components': [5, 10]}]

    We do not need to use get_grid directly, since GridSearchCV accepts
    a GridFactory directly.
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.datasets import make_classification
    >>> search = GridSearchCV(pipe, fac)
    >>> X, y = make_classification(random_state=0)
    >>> search.fit(X, y)
    GridSearchCV(estimator=Pipeline(steps=[('reduce', None),
                                           ('clf', LogisticRegression())]),
                 param_grid=GridFactory(...))
    """

    def __init__(self, estimator_grids=None):
        if estimator_grids is None:
            estimator_grids = {}
        self.estimator_grids = {**estimator_grids}

    @staticmethod
    def _update_grid(dest, src, prefix=None):
        # TODO: needs docs
        if src is None:
            return dest
        if prefix:
            src = [{prefix + k: v for k, v in d.items()} for d in src]
        out = []
        for d1, d2 in product(dest, src):
            out_d = d1.copy()
            out_d.update(d2)
            out.append(out_d)
        return out

    def _build_param_grid(self, estimator):
        grid = self.estimator_grids.get(estimator, {})
        if isinstance(grid, Mapping):
            grid = [grid]

        # handle estimator parameters having their own grids
        for param_name, value in estimator.get_params().items():
            if "__" not in param_name and hasattr(value, "get_params"):
                out = []
                value_grid = self._build_param_grid(value)
                for sub_grid in grid:
                    if param_name in sub_grid:
                        sub_grid = [sub_grid]
                    else:
                        sub_grid = self._update_grid(
                            [sub_grid], value_grid, param_name + "__"
                        )
                    out.extend(sub_grid)
                grid = out

        # handle grid values having their own grids
        out = []
        for out_d in grid:
            part = [out_d]
            for param_name, values in out_d.items():
                to_update = []
                no_sub_grid = []
                for v in values:
                    if hasattr(v, "get_params"):
                        sub_grid = self._build_param_grid(v)
                        if sub_grid is not None:
                            to_update.extend(
                                self._update_grid(
                                    [{param_name: [v]}], sub_grid, param_name + "__"
                                )
                            )
                            continue
                    no_sub_grid.append(v)

                if no_sub_grid:
                    to_update.append({param_name: no_sub_grid})

                part = self._update_grid(part, to_update)
            out.extend(part)

        if out == [{}]:
            return None

        return out

    def set(self, estimator, **grid):
        """Set the grid to search for the specified estimator

        Overwrites any previously set grid.

        Parameters
        ----------
        grid : dict (str -> list of values)
            Keyword arguments define the values to be searched for each specified
            parameter.

        Returns
        -------
        estimator
            Useful for chaining
        """
        self.estimator_grids[estimator] = grid
        return estimator

    def get_grid(self, estimator):
        """Determine the parameter grid for the given estimator

        Parameters
        ----------
        estimator : scikit-learn compatible estimator
        """
        out = self._build_param_grid(estimator)
        if out is None:
            return {}
        elif len(out) == 1:
            return out[0]
        return out

    def __repr__(self):
        clean_grids = dict(sorted(self.estimator_grids.items(), key=repr))
        return f"GridFactory({pformat(clean_grids)})"
