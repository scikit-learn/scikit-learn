
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


@dataclass
class _ColumnRule:
    keep: set
    other_label: str


class RareCategoryGrouper(BaseEstimator, TransformerMixin):
    """
    Group infrequent categories into a single label per selected column.

    Parameters
    ----------
    min_freq : int, default=10
        Minimum absolute count required for a category to be kept.
    min_prop : float or None, default=None
        Minimum proportion (0,1] of rows required to keep a category.
        If None, only `min_freq` applies. If set, a category is kept if it
        meets `min_freq` OR `min_prop`.
    columns : list of str or None, default=None
        Columns to process. If None, all object or category dtype columns
        will be processed.
    other_label : str, default="Other"
        Label used to replace rare/unseen categories.

    Attributes
    ----------
    category_maps_ : dict[str, _ColumnRule]
        Per-column rule with set of categories to keep and the label to use.

    n_features_in_ : int
        Number of features seen during fit.
    feature_names_in_ : list[str]
        Column names seen during fit.

    Examples
    --------
    See module-level examples.
    """

    def __init__(
        self,
        min_freq: int = 10,
        min_prop: Optional[float] = None,
        columns: Optional[List[str]] = None,
        other_label: str = "Other",
    ) -> None:
        self.min_freq = min_freq
        self.min_prop = min_prop
        self.columns = columns
        self.other_label = other_label

    def _select_columns(self, X: pd.DataFrame) -> List[str]:
        if self.columns is not None:
            return list(self.columns)
  
        return [c for c in X.columns if pd.api.types.is_object_dtype(X[c]) or pd.api.types.is_categorical_dtype(X[c])]

    def _kept_categories(self, s: pd.Series) -> set:
        counts = s.value_counts(dropna=False)
        keep_by_freq = counts[counts >= self.min_freq].index

        if self.min_prop is None:
            return set(keep_by_freq.tolist())

        prop = counts / len(s)
        keep_by_prop = prop[prop >= float(self.min_prop)].index
        return set(pd.Index(keep_by_freq).union(keep_by_prop).tolist())

    def fit(self, X: pd.DataFrame, y: Optional[Iterable] = None) -> "RareCategoryGrouper":
        if not isinstance(X, pd.DataFrame):
            raise TypeError("RareCategoryGrouper expects a pandas DataFrame as input.")

        cols = self._select_columns(X)
        category_maps: Dict[str, _ColumnRule] = {}

        for col in cols:
            kept = self._kept_categories(X[col])
            category_maps[col] = _ColumnRule(keep=kept, other_label=self.other_label)

        self.category_maps_ = category_maps
        self.n_features_in_ = X.shape[1]
        self.feature_names_in_ = list(X.columns)
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        if not isinstance(X, pd.DataFrame):
            raise TypeError("RareCategoryGrouper expects a pandas DataFrame as input.")

        if not hasattr(self, "category_maps_"):
            raise AttributeError("RareCategoryGrouper is not fitted yet. Call `fit` before `transform`.")

        X_out = X.copy()

        for col, rule in self.category_maps_.items():
            if col not in X_out.columns:
                continue

            s = X_out[col]

            if pd.api.types.is_categorical_dtype(s):
                if rule.other_label not in s.cat.categories:
                    s = s.cat.add_categories([rule.other_label])

            s = s.where(s.isin(rule.keep), other=rule.other_label)

            X_out[col] = s

        return X_out

    def get_kept_categories(self) -> Dict[str, List]:
        """
        Return the kept categories per column (after fitting).

        Returns
        -------
        dict[str, list]
        """
        if not hasattr(self, "category_maps_"):
            raise AttributeError("RareCategoryGrouper is not fitted yet.")
        return {col: sorted(list(rule.keep)) for col, rule in self.category_maps_.items()}


if __name__ == "__main__":
    df_demo = pd.DataFrame(
        {
            "city": ["SF", "SF", "SF", "NY", "LA", "LA", "SEA", "AUS", "SEA", "PDX"],
            "age": [25, 31, 29, 40, 38, 23, 27, 33, 36, 45],  # numeric column should be untouched
        }
    )
    print("Original:\n", df_demo, "\n")

    rcg = RareCategoryGrouper(min_freq=2, min_prop=None, columns=["city"], other_label="Other")
    rcg.fit(df_demo)
    print("Kept categories:", rcg.get_kept_categories(), "\n")

    print("Transformed:\n", rcg.transform(df_demo))
