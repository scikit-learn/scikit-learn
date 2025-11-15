"""Apriori algorithm for frequent itemsets and association rules."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import Counter
from itertools import chain, combinations
from typing import Dict, List

import numpy as np

try:
    import pandas as pd
except Exception:
    pd = None

from scipy.sparse import issparse

from sklearn.base import BaseEstimator, TransformerMixin, _fit_context
from sklearn.utils.validation import check_array


class Apriori(TransformerMixin, BaseEstimator):
    """
    Simple Apriori implementation for frequent itemsets and rules.
    """

    # Parameter constraints required by estimator checks
    _parameter_constraints = {
        "min_support": ["numeric"],
        "max_len": ["integer", "none"],
        "use_colnames": ["boolean"],
    }

    def __init__(self, *, min_support=0.5, max_len=None, use_colnames=True):
        self.min_support = min_support
        self.max_len = max_len
        self.use_colnames = use_colnames

    def _transactions_from_X(self, X) -> List[set]:
        """Convert X into a list of transaction sets."""
        # Case: DataFrame input (booleans or numeric one-hot)
        if pd is not None and isinstance(X, pd.DataFrame):
            valid_types = (
                (X.dtypes == bool).all()
                or np.issubdtype(X.values.dtype, np.integer)
                or np.issubdtype(X.values.dtype, np.floating)
            )
            if valid_types:
                # Use column names or auto names
                if self.use_colnames:
                    cols = list(X.columns)
                else:
                    cols = [f"item_{i}" for i in range(X.shape[1])]

                rows = []
                arr = X.to_numpy()

                # Convert rows into sets of present items
                for row in arr:
                    rows.append({cols[j] for j, v in enumerate(row) if bool(v)})
                return rows

        # Case: numpy 2D array
        if isinstance(X, np.ndarray) and X.ndim == 2:
            cols = [f"item_{i}" for i in range(X.shape[1])]
            rows = []
            for row in X:
                rows.append({cols[j] for j, v in enumerate(row) if bool(v)})
            return rows

        # Case: already list of iterables
        return [set(t) for t in X]

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        """Mine all frequent itemsets using Apriori."""
        # Validate common array-like inputs. For sparse matrices we intentionally
        # call check_array(..., accept_sparse=False) and allow it to raise so the
        # error message clearly states that sparse input is not supported.
        if pd is not None and isinstance(X, pd.DataFrame):
            check_array(X.to_numpy(), dtype=None, accept_sparse=False)
        elif isinstance(X, np.ndarray):
            check_array(X, dtype=None, accept_sparse=False)
        else:
            # If X is a scipy sparse matrix, raise the informative error.
            if issparse(X):
                check_array(X, dtype=None, accept_sparse=False)

        # Normalize transactions
        transactions = self._transactions_from_X(X)
        n_transactions = len(transactions)
        if n_transactions == 0:
            raise ValueError("No transactions provided")

        self.n_transactions_ = n_transactions

        # Count item frequencies (L1)
        item_counter = Counter(chain.from_iterable(transactions))
        L1 = {
            frozenset([item])
            for item, cnt in item_counter.items()
            if cnt / n_transactions >= self.min_support
        }

        # Support map for all itemsets
        support_map: Dict[frozenset, float] = {
            frozenset([item]): item_counter[item] / n_transactions
            for item in item_counter
        }

        frequent_itemsets = {}

        # Store frequent 1-itemsets
        for it in L1:
            frequent_itemsets[it] = support_map[it]

        current_L = L1
        k = 1

        # Apriori iterative expansion
        while current_L:
            k += 1
            # Stop if > max length
            if self.max_len is not None and k > self.max_len:
                break

            # Join step: generate candidates
            candidates = set()
            current_L_list = [tuple(sorted(s)) for s in current_L]

            for i in range(len(current_L_list)):
                for j in range(i + 1, len(current_L_list)):
                    union = frozenset(set(current_L_list[i]) | set(current_L_list[j]))
                    if len(union) == k:
                        # Prune step: all subsets must be frequent
                        subsets_ok = all(
                            frozenset(s) in current_L
                            for s in combinations(sorted(union), k - 1)
                        )
                        if subsets_ok:
                            candidates.add(union)

            if not candidates:
                break

            # Count support of candidates
            counts = Counter()
            for t in transactions:
                for c in candidates:
                    if c.issubset(t):
                        counts[c] += 1

            # Next level of frequent itemsets
            next_L = {
                c
                for c, cnt in counts.items()
                if cnt / n_transactions >= self.min_support
            }

            # Store supports
            for c in next_L:
                support = counts[c] / n_transactions
                support_map[c] = support
                frequent_itemsets[c] = support

            current_L = next_L

        # Save results
        self.frequent_itemsets_ = frequent_itemsets
        self.support_map_ = support_map
        return self

    def generate_rules(self, *, metric="confidence", min_threshold=0.7):
        """Generate association rules from discovered itemsets."""
        if not hasattr(self, "frequent_itemsets_"):
            raise ValueError("Fit before generating rules.")

        support_map = self.support_map_

        rules = []
        # Loop through all frequent itemsets
        for itemset, supp_itemset in support_map.items():
            if len(itemset) < 2:
                continue

            items = list(itemset)

            # Consider all possible antecedent splits
            for r in range(1, len(itemset)):
                for antecedent in combinations(items, r):
                    A = frozenset(antecedent)
                    B = itemset - A

                    supp_A = support_map.get(A)
                    supp_B = support_map.get(B)

                    if supp_A is None:
                        continue

                    # Compute rule metrics
                    confidence = supp_itemset / supp_A if supp_A > 0 else 0
                    lift = (
                        confidence / supp_B
                        if (supp_B is not None and supp_B > 0)
                        else np.inf
                    )
                    leverage = supp_itemset - (supp_A * (supp_B or 0))

                    metric_value = {
                        "confidence": confidence,
                        "lift": lift,
                        "leverage": leverage,
                    }[metric]

                    # Keep rule if metric threshold met
                    if metric_value >= min_threshold:
                        rules.append(
                            {
                                "antecedent": tuple(sorted(A)),
                                "consequent": tuple(sorted(B)),
                                "support": supp_itemset,
                                "confidence": confidence,
                                "lift": lift,
                                "leverage": leverage,
                            }
                        )

        # Return DataFrame if pandas available
        if pd is not None:
            return pd.DataFrame(rules)
        return rules
