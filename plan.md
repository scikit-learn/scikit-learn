# Plan: Interaction Constraints for `sklearn.tree` (Depth-First + Best Splitter)

## Scope
- Implement `interaction_cst` support for tree growth only when:
  - builder is `DepthFirstTreeBuilder` (`max_leaf_nodes is None`);
  - splitter is `"best"` (dense and sparse best splitter paths via `node_split_best`).
- Reuse the existing constant-feature caching pattern in `sklearn/tree/_splitter.pyx`.
- Add tests for validation, unsupported configurations, and behavioral correctness.

## High-Level Design
- Represent constraints with 2 CSR-like mappings:
  - feature -> groups:
    - `feature_to_groups_indptr` (shape `n_features + 1`)
    - `feature_to_groups_indices` (group ids)
  - group -> features:
    - `group_to_features_indptr` (shape `n_groups + 1`)
    - `group_to_features_indices` (feature ids)
- Add per-node state propagation for constraints in the depth-first stack:
  - number of active groups inherited from the parent path;
  - number of forbidden features inherited from the parent path.
- Do **not** store per-node arrays of active groups or forbidden features.
  - Store those arrays once in shared splitter buffers.
  - Store only scalar counters in node state (`n_active_interaction_groups`, `n_forbidden_features`) to read valid prefixes/suffixes of shared buffers.
  - This relies on depth-first growth ordering, where siblings are handled with restored invariants and counters.
- Keep forbidden features cached at the **end** of `features` (analogous to constant features cached at the beginning).

## Implementation Steps

1. Python-side API and validation (`sklearn/tree/_classes.py`)
- Add estimator parameter `interaction_cst` (default `None`) in `BaseDecisionTree` and concrete tree estimators.
- Extend `_parameter_constraints` and constructor signatures/docstrings.
- Add a private validator (similar spirit to HGB `_check_interaction_cst`) that:
  - accepts shorthands `"pairwise"` and `"no_interactions"`;
  - validates indices and containers;
  - canonicalizes to a list of feature sets;
  - appends unlisted features as singleton groups;
  - builds CSR-like arrays for `feature->groups` and `groups->features`.
- Enforce unsupported configurations with clear errors:
  - `splitter != "best"` when `interaction_cst is not None`;
  - `max_leaf_nodes is not None` when `interaction_cst is not None`.

2. Splitter interface extension (`sklearn/tree/_splitter.pxd`, `sklearn/tree/_splitter.pyx`)
- Extend `Splitter` ctor/signature to accept optional interaction CSR arrays.
- Store CSR arrays in splitter fields (typed memoryviews), plus flags (`with_interaction_cst`).
- Add preallocated working buffers:
  - active group ids cache;
  - forbidden features cache;
  - temporary marker arrays for fast membership/union assembly.
- These are shared buffers (single storage per splitter), not per-node allocations.
- Update `__reduce__`/serialization path to include new splitter ctor args.

3. Depth-first parent/stack state (`sklearn/tree/_tree.pxd`, `sklearn/tree/_tree.pyx`)
- Extend `ParentInfo` with:
  - `n_active_interaction_groups`;
  - `n_forbidden_features`.
- Extend `StackRecord` with the same scalar counters only (no per-node arrays) and initialize root values:
  - active groups = all groups;
  - forbidden features = 0.
- Propagate these fields to both children in `DepthFirstTreeBuilder`.
- Keep `BestFirstTreeBuilder` untouched (blocked by validation guard in step 1).

4. Best-split loop changes for forbidden-feature cache (`sklearn/tree/_splitter.pyx::node_split_best`)
- Integrate forbidden-feature boundaries into feature sampling loop so candidate draws exclude:
  - known constants (front of `features`);
  - known forbidden features (end of `features`).
- Maintain stable invariants for sibling compatibility (as done today for constants), now also for forbidden features.
- Explicitly maintain prefix/suffix invariants so scalar counters in `ParentInfo` are sufficient to restore/read shared buffers for both children.
- After selecting best split feature, compute child-node interaction state:
  - child active groups = parent active groups that contain `best_split.feature`;
  - child allowed features = union of features in child active groups;
  - child forbidden features = complement of allowed features.
- Ensure overlapping groups are fully supported with the same semantics as HGB/LightGBM-style propagation.
- Cache child forbidden features at the end of `features` and update `parent_record` counts so both children inherit the same state.
- Keep random splitter path unchanged (unsupported by validation).

5. Tests (`sklearn/tree/tests/`)
- Add tests for interaction constraint validation (shape/type/range, invalid content).
- Add tests for unsupported configurations:
  - `splitter="random"` + `interaction_cst` raises;
  - `max_leaf_nodes` + `interaction_cst` raises.
- Add behavioral tests (regression and classification variants where appropriate):
  - path-level structural check: for every root->node split sequence, used features must remain compatible with at least one constraint group;
  - numerical no-interaction test adapted from HGB (`[{0}, {1}]` should suppress interaction term);
  - overlapping-groups test to ensure descendant allowed features are correctly narrowed.
- Update/extend splitter serialization test if splitter signature changes require it.

6. Verification pass
- Run targeted tests first:
  - `pytest sklearn/tree/tests/test_tree.py -k interaction_cst`
  - (or dedicated new test file if added)
- Run both dense and sparse variants in interaction tests where applicable.
- Run monotonic-related tree tests to catch regressions in shared split logic:
  - `pytest sklearn/tree/tests/test_monotonic_tree.py`

## Confirmed Decisions
1. Overlapping interaction groups are supported.
2. Sparse input support is included for `splitter="best"` because the core logic is shared in `node_split_best` for dense/sparse partitioners.
3. String shorthands are supported: `"pairwise"` and `"no_interactions"`.
