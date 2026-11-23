- The parameter `multiclass_colors` was deprecated in favour of `target_colors` in
  :class:`inspection.DecisionBoundaryDisplay`. The attribute `multiclass_colors_` was
  also renamed to `target_colors_`. Now they can be used for binary problems as well
  without causing confusion (which will be added in a follow-up PR).
  By :user:`Anne Beyer <AnneBeyer>`.
