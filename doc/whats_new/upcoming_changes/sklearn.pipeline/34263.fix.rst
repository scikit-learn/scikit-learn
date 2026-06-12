- The default value for the `transform_input` parameter of :class:`Pipeline` was changed
  from `None` to `("X_val",)` so that the validation set, when passed to `fit`, is
  always transformed alongside `X`, to prevent easy to miss mistakes.
  By :user:`Jérémie du Boisberranger <jeremiedbb>`.
