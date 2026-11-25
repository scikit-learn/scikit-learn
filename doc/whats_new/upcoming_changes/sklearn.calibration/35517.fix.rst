- :class:`calibration.CalibratedClassifierCV` now assigns sigmoid and isotonic
  probabilities to both trained classes when a fold in a multiclass ensemble
  trains on only two classes. Unseen classes still receive zero probability.
  By :user:`Nishan Kashyap <Y2XNishan>`.
