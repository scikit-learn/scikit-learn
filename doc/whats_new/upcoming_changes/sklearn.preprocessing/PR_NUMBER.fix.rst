- :meth:`preprocessing.OrdinalEncoder.inverse_transform` now correctly decodes
  features that have both infrequent categories and missing values. Previously,
  infrequent and missing values could be swapped or raise an `IndexError`.
  By :user:`digit50 <digit50>`.
