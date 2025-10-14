
### RareCategoryGrouper
-------------------
A tiny, self-contained example transformer that groups infrequent categories
into a single label (e.g., "Other"). Useful to reduce cardinality before
OneHotEncoding.


Usage
-----
>>> import pandas as pd
>>> from rare_category_grouper import RareCategoryGrouper
>>> df = pd.DataFrame({"city": ["SF","SF","NY","LA","LA","LA","SEA"]})
>>> rcg = RareCategoryGrouper(min_freq=2, columns=["city"], other_label="Other")
>>> rcg.fit(df)
RareCategoryGrouper(...)
>>> rcg.transform(df)
   city
0    SF
1    SF
2    Other
3     LA
4     LA
5     LA
6   Other

Notes
-----
- Input must be a pandas DataFrame. (Kept simple for teaching.)
- You can choose frequency threshold either by absolute count (min_freq)
  or by proportion (min_prop). If both are set, they are combined: keep
  a category if it meets *either* threshold (logical OR).
- Unseen categories at transform time are mapped to `other_label`.

