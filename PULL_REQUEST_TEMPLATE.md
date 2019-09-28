<!--
Thanks for contributing a pull request! Please ensure you have taken a look at
the contribution guidelines: https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md#pull-request-checklist
-->

#### Reference Issues/PRs
<!--
Example: Fixes #1234. See also #3456.
Please use keywords (e.g., Fixes) to create link to the issues or pull requests
you resolved, so that they will automatically be closed when your pull request
is merged. See https://github.com/blog/1506-closing-issues-via-pull-requests
-->


#### What does this implement/fix? Explain your changes.

This change adds the functionality of handling new values in the OrdinalEncoder class.
It implements the handle_unknown variable on construction as the OneHotEncoder does.
Returns -1 when a new class is seen by the encoder.

#### Any other comments?

Example:

```python
import pandas as pd

df = pd.DataFrame(['ola', 'k', 'ase'])
df_test = pd.DataFrame(['ase', 'k', 'ase'])
df_test_unknown = pd.DataFrame(['tu', 'k', 'ases'])

from sklearn.preprocessing import OrdinalEncoder

encoder = OrdinalEncoder(handle_unknown='ignore')
encoder.fit(df)

print(f'Good shape: {encoder.transform(df_test)}')
print(f'Bad shape: {encoder.transform(df_test_unknown)}')

```


<!--
Please be aware that we are a loose team of volunteers so patience is
necessary; assistance handling other issues is very welcome. We value
all user contributions, no matter how minor they are. If we are slow to
review, either the pull request needs some benchmarking, tinkering,
convincing, etc. or more likely the reviewers are simply busy. In either
case, we ask for your understanding during the review process.
For more information, see our FAQ on this topic:
http://scikit-learn.org/dev/faq.html#why-is-my-pull-request-not-getting-any-attention.

Thanks for contributing!
-->
