---
name: Bug report
about: Create a report to help us reproduce and correct the bug
title: ''
labels: 'Bug: triage'
assignees: ''

---

<!--
Before submitting a bug, please make sure the issue hasn't been already
addressed by searching through the past issues.
-->

#### Describe the bug
<!--
A clear and concise description of what the bug is.
-->

#### Steps/Code to Reproduce
<!--
Please add a minimal example that we can reproduce the error by running the
code. Be as succinct as possible, do not depend on external data. In short, we
are going to copy-paste your code and we expect to get the same
result as you.

Example:
```python
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.decomposition import LatentDirichletAllocation
docs = ["Help I have a bug" for i in range(1000)]
vectorizer = CountVectorizer(input=docs, analyzer='word')
lda_features = vectorizer.fit_transform(docs)
lda_model = LatentDirichletAllocation(
    n_topics=10,
    learning_method='online',
    evaluate_every=10,
    n_jobs=4,
)
model = lda_model.fit(lda_features)
```
If the code is too long, feel free to put it in a public gist and link
it in the issue: https://gist.github.com
-->

```
Sample code to reproduce the problem
```

#### Expected Results
<!-- Example: No error is thrown. Please paste or describe the expected results.-->

#### Actual Results
<!-- Please paste or specifically describe the actual output or traceback. -->

#### Versions
<!--
Please run the following snippet and paste the output below.
For scikit-learn >= 0.20:
import sklearn; sklearn.show_versions()
For scikit-learn < 0.20:
import platform; print(platform.platform())
import sys; print("Python", sys.version)
import numpy; print("NumPy", numpy.__version__)
import scipy; print("SciPy", scipy.__version__)
import sklearn; print("Scikit-Learn", sklearn.__version__)
import imblearn; print("Imbalanced-Learn", imblearn.__version__)
-->


<!-- Thanks for contributing! -->
