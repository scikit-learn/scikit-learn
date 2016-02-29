Please try to adhere to the guidelines below as much as possible when
submitting your issue.
- Verify that your issue is not being currently addressed by other
[issues](https://github.com/scikit-learn/scikit-learn/issues?q=)
or [pull requests](https://github.com/scikit-learn/scikit-learn/pulls?q=).
- If your issue is a usage question or does not potentially require
changes to the codebase to be solved, then
[StackOverflow](http://stackoverflow.com/questions/tagged/scikit-learn)
(using the`[scikit-learn]` tag) or our
[mailing list](https://lists.sourceforge.net/lists/listinfo/scikit-learn-general)
may be a better place to bring it up. For more information, see
[User Questions](http://scikit-learn.org/stable/support.html#user-questions).

If you are submitting a bug issue:
- Please include your operating system type and version number, as well
as your Python, scikit-learn, numpy, and scipy versions. This information
can be found by runnning the following code snippet:
```
import platform; print(platform.platform())
import sys; print("Python", sys.version)
import numpy; print("NumPy", numpy.__version__)
import scipy; print("SciPy", scipy.__version__)
import sklearn; print("Scikit-Learn", sklearn.__version__)
```
- Please be specific about what estimators and/or functions are involved
and the shape of the data, as appropriate; please include a
[reproducible](http://stackoverflow.com/help/mcve) code snippet
or link to a [gist](https://gist.github.com). If an exception is raised,
please provide the traceback.
- Please ensure all code snippets and error messages are formatted in
appropriate code blocks.
See ["Creating and highlighting code blocks"](https://help.github.com/articles/creating-and-highlighting-code-blocks).

If you are submitting an algorithm or feature request:
- Please verify that the algorithm fulfills our
[new algorithm requirements](http://scikit-learn.org/stable/faq.html#can-i-add-this-new-algorithm-that-i-or-someone-else-just-published).

Thanks for contributing! Please delete these guidelines before submitting
your issue.