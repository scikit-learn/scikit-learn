# Changelog instructions

This directory (`doc/whats_new/upcoming_changes`) contains "news fragments",
which are short files that contain a small **ReST**-formatted text that will be
added to the next release changelog.

Each file should be named like `<PULL REQUEST>.<TYPE>.rst`, where
`<PULL REQUEST>` is a pull request number, and `<TYPE>` is one of:

* `major-feature`
* `feature`
* `efficiency`
* `enhancement`
* `fix`
* `api`
* `other` (see [](#custom-top-level-folder))

See [this](https://github.com/scikit-learn/scikit-learn/blob/main/doc/whats_new/changelog_legend.inc)
for more details about the meaning of each type.

This file needs to be added to the right folder like `sklearn.linear_model` or
`sklearn.tree` depending on which part of scikit-learn your PR changes. There
are also a few folders for some topics like `array-api`, `metadata-routing` or `security`.

In almost all cases, your fragment should be formatted as a bullet point.

For example, `28268.feature.rst` would be added to the `sklearn.ensemble`
folder with the following content::

```rst
- :class:`ensemble.ExtraTreesClassifier` and :class:`ensemble.ExtraTreesRegressor`
  now supports missing values in the data matrix `X`. Missing-values are
  handled by randomly moving all of the samples to the left, or right child
  node as the tree is traversed.
  By :user:`Adam Li <adam2392>`
```

If you are unsure how to name the news fragment or which folder to use, don't
hesitate to ask in your pull request!

You can install [`towncrier`](https://github.com/twisted/towncrier) and run
`towncrier create` to help you create a news fragment. You can also run
`towncrier build --draft --version <version_number>` if
you want to get a preview of how your change will look in the final release
notes.


## `custom-top-level` folder

The `custom-top-level` folder is for changes for which there is no good
folder and are somewhat one-off topics. Type `other` is mostly meant to be used
in the `custom-top-level` section.
