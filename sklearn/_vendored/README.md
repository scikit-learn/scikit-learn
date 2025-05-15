The contents of this folder is managed by the `vendoring` tool.

Install the developer branch of https://github.com/pradyunsg/vendoring/ or
version 1.3+ once released. This is needed to in particular to keep the SBOM
file up-to-date.

Bump the version numbers in the `sklearn/_vendored/vendoring-requirements.txt`
to the latest published versions using (**to be run from the top-level project
folder**):

```
python -m vendoring update
```

Actually synchronized the vendored code under `sklearn/_vendored/` to match the
versions configured in the requirements file using:

```
python -m vendoring sync
```

In addition to this synchronizing the code, the vendoring tool maintains the
contents of `sklearn/_vendored/bom.cdx.json` up to date to make the vendored
dependencies explicit for software composition analysis (SCA) tools such as
such as [syft](https://github.com/anchore/syft).

Some files to ignore (such as `conftest.py` or this  for instance) are
configured under the `tool.vendoring` section of the top level
`pyproject.toml`.
