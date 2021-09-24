#!/usr/bin/env bash
set -x
set -e

# Decide what kind of documentation build to run, and run it.
#
# If the last commit message has a "[doc skip]" marker, do not build
# the doc. On the contrary if a "[doc build]" marker is found, build the doc
# instead of relying on the subsequent rules.
#
# We always build the documentation for jobs that are not related to a specific
# PR (e.g. a merge to main or a maintenance branch).
#
# If this is a PR, do a full build if there are some files in this PR that are
# under the "doc/" or "examples/" folders, otherwise perform a quick build.
#
# If the inspection of the current commit fails for any reason, the default
# behavior is to quick build the documentation.

get_build_type() {
    if [ -z "$GITHUB_SHA" ]
    then
        echo SKIP: undefined GITHUB_SHA
        return
    fi
    commit_msg=$(git log --format=%B -n 1 $GITHUB_SHA)
    if [ -z "$commit_msg" ]
    then
        echo QUICK BUILD: failed to inspect commit $GITHUB_SHA
        return
    fi
    if [[ "$commit_msg" =~ \[doc\ skip\] ]]
    then
        echo SKIP: [doc skip] marker found
        return
    fi
    if [[ "$commit_msg" =~ \[doc\ quick\] ]]
    then
        echo QUICK: [doc quick] marker found
        return
    fi
    if [[ "$commit_msg" =~ \[doc\ build\] ]]
    then
        echo BUILD: [doc build] marker found
        return
    fi
    if [[ "$GITHUB_EVENT_NAME" != pull_request ]]
    then
        echo BUILD: not a pull request
        return
    fi
    git_range="origin/main...$GITHUB_SHA"
    git fetch origin main >&2 || (echo QUICK BUILD: failed to get changed filenames for $git_range; return)
    filenames=$(git diff --name-only $git_range)
    if [ -z "$filenames" ]
    then
        echo QUICK BUILD: no changed filenames for $git_range
        return
    fi
    changed_examples=$(echo "$filenames" | grep -E "^examples/(.*/)*plot_")

    # The following is used to extract the list of filenames of example python
    # files that sphinx-gallery needs to run to generate png files used as
    # figures or images in the .rst files  from the documentation.
    # If the contributor changes a .rst file in a PR we need to run all
    # the examples mentioned in that file to get sphinx build the
    # documentation without generating spurious warnings related to missing
    # png files.

    if [[ -n "$filenames" ]]
    then
        # get rst files
        rst_files="$(echo "$filenames" | grep -E "rst$")"

        # get lines with figure or images
        img_fig_lines="$(echo "$rst_files" | xargs grep -shE "(figure|image)::")"

        # get only auto_examples
        auto_example_files="$(echo "$img_fig_lines" | grep auto_examples | awk -F "/" '{print $NF}')"

        # remove "sphx_glr_" from path and accept replace _(\d\d\d|thumb).png with .py
        scripts_names="$(echo "$auto_example_files" | sed 's/sphx_glr_//' | sed -E 's/_([[:digit:]][[:digit:]][[:digit:]]|thumb).png/.py/')"

        # get unique values
        examples_in_rst="$(echo "$scripts_names" | uniq )"
    fi

    # executed only if there are examples in the modified rst files
    if [[ -n "$examples_in_rst" ]]
    then
        if [[ -n "$changed_examples" ]]
        then
            changed_examples="$changed_examples|$examples_in_rst"
        else
            changed_examples="$examples_in_rst"
        fi
    fi

    if [[ -n "$changed_examples" ]]
    then
        echo BUILD: detected examples/ filename modified in $git_range: $changed_examples
        pattern=$(echo "$changed_examples" | paste -sd '|')
        # pattern for examples to run is the last line of output
        echo "$pattern"
        return
    fi
    echo QUICK BUILD: no examples/ filename modified in $git_range:
    echo "$filenames"
}

build_type=$(get_build_type)
if [[ "$build_type" =~ ^SKIP ]]
then
    exit 0
fi

if [[ "$GITHUB_REF" =~ ^main$|^[0-9]+\.[0-9]+\.X$ && "$GITHUB_EVENT_NAME" != pull_request ]]
then
    # ZIP linked into HTML
    make_args=dist
elif [[ "$build_type" =~ ^QUICK ]]
then
    make_args=html-noplot
elif [[ "$build_type" =~ ^'BUILD: detected examples' ]]
then
    # pattern for examples to run is the last line of output
    pattern=$(echo "$build_type" | tail -n 1)
    make_args="html EXAMPLES_PATTERN=$pattern"
else
    make_args=html
fi

make_args="SPHINXOPTS=-T $make_args"  # show full traceback on exception

# Installing required system packages to support the rendering of math
# notation in the HTML documentation and to optimize the image files
sudo -E apt-get -yq update --allow-releaseinfo-change
sudo -E apt-get -yq --no-install-suggests --no-install-recommends \
    install dvipng gsfonts ccache zip optipng

# imports get_dep
source build_tools/shared.sh

pip install "$(get_dep numpy $NUMPY_VERSION)" \
            "$(get_dep scipy $SCIPY_VERSION)" \
            "$(get_dep cython $CYTHON_VERSION)" \
            "$(get_dep matplotlib $MATPLOTLIB_VERSION)" \
            "$(get_dep sphinx $SPHINX_VERSION)" \
            "$(get_dep pandas $PANDAS_VERSION)" \
            "$(get_dep scikit-image $SCIKIT_IMAGE_VERSION)"  \
            "$(get_dep sphinx-gallery $SPHINX_GALLERY_VERSION)" \
            "$(get_dep numpydoc $NUMPYDOC_VERSION)" \
            "$(get_dep sphinx-prompt $SPHINX_PROMPT_VERSION)" \
            "$(get_dep sphinxext-opengraph $SPHINXEXT_OPENGRAPH_VERSION)" \
            joblib memory_profiler seaborn pillow pytest coverage

# Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
# workers with 2 cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=3
python setup.py develop

export OMP_NUM_THREADS=1

if [[ "$GITHUB_SHA" =~ ^main$ && "$GITHUB_EVENT_NAME" == pull_request ]]
then
    # List available documentation versions if on main
    python build_tools/circle/list_versions.py > doc/versions.rst
fi

# The pipefail is requested to propagate exit code
cd doc && make $make_args

# Insert the version warning for deployment
find _build/html/stable -name "*.html" | xargs sed -i '/<\/body>/ i \
\    <script src="https://scikit-learn.org/versionwarning.js"></script>'
