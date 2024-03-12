#!/usr/bin/env bash
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

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

if [ -n "$GITHUB_ACTION" ]
then
    # Map the variables from Github Action to CircleCI
    CIRCLE_SHA1=$(git log -1 --pretty=format:%H)

    CIRCLE_JOB=$GITHUB_JOB

    if [ "$GITHUB_EVENT_NAME" == "pull_request" ]
    then
        CIRCLE_BRANCH=$GITHUB_HEAD_REF
        CI_PULL_REQUEST=true
    else
        CIRCLE_BRANCH=$GITHUB_REF_NAME
    fi
fi

get_build_type() {
    if [ -z "$CIRCLE_SHA1" ]
    then
        echo SKIP: undefined CIRCLE_SHA1
        return
    fi
    commit_msg=$(git log --format=%B -n 1 $CIRCLE_SHA1)
    if [ -z "$commit_msg" ]
    then
        echo QUICK BUILD: failed to inspect commit $CIRCLE_SHA1
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
    if [ -z "$CI_PULL_REQUEST" ]
    then
        echo BUILD: not a pull request
        return
    fi
    git_range="origin/main...$CIRCLE_SHA1"
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

if [[ "$CIRCLE_BRANCH" =~ ^main$|^[0-9]+\.[0-9]+\.X$ && -z "$CI_PULL_REQUEST" ]]
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

# Installing required system packages to support the rendering of math
# notation in the HTML documentation and to optimize the image files
sudo -E apt-get -yq update --allow-releaseinfo-change
sudo -E apt-get -yq --no-install-suggests --no-install-recommends \
    install dvipng gsfonts ccache zip optipng

# deactivate circleci virtualenv and setup a conda env instead
if [[ `type -t deactivate` ]]; then
  deactivate
fi

MAMBAFORGE_PATH=$HOME/mambaforge
# Install dependencies with mamba
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    -O mambaforge.sh
chmod +x mambaforge.sh && ./mambaforge.sh -b -p $MAMBAFORGE_PATH
export PATH="/usr/lib/ccache:$MAMBAFORGE_PATH/bin:$PATH"

ccache -M 512M
export CCACHE_COMPRESS=1

# pin conda-lock to latest released version (needs manual update from time to time)
mamba install "$(get_dep conda-lock min)" -y

conda-lock install --log-level DEBUG --name $CONDA_ENV_NAME $LOCK_FILE
source activate $CONDA_ENV_NAME

show_installed_libraries

# Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
# workers with 2 cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=3
pip install -e . --no-build-isolation

echo "ccache build summary:"
ccache -s

export OMP_NUM_THREADS=1

if [[ "$CIRCLE_BRANCH" =~ ^main$ && -z "$CI_PULL_REQUEST" ]]
then
    # List available documentation versions if on main
    python build_tools/circle/list_versions.py > doc/versions.rst
fi


# The pipefail is requested to propagate exit code
set -o pipefail && cd doc && make $make_args 2>&1 | tee ~/log.txt

# Insert the version warning for deployment
find _build/html/stable -name "*.html" | xargs sed -i '/<\/body>/ i \
\    <script src="https://scikit-learn.org/versionwarning.js"></script>'

cd -
set +o pipefail

affected_doc_paths() {
    files=$(git diff --name-only origin/main...$CIRCLE_SHA1)
    echo "$files" | grep ^doc/.*\.rst | sed 's/^doc\/\(.*\)\.rst$/\1.html/'
    echo "$files" | grep ^examples/.*.py | sed 's/^\(.*\)\.py$/auto_\1.html/'
    sklearn_files=$(echo "$files" | grep '^sklearn/')
    if [ -n "$sklearn_files" ]
    then
        grep -hlR -f<(echo "$sklearn_files" | sed 's/^/scikit-learn\/blob\/[a-z0-9]*\//') doc/_build/html/stable/modules/generated | cut -d/ -f5-
    fi
}

affected_doc_warnings() {
    files=$(git diff --name-only origin/main...$CIRCLE_SHA1)
    # Look for sphinx warnings only in files affected by the PR
    if [ -n "$files" ]
    then
        for af in ${files[@]}
        do
          warn+=`grep WARNING ~/log.txt | grep $af`
        done
    fi
    echo "$warn"
}

if [ -n "$CI_PULL_REQUEST" ]
then
    echo "The following documentation warnings may have been generated by PR #$CI_PULL_REQUEST:"
    warnings=$(affected_doc_warnings)
    if [ -z "$warnings" ]
    then
        warnings="/home/circleci/project/ no warnings"
    fi
    echo "$warnings"

    echo "The following documentation files may have been changed by PR #$CI_PULL_REQUEST:"
    affected=$(affected_doc_paths)
    echo "$affected"
    (
    echo '<html><body><ul>'
    echo "$affected" | sed 's|.*|<li><a href="&">&</a> [<a href="https://scikit-learn.org/dev/&">dev</a>, <a href="https://scikit-learn.org/stable/&">stable</a>]</li>|'
    echo '</ul><p>General: <a href="index.html">Home</a> | <a href="modules/classes.html">API Reference</a> | <a href="auto_examples/index.html">Examples</a></p>'
    echo '<strong>Sphinx Warnings in affected files</strong><ul>'
    echo "$warnings" | sed 's/\/home\/circleci\/project\//<li>/g'
    echo '</ul></body></html>'
    ) > 'doc/_build/html/stable/_changed.html'

    if [ "$warnings" != "/home/circleci/project/ no warnings" ]
    then
        echo "Sphinx generated warnings when building the documentation related to files modified in this PR."
        echo "Please check doc/_build/html/stable/_changed.html"
        exit 1
    fi
fi
