set -xe

PLATFORM_ID=$1
# Set environment variables to make our wheel build easier to reproduce byte
# for byte from source. See https://reproducible-builds.org/. The long term
# motivation would be to be able to detect supply chain attacks.
#
# In particular we set SOURCE_DATE_EPOCH to the commit date of the last commit.
#
# XXX: setting those environment variables is not enough. See the following
# issue for more details on what remains to do:
# https://github.com/scikit-learn/scikit-learn/issues/28151
echo "SOURCE_DATE_EPOCH=$(git log -1 --pretty=%ct)" >> "$GITHUB_ENV"
# TODO PYTHONHASHSEED needs to be passed into the container in the Linux case.
# SOURCE_DATE_EPOCH is always passed in. Maybe use both in CIBW_ENVIRONMENT for
# simplicity?
echo PYTHONHASHSEED=0 >> "$GITHUB_ENV"

# OpenMP is not present on macOS by default
if [[ $(uname) == "Darwin" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos.

    # TODO this is the problem here, CIBW_BUILD does not exist ...
    if [[ "$PLATFORM_ID" == macosx_arm64 ]]; then
        # SciPy requires 12.0 on arm to prevent kernel panics
        # https://github.com/scipy/scipy/issues/14688
        # We use the same deployment target to match SciPy.
        export MACOSX_DEPLOYMENT_TARGET=12.0
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-arm64/llvm-openmp-11.1.0-hf3c4609_1.tar.bz2"
    else
        export MACOSX_DEPLOYMENT_TARGET=10.9
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-64/llvm-openmp-11.1.0-hda6cdc1_1.tar.bz2"
    fi

    conda create -n build $OPENMP_URL
    PREFIX="$HOME/miniconda3/envs/build"

    echo CC=/usr/bin/clang >> "$GITHUB_ENV"
    echo CXX=/usr/bin/clang++ >> "$GITHUB_ENV"
    echo "CPPFLAGS=$CPPFLAGS -Xpreprocessor -fopenmp" >> "$GITHUB_ENV"
    echo "CFLAGS=$CFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo "CXXFLAGS=$CXXFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo "LDFLAGS=$LDFLAGS -Wl,-rpath,$PREFIX/lib -L$PREFIX/lib -lomp" >> "$GITHUB_ENV"
fi
