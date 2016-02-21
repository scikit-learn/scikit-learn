set -x
set -e

# Introspect the commit to know whether or not we should skip building the
# documentation: a pull request that does not change any file in doc/ or
# examples/ folder should be skipped unless the "[doc: build]" is found the
# commit message.
BUILD_DOC=`python build_tools/circle/check_build_doc.py`
echo -e $BUILD_DOC
if [[ $BUILD_DOC == "SKIP:"* ]]; then
    touch ~/log.txt  # the "test" segment needs that file
    exit 0
fi

# Installing required system packages to support the rendering of match
# notation in the HTML documentation
sudo -E apt-get -yq update
sudo -E apt-get -yq remove texlive-binaries --purge
sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes \
    install dvipng texlive-latex-base texlive-latex-extra

# deactivate circleci virtualenv and setup a miniconda env instead
if [[ `type -t deactivate` ]]; then
  deactivate
fi

# Install dependencies with miniconda
pushd .
cd
mkdir -p download
cd download
echo "Cached in $HOME/download :"
ls -l
if [[ ! -f miniconda.sh ]]
then
   wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
   -O miniconda.sh
fi
chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda
cd ..
export PATH="$HOME/miniconda/bin:$PATH"
conda update --yes --quiet conda
popd

# Configure the conda environment and put it in the path using the
# provided versions
conda create -n testenv --yes --quiet python numpy scipy \
  cython nose coverage matplotlib sphinx pillow
source /home/ubuntu/miniconda/envs/testenv/bin/activate testenv

# Build and install scikit-learn in dev mode
python setup.py develop

# The pipefail is requested to propagate exit code
set -o pipefail && cd doc && make html 2>&1 | tee ~/log.txt
