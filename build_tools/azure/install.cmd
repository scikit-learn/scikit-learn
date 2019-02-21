set PIP_INSTALL=pip install -q

python -m pip install -U pip
python --version
pip --version

@rem Install the build and runtime dependencies of the project.
python setup.py bdist_wheel bdist_wininst -b doc\logos\scikit-learn-logo.bmp

@rem Install the generated wheel package to test it
pip install --pre --no-index --find-links dist\ scikit-learn

if %errorlevel% neq 0 exit /b %errorlevel%
