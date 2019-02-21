@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/setup_conda_environment.cmd
@rem The cmd /C hack circumvents a regression where conda installs a conda.bat
@rem script in non-root environments.
set CONDA_INSTALL=cmd /C conda install -q -y
set PIP_INSTALL=pip install -q

@echo on

@rem Deactivate any environment
call deactivate
@rem Display root environment (for debugging)
conda list
@rem Clean up any left-over from a previous build
conda remove --all -q -y -n %VIRTUALENV%
conda env create -n %VIRTUALENV% python=%CONDA_PY% -y

call activate %VIRTUALENV%
python -m pip install -U pip
python --version
pip --version

pip install --timeout=60 --trusted-host 28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com -r build_tools\appveyor\requirements.txt

@rem Install the build and runtime dependencies of the project.
python setup.py bdist_wheel bdist_wininst -b doc\logos\scikit-learn-logo.bmp

@rem Install the generated wheel package to test it
pip install --pre --no-index --find-links dist\ scikit-learn

if %errorlevel% neq 0 exit /b %errorlevel%
