@echo on

@rem Only 64 bit uses conda and uses a python newer than 3.5
IF "%PYTHON_ARCH%"=="64" (
    call activate %VIRTUALENV%
)

mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

if "%PYTEST_XDIST%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% -n2
)

if "%CHECK_WARNINGS%" == "true" (
    REM numpy's 1.19.0's tostring() deprecation is ignored until scipy and joblib removes its usage
    set PYTEST_ARGS=%PYTEST_ARGS% -Werror::DeprecationWarning -Werror::FutureWarning -Wignore:tostring:DeprecationWarning
)

if "%COVERAGE%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% --cov sklearn
)

pytest --junitxml=%JUNITXML% --showlocals --durations=20 %PYTEST_ARGS% --pyargs sklearn
