@echo on

@rem Only 64 bit uses conda
IF "%PYTHON_ARCH%"=="64" (
    call activate %VIRTUALENV%
)

mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

if "%CHECK_WARNINGS%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% -Werror::DeprecationWarning -Werror::FutureWarning
)

if "%COVERAGE%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% --cov sklearn
)

pytest --junitxml=%JUNITXML% --showlocals --durations=20 %PYTEST_ARGS% --pyargs sklearn
