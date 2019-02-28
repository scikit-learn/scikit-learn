call activate %VIRTUALENV%

mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

if "%CHECK_WARNINGS%" == "true" (
    pytest --junitxml=%JUNITXML% --showlocals --durations=20 -Werror::DeprecationWarning -Werror::FutureWarning --pyargs sklearn
) else (
    pytest --junitxml=%JUNITXML% --showlocals --durations=20 --pyargs sklearn
)
