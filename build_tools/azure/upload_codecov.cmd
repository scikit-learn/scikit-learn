@echo on

@rem Only 64 bit uses conda
IF "%PYTHON_ARCH%"=="64" (
    call activate %VIRTUALENV%
)

copy %TMP_FOLDER%\.coverage %BUILD_REPOSITORY_LOCALPATH%

codecov --root %BUILD_REPOSITORY_LOCALPATH% -t %CODECOV_TOKEN%
