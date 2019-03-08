@echo on

call activate %VIRTUALENV%

copy %TMP_FOLDER%\.coverage %BUILD_REPOSITORY_LOCALPATH%

codecov --root %BUILD_REPOSITORY_LOCALPATH% -t %CODECOV_TOKEN%
