@ECHO OFF

REM Command file for Sphinx documentation

set SPHINXOPTS = -v

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set BUILDDIR=_build
set ALLSPHINXOPTS=-d %BUILDDIR%/doctrees %SPHINXOPTS% .

if "%1" == "clean" (
	for /d %%i in (%BUILDDIR%\*) do rmdir /q /s %%i
	del /q /s %BUILDDIR%\*
	if exist auto_examples rd /q /s auto_examples
	if exist auto_plotly_examples rd /q /s auto_plotly_examples
	if exist auto_pyvista_examples rd /q /s auto_pyvista_examples
	if exist tutorials rd /q /s tutorials
	if exist gen_modules rd /q /s gen_modules
	goto end
)


if "%1" == "html" (
	%SPHINXBUILD% -b html %ALLSPHINXOPTS% %BUILDDIR%/html
	if errorlevel 1 exit /b 1
	echo.
	echo.Build finished. The HTML pages are in %BUILDDIR%/html.
	goto end
)


if "%1" == "latexpdf" (
	%SPHINXBUILD% -b latex %ALLSPHINXOPTS% %BUILDDIR%/latex
	cd %BUILDDIR%/latex
	make all-pdf
	cd %BUILDDIR%/..
	echo.
	echo.Build finished; the PDF files are in %BUILDDIR%/latex.
	goto end
)


:end
