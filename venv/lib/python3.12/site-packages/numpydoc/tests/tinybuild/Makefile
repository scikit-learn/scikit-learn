all: clean html show

clean:
	rm -rf _build/*
	rm -rf generated/

html:
	sphinx-build -nWT --keep-going -b html -d _build/doctrees . _build/html

show:
	@python -c "import webbrowser; webbrowser.open_new_tab('file://$(PWD)/_build/html/index.html')"
