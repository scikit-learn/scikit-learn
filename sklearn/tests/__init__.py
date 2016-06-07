'''
The __init__.py file makes Python treat directories containing it as module/package.
If test folders have a __init__.py then it is installed in site-packages so when you run nosetests sklearn 
outside of the source tree after doing python setup.py install you will find the tests that have been installed in site-packages.
Omiting the __init__.py before setup.py installation don't install the test!
'''
