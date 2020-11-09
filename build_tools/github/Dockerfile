ARG PYTHON_VERSION

# Get the Python version of the base image from a build argument
FROM winamd64/python:$PYTHON_VERSION-windowsservercore

# Copy the Windows wheel to a temporary directory
ADD sklearn/setup.py /tmp/
CMD ls -l /tmp/

RUN python -c "import sklearn" 
