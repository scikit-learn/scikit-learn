FROM ubuntu

RUN apt-get update
RUN apt-get install -y -q \
        python-dev \
        python-setuptools \ 
        make \
        pkg-config \
        libfreetype6 \
        libfreetype6-dev \
        libpng12-dev \
        g++ \
        ipython
RUN easy_install pip
RUN pip install \ 
    'scipy>=0.9' \
    'nose>=1.1.2' \
    'matplotlib>=1.1.1' \
    'numpy>=1.6.1' \
    'cython>=0.23'
