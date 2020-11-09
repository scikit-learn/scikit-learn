# Python version for the base image
ARG PYTHON_VERSION

FROM winamd64/python:$PYTHON_VERSION-windowsservercore

# Run the container in user mode
ARG USERNAME=scikit-learn
ENV HOME /home/$USERNAME
RUN useradd -ms /bin/bash $USERNAME
USER $USERNAME

# Move to the workspace
RUN mkdir -p $HOME/workspace
WORKDIR $HOME/workspace

CMD latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf $FILE
