#!/bin/bash

set -e

"${SHELL}" <(curl -Ls micro.mamba.pm/install.sh) < /dev/null

# >>> mamba initialize >>>
# !! Contents within this block are managed by 'micromamba shell init' !!
export MAMBA_EXE='/home/vscode/.local/bin/micromamba';
export MAMBA_ROOT_PREFIX='/home/vscode/micromamba';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias micromamba="$MAMBA_EXE"  # Fallback on help from micromamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

micromamba shell init --shell bash

micromamba env create -f build_tools/circle/doc_environment.yml -n sklearn-dev --yes
# Auto-activate sklearn-dev in terminal
echo "micromamba activate sklearn-dev" >> $HOME/.bashrc

# Enables users to activate environment without having to specify the full path
# echo "envs_dirs:
#  - /home/codespace/micromamba/envs" > /opt/conda/.condarc
