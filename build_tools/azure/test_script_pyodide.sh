#!/bin/bash

set -e

node build_tools/azure/pytest-pyodide.js --pyargs sklearn
