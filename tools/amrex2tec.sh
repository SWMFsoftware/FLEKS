#!/bin/bash

SCRIPT_DIR=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)
python $SCRIPT_DIR/amrex2tec.py "$@"
