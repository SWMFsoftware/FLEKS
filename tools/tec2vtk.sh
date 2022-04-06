#!/bin/bash

SCRIPT_DIR=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)
pvpython $SCRIPT_DIR/converter.py "$@"
