#!/bin/bash

SCRIPT_DIR=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

for f in "$@"
do
stdbuf -oL pvpython $SCRIPT_DIR/converter.py $f
done
