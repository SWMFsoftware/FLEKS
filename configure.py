#!/usr/bin/env python

from __future__ import print_function

import sys
import os

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for configure.py")

import argparse

def configure(argv):
    argv[0] = "configure" # So the help message print it

    f = open("Makefile.def","w")
    f.write("INSTALL_DIR = " + os.getcwd() + "\n")
    f.write("\n")

    f.close()

if __name__ == "__main__":
    configure(sys.argv)
