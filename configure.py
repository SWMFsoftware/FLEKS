#!/usr/bin/env python

from __future__ import print_function

import sys
import os

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for configure.py")

import argparse
import subprocess

def configure(argv):
    argv[0] = "configure" # So the help message print it

    # Check if we are inside SWMF or standalone
    if not os.path.exists("../../share") and not os.path.exists("share"):
        print("--- Cloning SWMFsoftware/share and util ---")
        GITDIR = "git@github.com:SWMFsoftware"
        subprocess.call(["git", "clone", GITDIR + "/share"])
        subprocess.call(["git", "clone", GITDIR + "/util"])

    f = open("Makefile.def","w")
    f.write("INSTALL_DIR = " + os.getcwd() + "\n")
    f.write("MYDIR = " + os.getcwd() + "\n")
    f.write("BINDIR = " + os.getcwd() + "/bin\n")
    f.write("LIBDIR = " + os.getcwd() + "/lib\n")
    
    if os.path.exists("../../share"):
        f.write("SHAREDIR = " + os.path.abspath("../../share") + "/Library/src\n")
        f.write("SCRIPTDIR = " + os.path.abspath("../../share") + "/Scripts\n")
    else:
        f.write("SHAREDIR = " + os.getcwd() + "/share/Library/src\n")
        f.write("SCRIPTDIR = " + os.getcwd() + "/share/Scripts\n")

    f.write("\n")
    f.close()
    
    print("Created Makefile.def")

if __name__ == "__main__":
    configure(sys.argv)
