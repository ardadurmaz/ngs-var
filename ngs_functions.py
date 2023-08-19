#!/usr/bin/python
######################################################
import sys
import os
import argparse
import os.path
import re
import glob
######################################################

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)
import ngs_classes

def print_log(inputs, text):
    with open(inputs.log, 'a') as f:
        f.write("\n" + text)

def run_command(inputs, command, error_msg = "[ERROR] Command failed"):
    if inputs.verbose:
        print_log(inputs, f"[INFO:COMMAND] {command}")
    if not inputs.dry:
        ret_val = os.system(command)
        if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept(error_msg)
