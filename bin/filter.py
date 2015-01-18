#!/usr/bin/env python
from tomopy.local.filterseism import seism_filter
from tomopy.local.utility import read_option, get_gnsrc
import sys


if __name__ == '__main__':

    help_string = """
This program perform a filter to waveform.
The source can be specified, and will be set corresponding to
iteration if not. A working path can be specified.

--help              help information
--path=data         default as the current directory

for instance:
filter.py
"""
    option_dict = read_option(sys, help_string, 1, 3)

    working_path = option_dict.setdefault('path', '.')


    seism_filter(working_path)
