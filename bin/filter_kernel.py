#!/usr/bin/env python
from tomopy.local.param import Fd2dParam, associate_blocks
from tomopy.local.utility import read_option
from tomopy.local.filterkernel import kernel_copy_file, kernel_export, filter2d
import subprocess
import sys

if __name__ == '__main__':

    help_string = """
This program is to filter kernel.

--help          show help information
--path=data     give the working path
--mf=5          specify the parameter of median filter
--gf=1.0        specify the parameter of gaussian filter
--plot=yes      plot the kernel and don't save the result
--comp=ka       specify the component of kernel to plot

for instance:
filter_kernel.py --mf=11
filter_kernel.py --gf=1.0
filter_kernel.py --mf==11 --gf=1.0
"""
    option_dict = read_option(sys, help_string, 1, 6)

    working_path = option_dict.setdefault('path', '.')
    plot_flag = option_dict.setdefault('plot', 'no')
    mf = int(option_dict.setdefault('mf', 11))
    gf = float(option_dict.setdefault('gf', 1.0))
    comp = option_dict.setdefault('comp', 'ka')

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim

    pnm_kernel = working_path + '/' + 'kernel'
    pnm_tmp = working_path + '/' + 'temp'
    ker_comp = ('ka', 'kb')
    ker = {x : associate_blocks(pnm_kernel, 'kernel', x, dim1, dim2) for x in ker_comp}

    ker_filter = { x : filter2d(ker[x], mf, gf) for x in ker_comp }

    kernel_copy_file(pnm_kernel, pnm_tmp, dim1, dim2)
    kernel_export(pnm_kernel, ker_filter, dim1, dim2)

    if plot_flag == 'yes':
        if comp == 'ka':
            subprocess.call("plot_kernel.py", shell=True)
        elif comp == 'kb':
            subprocess.call("plot_kernel.py --comp=kb", shell=True)
        kernel_copy_file(pnm_tmp, pnm_kernel, dim1, dim2)




