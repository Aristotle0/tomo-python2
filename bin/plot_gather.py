#!/usr/bin/env python
from tomopy.local.param import Fd2dParam, get_numpt, get_sta_coord
from tomopy.local.ioseism import read_seism, gather_seism
from tomopy.local.utility import read_option
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':

    help_string = """
This program plots a gather of seismograms.
Directory and component need to be specified.
The source can be specified, and will be set corresponding to
iteration if not.

--help              help information
--src=12            specify source number
--dir1=seism_syn    specify the directory that seismograms are saved, syn_seism by default
--dir2=seism_obs    specify the directory that seismograms are saved, plot only dir1 if not
--path=data         the current directory by default
--comp=x            specify which compontent to display, x by default

for instance:
plot_gather --src=200 --dir1=seism_syn --comp=x
plot_gather --src=200 --dir1=seism_syn --dir2=seism_obs --comp=x
"""
    option_dict = read_option(sys, help_string, 2, 6)

    working_path = option_dict.setdefault('path', '.')
    nsrc = int(option_dict.setdefault('src', 1))
    pnm_seism1 = option_dict.setdefault('dir1', 'seism_syn')
    pnm_seism2 = option_dict.setdefault('dir2', None)
    comp = option_dict.setdefault('comp', 'x')

    if comp == 'x':
        id_comp = 0
    elif comp == 'z':
        id_comp = 1

    if pnm_seism2:
        is_two_flag = True
    else:
        is_two_flag = False

    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nt = para.nt
    stept = para.stept
    seism_gather1, coord_gather = gather_seism(working_path, pnm_seism1, id_comp, dim1, dim2, nsrc, nt)
    if is_two_flag:
        seism_gather2, coord_gather = gather_seism(working_path, pnm_seism2, id_comp, dim1, dim2, nsrc, nt)
    time = np.arange(0, nt*stept, stept)

    smax = np.amax(np.abs(seism_gather1))
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    for nsta in range(coord_gather.shape[0]):
        p1, = ax1.plot(time, seism_gather1[nsta, :]/smax + coord_gather[nsta], 'b')
        if is_two_flag:
            p2, = ax1.plot(time, seism_gather2[nsta, :]/smax + coord_gather[nsta], 'r')

    ax1.set_xlabel("time (s)")
    ax1.set_ylabel("location of receivers (km)")
    ax1.set_xlim([time.min(), time.max()])
    if is_two_flag:
        ax1.legend((p1, p2), (pnm_seism1, pnm_seism2))
    else:
        ax1.legend((p1,), (pnm_seism1,))
    plt.title('the %s component of seismograms for the source %i' % (comp, nsrc))
    plt.show()
