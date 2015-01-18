from tomopy.local.utility import read_option, get_gnsrc
from tomopy.local.misfitseism import get_adjs
import sys


if __name__ == '__main__':
    help_string = """
This program is to compute the misfits and adjoint sources.
The type of misfit functionals can be specified.
--help              show help information
--type=l2           the specified misfit funcitonal, L2 by default
--path=data         the current directory by default
--comp=x            specify the component of waveform to calculate the misfits, x by default
--fmin=10           minimium frequency for CWT
--fmax=50           maximium frequency for CWT
--nf=50             number of frequency sampling

type of misfit functionals:

* --type=l2         L2 waveform difference
* --type=wavelet    Continuous Wavelet transform
* --type=cc         cross-correlation (undifined)

for instance:
misfit.py --help
misfit.py --type=l2
misfit.py --type=wavelet --fmin=20 --fmax=50 --nf=50
"""
    option_dict = read_option(sys, help_string, 1, 8)

    working_path = option_dict.setdefault('path', '.')
    mtype = option_dict.setdefault('type', 'l2')
    comp = option_dict.setdefault('comp', 'x')
    fmin = float(option_dict.setdefault('fmin', 0))
    fmax = float(option_dict.setdefault('fmax', 0))
    nf = float(option_dict.setdefault('nf', 0))

    get_adjs(working_path, mtype, comp, fmin, fmax, nf)





