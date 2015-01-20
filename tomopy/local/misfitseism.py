from ..signal.misfit_functionals import wavemisf, waveletmisf, ccmisf
from .param import Fd2dParam, get_numpt
from .ioseism import read_seism
import numpy as np
from netCDF4 import Dataset

def write_misf(misf, adjs, folder, gnsrc, n_i, n_k, nt, num_pt, nd):
    """ Write misfits

    Parameters
    ----------
    misf : ndarray
        (3, num_pt) misfit values
    ajds : ndarray
        (nt, 3, num_pt) adjoint source function
    folder : string
        Folder name where the file we want is
    gnsrc : int
        No. of the seismic source
    n_i, n_k : int
        Id of the MPI block
    nt : int
        Number of time samplings
    num_pt : int
        Number of receivers
    nd : int
        the component id
    """
    str1 = (folder, gnsrc, n_i, n_k)
    filenm = '%s/misf_s%03i_mpi%02i%02i.nc' % str1
    fnc = Dataset(filenm, 'w', format='NETCDF3_CLASSIC')
    fnc.createDimension('geo', 3)
    fnc.createDimension('num_pt', num_pt)
    fnc.createDimension('time', nt)
    misfid = fnc.createVariable('misf', 'f', ('geo', 'num_pt',))
    misfid[:, :] = misf
    adjsid = fnc.createVariable('adjs', 'f', ('time', 'geo', 'num_pt'))
    adjsid[:, :, :] = adjs
    fnc.close()



def get_adjs(working_path, mtype, comp, fmin, fmax, nf):
    """ calculate the misfits and the corresponding adjoint source.

    Parameters
    ----------
    working_path : string
        the current working path
    mtype : string
        l2 for L2 waveform difference misfits
        cc for cross-correlation misfits
        wavelet for CWT misfits
    comp : string
        the component of seismograms
    fmin : float
        minimium frequency for CWT
    fmax : float
        maximium frequency for CWT
    nf : int
        number of frequency sampling
    """


    para = Fd2dParam(working_path)
    dim1, dim2 = para.dim
    nt = para.nt
    num_src = para.number_of_moment_source
    stept = para.stept

    if mtype == 'l2':
        misf_func = lambda x, y : wavemisf(x, y, stept)
    elif mtype == 'wavelet':
        misf_func = lambda x, y : waveletmisf(x, y, fmin, fmax, nf)
    elif mtype == 'cc':
        misf_func = lambda x, y : ccmisf(x, y, stept)

    if comp == 'x':
        id_comp = 0
        id_misf = 0
    elif comp == 'z':
        id_comp = 1
        id_misf = 2

    pnm_seism = (para.pnm_syn_filter, para.pnm_obs_filter)
    pnm_seism = [working_path + '/' + x for x in pnm_seism]
    pnm_misf = working_path + '/' + 'misf'

    for nsrc in range(1, num_src+1):
        for n_i in range(dim1):
            for n_k in range(dim2):
                num_pt = get_numpt(n_i, n_k, working_path)
                if num_pt == 0:
                    continue
                seism_syn, time = read_seism(pnm_seism[0], nsrc, n_i, n_k, nt, num_pt)
                seism_obs, time = read_seism(pnm_seism[1], nsrc, n_i, n_k, nt, num_pt)
                misf = np.zeros([3, num_pt], dtype=np.float32)
                adjs = np.zeros([nt, 3, num_pt], dtype=np.float32, order='F')
                for nsta in range(num_pt):
                    misf[id_misf, nsta], adjs[:, id_misf, nsta] = \
                            misf_func(seism_obs[id_comp, :, nsta], seism_syn[id_comp, :, nsta])
                write_misf(misf, adjs, pnm_misf, nsrc, n_i, n_k, nt, num_pt, id_misf)


