from netCDF4 import Dataset
import shutil
from scipy.signal import medfilt2d
from scipy.ndimage import gaussian_filter


def kernel_export(pnm_tmp, ker_filter, dim1, dim2):
    """ export kernel as nc files in temp directory
    """
    for kname, kvalue in ker_filter.items():
        ni1 = ni2 = 0
        for n_i in range(dim1):
            nk1 = nk2 = 0
            for n_k in range(dim2):
                fnm = "%s/kernel_mpi%02i%02i.nc" % (pnm_tmp, n_i, n_k)
                fnc = Dataset(fnm, 'a', format='NETCDF3_CLASSIC')
                li = len(fnc.dimensions['I'])
                lk = len(fnc.dimensions['K'])
                if n_k == 0:
                    ni2 += li
                nk2 += lk
                block = fnc.variables[kname]
                block[:, :] = kvalue[nk1:nk2, ni1:ni2]
                fnc.close()
                nk1 += lk
            ni1 += li


def kernel_copy_file(pnm_tmp, pnm_kernel, dim1, dim2):
    """ copy kernel files in temp/ to kernel/
    """
    for n_i in range(dim1):
        for n_k in range(dim2):
            fnm_tmp = "%s/kernel_mpi%02i%02i.nc" % (pnm_tmp, n_i, n_k)
            fnm_kernel = "%s/kernel_mpi%02i%02i.nc" % (pnm_kernel, n_i, n_k)
            shutil.copyfile(fnm_tmp, fnm_kernel)

def filter2d(ker, mf, gf):
    """ perform median filter and gaussian filter on kernel
    """
    ret = ker[:, :]
    if mf:
        ret = medfilt2d(ret, mf)

    if gf:
        ret = gaussian_filter(ret, gf)
    return ret

