import shutil
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata

__all__ = ['import_inv', 'get_coord', 'exp_media', 'myinterp2d',
    'media_extend', 'media_exchange']

def import_inv(path):
    """ Import inversion results
    """
    fnc = Dataset(path + '/' + 'inversion.nc', 'r')
    lam = fnc.variables['lambda'][:, :]
    mu = fnc.variables['mu'][:, :]
    x = fnc.variables['coord_x'][:, :]
    z = fnc.variables['coord_z'][:, :]
    fnc.close()
    return (lam, mu, x, z)

def get_coord(path, n_i, n_k):
    """ Get the coordinate of grid
    """
    fnm = "%s/coord_mpi%02i%02i.nc" % (path, n_i, n_k)
    fnc = Dataset(fnm, 'r')
    x = fnc.variables['x'][:, :]
    z = fnc.variables['z'][:, :]
    return x, z


def exp_media(dst, src, lam, mu, n_i, n_k):
    """ update media input files after upsamplint and backup old files
    """
    fsrc = "%s/media_mpi%02i%02i.nc" % (src, n_i, n_k)
    fdst = "%s/media_mpi%02i%02i.nc" % (dst, n_i, n_k)
    shutil.copyfile(fsrc, fdst)
    fnc = Dataset(fsrc, 'a', format='NETCDF3_CLASSIC')
    lamid = fnc.variables['lambda']
    muid = fnc.variables['mu']
    lamid[:, :] = np.exp(lam)
    muid[:, :] = np.exp(mu)
    fnc.close()

def myinterp2d(x, y, z, xnew, ynew, method='linear'):
    """ make 2d interpolation according to the axis
    """
    x = np.ravel(x)
    y = np.ravel(y)
    z = np.ravel(z)
    znew = griddata((x, y), z, (xnew, ynew), method=method, fill_value=0.)
    return znew

def extend_equal(w, ni, nk):
    LenFD = 3
    ni1 = LenFD
    ni2 = ni + LenFD -1
    nk1 = LenFD
    nk2 = nk + LenFD -1
    nx1 = 0
    nx2 = ni + LenFD*2 -1
    nz1 = 0
    nz2 = nk + LenFD*2 -1
    for k in range(nz1, nz2+1):
        for n in range(1, LenFD+1):
            w[k, ni1-n] = w[k, ni1]
            w[k, ni2+n] = w[k, ni2]
    for i in range(nx1, nx2+1):
        for n in range(1, LenFD+1):
            w[nk1-n, i] = w[nk1, i]
            w[nk2+n, i] = w[nk2, i]
    return w

def media_extend(dim1, dim2, ni, nk):
    for i in range(dim1):
        for k in range(dim2):
            fnm = "input/media_mpi%02i%02i.nc" % (i, k)
            f = Dataset(fnm, 'a')
            lam = f.variables['lambda']
            mu = f.variables['mu']
            lam[:, :] = extend_equal(lam[:, :], ni, nk)
            mu[:, :] = extend_equal(mu[:, :], ni, nk)
            f.close()


def media_exchange(dim1, dim2, ni, nk):
    LenFD = 3
    ni1 = LenFD
    ni2 = ni + LenFD -1
    nk1 = LenFD
    nk2 = nk + LenFD -1
    nx1 = 0
    nx2 = ni + LenFD*2 -1
    nz1 = 0
    nz2 = nk + LenFD*2 -1

    x1_src = [nz1, nz2+1, ni1, ni1+LenFD]
    x2_src = [nz1, nz2+1, ni2-LenFD+1, ni2+1]
    z1_src = [nk1, nk1+LenFD, nx1, nx2+1]
    z2_src = [nk2-LenFD+1, nk2+1, nx1, nx2+1]

    x2_dst = [nz1, nz2+1, ni2+1, nx2+1]
    x1_dst = [nz1, nz2+1, nx1, ni1]
    z2_dst = [nk2+1, nz2+1, nx1, nx2+1]
    z1_dst = [nz1, nk1, nx1, nx2+1]
    for i in range(dim1):
        for k in range(dim2):
            fnm = "input/media_mpi%02i%02i.nc" % (i, k)
            f = Dataset(fnm, 'r')
            lam_src = f.variables['lambda'][:, :]
            mu_src = f.variables['mu'][:, :]
            f.close()
            # transmit left
            if i>0:
                exchange(i-1, k, lam_src, mu_src,
                    x1_src, x2_dst)
            # transmit right
            if i<dim1-1:
                exchange(i+1, k, lam_src, mu_src,
                    x2_src, x1_dst)
            # transmit down
            if k>0:
                exchange(i, k-1, lam_src, mu_src,
                    z1_src, z2_dst)
            if k<dim2-1:
                exchange(i, k+1, lam_src, mu_src,
                    z2_src, z1_dst)

def exchange(i, k, lam_src, mu_src, sl, dl):
    fnm = "input/media_mpi%02i%02i.nc" % (i, k)
    f = Dataset(fnm, 'a')
    lam_dst = f.variables['lambda']
    mu_dst = f.variables['mu']
    lam_dst[dl[0]:dl[1], dl[2]:dl[3]] = (
        lam_src[sl[0]:sl[1], sl[2]:sl[3]])
    mu_dst[dl[0]:dl[1], dl[2]:dl[3]] = (
        mu_src[sl[0]:sl[1], sl[2]:sl[3]])
    f.close()

