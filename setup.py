from setuptools import setup, find_packages
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
setup(
    name = "tomopy",
    version = "0.0.1",
    packages = ['tomopy'],

    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'netCDF4',
    ],
    scripts=[
        'bin/filter_kernel.py',
        'bin/filter_seism.py',
        'bin/misfit.py',
        'bin/plot_gather.py',
        'bin/plot_kernel.py',
        'bin/plot_media.py',
        'bin/upsampling.py',
    ],

    author = "Lei Pan",
    author_email = "panlei7@gmail.com",
    description = "An additional package for tomo7",
    long_description = read('README.md'),
    keywords = "tomo7, FWI, model",

)
