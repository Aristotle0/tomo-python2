from setuptools import setup, find_packages
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
setup(
    name = "tomopy",
    version = "0.0.1",
    packages = find_packages(tomopy),

    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'netCDF4',
    ]
    scripts=[
        'bin/filter.py',
        'bin/misfit.py',
        'bin/plot_gather.py',
        'bin/plot_kernel.py',
        'bin/plot_media.py',
        'bin/upsampling.py',
    ]

    author = "Lei Pan",
    author_email = "panlei7@gmail.com",
    description = "An additional package for tomo7",
    long_description = read('README'),
    keywords = "tomo7, FWI, model",

)
