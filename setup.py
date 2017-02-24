#!/usr/bin/env python

"""
setup.py for pyReefCore
"""
from numpy.distutils.core import setup, Extension

ext_modules = []

setup(
    name="pyReefCore",
    version="0.1",
    author="Tristan Salles",
    author_email="",
    description=("Synthetic carbonate platform core model"),
    long_description=open('README.md').read(),
    classifiers=[
        "Development Status :: 1 - Alpha",
    ],
    packages=['pyReefCore'],
    ext_package='pyReefCore',
    ext_modules=ext_modules,
    scripts=[],
)
