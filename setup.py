## -*- encoding: utf-8 -*-
import os
import sys
from setuptools import setup
from codecs import open

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding="utf-8") as f:
        return f.read()

setup(
    name="cvolume",
    version=readfile("VERSION").strip(), # the VERSION file is shared with the documentation
    description="Completed volumes of strata of quadratic differentials with odd zeros",
    long_description=readfile("README.rst"), # get the long description from the README
    url="",
    author="Vincent Delecroix, Eduard Duryev",
    author_email="edwardduriev@gmail.com",
    license="GPLv2+",
    classifiers=[
      "Development Status :: 4 - Beta",
      "Intended Audience :: Science/Research",
      "Topic :: Software Development :: Build Tools",
      "Topic :: Scientific/Engineering :: Mathematics",
      "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
      "Programming Language :: Python :: 3.6",
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = "SageMath geometry moduli space curve differential",
    packages = ["cvolume"],
)
