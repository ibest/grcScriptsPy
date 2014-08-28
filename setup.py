#! /usr/bin/env python

# Copyright 2014, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
""" Setup for grcScriptsPy project. """

import os
import sys
from os import listdir as os_listdir
from os.path import join as path_join

from setuptools import setup
from setuptools import Extension
from distutils.sysconfig import get_config_vars
    
# Setup versioneer
import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'grcScriptsPy/_version.py'
versioneer.versionfile_build = 'grcScriptsPy/_version.py'
versioneer.tag_prefix = ''  # tags are like 1.2.0
versioneer.parentdir_prefix = '.'

CMDCLASS = versioneer.get_cmdclass()

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
OPT = get_config_vars('OPT')[0]
os.environ['OPT'] = " ".join(
    flag for flag in OPT.split() if flag != '-Wstrict-prototypes'
)

BUILD_DEPENDS = []
#BUILD_DEPENDS.extend(path_join("src", bn + ".h") for bn in [
#    "editdist"])

SOURCES = []

#SOURCES.extend(path_join("src", bn + ".c") for bn in [
#   "editdist"])
#    ])
SOURCES.extend(path_join("src", bn + ".cc") for bn in [
   "_grcScriptsPymodule", "editdist"])

EXTRA_COMPILE_ARGS = ['-O3']

#if sys.platform == 'darwin':
#    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64'])  # force 64bit only builds

EXTENSION_MOD_DICT = \
    {
        "sources": SOURCES,
        "depends": BUILD_DEPENDS,
        "extra_compile_args": EXTRA_COMPILE_ARGS,
        "language": "c++",
        "define_macros": [("VERSION", versioneer.get_version()), ],
    }

EXTENSION_MOD = Extension("_grcScripts",  # pylint: disable=W0142
                          ** EXTENSION_MOD_DICT)

SCRIPTS = []
SCRIPTS.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

INSTALL_REQUIRES = []

try:
    import argparse
    del argparse
except ImportError:
    INSTALL_REQUIRES.append("argparse >= 1.2.1")

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

SETUP_METADATA = \
    {
        "name": "grcScriptsPy",
        "version": versioneer.get_version(),
        "description": 'University of Idaho, Genomics Resources Core Python Bioinformatic Utilities',
        "long_description": read('README.md'),
        "author": 'Sam Hunter and Matthew L. Settles',
        "author_email": 'msettles@uidaho.edu',
        "url": 'https://github.com/ibest/grcScriptsPy',
        "packages": ['grcScriptsPy'],
        "install_requires": INSTALL_REQUIRES,
        "extras_require": {'docs': ['sphinx', 'sphinxcontrib-autoprogram'],
                           'tests': ['nose >= 1.0']},
        "scripts": SCRIPTS,
        "ext_modules": [EXTENSION_MOD, ],
#        "include_package_data": True,
        "zip_safe": False,
        "classifiers": [
            "Development Status :: 1 - Planning",
            "Environment :: Console",
            "Environment :: MacOS X",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: Apache Software License",
            "Natural Language :: English",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS :: MacOS X",
            "Programming Language :: C",
            "Programming Language :: C++",
            "Programming Language :: Python :: 2.7",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    }

# pylint: disable=W0142
setup(cmdclass=CMDCLASS,
    **SETUP_METADATA)
