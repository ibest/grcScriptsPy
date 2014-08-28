#!/usr/bin/env python

# This file is part of grcScriptsPy, http://github.com/ibest/grcScriptsPy/
# Copyright 2014, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License; see LICENSE.txt.
# Contact: msettles@uidaho.edu



def print_hello_world_py():
    print("Hello World, Py Version")


#### Misc functions

import sys, os, errno
from subprocess import Popen, PIPE
import glob
import shlex

'''
Gzip utilities, run gzip in a subprocess
'''

### Add try Popen except try gzip.open except
def sp_gzip_read(file):
    p = Popen(shlex.split('gzip --decompress --to-stdout') + [file], stdout = PIPE, stderr = PIPE, bufsize=-1)
## gzip.open(read2, 'rb')
    return p.stdout

def sp_gzip_write(file):
    p = Popen('gzip > ' + file,stdin=PIPE,shell=True)
    return p.stdin


def infer_read_file_name(baseread, seakread):
    ''' Find other read filenames (ex. R1, R2, R3, R4) in the directory based on Read 1 filename '''
    basename = os.path.basename(baseread)
    path = os.path.dirname(os.path.realpath(baseread))
    testname = glob.glob(path + '/*' + os.path.splitext(baseread)[1])
    count = 0
    pos = -1
    read = []
    for name in testname:
        count = 0
        if os.path.basename(name) == basename:  ## ignore the same file
            continue
        elif len(os.path.basename(name)) != len(basename): ## must be the same length
            continue
        else:
            for i, (ch1, ch2) in enumerate(zip(os.path.basename(name), basename)): ## calculate the hamming distance
                if ch1 != ch2 and ch2 == '1' and ch1 == seakread:
                    count += 1
                    pos = i
            if count == 1:
                read.append(path + '/' + basename[0:pos] + seakread + basename[pos+1:])
                continue
    if len(read) == 1:
        return read[0]
    else:
        raise Exception("Error inferring read " + seakread + " from read 1, found " + str(len(read)) + " suitable matches.")


def make_sure_path_exists(path):
    """
    Try and create a path, if not error
    """
    if path != '':
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def expand_path(list):
    newlist = []
    for file in list:
        newlist.append(os.path.realpath(file))
    return newlist

def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])

iupacdict = {
    'A':['A'],
    'C':['C'],
    'G':['G'],
    'T':['T'],
    'M':['A','C'],
    'R':['A','G'],
    'W':['A','T'],
    'S':['C','G'],
    'Y':['C','T'],
    'K':['G','T'],
    'V':['A','C','G'],
    'H':['A','C','T'],
    'D':['A','G','T'],
    'B':['C','G','T'],
    'X':['A','C','G','T'],
    'N':['A','C','G','T']}


def expand_iupac(seq):
    res = ['']
    for m in seq:
        newres = []
        for i in res:
            for j in iupacdict.get(m):
                newres.append(i + j)
        res = newres
    return res


#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
# FROM file.py
'''
File handling/checking utilities for command-line scripts.
'''


def check_space(in_files):
    """
    Estimate size of input files passed, then calculate
    disk space available. Exit if insufficient disk space,
    """

    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    in_file = in_files[0]

    dir_path = os.path.dirname(os.path.realpath(in_file))
    target = os.statvfs(dir_path)
    free_space = target.f_frsize * target.f_bavail

    # Check input file array, remove corrupt files
    valid_files = [f for f in in_files if os.path.isfile(f)]

    # Get input file size as worst case estimate of
    # output file size
    file_sizes = [os.stat(f).st_size for f in valid_files]
    total_size = reduce(lambda f1, f2: f1 + f2, file_sizes)

    size_diff = total_size - free_space
    if size_diff > 0:
        print >>sys.stderr, "ERROR: Not enough free space on disk " \
                            "need at least %s more" % str(size_diff)
        sys.exit(1)


def check_space_for_hashtable(hash_size):
    """
    Check we have enough size to write a hash table
    """
    cwd = os.getcwd()
    dir_path = os.path.dirname(os.path.realpath(cwd))
    target = os.statvfs(dir_path)
    free_space = target.f_frsize * target.f_bavail

    size_diff = hash_size - free_space
    if size_diff > 0:
        print >>sys.stderr, "ERROR: Not enough free space on disk, " \
                            "need at least %s more," % str(size_diff)
        sys.exit(1)

