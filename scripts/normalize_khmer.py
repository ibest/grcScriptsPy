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

import os
import sys


profile = False



def main():
    """
    main function
    """
    lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
    if lib_path not in sys.path:
        sys.path.insert(0, lib_path)
    print "Hello World"
    

if __name__ == '__main__':
    main()
 
