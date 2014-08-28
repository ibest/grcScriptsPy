/*
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
*/

#ifndef EDITDIST_H
#define EDITDIST_H



#if defined(_MSC_VER)
typedef unsigned __int8 u_int8_t;
#endif

#ifndef MIN
# define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct _Tuple {
    int dist;
    int pos;
} Tuple;

Tuple bounded_editdist(const char *a, int alen, const char *b, int blen, int k, int m);

int edit_distance(const char *a, int alen, const char *b, int blen);

int hammingdist(const char *a, int alen, const char *b, int blen);

//EDITDIST_H
#endif  
