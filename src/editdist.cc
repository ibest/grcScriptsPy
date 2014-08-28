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

/* Original Copyright for Levenstein Dist function edit_dist
 * Copyright (c) 2006 Damien Miller <djm@mindrot.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/* $Id: editdist.c,v 1.5 2007/05/03 23:36:36 djm Exp $ */
/* $Id: editdist.c,v 0.4 2013/12/31 mls $ */
/* $Id: editdist.c,v 0.4 2014/4/17 mls added hamming distance */
/* $Id: editdist.c,v 0.5 2014/7/14 mls  moved editdist over to grcScriptsPy*/

#include "editdist.hh"

#include <stdlib.h>
#include <string.h>
/* 
compute the Levenstein distance between a (primer) and b (sequence read)
pegged to the 5' end and bounded by edit distance k
*/
Tuple
bounded_editdist(const char *a, int alen, const char *b, int blen, int k, int m)
{
    // a is primer, b is seq, k is max error and m is end matches
    int i, j;
    int *current, *previous, *tmpl, add, del, chg;
    Tuple val = { k+1, alen };
    /* a (primer) should always be < b (read) + k */
//  if (alen > (blen-k)) {
    if (alen > (blen)) {
        return (val);
    }

    if (alen == 0){
        val.dist = -2;
        return (val);
    }
    if ((previous = (int *)calloc(alen + 1, sizeof(int))) == NULL) {
        free(previous);
        val.dist = -1;
        return (val);
    }
    if ((current = (int *)calloc(alen + 1, sizeof(int))) == NULL) {
        free(current);
        val.dist = -1;
        return (val);
    }

    for (i = 0; i < alen - m + 1; i++)
        previous[i] = i;

    for (i = 1; i < alen - m + k + 1 ; i++) {
        if (i > 1) {
            memset(previous, 0, (alen - m + 1) * sizeof(*previous));
            tmpl = previous;
            previous = current;
            current = tmpl;
        }
        current[0] = i;
        for (j = 1; j < alen -m + 1; j++) {
            add = previous[j] + 1;
            del = current[j - 1] + 1;
            chg = previous[j - 1];
            if (a[j - 1] != b[i - 1])
                chg++;
            current[j] = MIN(add, del);
            current[j] = MIN(current[j], chg);
        }
        if (current[alen - m] <= val.dist ){
            val.dist = current[alen - m ]; // bottom right node, global alignment
            val.pos = i+m;
        }
    }

    for (i = 1; i <= m; i++){
        if (a[alen - i] != b[val.pos - i]){
            val.dist = val.dist+100; // penelty of 100 for not meeting end match criteria
            break;          
        }       
    }
    free(previous);
    free(current);

    return (val);
}

/*
Compute the Levenstein distance between a and b 
*/
int
edit_distance(const char *a, int alen, const char *b, int blen)
{
    int tmplen, i, j;
    const char *tmp;
    int *current, *previous, *tmpl, add, del, chg, r;

    /* Swap to reduce worst-case memory requirement */
    if (alen > blen) {
        tmp = a;
        a = b;
        b = tmp;
        tmplen = alen;
        alen = blen;
        blen = tmplen;
    }

    if (alen == 0)
        return (blen);

    if ((previous = (int *)calloc(alen + 1, sizeof(int))) == NULL) {
        free(previous);
        return (-1);
    }
    if ((current = (int *)calloc(alen + 1, sizeof(int))) == NULL) {
        free(current);
        return (-1);
    }
    for (i = 0; i < alen + 1; i++)
        previous[i] = i;

    for (i = 1; i < blen + 1; i++) {
        if (i > 1) {
            memset(previous, 0, (alen + 1) * sizeof(*previous));
            tmpl = previous;
            previous = current;
            current = tmpl;
        }
        current[0] = i;
        for (j = 1; j < alen + 1; j++) {
            add = previous[j] + 1;
            del = current[j - 1] + 1;
            chg = previous[j - 1];
            if (a[j - 1] != b[i - 1])
                chg++;
            current[j] = MIN(add, del);
            current[j] = MIN(current[j], chg);
        }
    }
    r = current[alen];
    free(previous);
    free(current);
    return (r);
}

/*
Compute the Levenstein distance between a and b 
*/
int
hammingdist(const char *a, int alen, const char *b, int blen)
{
    int i, diff;

    // a and b should be same length, otherwise return error
    if (alen != blen){
        return (-2);
    }
    // length should not be 0, otherwise return error
    if (alen == 0)
        return (-2);

    diff = 0;
    for (i = 0; i < alen; i++){
        if (a[i] != b[i])
            diff++;
    }
    return (diff);
}

