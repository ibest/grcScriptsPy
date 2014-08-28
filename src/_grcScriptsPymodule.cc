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

/*
Python C/C++ interface functions
*/

#include <Python.h>

#include "editdist.hh"


static PyObject *
print_hello_world( PyObject * self, PyObject * args )
{
    printf("Hello World, C version\n");
    return Py_BuildValue("");
}

PyDoc_STRVAR(hello_world_doc,
"hello_wold()\n\
    Displays hello world on the screen\n");


// Python C interface functions for C functions in editdist
// Calculate Hamming distance, Levenshtein's edit distance, and a edge bounded Levenshtein's edit distance.

PyDoc_STRVAR(bounded_editdist_distance_doc,
"bounded_edit_distance(a, b, k, m) -> int, int\n\
    Calculates the bounded Levenshtein's edit distance between strings \"a\" and \"b\" with bound \"k\" and \"m\" matching bases at end\n");

static PyObject *
bounded_editdist_distance(PyObject *self, PyObject *args)
{
    Tuple r;
    char *a, *b;
    int alen, blen, k, m;

    if (!PyArg_ParseTuple(args, "s#s#ii", &a, &alen, &b, &blen, &k, &m))
                return NULL;
    r = bounded_editdist(a, alen, b, blen, k, m);
    if (r.dist== -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    if (r.dist== -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return Py_BuildValue("ii", r.dist, r.pos);
}

PyDoc_STRVAR(editdist_distance_doc,
"distance(a, b) -> int\n\
    Calculates Levenshtein's edit distance between strings \"a\" and \"b\"\n");

static PyObject *
editdist_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;
    r = edit_distance(a, alen, b, blen);
    if (r == -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    return PyInt_FromLong(r);
}

PyDoc_STRVAR(hammingdist_distance_doc,
"distance(a, b) -> int\n\
    Calculates hamming distance between two equal length strings \"a\" and \"b\"\n");

static PyObject *
hammingdist_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;
    r = hammingdist(a, alen, b, blen);
    if (r == -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return PyInt_FromLong(r);
}
// END interface functions for editdist

// grcScripts_methods
// Module declarations

static PyMethodDef grcScripts_methods[] = {
    {   
        "hello_world", (PyCFunction)print_hello_world,
        METH_NOARGS,    hello_world_doc       
    },
    {   "edit_distance", (PyCFunction)editdist_distance,
        METH_VARARGS,   editdist_distance_doc
    },
    {   "bounded_distance", (PyCFunction)bounded_editdist_distance,
        METH_VARARGS,   bounded_editdist_distance_doc
    },
    {   "hamming_distance", (PyCFunction)hammingdist_distance,
        METH_VARARGS,    hammingdist_distance_doc
    },
    { NULL, NULL, 0, NULL }  /* sentinel */
};


PyDoc_STRVAR(module_doc, "Interface for the grcScriptsPy module low-level extensions\n");


PyMODINIT_FUNC
init_grcScripts(void)
{
    PyObject *m;

    m = Py_InitModule3("_grcScripts", grcScripts_methods, module_doc);
    if (m == NULL)
        return;
}

