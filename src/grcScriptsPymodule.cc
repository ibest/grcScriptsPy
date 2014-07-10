/*
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
*/

include "Python.h"

using namespace grcScripts;


//
// Function necessary for Python loading:
//

extern "C" {
    void init_grcScripts();
}

void
print_hello_world( PyObject * self, PyObject * args )
{
    printf("Hello World\n");
    return;
}

PyDoc_STRVAR(hello_world_doc,
"hello_wold()\n\
    Displays hello world on the screen\n");

static PyMethodDef grcScripts_methods[] = {
    {   
        "hello_world", (PyCFunction)print_hello_wold,
        METH_NOARGS,    hello_world_doc       
    },
    { NULL, NULL, 0, NULL }  /* sentinel */
};


PyDoc_STRVAR(module_doc, "Interface for the grcScriptsPy module low-level extensions\n");


PyMODINIT_FUNC
init_grcScripts(void)
{
    PyObject *m;

    m = Py_InitModule3("grcScripts", grcScripts_methods, module_doc);
}

