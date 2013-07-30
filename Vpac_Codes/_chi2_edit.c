#include <Python.h>
#include <numpy/arrayobject.h>
#include "chi2.h"



static char module_docstring[] =
	"This module provides an interface for calculating chi-squared using C.";
static char chi2_docstring[] =
	"Calculate the chi-squared of some data given a model";
static PyObject *chi2_chi2(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
	{"chi2",chi2_chi2,METH_VARARGS,chi2_docstring},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_chi2(void)
{
	PyObject *m = Py_InitModule3("_chi2",module_methods, module_docstring);
	if(m == NULL )
		return;

	/*load `numpy' functionality*/

	import_array();



}

static PyObject *chi2_chi2(PyObject *self, PyObject *args)
{

PyObject *params_obj; 
int lum_sample;

 /* Parse the input tuple */
 if (!PyArg_ParseTuple(args, "id",&lum_sample,&params_obj))
          return NULL;

PyObject *params_array = PyArray_FROM_OTF(params_obj,NPY_DOUBLE,NPY_IN_ARRAY);

double *params = (double *)PyArray_DATA(params_array);

/* Call the external C function to compute the chi-squared. */
     double value = chi2(lum_sample,params);

Py_DECREF(params_array);

if (value < 0.0) {
         PyErr_SetString(PyExc_RuntimeError,
                     "Chi-squared returned an impossible value.");
         return NULL;
     }

/* Build the output tuple */
     PyObject *ret = Py_BuildValue("d", value);
     return ret;

}
