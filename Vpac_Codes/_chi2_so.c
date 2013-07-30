#include <Python.h>
#include <numpy/arrayobject.h>
#include "chi2_so.h"



static char module_docstring[] =
	"This module provides an interface for calculating chi-squared using C.";
static char chi2_docstring[] =
	"Calculate the chi-squared of some data given a model";
static PyObject *chi2_so_chi2_so(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
	{"chi2_so",chi2_so_chi2_so,METH_VARARGS,chi2_docstring},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_chi2_so(void)
{
	PyObject *m = Py_InitModule3("_chi2_so",module_methods, module_docstring);
	if(m == NULL )
		return;

	/*load `numpy' functionality*/

	import_array();



}

static PyObject *chi2_so_chi2_so(PyObject *self, PyObject *args)
{

double siglogM,logM0,logM1,alpha,gamma,fgal;
int lum_sample,box;

 /* Parse the input tuple */
 if (!PyArg_ParseTuple(args, "iidddddd",&lum_sample,&box,&siglogM,&logM0,&logM1,&alpha,&gamma,&fgal))
          return NULL;

/* Call the external C function to compute the chi-squared. */
     double value = chi2_so(lum_sample,box,siglogM,logM0,logM1,alpha,gamma,fgal);

if (value < 0.0) {
         PyErr_SetString(PyExc_RuntimeError,
                     "Chi-squared returned an impossible value.");
         return NULL;
     }

/* Build the output tuple */
     PyObject *ret = Py_BuildValue("d", value);
     return ret;

}
