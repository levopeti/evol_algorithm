////
//// Created by biot on 2018.05.23..
////
//#include <iostream>
//#include <math.h>
//#include <float.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
//#include <vector>
//#include <fstream>
//#include "stdafx.h"
//#include <Python.h>
//#include <stdio.h>
//#include <stdlib.h>
//
//using namespace std;
//
//void add(int a, int b)
//{
//    cout<<a+b<<endl;
//}
//
//
//unsigned int add_p(long a, long b)
//{
//
//    Py_Initialize();
//   // PyRun_SimpleString("import sys");
//   // PyRun_SimpleString("sys.path.append(\".\")");
//    //PyRun_SimpleString("sys.path.append(\"/home/biot/projects/szakdolgozat/Evolutionary_algorithm/\")");
//   // PyRun_SimpleString("import add");
//   // Py_Finalize();
//
//    PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pValue1, *pValue2;
//
//    pName = PyUnicode_FromString("/home/biot/projects/szakdolgozat/Evolutionary_algorithm/add.py");
//    pModule = PyImport_Import(pName);
//    if (pModule == nullptr)
//    {
//        PyErr_Print();
//        std::exit(1);
//    }
//    pDict = PyModule_GetDict(pModule);
//    pFunc = PyDict_GetItemString(pDict, "add_in_py");
//    pArgs = PyTuple_New(2);
//    pValue1 = PyLong_FromLong(a);
//    pValue2 = PyLong_FromLong(b);
//    PyTuple_SetItem(pArgs, 0, pValue1);
//    PyTuple_SetItem(pArgs, 1, pValue1);
//
//    PyObject *pResult = PyObject_CallObject(pFunc, pArgs);
//    long result = PyLong_AsLong(pResult);
//
//    Py_Finalize();
//
//    //unsigned int result = system("/usr/bin/python /home/biot/projects/szakdolgozat/Python_for_EA/CNN/cnn.py");
//
//    return result;
//
//
//}