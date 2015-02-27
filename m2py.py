#!/usr/bin/env python
''' Debug & Test support for matplot to python conversion.
'''
import os
import numpy as np
from scipy.io import loadmat


def dmpdat(s, e):
    """ Dump a data structure with its name & shape.
    Params:
    -------
    s: str. The name of the structure
    e: expression. An expression to dump. Implicitly assumes e is
    array_like
    """
    print("%s:" % s)
    print(e)
    print("%s.shape:" % s)
    print(e.shape)
    print("%s.dtype:" % s)
    print(e.dtype)
    print("-------------------------------------------")


def hbrk(msg=None):
    if msg is not None:
        print(msg)
    exit(-1)


def brk(s, e):
    """ Used for debugging, just break the script, dumping data.
    """
    dmpdat(s, e)
    exit(-1)


def chkdat(t, s, e, rtol=1e-05, atol=1e-08):
    """ Check this matrix against data dumped by octave, with
    given tolerance
    """
    mat = loadmat(os.path.join('check_data', t, s) + '.mat')['ex']
    is_equal = np.allclose(e, mat, rtol=rtol, atol=atol)
    #is_equal =  np.array_equal(e, mat)
    print("%s:%s:iEqual=%d" % (t, s, is_equal))
    if not is_equal:
        dmpdat(s + '<python>', e)
        dmpdat(s + '<matlab>', mat)
        np.savetxt(os.path.join("check_data", t, s) + '_python_err', e)
        np.savetxt(os.path.join("check_data", t, s) + '_matlab_err', mat)
        print("FAILED check on expr: %s, signal: %s" % (s, t))
        #hbrk()
    return is_equal
