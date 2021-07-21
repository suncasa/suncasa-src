import ctypes
from numpy.ctypeslib import ndpointer

def initGET_MW(libname):
    _intp=ndpointer(dtype=ctypes.c_int32, flags='F')
    _doublep=ndpointer(dtype=ctypes.c_double, flags='F')
    
    libc_mw=ctypes.CDLL(libname)
    mwfunc=libc_mw.pyGET_MW
    mwfunc.argtypes=[_intp, _doublep, _doublep, _doublep, _doublep, _doublep, _doublep]
    mwfunc.restype=ctypes.c_int

    return mwfunc

def initGET_MW_SLICE(libname):
    _intp=ndpointer(dtype=ctypes.c_int32, flags='F')
    _doublep=ndpointer(dtype=ctypes.c_double, flags='F')
    
    libc_mw=ctypes.CDLL(libname)
    mwfunc=libc_mw.pyGET_MW_SLICE
    mwfunc.argtypes=[_intp, _doublep, _doublep, _doublep, _doublep, _doublep, _doublep]
    mwfunc.restype=ctypes.c_int

    return mwfunc
