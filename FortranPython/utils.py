from ctypes import Structure, c_double, POINTER, CFUNCTYPE  

getFunVal_proc = CFUNCTYPE(c_double, c_double)

class CDoubleComplex(Structure):

    _fields_ = [("real", c_double), ("imag", c_double)]

    @property 
    def value(self):
        return self.real + 1j * self.imag 

def pointer_cmplx_array(cmplx_array):
    return cmplx_array.ctypes.data_as(POINTER(CDoubleComplex))

def pointer_array(array):
    return array.ctypes.data_as(POINTER(c_double))
