import numpy as np
import ctypes as ct 
import sys  
import os.path as osp

from .utils import pointer_array, pointer_cmplx_array, getFunVal_proc

class FortranWrapper: 

    is_linux = 'linux' in sys.platform 
    dll_name = 'my_fortran_library'

    def __init__(self):

        self.dll_path = osp.abspath(osp.join(sys.exec_prefix, 'DLLs', self.dll_name))

        if self.is_linux:
            self.dll_path += '.so'
        else:
            self.dll_path += '.dll'

        try:
            self.fort_lib = ct.CDLL(self.dll_path) if self.is_linux else ct.WinDLL(self.dll_path)
            print(f"Fortran Library Object '{self.fort_lib._name}' is loaded in memory")
            print()
        except Exception as e:
            print(str(e))
            print()


        sys.stdout.flush()

    def integrate(self, f, r0, r1, i=10000):
        """Wrapper for Integrate fortran subroutine 

        Args:
            f (function): callable function to numerically integrate
            r0 (float): lower extreme of integration
            r1 (float): upper extrene of integration
            i (int, optional): Number of iterations. Defaults to 10000.

        Returns:
            float: numeric result of the definite integral
        """
        Integrate = self.fort_lib.Integrate 
        Integrate.restype = None 
        getFun_ptr = getFunVal_proc(f)
        x0 = ct.c_double(r0)
        x1 = ct.c_double(r1)
        outputReal = ct.c_double(.0)
        ii = ct.c_int(i)
        Integrate(getFun_ptr, ct.byref(x0), ct.byref(x1), ct.byref(ii), ct.byref(outputReal))
        return outputReal.value 
    
    def getDeterminant(self, A):
        """Wrapper for GetDeterminant fortran subroutine 

        Args:
            A (ndarray): square matrix for determinant calculation

        Returns:
            float: determinant result
        """
        GetDeterminant = self.fort_lib.GetDeterminant 
        GetDeterminant.restype = None 
        a = np.array(A, order='F')
        shape = a.shape 
        assert shape[0] == shape[1], 'A must be Square Matrix'
        n = ct.c_int(shape[0])
        det = ct.c_double()
        GetDeterminant(ct.byref(n), pointer_array(a), ct.byref(det))
        return det.value 

    def fastFourierTransform(self, signal):
        """Wrapper for FastFourierTransform fortran subroutine 

        Args:
            signal (complex array): complex input signal to transform

        Returns:
            complex array: fourier transform of input signal
        """
        FastFourierTransform = self.fort_lib.FastFourierTransform
        FastFourierTransform.restype = None 
        f = np.array(signal, order='F')
        n = ct.c_int(len(signal))
        FastFourierTransform(ct.byref(n), pointer_cmplx_array(f))
        return f 
    
    def jacobiDiagonalization(self, A):
        """Wrapper for JacobiDiagonalization fortran subroutine 

        Args:
            A (ndarray): _description_

        Returns:
           list(float): eigenvalues 
           ndarray: eigenvectors
        """
        JacobiDiagonalization = self.fort_lib.JacobiDiagonalization
        JacobiDiagonalization.restype = None 
        shape = np.array(A).shape 
        assert shape[0] == shape[1], 'Matrix A must be Square Matrix'
        n = ct.c_int(shape[0])
        a = np.array(A, order='F')
        w = np.zeros(shape[0], order='F')
        JacobiDiagonalization(ct.byref(n), pointer_array(a), pointer_array(w))
        eigenvectors = a
        eigenvalues = w
        return eigenvalues, eigenvectors 
    
    def leastSquareFit(self, X, Y):
        """_summary_

        Args:
            X (_type_): _description_
            Y (_type_): _description_

        Returns:
            _type_: _description_
        """
        LeastSquareFit = self.fort_lib.LeastSquareFit 
        LeastSquareFit.restype = None 
        x = np.array(X, order='F')
        y = np.array(Y, order='F')
        n1 = x.shape[0]
        n2 = y.shape[0]
        assert n1 == n2, 'X and Y must have the same length'
        n = ct.c_int(n1)
        m = ct.c_double()
        q = ct.c_double()
        erm = ct.c_double()
        erq = ct.c_double()
        r2 = ct.c_double()
        LeastSquareFit(ct.byref(n), pointer_array(x), pointer_array(y), ct.byref(m), ct.byref(q), ct.byref(r2), ct.byref(erm), ct.byref(erq))
        return m.value, q.value, erm.value, erq.value, r2.value 

    def gaussJordan(self, A, B):
        GaussJordan = self.fort_lib.GaussJordan 
        GaussJordan.restype = None 
        a = np.array(A, order='F')
        b = np.array(B, order='F')
        na = a.shape 
        assert na[0] == na[1], 'A must be a Square Matrix'
        nb = b.shape 
        assert na[1] == nb[0], 'A 2nd dimension and B 1st dimension must be equal'
        n = ct.c_int(nb[0])
        m = ct.c_int(nb[1])
        GaussJordan(ct.byref(n), ct.byref(m), pointer_array(a), pointer_array(b))
        invA, sol = a, b
        return invA, sol  



