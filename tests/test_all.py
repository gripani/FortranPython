import numpy as np
from scipy.integrate import quad 
from scipy.fft import fft 

from FortranPython.fpython import wrapper as flib

def test_determinant(get_A_for_det):
    fort_det = flib.getDeterminant(get_A_for_det)
    numpy_det = np.linalg.det(get_A_for_det)
    assert np.isclose(fort_det, numpy_det, rtol=1e-05, atol=1e-08, equal_nan=False)

def test_diagonalization(get_A_for_diag):
    fort_eigvals, fort_eigvecs = flib.jacobiDiagonalization(get_A_for_diag)
    idx = fort_eigvals.argsort()[::-1]
    fort_eigvals, fort_eigvecs = fort_eigvals[idx], fort_eigvecs[:, idx] 
    numpy_eigvals, numpy_eigvecs = np.linalg.eig(get_A_for_diag)
    idx = numpy_eigvals.argsort()[::-1]
    numpy_eigvals, numpy_eigvecs = numpy_eigvals[idx], numpy_eigvecs[:, idx]
    for i in range(len(fort_eigvals)):
        assert np.isclose(fort_eigvals[i], numpy_eigvals[i],  rtol=1e-05, atol=1e-08, equal_nan=False)
        try:
            fort_eigvec_i = fort_eigvecs[:, i]
            numpy_eigvec_i = numpy_eigvecs[:, i]
            assert np.allclose(fort_eigvec_i, numpy_eigvec_i,  rtol=1e-05, atol=1e-08, equal_nan=False)
        except:
            fort_eigvec_i = fort_eigvecs[:, i]
            numpy_eigvec_i = -numpy_eigvecs[:, i]
            assert np.allclose(fort_eigvec_i, numpy_eigvec_i,  rtol=1e-05, atol=1e-08, equal_nan=False) 

def test_elimination(get_Ab_for_GaussJordan):
    a, b = get_Ab_for_GaussJordan
    fort_invA, fort_sol = flib.gaussJordan(a, b)
    numpy_sol = np.linalg.solve(a, b)
    numpy_invA = np.linalg.inv(a)
    assert np.allclose(fort_sol, numpy_sol, rtol=1e-05, atol=1e-08, equal_nan=False)
    assert np.allclose(fort_invA, numpy_invA, rtol=1e-05, atol=1e-08, equal_nan=False)
    
def test_integrate(get_function_for_callback):
    extr, getFunction = get_function_for_callback
    x0, x1 = extr 
    f_res = flib.integrate(getFunction, x0, x1)
    sp_res, _ = quad(getFunction, x0, x1)
    print(f_res, sp_res)
    assert np.isclose(f_res, sp_res, rtol=1e-05, atol=1e-08, equal_nan=False)

def test_fourier_transform(get_signal_for_fft):
    signal = get_signal_for_fft 
    n = len(signal)
    fort_fft = flib.fastFourierTransform(signal)
    sp_fft = fft(signal)
    assert np.allclose(abs(fort_fft[:n//2]), abs(sp_fft[:n//2]), rtol=1e-05, atol=1e-08, equal_nan=False)