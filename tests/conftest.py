import pytest 
from numpy import linspace, sqrt, cos, sin, pi 

@pytest.fixture 
def get_A_for_det():
    a = [  
        [10., .0, -3.],
        [-2.,-4., 1.],
        [3., .0, 2.]
    ]
    return a 

@pytest.fixture 
def get_A_for_diag():
    a = [ 
        [1., sqrt(2.), 2.],
        [sqrt(2), 3., sqrt(2)],
        [2., sqrt(2.), 1.]
    ]
    return a 

@pytest.fixture 
def get_Ab_for_GaussJordan():
    a = [
        [ 2., 1.,-1., 2.],
        [ 4., 5.,-3., 6.],
        [-2., 5.,-2., 6.],
        [ 4.,11.,-4., 8.]
    ]
    b = [
        [ 5.,10.],
        [ 9.,18.],
        [ 4., 8.],
        [ 2., 4.]
    ]
    return a, b 

@pytest.fixture
def get_function_for_callback():
    getFunction = lambda x: x * cos(10.*x**2) / (x**2 + 1.)
    extr = [.0, pi]
    return extr, getFunction 

@pytest.fixture 
def get_signal_for_fft():
    N = 2048 
    T = 1. / 720. 
    t = linspace(.0, N*T, N)
    f1 = 50.*2.*pi 
    f2 = 80.*2.*pi 
    signal = sin(f1*t) + .5*sin(f2*t) + 0j*t
    return signal 
