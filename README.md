Usage:
```python
import numpy as np 
from FortranPython.fpython import wrapper as flib

#Least square fit 
N = int(8e3)

noise = np.random.randn(N)
X = [1.*i for i in range(N)]
Y = [(2.*X[i] + noise[i]) for i in range(N)]
m, q, erm, erq, r2 = flib.leastSquareFit(X, Y)

print(f'm, q values = {m:.3f}, {q:.3f}')
print(f'with errors = {erm:.3f}, {erq:.3f}')
print(f'coefficient of determination r^2 = {r2:.3f}')
print('')

#Determinant
A = [ 
     [10., .0, -3.],
     [-2.,-4., 1.],
     [3., .0, 2.]
]
det = flib.getDeterminant(A)

print(f'matrix A = \n{np.array(A)}\n')
print(f'determinant = {det}')
print('')

#Jacobi diagonalization
A = [[1., np.sqrt(2.), 2.],
     [np.sqrt(2), 3., np.sqrt(2)],
     [2., np.sqrt(2.), 1.]]

eigvals, eigvecs = flib.jacobiDiagonalization(A)

idx = eigvals.argsort()[::-1]
eigvals, eigvecs = eigvals[idx], eigvecs[:, idx] 

print(f'matrix A = \n{np.array(A)}\n')
print('eigvecs, eigvals')
for eigval, eigvec in zip(eigvals, eigvecs):
    print(eigvec, eigval)
print('')

#Gauss Jordan elimination
a = [[ 2., 1.,-1., 2.],
     [ 4., 5.,-3., 6.],
     [-2., 5.,-2., 6.],
     [ 4.,11.,-4., 8.]]
b = [[ 5.,10.],
     [ 9.,18.],
     [ 4., 8.],
     [ 2., 4.]]
inv, sol = flib.gaussJordan(a, b)

print(f'coefficient matrix = \n{np.array(a)}')
print(f'RHS vector = \n{np.array(b)}')
print('')
print(f'solution x = \n{sol}')
print(f'inverse of coefficient matrix = \n{inv}')

#Callback function integrate
res = flib.integrate(lambda x: x * np.cos(10.*x**2) / (x**2 + 1.), .0, np.pi)

print(f'result of the numeric integral = {res:.7f}')
print('')

#Fast Fourier Transform
from numpy import linspace, sin, pi
  
N = 2048 
T = 1. / 720.
time = linspace(.0, N*T, N)
f1 = 50.*2.*pi
f2 = 80.*2.*pi
signal = sin(f1*time) + .5*sin(f2*time) + 0j*time
ft = flib.fastFourierTransform(signal)

import matplotlib.pyplot as plt 
print('Plotting...')
freq = linspace(.0, 1.//T, N)
_, axes = plt.subplots(2)
axes[0].plot(time, signal.real);
axes[0].plot(time, signal.imag);
axes[1].plot(freq[:N//2], 2./N * abs(ft[:N//2]))
plt.show()
```