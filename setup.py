from distutils.core import setup
from sys import argv
from utils import build_fortran 

if __name__ == '__main__':

    if argv[1] == 'build':
        build_fortran()

    setup(
        name='FortranPython', 
        version='0.0.1', 
        packages=['FortranPython'], 
        author='Giorgio Ripani', 
        author_email='g.ripani93@gmail.com'
    )