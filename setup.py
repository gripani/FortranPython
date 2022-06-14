from distutils.core import setup
from sys import argv
from utils import build_fortran 
from pathlib import Path 

if __name__ == '__main__':

    if argv[1] == 'build':
        build_fortran()

    pwd = Path(__file__).parent 
    long_description = (pwd / "README.md").read_text()

    setup(
        name='FortranPython', 
        version='0.0.1', 
        packages=['FortranPython', 'tests'], 
        author='Giorgio Ripani', 
        author_email='g.ripani93@gmail.com',
        long_description=long_description,
    )