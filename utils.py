from os import chdir, system
import os.path as osp 
from sys import platform, exec_prefix 
from shutil import copy2
import json

def build_fortran():
   
    with open('config.json', 'r') as f:
        config = json.load(f)
    fortran_compiler = config['fortran_compiler']
    lib_name = config['lib_name']
    lib_dir = config['lib_dir']
    ext = '.so' if 'linux' in platform else '.dll'
    chdir(lib_dir)
    print(f'compiling {lib_name}...')
    system(f'"{fortran_compiler} -fPIC -shared {lib_name}.f -o {lib_name+ext} -ffixed-line-length-none"')
    print('compiled')
    dest =  osp.join(exec_prefix, 'DLLs')
    print(f'copying {lib_name+ext} in {dest}...')
    print('copied')
    print('')
    copy2(lib_name+ext, dest)
    chdir('..')