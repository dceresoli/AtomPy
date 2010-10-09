from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

files = ["Atomic/*", "Atomic/clib/*", "Doc/*"]
ext_modules = [Extension('Atomic/clib/shoot', ['Atomic/clib/shoot.pyx'])]

setup(
    name = "AtomPy",
    version = "0.1",
    description = "Atomic calculations in python",
    author = "Davide Ceresoli",
    author_email = "dceresoli@gmail.com",
    url = "",
    packages = ['Atomic', 'Atomic/clib'],
    scripts = ['run_atom.py', 'plot_atom.py'],
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext},
) 


