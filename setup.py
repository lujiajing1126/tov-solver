from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
import glob

module = Extension('tov',
                    sources = glob.glob("./src/*.c"),
                    libraries = ['gsl'])

setup (name = 'TOVSolver',
       version = '1.0.0',
       description = 'This is a high performance python-C extension for solving TOV equation',
       author = 'Megrez Lu',
       author_email = 'lujiajing1126@gmail.com',
       include_dirs = [*get_numpy_include_dirs()],
       ext_modules = [module])