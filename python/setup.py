"""
fast5.setup.py
(c) 2016: Matei David, Ontario Institute for Cancer Research
MIT License
"""

import os
import re
import pkg_resources
import sys
from setuptools import setup, Extension

exec(open('fast5/version.py').read())

# check HDF5 include and lib dirs
hdf5_dir = os.environ.get('HDF5_DIR', '/usr')
hdf5_include_dir = os.environ.get('HDF5_INCLUDE_DIR', os.path.join(hdf5_dir, 'include'))
hdf5_lib_dir = os.environ.get('HDF5_LIB_DIR', os.path.join(hdf5_dir, 'lib'))
hdf5_lib = os.environ.get('HDF5_LIB', 'hdf5')
if not os.path.isfile(os.path.join(hdf5_include_dir, 'H5pubconf.h')):
    sys.exit(hdf5_include_dir + ': could not find HDF5 header files; use HDF5_DIR or HDF5_INCLUDE_DIR')
if (not os.path.isfile(os.path.join(hdf5_lib_dir, 'lib' + hdf5_lib + '.so'))
    and not os.path.isfile(os.path.join(hdf5_lib_dir, 'lib' + hdf5_lib + '.a'))):
    sys.exit(hdf5_lib_dir + ': could not find HDF5 library file; use HDF5_DIR or HDF5_LIB_DIR/HDF5_LIB')

# check Boost.Python include and lib dirs
boost_dir = os.environ.get('BOOST_DIR', '/usr')
boost_include_dir = os.environ.get('BOOST_INCLUDE_DIR', os.path.join(boost_dir, 'include'))
boost_lib_dir = os.environ.get('BOOST_LIB_DIR', os.path.join(boost_dir, 'lib'))
boost_python_lib = os.environ.get('BOOST_PYTHON_LIB', 'boost_python')
if not os.path.isfile(os.path.join(boost_include_dir, 'boost', 'python.hpp')):
    sys.exit(boost_include_dir + ': could not find Boost Python header files; use BOOST_DIR or BOOST_INCLUDE_DIR')
if (not os.path.isfile(os.path.join(boost_lib_dir, 'lib' + boost_python_lib + '.so'))
    and not os.path.isfile(os.path.join(boost_lib_dir, 'lib' + boost_python_lib + '.a'))):
    sys.exit(boost_lib_dir + ': could not find Boost Python library file; use BOOST_DIR or BOOST_LIB_DIR/BOOST_PYTHON_LIB')

fast5_dir = os.environ.get('FAST5_DIR', os.path.join('..', 'src'))

extra_compile_args = [
    '-std=c++11',
    '-Wall', '-Wextra', '-Wpedantic',
]
# don't indiscriminately add /usr/include to work around bug:
# https://lists.fedoraproject.org/archives/list/devel@lists.fedoraproject.org/thread/Q5SWCUUMWQ4EMS7CU2CBOZHV3WZYOOTT/
for d in [hdf5_include_dir, boost_include_dir]:
    if d != '/usr/include':
        extra_compile_args += ['-isystem', d]

#extra_compile_args += ['-O0', '-g3', '-ggdb', '-fno-eliminate-unused-debug-types', '-v']
extra_link_args = []
#extra_link_args += ['-v']

extensions = [
    Extension(
        'fast5.fast5',
        include_dirs=[
            fast5_dir,
        ],
        sources=[
            os.path.join('fast5', 'source', 'fast5.cpp'),
        ],
        depends=[
            os.path.join(fast5_dir, fn)
            for fn in ['fast5.hpp', 'hdf5_tools.hpp']
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        library_dirs=[
            hdf5_lib_dir,
            boost_lib_dir,
        ],
        runtime_library_dirs=[
            hdf5_lib_dir,
            boost_lib_dir,
        ],
        libraries=[
            hdf5_lib,
            boost_python_lib,
        ],
    ),
]

setup(
    name='fast5',
    description='Fast5 file interface.',
    version=__version__,
    #long_description=open('README').read(),
    author='Matei David, Ontario Institute for Cancer Research',
    author_email='matei.david at oicr.on.ca',
    license='MIT',
    url='https://github.com/mateidavid/fast5',
    packages=['fast5'],
    exclude_package_data={
        '': ['*.c', '*.cpp', '*.h', '*.hpp'],
    },
    ext_modules=extensions,
    scripts=[],
)
