try:
    from setuptools import setup, Extension
    from setuptools.command.install import install

except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

import sys

include_dirs = []
# import numpy as np

# creat dummy closures for metadata reading on first parse of setup.py
# that way it picks up the requirements and installs them, then can use them
# for the install.
try:
    import numpy as np
    include_dirs = ['include/', 'src/','python/', 'slow5lib/include/', 'slow5lib/thirdparty/streamvbyte/include/', np.get_include()]
except ImportError:
    include_dirs = ['include/', 'src/','python/', 'slow5lib/include/', 'slow5lib/thirdparty/streamvbyte/include/']
    def np(*args, ** kwargs ):
        import numpy as np
        return np(*args, ** kwargs)

try:
    from Cython.Build import build_ext
except ImportError:
    def build_ext(*args, ** kwargs ):
        from Cython.Build import build_ext
        return build_ext(*args, ** kwargs)

# from Cython.Build import build_ext

#adapted from https://github.com/lh3/minimap2/blob/master/setup.py

sources=['python/pysquigulator.pyx',
         'src/format.c',
         'src/genread.c',
         'src/gensig.c',
         'src/methmodel.c',
         'src/misc.c',
         'src/model.c',
         'src/ref.c',
         'src/sim.c',
         'src/sq_api.c',
         'src/thread.c',
         'slow5lib/src/slow5.c',
         'slow5lib/src/slow5_press.c',
         'slow5lib/src/slow5_misc.c',
         'slow5lib/src/slow5_idx.c',
         'slow5lib/thirdparty/streamvbyte/src/streamvbyte_zigzag.c',
         'slow5lib/thirdparty/streamvbyte/src/streamvbyte_decode.c',
         'slow5lib/thirdparty/streamvbyte/src/streamvbyte_encode.c']
depends=['python/pysquigulator.pxd',
         'python/pysquigulator.h',
         'include/squigulator.h',
         'src/error.h',
         'src/format.h',
         'src/khash.h',
         'src/kseq.h',
         'src/ksort.h',
         'src/misc.h',
         'src/model.h',
         'src/rand.h',
         'src/ref.h',
         'src/seq.h',
         'src/sq.h',
         'src/str.h',
         'src/version.h',
         'slow5lib/include/slow5/slow5.h',
         'slow5lib/include/slow5/slow5_defs.h',
         'slow5lib/include/slow5/slow5_error.h',
         'slow5lib/include/slow5/slow5_press.h',
         'slow5lib/include/slow5/klib/khash.h',
         'slow5lib/include/slow5/klib/kvec.h',
         'slow5lib/src/slow5_extra.h',
         'slow5lib/src/slow5_idx.h',
         'slow5lib/src/slow5_misc.h',
         'slow5lib/src/slow5_byte.h',
         'slow5lib/thirdparty/streamvbyte/include/streamvbyte.h',
         'slow5lib/thirdparty/streamvbyte/include/streamvbyte_zigzag.h']
extra_compile_args = ['-g', '-Wall', '-O2', '-std=c99']


libraries = ['m', 'z']
library_dirs = ['.']

extensions = [Extension('pysquigulator',
                  sources = sources,
                  depends = depends,
                  extra_compile_args = extra_compile_args,
                  libraries = libraries,
                  include_dirs = include_dirs,
                  library_dirs = library_dirs,
                  language = 'c' )]

def readme():
	with open('python/README.md') as f:
		return f.read()


setup(
    name = 'pysquigulator',
    version='0.5.0-dev',
    url = 'https://github.com/hasindu2008/squigulator',
    description='pysquigulator python bindings',
    long_description=readme(),
    long_description_content_type='text/markdown',
    author='Hasindu Gamaarachchi',
    author_email='hasindu@garvan.org.au',
    maintainer='Hasindu Gamaarachchi',
    maintainer_email='hasindu@garvan.org.au',
    license = 'MIT',
    keywords = ['nanopore', 'signal'],
    ext_modules=extensions,
    cmdclass= {'build_ext': build_ext},
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    python_requires='>=3.4.3',
    install_requires=["numpy"],
    setup_requires=["Cython", "numpy"]
    )
