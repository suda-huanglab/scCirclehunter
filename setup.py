#!/usr/bin/env python3

from setuptools import setup, Extension
import os

import numpy as np
import pysam

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


extensions = [
    Extension(
        "circlehunter2.tagTrack", ["circlehunter2/tagTrack.pyx"],
        include_dirs=[np.get_include()],
        # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    ),
    Extension(
        "circlehunter2.tagTrackCollection",
        ["circlehunter2/tagTrackCollection.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.pairedTagTrack",
        ["circlehunter2/pairedTagTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.pairedTagTrackCollection",
        ["circlehunter2/pairedTagTrackCollection.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.utils", ["circlehunter2/utils.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.bedGraphTrack", ["circlehunter2/bedGraphTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.peakTrack", ["circlehunter2/peakTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.pairedPeakTrack",
        ["circlehunter2/pairedPeakTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.pairedBedGraphTrack",
        ["circlehunter2/pairedBedGraphTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.bedTrack",
        ["circlehunter2/bedTrack.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "circlehunter2.bamReader",
        ["circlehunter2/bamReader.pyx"],
        include_dirs=[np.get_include()] + pysam.get_include()
    )
]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

with open("requirements.txt") as fp:
    install_requires = fp.read().strip().split("\n")

with open("requirements-dev.txt") as fp:
    dev_requires = fp.read().strip().split("\n")

setup(
    name="circlehunter2",
    version="0.1.2",
    packages=["circlehunter2"],
    ext_modules=extensions,
    install_requires=install_requires,
    extras_require={
        "dev": dev_requires,
    },
    entry_points={
        'console_scripts': [
            'circlehunter2 = circlehunter2.__main__:main',
        ]
    }
)
