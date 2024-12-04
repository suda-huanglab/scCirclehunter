# cython: language_level=3

from numpy cimport uint32_t, float64_t, ndarray

cpdef tuple make_windows(dict chrom_sizes, long size)

cdef tuple _pileup(ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends)

cpdef tuple pileup(ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends)

cdef float64_t _logp(float64_t x, float64_t mu)

cpdef float64_t logp(float64_t x, float64_t mu)
