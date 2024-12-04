# cython: language_level=3

from numpy cimport uint32_t, float64_t
from numpy cimport ndarray
from cpython cimport bool


cdef class BedGraphTrack:
    """
    Class representing values associated with a genome region like bedGraph file.
    """
    cdef:
        object _chroms
        dict _tmp_positions
        dict _tmp_values
        dict _counts
        uint32_t[:] _positions
        float64_t[:] _values
        float64_t _default
        bool _closed

    cpdef add_chromosome(self,
                         str chrom,
                         ndarray[uint32_t, ndim=1] positions,
                         ndarray[float64_t, ndim=1] values)

    cpdef close(self)

    cpdef coverage(self, str chrom, uint32_t start, uint32_t end, uint32_t cutoff)

    cpdef BedGraphTrack new_values(self, ndarray[float64_t, ndim=1] values)

    cpdef BedGraphTrack new_value(self, float64_t value)

    cpdef float64_t mean(self)
