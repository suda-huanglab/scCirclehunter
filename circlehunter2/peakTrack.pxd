# cython: language_level=3

from numpy cimport uint8_t, uint32_t, float64_t, ndarray
from cpython cimport bool


cdef class PeakTrack:
    """
    Class representing Peaks.
    """

    cdef:
        uint32_t _buffer_size

        object _chroms
        dict _pointers
        dict _sizes
        dict _tmp_starts
        dict _tmp_ends
        dict _tmp_max_logp
        dict _tmp_max_fc
        dict _tmp_mean_logp
        dict _tmp_mean_fc
        dict _counts

        uint32_t[:] _starts
        uint32_t[:] _ends
        float64_t[:] _max_logp
        float64_t[:] _max_fc
        float64_t[:] _mean_logp
        float64_t[:] _mean_fc

        bool _closed

    cpdef void add_peak(self, str chrom, uint32_t start, uint32_t end,
                        float64_t max_logp, float64_t max_fc,
                        float64_t mean_logp, float64_t mean_fc)

    cpdef close(self)

    cpdef ndarray[uint32_t, ndim=1] count_overlaps(self, PeakTrack other)

    cpdef float64_t coverage(self, str chrom, uint32_t start, uint32_t end)

    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep)

    cpdef void to_bed(self, object writable, str name_prefix)
