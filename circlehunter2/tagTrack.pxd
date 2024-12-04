# cython: language_level=3

from numpy cimport uint32_t, uint64_t, ndarray
from cpython cimport bool

cdef class TagTrack:
    cdef:
        uint32_t buffer_size
        dict _chroms
        dict _starts
        dict _ends
        dict _reverses
        dict _inserts
        dict _pointers
        dict _sizes
        uint64_t _total
        uint64_t _total_length
        bool _closed

    cpdef void add_tag(self, str chrom, uint32_t start, uint32_t end, bint reverse, uint32_t insert)

    cpdef  void add_chromosome(
            self, str chrom, ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends,
            ndarray[uint32_t, ndim=1] reverses, ndarray[uint32_t, ndim=1] inserts)

    cpdef void add_track(self, TagTrack track)

    cpdef void close(self)

    cpdef set_chromosome_sizes(self, dict chrom_sizes)

    cpdef TagTrack filter(self, bool reverse, uint32_t min_insert)

    cpdef void extend(self, uint32_t left, uint32_t right)

    cpdef tuple make_windows(self, window_size=?)

    cpdef ndarray[uint32_t, ndim=1] count_windows(self, uint32_t window_size=?)

    cpdef uint32_t count(self, str chrom, uint32_t start, uint32_t end)

    cpdef void to_bed(self, object writable, str name_prefix)
