# cython: language_level=3

from numpy cimport uint32_t
from cpython cimport bool


cdef class BedTrack:
    cdef:
        uint32_t _buffer_size

        dict _chrom_indices
        dict _index_chroms

        uint32_t[:] _chrom_sizes
        uint32_t[:] _left_limits
        uint32_t[:] _right_limits

        uint32_t[:] _lefts
        uint32_t[:] _rights

        uint32_t _pointer
        uint32_t _size

        bool _closed

    cpdef uint32_t get_chromosome_size(self, str chrom)

    cdef uint32_t _locate_position_chromosome(self, uint32_t position)

    cpdef void add_region(self, str chrom, uint32_t start, uint32_t end)

    cpdef void close(self)

    cpdef set_chromosome_sizes(self, dict chrom_sizes)

    cpdef void extend(self, uint32_t left, uint32_t right)

    cpdef bint overlap(self, str chrom, uint32_t position, uint32_t extend=?)
