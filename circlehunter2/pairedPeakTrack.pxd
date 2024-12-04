# cython: language_level=3

from numpy cimport uint8_t, uint32_t, ndarray
from cpython cimport bool


cdef class PairedPeakTrack:
    """
    Class representing Peaks.
    """
    cdef:
        uint32_t _buffer_size

        dict _chrom_indices
        dict _index_chroms

        uint32_t[:] _chrom_sizes
        uint32_t[:] _left_limits
        uint32_t[:] _right_limits

        uint32_t[:] _lefts1
        uint32_t[:] _rights1
        uint8_t[:] _reverses1
        uint32_t[:] _lefts2
        uint32_t[:] _rights2
        uint8_t[:] _reverses2
        uint32_t[:] _counts

        uint32_t _pointer
        uint32_t _size

        bool _closed

    cpdef uint32_t get_chromosome_size(self, str chrom)

    cdef uint32_t _locate_position_chromosome(self, uint32_t position)

    cpdef void add_peak_pair(self,
                        str chrom1, uint32_t start1, uint32_t end1, bint reverses1,
                        str chrom2, uint32_t start2, uint32_t end2, bint reverses2,
                        uint32_t count)

    cpdef void close(self)

    cpdef void keep_bi_peaks(self)

    cpdef duplicated(self)

    cpdef calculate_gc_percent(self, str fasta_file)

    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep)

    cpdef dict to_graph_data(self)

    cpdef void to_bedpe(self, object writable)

    cpdef void to_bed(self, object writable)
