# cython: language_level=3

from numpy cimport uint32_t, uint8_t, ndarray
from cpython cimport bool

from circlehunter2.bedTrack cimport BedTrack
from circlehunter2.pairedPeakTrack cimport PairedPeakTrack


cdef class PairedTagTrack:
    """
    Class representing PE reads alignment.
    """

    cdef:
        uint32_t _buffer_size

        dict _chroms
        dict _shifts

        dict _data

        uint32_t _pointer
        uint32_t _size

        bool _closed

    cpdef void add_tag_pair(self,
                            str chrom1, uint32_t start1, bint reverse1,
                            str chrom2, uint32_t start2, bint reverse2)

    cpdef void add_track(self, PairedTagTrack track)

    cpdef void close(self)

    cpdef set_chromosome_sizes(self, dict chrom_sizes)

    cpdef overlap(self, BedTrack track)

    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep)

    cdef uint32_t[:] _merge_into_regions(self, uint32_t distance)

    cpdef PairedPeakTrack call_peaks(self, uint32_t cutoff, uint32_t distance)

    cpdef void to_bedpe(self, object writable, uint32_t extend=?)
