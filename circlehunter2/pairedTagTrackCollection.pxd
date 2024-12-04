# cython: language_level=3

from numpy cimport uint32_t
from cpython cimport bool

from circlehunter2.pairedTagTrack cimport PairedTagTrack


cdef class PairedTagTrackCollection:
    cdef:
        uint32_t _buffer_size

        dict _chroms
        dict _shifts

        dict _tracks

        uint32_t _total

        bool _closed

    cpdef void add_tag_pair(self, str barcode,
                            str chrom1, uint32_t start1, bint reverse1,
                            str chrom2, uint32_t start2, bint reverse2)

    cpdef void add_track(self, str barcode, PairedTagTrack track)

    cpdef void add_collection(self, PairedTagTrackCollection collection)

    cpdef void close(self)

    cpdef set_chromosome_sizes(self, dict chrom_sizes)

    cpdef void set_barcodes(self, list barcodes)

    cpdef PairedTagTrack pool_barcodes(self, list barcodes)
