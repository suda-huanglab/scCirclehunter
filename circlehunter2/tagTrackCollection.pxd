# cython: language_level=3

from numpy cimport uint32_t, uint64_t
from cpython cimport bool

from circlehunter2.tagTrack cimport TagTrack


cdef class TagTrackCollection:
    cdef:
        uint32_t buffer_size
        dict _chroms
        dict _tracks
        uint64_t _total
        uint64_t _total_length
        bool _closed

    cpdef void add_tag(
            self, str barcode, str chrom, uint32_t start,
            uint32_t end, bint reverse, uint32_t insert)

    cpdef void add_track(self, str barcode, TagTrack track)

    cpdef void add_collection(self, TagTrackCollection collection)

    cpdef void close(self)

    cpdef set_chromosome_sizes(self, dict chrom_sizes)

    cpdef void set_barcodes(self, list barcodes)

    cpdef tuple make_windows(self, window_size=?)

    cpdef count_windows(self, window_size=?)

    cpdef count(self, str chrom, uint32_t start, uint32_t end)

    cpdef TagTrackCollection filter(self, bool reverse, uint32_t min_insert)

    cpdef pool_barcodes(self, list barcodes)

    cpdef void extend(self, uint32_t left, uint32_t right)
