# cython: language_level=3

import numpy as np

from numpy cimport uint32_t
cimport cython

from natsort import natsorted


cdef class BedTrack:
    """
    Class representing Chromosome regions.
    """

    def __init__(self, dict chrom_sizes, uint32_t buffer_size=256):
        assert buffer_size > 0

        self._buffer_size = buffer_size

        chroms = natsorted(chrom_sizes.keys())
        chroms_len = len(chroms)

        # init chromosome information
        self._chrom_indices = dict()
        self._index_chroms = dict()
        self._chrom_sizes = np.empty(chroms_len, dtype=np.uint32)
        self._left_limits = np.empty(chroms_len, dtype=np.uint32)
        self._right_limits = np.empty(chroms_len, dtype=np.uint32)

        shift = 0
        for i in range(chroms_len):
            chrom = chroms[i]

            self._chrom_indices[chrom] = i
            self._index_chroms[i] = chrom

            size = chrom_sizes[chrom]
            self._chrom_sizes[i] = size
            self._left_limits[i] = shift
            self._right_limits[i] = shift + size
            shift += size

        self._lefts = np.empty(self._buffer_size, dtype=np.uint32)
        self._rights = np.empty(self._buffer_size, dtype=np.uint32)

        self._pointer = 0
        self._size = self._buffer_size

        self._closed = False


    def get_chromosomes(self):
        """
        Chromosomes stored in this track.
        """
        return tuple(natsorted(self._chrom_indices.keys()))


    cpdef uint32_t get_chromosome_size(self, str chrom):
        """
        Return the size of the given chromosome
        """
        cdef uint32_t i = 0

        if chrom not in self._chrom_indices:
            raise KeyError(f"No such chromosome: {chrom}")

        i = self._chrom_indices[chrom]
        return self._chrom_sizes[i]


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef uint32_t _locate_position_chromosome(self, uint32_t position):
        """
        Find the chromosome index of the given position.
        """
        cdef long i = 0

        for i in range(self._right_limits.size):
            if position < self._right_limits[i]:
                break

        return i


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void add_region(self, str chrom, uint32_t start, uint32_t end):
        """
        Add a region into this track.
        """
        cdef:
            long j

            uint32_t shift

        assert not self._closed, "track closed!"

        assert chrom in self._chrom_indices, f"no such chromosome: {chrom}"
        assert start < end, "start >= end"
        assert end <= self.get_chromosome_size(chrom), f"end {end} out of chromosome bound: {self.get_chromosome_size(chrom)}"

        if self._pointer >= self._size:  # buffer size increase needed
            new_size = self._pointer + self._buffer_size
            self._lefts = np.resize(self._lefts, new_size)
            self._rights = np.resize(self._rights, new_size)
            self._size = new_size

        # find shift
        j = self._chrom_indices[chrom]
        shift = self._left_limits[j]

        # add tag
        self._lefts[self._pointer] = start + shift
        self._rights[self._pointer] = end + shift

        # next pointer
        self._pointer += 1


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void close(self):
        """
        Close this track.
        Indicating that there would not be any peak added into this track.
        All buffer will be released.
        """
        cdef:
            long new_size

        new_size = self._pointer
        # release memory
        self._lefts = np.resize(self._lefts, new_size)
        self._rights = np.resize(self._rights, new_size)
        self._size = new_size

        # sort by position
        indices = np.lexsort([self._lefts, self._rights])
        self._lefts = self._lefts.base[indices]
        self._rights = self._rights.base[indices]

        self._closed = True


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef set_chromosome_sizes(self, dict chrom_sizes):
        """
        Set the chromosome size of this track.
        Note: Regions out of this chromosome list or chromosome regions is filter out.
        """
        cdef:
            long new_chroms_len
            dict new_chrom_indices
            dict new_index_chroms

            uint32_t[:] new_chrom_sizes
            uint32_t[:] new_left_limits
            uint32_t[:] new_right_limits

            uint32_t shift
            uint32_t size

            uint32_t[:] new_lefts
            uint32_t[:] new_rights

            uint32_t new_size
            uint32_t new_shift

            uint32_t new_left
            uint32_t new_right

            long i
            long j
            long k

        assert self._closed, "Track not closed!"

        # calculate new shift
        new_chroms = natsorted(chrom_sizes.keys())
        new_chroms_len = len(new_chroms)
        new_chrom_indices = dict()
        new_index_chroms = dict()
        new_chrom_sizes = np.empty(new_chroms_len, dtype=np.uint32)
        new_left_limits = np.empty(new_chroms_len, dtype=np.uint32)
        new_right_limits = np.empty(new_chroms_len, dtype=np.uint32)

        shift = 0
        for i in range(new_chroms_len):
            chrom = new_chroms[i]

            new_chrom_indices[chrom] = i
            new_index_chroms[i] = chrom

            size = chrom_sizes[chrom]
            new_chrom_sizes[i] = size
            new_left_limits[i] = shift
            new_right_limits[i] = shift + size
            shift += size

        # new data
        new_lefts = np.empty(self._size, dtype=np.uint32)
        new_rights = np.empty(self._size, dtype=np.uint32)

        new_size = 0
        for i in range(self._size):
            # find old chrom and filter
            j = self._locate_position_chromosome(self._lefts[i])
            new_chrom = self._index_chroms[j]
            if new_chrom not in new_chrom_indices:
                continue

            # shift to new position
            k = new_chrom_indices[new_chrom]
            new_shift = new_left_limits[k]

            # check if left out of chromosome bound
            new_left = self._lefts[i] - self._left_limits[j] + new_shift
            if new_left < new_left_limits[k] or new_left >= new_right_limits[k]:
                continue

            new_right = self._rights[i] - self._left_limits[j] + new_shift
            # check if right out of chromosome bound
            if new_right >= new_right_limits[k]:
                continue

            new_lefts[new_size] = new_left
            new_rights[new_size] = new_right

            new_size += 1

        # release unused memory
        self._lefts = np.resize(new_lefts, new_size)
        self._rights = np.resize(new_rights, new_size)
        self._size = new_size

        self._chrom_indices = new_chrom_indices
        self._index_chroms = new_index_chroms
        self._chrom_sizes = new_chrom_sizes
        self._left_limits = new_left_limits
        self._right_limits = new_right_limits


    def get_data(self):
        """
        Return a 2D array containing all data of tags.
        """
        assert self._closed

        return np.rec.fromarrays(
            [self._lefts.base, self._rights.base], names=["lefts", "rights"]
        )


    cpdef void extend(self, uint32_t left, uint32_t right):
        """
        Extend regions in this track.
        """
        cdef:
            long i
            uint32_t j

            uint32_t new_left
            uint32_t new_right

        # loop through regions
        for i in range(self._size):
            j = self._locate_position_chromosome(self._lefts[i])
            new_left = self._lefts[i] - left
            if self._lefts[i] <= self._left_limits[j] + left:  # not out of bound
                new_left = self._left_limits[j]
            new_right = self._rights[i] + right
            if self._rights[i] >= self._right_limits[j] - right:  # not out of bound
                new_right = self._right_limits[j]
            self._lefts[i] = new_left
            self._rights[i] = new_right


    cpdef bint overlap(self, str chrom, uint32_t position, uint32_t extend=0):
        cdef:
            long i
            long j

            uint32_t shift
            uint32_t left
            uint32_t right

        assert self._closed, "track not closed!"

        assert chrom in self._chrom_indices, f"no such chromosome: {chrom}"

        j = self._chrom_indices[chrom]
        shift = self._left_limits[j]

        left = position + shift - extend
        right = position + shift + extend

        if left <= self._left_limits[j] + extend:
            left = self._left_limits[j]
        if right >= self._right_limits[j] - extend:
            right = self._right_limits[j]

        for i in range(self._size):
            if right < self._lefts[i]:
                continue
            if left >= self._rights[i]:
                break
            else:
                return 1
        return 0


    @classmethod
    def from_bed(cls, dict chrom_sizes, str filename, uint32_t buffer_size=256):
        track = BedTrack(chrom_sizes=chrom_sizes, buffer_size=buffer_size)
        with open(filename, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                chrom, start, end = line.strip().split("\t")
                track.add_region(chrom, int(start), int(end))
        track.close()
        return track

