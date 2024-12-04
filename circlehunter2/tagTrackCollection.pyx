# cython: language_level=3

from scipy.sparse import coo_matrix
from natsort import natsorted

from numpy cimport uint32_t
from cpython cimport bool
import numpy as np
cimport numpy as np
cimport cython

from circlehunter2.tagTrack cimport TagTrack
from circlehunter2.bedGraphTrack import BedGraphTrack
from circlehunter2.utils cimport _pileup as pileup_as_depth

from circlehunter2.utils import make_windows

cdef class TagTrackCollection:
    """
    Class representing all tags from multi samples or cells.
    """
    def __init__(self, uint32_t buffer_size=256):
        self.buffer_size = buffer_size
        self._tracks = dict()
        self._chroms = dict()
        self._total = 0
        self._total_length = 0
        self._closed = False


    cpdef void add_tag(
            self, str barcode, str chrom, uint32_t start,
            uint32_t end, bint reverse, uint32_t insert):
        """
        Add tag to this collection
        """
        assert not self._closed

        if barcode not in self._tracks:
            self._tracks[barcode] = TagTrack(buffer_size=self.buffer_size)
        self._tracks[barcode].add_tag(chrom, start, end, reverse, insert)


    cpdef void add_track(self, str barcode, TagTrack track):
        """
        Add track to this collection
        """
        assert not self._closed

        if barcode not in self._tracks:  # new barcode
            self._tracks[barcode] = TagTrack(buffer_size=self.buffer_size)
        # add track
        self._tracks[barcode].add_track(track)


    cpdef void add_collection(self, TagTrackCollection collection):
        """
        Add collection to this collection
        """
        assert not self._closed

        barcodes = collection.get_barcodes()
        # loop through barcodes
        for barcode in barcodes:
            track = collection.get_track(barcode)
            self.add_track(barcode, track)


    cpdef void close(self):
        """
        Close this collection.
        Indicating that there would not be any tag added into this collection.
        All buffer will be released.
        """
        cdef:
            str barcode

        assert not self._closed

        # close all tracks
        for track in self._tracks.values():
            track.close()

            # chrom length
            chroms = track.get_chromosomes()
            for chrom in chroms:
                length = track.get_chromosome_size(chrom)
                self._chroms[chrom] = max(length, self._chroms.get(chrom, 0))

            # stats
            self._total += track.total
            self._total_length += track.total_length

        # mark as closed
        self._closed = True

        # this have to set after closed!
        self.set_chromosome_sizes(self._chroms)


    cpdef set_chromosome_sizes(self, dict chrom_sizes):
        """
        Set the chromosome size of this collection.
        Note: Tags out of this chromosome list or chromosome regions is filter out.
        """
        assert self._closed, "Track not closed!"

        self._total = 0
        self._total_length = 0

        barcodes = list(self._tracks.keys())

        # loop through barcodes
        for barcode in barcodes:
            self._tracks[barcode].set_chromosome_sizes(chrom_sizes)
            self._total += self._tracks[barcode].total
            self._total_length += self._tracks[barcode].total_length

        # set size of this collection
        chroms = list(self._chroms.keys())
        for chrom in chroms:
            if chrom in chrom_sizes:
                size = chrom_sizes[chrom]
                self._chroms[chrom] = size
            else:
                self._chroms.pop(chrom)


    cpdef void set_barcodes(self, list barcodes):
        """
        Set the barcodes of this collection.
        Note: Barcodes not in the input list is filter out.
        """
        assert self._closed

        input_barcodes = set(barcodes)

        exist_barcodes = list(self._tracks.keys())
        # loop through barcodes
        for barcode in exist_barcodes:
            if barcode not in input_barcodes:
                # filter out barcodes
                track = self._tracks.pop(barcode)
                # re-calculate stats
                self._total -= track.total
                self._total_length -= track.total_length
                # release memory
                track = None


    def get_chromosomes(self):
        """
        Chromosomes stored in this collection.
        """
        assert self._closed

        return tuple(self._chroms.keys())


    def get_barcodes(self):
        """
        Chromosomes stored in this collection.
        """
        assert self._closed

        return tuple(natsorted(self._tracks.keys()))


    @property
    def total(self):
        """
        Tag counts in this collection
        """
        assert self._closed, "track not closed!"

        return self._total


    @property
    def total_length(self):
        """
        Total tag length in this collection
        """
        assert self._closed, "track not closed!"

        return self._total_length


    def get_chromosome_size(self, str chrom):
        """
        Return the length of the given chromosome
        """
        assert self._closed, "track not closed!"

        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")

        return self._chroms[chrom]


    def get_track(self, str barcode):
        """
        Return the track associated with given barcode
        """
        assert self._closed, "track not closed!"

        if barcode not in self._tracks:
            raise KeyError(f"No such track associated with barcode: {barcode}")

        return self._tracks[barcode]


    cpdef void extend(self, uint32_t left, uint32_t right):
        """
        Extend tags stored in this track.
        """
        cdef:
            str  barcode
            TagTrack track

        assert self._closed, "track not closed!"

        self._total_length = 0

        # loop through barcodes
        barcodes = list(self._tracks.keys())

        # loop through barcodes
        for barcode in barcodes:
            track = self.get_track(barcode)
            track.extend(left, right)
            self._total_length += track.total_length


    cpdef tuple make_windows(self, window_size=100_000):
        """
        Return the windows of chromosomes by the given windows size
        """
        return make_windows(self._chroms, window_size)


    cpdef count_windows(self, window_size=100_000):
        """
        Count tag fall in windows make by the given window size.
        This func return a sparse matrix with windows as rows and barcode as columns.
        Can be used to generate AnnData object.
        """
        cdef:
            tuple chrom_windows
            tuple chrom_window
            dict window_indices
            tuple barcodes
            str barcode
            dict barcode_indices
            long i
            long j
            uint32_t[:] counts
            uint32_t[:] indices
            list I, J, V

        chrom_windows = self.make_windows(window_size)
        window_indices = dict()
        i = 0
        for chrom_window in iter(chrom_windows):
            window_indices[chrom_window] = i
            i += 1

        I = list()
        J = list()
        V = list()

        barcodes = self.get_barcodes()
        i = 0
        for barcode in iter(barcodes):
            track = self.get_track(barcode)

            counts = track.count_windows(window_size)
            V.append(counts)

            I.append(np.full_like(counts, i))
            i += 1

            indices = np.empty_like(counts)
            chrom_windows = track.make_windows(window_size)
            j = 0
            for chrom_window in iter(chrom_windows):
                indices[j] = window_indices[chrom_window]
                j += 1
            J.append(indices.base)

        matrix = (np.hstack(V), (np.hstack(I), np.hstack(J)))

        return coo_matrix(matrix, shape=(i, len(window_indices))).tocsr()


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef count(self, str chrom, uint32_t start, uint32_t end):
        cdef:
            long i
            TagTrack track
            uint32_t[:] counts

        assert self._closed
        assert chrom in self._chroms
        assert start < end

        barcodes = self.get_barcodes()
        counts = np.empty(len(barcodes), dtype=np.uint32)
        for i in range(len(barcodes)):
            track = self.get_track(barcodes[i])
            counts[i] = track.count(chrom, start, end)

        return counts.base


    cpdef pool_barcodes(self, list barcodes):
        """
        Pool tracks within the given barcodes and pileup as BedGraphTrack
        """
        assert self._closed, "track not closed!"

        # barcodes intersect
        pool = set(barcodes) & set(self.get_barcodes())

        # prepare bedGraph
        pileup = BedGraphTrack(chrom_sizes=self._chroms)

        # loop through chromosomes
        chroms = self.get_chromosomes()
        for chrom in chroms:
            # collect positions from tracks
            positions_list = list()
            for barcode in pool:
                track = self.get_track(barcode)
                try:
                    data = track.get_data(chrom)
                    positions_list.append(data[:2])
                except KeyError:  # chrom not exist in this track
                    continue

            if len(positions_list) > 0:  # no such chromosome
                # concat positions
                positions = np.hstack(positions_list)
                # pileup and add to bedGraph
                pos, val = pileup_as_depth(positions[0], positions[1])
                pileup.add_chromosome(chrom, pos, val)

        pileup.close()

        return pileup


    cpdef TagTrackCollection filter(self, bool reverse, uint32_t min_insert):
        """
        Filter tags by strand and insert size.
        Note: This func return a new instance
        """
        assert self._closed, "track not closed!"

        collection = TagTrackCollection(buffer_size=self.buffer_size)

        # loop through barcodes
        barcodes = self.get_barcodes()
        for barcode in barcodes:
            track = self.get_track(barcode)
            filtered = track.filter(reverse, min_insert)
            collection.add_track(barcode, filtered)

        collection.close()

        # reset chromosome sizes
        chroms = self.get_chromosomes()
        chrom_sizes = {
            chrom: self.get_chromosome_size(chrom) for chrom in chroms
        }
        collection.set_chromosome_sizes(chrom_sizes)

        return collection
