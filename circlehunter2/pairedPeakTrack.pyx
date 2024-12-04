# cython: language_level=3

from pyfaidx import Fasta
import numpy as np

from numpy cimport uint32_t, float64_t, ndarray
cimport cython

from natsort import natsorted


cdef class PairedPeakTrack:
    """
    Class representing Peaks.
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

        self._lefts1 = np.empty(self._buffer_size, dtype=np.uint32)
        self._rights1 = np.empty(self._buffer_size, dtype=np.uint32)
        self._reverses1 = np.empty(self._buffer_size, dtype=np.uint8)
        self._lefts2 = np.empty(self._buffer_size, dtype=np.uint32)
        self._rights2 = np.empty(self._buffer_size, dtype=np.uint32)
        self._reverses2 = np.empty(self._buffer_size, dtype=np.uint8)
        self._counts = np.empty(self._buffer_size, dtype=np.uint32)

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
        cdef long i = 0

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
    cpdef void add_peak_pair(self,
                        str chrom1, uint32_t start1, uint32_t end1, bint reverses1,
                        str chrom2, uint32_t start2, uint32_t end2, bint reverses2,
                        uint32_t count):
        """
        Add a paired peak into this track.
        """
        cdef:
            long j
            long k

            uint32_t shift1
            uint32_t shift2

        assert not self._closed, "track closed!"

        assert chrom1 in self._chrom_indices, f"no such chromosome: {chrom1}"
        assert start1 < end1, "start1 >= end1"
        assert end1 <= self.get_chromosome_size(chrom1), f"end1 {end1} out of chromosome bound: {self.get_chromosome_size(chrom1)}"

        assert chrom2 in self._chrom_indices, f"no such chromosome: {chrom2}"
        assert start2 < end2, "start2 >= end2"
        assert end2 <= self.get_chromosome_size(chrom2), f"end2 {end2} out of chromosome bound: {self.get_chromosome_size(chrom2)}"

        if self._pointer >= self._size:  # buffer size increase needed
            new_size = self._pointer + self._buffer_size
            self._lefts1 = np.resize(self._lefts1, new_size)
            self._rights1 = np.resize(self._rights1, new_size)
            self._reverses1 = np.resize(self._reverses1, new_size)
            self._lefts2 = np.resize(self._lefts2, new_size)
            self._rights2 = np.resize(self._rights2, new_size)
            self._reverses2 = np.resize(self._reverses2, new_size)
            self._counts = np.resize(self._counts, new_size)
            self._size = new_size

        # find shift
        j = self._chrom_indices[chrom1]
        shift1 = self._left_limits[j]
        k = self._chrom_indices[chrom2]
        shift2 = self._left_limits[k]

        # add tag
        self._lefts1[self._pointer] = start1 + shift1
        self._rights1[self._pointer] = end1 + shift1
        self._reverses1[self._pointer] = reverses1
        self._lefts2[self._pointer] = start2 + shift2
        self._rights2[self._pointer] = end2 + shift2
        self._reverses2[self._pointer] = reverses2
        self._counts[self._pointer] = count

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
        self._lefts1 = np.resize(self._lefts1, new_size)
        self._rights1 = np.resize(self._rights1, new_size)
        self._reverses1 = np.resize(self._reverses1, new_size)
        self._lefts2 = np.resize(self._lefts2, new_size)
        self._rights2 = np.resize(self._rights2, new_size)
        self._reverses2 = np.resize(self._reverses2, new_size)
        self._counts = np.resize(self._counts, new_size)
        self._size = new_size

        self._closed = True


    @property
    def total(self):
        """
        Return peak count of this track.
        """
        assert self._closed

        return self._size


    def get_data(self):
        """
        Return a 2D array containing all data of tags.
        """
        assert self._closed

        return np.rec.fromarrays([
            self._lefts1.base, self._rights1.base, self._reverses1.base,
            self._lefts2.base, self._rights2.base, self._reverses2.base,
            self._counts.base
        ], names=["lefts1", "rights1", "reverses1", "lefts2", "rights2", "reverses2", "counts"])


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void keep_bi_peaks(self):
        """
        Look through peaks than only keep thoughts than are appears in both forward and reverse.
        """
        cdef:
            uint32_t[:] new_lefts1
            uint32_t[:] new_rights1
            uint8_t[:] new_reverses1
            uint32_t[:] new_lefts2
            uint32_t[:] new_rights2
            uint8_t[:] new_reverses2
            uint32_t[:] new_counts

            uint8_t[:] visited

            uint32_t new_size
            uint32_t overlap

            long i
            long j

        assert self._closed

        new_lefts1 = np.empty(self._size, dtype=np.uint32)
        new_rights1 = np.empty(self._size, dtype=np.uint32)
        new_reverses1 = np.empty(self._size, dtype=np.uint8)
        new_lefts2 = np.empty(self._size, dtype=np.uint32)
        new_rights2 = np.empty(self._size, dtype=np.uint32)
        new_reverses2 = np.empty(self._size, dtype=np.uint8)
        new_counts = np.empty(self._size, dtype=np.uint32)

        visited = np.zeros(self._size, dtype=np.uint8)

        new_size = 0
        for i in range(self._size):
            if visited[i] > 0:
                continue

            # new peak pair
            new_lefts1[new_size] = self._lefts1[i]
            new_rights1[new_size] = self._rights1[i]
            new_reverses1[new_size] = self._reverses1[i]
            new_lefts2[new_size] = self._lefts2[i]
            new_rights2[new_size] = self._rights2[i]
            new_reverses2[new_size] = self._reverses2[i]
            new_counts[new_size] = self._counts[i]

            # check if there is a bi peaks
            overlap = 0
            for j in range(i + 1, self._size):
                if visited[j] > 0:
                    continue
                # check if peaks overlap
                if new_reverses1[new_size] != self._reverses2[j]:
                    continue
                if new_reverses2[new_size] != self._reverses1[j]:
                    continue
                if new_lefts1[new_size] >= self._rights2[j]:
                    continue
                if new_rights1[new_size] <= self._lefts2[j]:
                    continue
                if new_lefts2[new_size] >= self._rights1[j]:
                    continue
                if new_rights2[new_size] <= self._lefts1[j]:
                    continue
                # merge peaks
                new_lefts1[new_size] = min(self._lefts2[j], new_lefts1[new_size])
                new_rights1[new_size] = max(self._rights2[j], new_rights1[new_size])
                new_reverses1[new_size] = self._reverses1[i]
                new_lefts2[new_size] = min(self._lefts1[j], new_lefts2[new_size])
                new_rights2[new_size] = max(self._rights1[j], new_rights2[new_size])
                new_reverses2[new_size] = self._reverses2[i]
                new_counts[new_size] = max(self._counts[j], new_counts[new_size])
                visited[j] = 1
                overlap += 1
            if overlap > 0:
                new_size += 1
        # assign data
        self._lefts1 = np.resize(new_lefts1, new_size)
        self._rights1 = np.resize(new_rights1, new_size)
        self._reverses1 = np.resize(new_reverses1, new_size)
        self._lefts2 = np.resize(new_lefts2, new_size)
        self._rights2 = np.resize(new_rights2, new_size)
        self._reverses2 = np.resize(new_reverses2, new_size)
        self._counts = np.resize(new_counts, new_size)
        self._size = new_size


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef duplicated(self):
        """
        Return a 1D array to indicate whether the peaks is duplicated.
        """
        cdef:
            uint8_t[:] duplicated

            uint32_t left1
            uint32_t right1
            uint8_t reverse1
            uint32_t left2
            uint32_t right2
            uint8_t reverse2

            long i
            long j
            long k

        assert self._closed

        duplicated = np.zeros(self._size, np.uint8)

        for i in range(self._size):
            if duplicated[i] > 0:
                continue
            k = i
            for j in range(i + 1, self._size):
                if duplicated[j] > 0:
                    continue

                # 1 vs 1, 2 vs 2
                left1 = self._lefts1[k] - self._lefts1[j]
                right1 = self._rights1[k] - self._rights1[j]
                reverse1 = self._rights1[k] - self._rights1[j]
                left2 = self._lefts2[k] - self._lefts2[j]
                right2 = self._rights2[k] - self._rights2[j]
                reverse2 = self._rights2[k] - self._rights2[j]

                if (left1 + right1 + left2 + right2) == 0 and (reverse1 + reverse2) == 0:
                    if self._counts[j] > self._counts[k]:
                        duplicated[k] = 1
                        k = j
                    else:
                        duplicated[j] = 1
                    continue

                # 1vs 2, 2 vs 1
                left1 = self._lefts1[k] - self._lefts2[j]
                right1 = self._rights1[k] - self._rights2[j]
                reverse1 = self._rights1[k] - self._rights2[j]
                left2 = self._lefts2[k] - self._lefts1[j]
                right2 = self._rights2[k] - self._rights1[j]
                reverse2 = self._rights2[k] - self._rights1[j]

                if (left1 + right1 + left2 + right2) == 0 and (reverse1 + reverse2) == 0:
                    if self._counts[j] > self._counts[k]:
                        duplicated[k] = 1
                        k = j
                    else:
                        duplicated[j] = 1
                    continue

        return duplicated.base


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef calculate_gc_percent(self, str fasta_file):
        """
        Return a 2D array containing GC content percent of each peak.
        """
        cdef:
            float64_t[:] gc_pct1
            float64_t[:] gc_pct2

            str chrom
            uint32_t start
            uint32_t end

            float64_t gc1
            float64_t gc2

            long i
            uint32_t j

        assert self._closed

        fasta = Fasta(fasta_file, one_based_attributes=False)

        gc_pct1 = np.empty(self._size, dtype=np.float64)
        gc_pct2 = np.empty(self._size, dtype=np.float64)

        for i in range(self._size):
            j = self._locate_position_chromosome(self._lefts1[i])
            chrom = self._index_chroms[j]
            start = self._lefts1[i] - self._left_limits[j]
            end = self._rights1[i] - self._left_limits[j]
            gc1 = fasta[chrom][start:end].gc

            j = self._locate_position_chromosome(self._lefts2[i])
            chrom = self._index_chroms[j]
            start = self._lefts2[i] - self._left_limits[j]
            end = self._rights2[i] - self._left_limits[j]
            gc2 = fasta[chrom][start:end].gc

            gc_pct1[i] = gc1
            gc_pct2[i] = gc2

        return np.rec.fromarrays([gc_pct1.base, gc_pct2.base], names=["gc_pct1", "gc_pct2"])


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep):
        """
        Filter peaks according to the given keep array.
        """
        cdef:
            uint32_t new_size
            uint32_t[:] new_lefts1
            uint32_t[:] new_rights1
            uint8_t[:] new_reverses1
            uint32_t[:] new_lefts2
            uint32_t[:] new_rights2
            uint8_t[:] new_reverses2
            uint32_t[:] new_counts

            uint8_t[:] keep_view

            long i
            long j


        assert self._closed
        assert keep.size == self._size

        # new data
        new_size = np.sum(keep)
        new_lefts1 = np.empty(new_size, dtype=np.uint32)
        new_rights1 = np.empty(new_size, dtype=np.uint32)
        new_reverses1 = np.empty(new_size, dtype=np.uint8)
        new_lefts2 = np.empty(new_size, dtype=np.uint32)
        new_rights2 = np.empty(new_size, dtype=np.uint32)
        new_reverses2 = np.empty(new_size, dtype=np.uint8)
        new_counts = np.empty(new_size, dtype=np.uint32)

        keep_view = keep
        # move data to new_data
        i = j = 0
        for i in range(self._size):
            if keep_view[i] > 0:
                new_lefts1[j] = self._lefts1[i]
                new_rights1[j] = self._rights1[i]
                new_reverses1[j] = self._reverses1[i]
                new_lefts2[j] = self._lefts2[i]
                new_rights2[j] = self._rights2[i]
                new_reverses2[j] = self._reverses2[i]
                new_counts[j] = self._counts[i]
                j += 1
        # assign data
        self._lefts1 = new_lefts1
        self._rights1 = new_rights1
        self._reverses1 = new_reverses1
        self._lefts2 = new_lefts2
        self._rights2 = new_rights2
        self._reverses2 = new_reverses2
        self._counts = new_counts
        self._size = new_size


    cpdef dict to_graph_data(self):
        """
        Generate a breakpoint graph using paired peaks store in this track.
        """
        cdef:
            str chrom1
            uint32_t start1
            uint32_t end1
            bool reverse1

            str chrom2
            uint32_t start2
            uint32_t end2
            bool reverse2

            uint32_t count

            long i
            uint32_t j
            uint32_t k

        assert self._closed

        graph = dict()

        # loop through tags
        for i in range(self._size):
            # extract data
            j = self._locate_position_chromosome(self._lefts1[i])
            chrom1 = self._index_chroms[j]
            start1 = self._lefts1[i] - self._left_limits[j]
            end1 = self._rights1[i] - self._left_limits[j]
            reverse1 = self._reverses1[i] > 0

            k = self._locate_position_chromosome(self._lefts2[i])
            chrom2 = self._index_chroms[k]
            start2 = self._lefts2[i] - self._left_limits[k]
            end2 = self._rights2[i] - self._left_limits[k]
            reverse2 = self._reverses2[i] > 0

            count = self._counts[i]

            graph[(chrom1, start1, end1, reverse1)] = {
                (chrom2, start2, end2, reverse2): {"DISCORDANT": {"count": count}}
            }

        return graph


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void to_bedpe(self, object writable):
        """
        Write peaks to a BEDPE format file.
        Columns: chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2
        """
        cdef:
            str chrom1
            uint32_t start1
            uint32_t end1
            str strand1

            str chrom2
            uint32_t start2
            uint32_t end2
            str strand2

            uint32_t score

            long order

            long i
            uint32_t j
            uint32_t k

        assert self._closed

        writable.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n")

        order = 1
        # loop through tags
        for i in range(self._size):

            # extract data
            j = self._locate_position_chromosome(self._lefts1[i])
            chrom1 = self._index_chroms[j]
            start1 = self._lefts1[i] - self._left_limits[j]
            end1 = self._rights1[i] - self._left_limits[j]
            strand1 = "-" if self._reverses1[i] else "+"

            k = self._locate_position_chromosome(self._lefts2[i])
            chrom2 = self._index_chroms[k]
            start2 = self._lefts2[i] - self._left_limits[k]
            end2 = self._rights2[i] - self._left_limits[k]
            strand2 = "-" if self._reverses2[i] else "+"

            name = f"paired_peak_{order}"
            score = self._counts[i]

            writable.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{strand1}\t{strand2}\n")

            order += 1
            # end peak
        # end for loop


    @staticmethod
    def from_bedpe(chrom_sizes, readable, buffer_size=256):
        track = PairedPeakTrack(chrom_sizes=chrom_sizes, buffer_size=buffer_size)
        while True:
            line = readable.readline()
            if not line:
                break
            if line.startswith("#"):
                continue
            chrom1, start1, end1, chrom2, start2, end2, name, count, reverse1, reverse2 = line.strip().split("\t")
            start1, end1, start2, end2, count = map(int, [start1, end1, start2, end2, count])
            reverse1, reverse2 = map(lambda x: x == "-", [reverse1, reverse2])
            track.add_peak_pair(chrom1, start1, end1, reverse1, chrom2, start2, end2, reverse2, count)
        track.close()
        return track


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void to_bed(self, object writable):
        """
        Write peaks to a BED format file.
        Columns: chrom1, start1, end1, name-1, score, strand1
                 chrom2, start2, end2, name-2, score, strand2
        Might use to visualize the enrich region.
        """
        cdef:
            str chrom1
            uint32_t start1
            uint32_t end1
            str strand1

            str chrom2
            uint32_t start2
            uint32_t end2
            str strand2

            str name1
            str name2

            uint32_t score

            long order

            long i
            uint32_t j
            uint32_t k

        assert self._closed

        writable.write("#chrom\tstart\tend\tname\tscore\tstrand\n")

        order = 1
        # loop through tags
        for i in range(self._size):
            # extract data
            j = self._locate_position_chromosome(self._lefts1[i])
            chrom1 = self._index_chroms[j]
            start1 = self._lefts1[i] - self._left_limits[j]
            end1 = self._rights1[i] - self._left_limits[j]
            strand1 = "-" if self._reverses1[i] else "+"

            k = self._locate_position_chromosome(self._lefts2[i])
            chrom2 = self._index_chroms[k]
            start2 = self._lefts2[i] - self._left_limits[k]
            end2 = self._rights2[i] - self._left_limits[k]
            strand2 = "-" if self._reverses2[i] else "+"

            name1 = f"paired_peak_{order}_1"
            name2 = f"paired_peak_{order}_2"
            score = self._counts[i]

            writable.write((
                f"{chrom1}\t{start1}\t{end1}\t{name1}\t{score}\t{strand1}\n"
                f"{chrom2}\t{start2}\t{end2}\t{name2}\t{score}\t{strand2}\n"
            ))

            order += 1
            # end peak
        # end for loop