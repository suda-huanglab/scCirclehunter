# cython: language_level=3

from numpy cimport uint8_t, uint32_t, uint64_t, ndarray
from cpython cimport bool
import numpy as np
cimport numpy as np
cimport cython

from circlehunter2.utils import nasort_key
from circlehunter2.utils import make_windows


cdef class TagTrack:
    """
    Class representing all tags from one sample or cell.
    """
    def __init__(self, uint32_t buffer_size=256):
        self.buffer_size = buffer_size
        self._chroms = dict()
        self._starts = dict()
        self._ends = dict()
        self._reverses = dict()
        self._inserts = dict()
        self._pointers = dict()
        self._sizes = dict()
        self._total = 0
        self._total_length = 0
        self._closed = False


    # optimize performance by turn off these checks
    @cython.nonecheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void add_tag(self, str chrom, uint32_t start, uint32_t end, bint reverse, uint32_t insert):
        """
        Add one tag to this track.
        """
        cdef long pointer

        assert not self._closed
        assert start < end

        if chrom not in self._sizes:  # new chrom
            pointer = 0
            self._starts[chrom] = np.empty(self.buffer_size, dtype=np.uint32)
            self._ends[chrom] = np.empty(self.buffer_size, dtype=np.uint32)
            self._reverses[chrom] = np.empty(self.buffer_size, dtype=np.bool_)
            self._inserts[chrom] = np.empty(self.buffer_size, dtype=np.uint32)
            self._sizes[chrom] = self.buffer_size
        else:
            pointer = self._pointers[chrom] + 1
            if pointer >= self._sizes[chrom]:  # buffer size increase needed
                self._starts[chrom].resize(pointer + self.buffer_size)
                self._ends[chrom].resize(pointer + self.buffer_size)
                self._reverses[chrom].resize(pointer + self.buffer_size)
                self._inserts[chrom].resize(pointer + self.buffer_size)
                self._sizes[chrom] += self.buffer_size

        # add tag
        self._starts[chrom][pointer] = start
        self._ends[chrom][pointer] = end
        self._reverses[chrom][pointer] = reverse
        self._inserts[chrom][pointer] = insert
        self._pointers[chrom] = pointer

        # add stats
        self._total += 1
        self._total_length += <uint64_t> (end - start)


    cpdef  void add_chromosome(
            self, str chrom, ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends,
            ndarray[uint32_t, ndim=1] reverses, ndarray[uint32_t, ndim=1] inserts):
        """
        Add a chromosome to this track.
        """
        assert not self._closed, "track closed"
        assert starts.shape[0] == ends.shape[0], "shape not equal"
        assert starts.shape[0] == reverses.shape[0], "shape not equal"
        assert starts.shape[0] == inserts.shape[0], "shape not equal"

        size = starts.shape[0]

        if chrom not in self._sizes:  # new chrom
            pointer = 0
            self._starts[chrom] = np.empty(size, dtype=np.uint32)
            self._ends[chrom] = np.empty(size, dtype=np.uint32)
            self._reverses[chrom] = np.empty(size, dtype=np.bool_)
            self._inserts[chrom] = np.empty(size, dtype=np.uint32)
            self._sizes[chrom] = size
        else:
            pointer = self._pointers[chrom] + 1
            if pointer + size >= self._sizes[chrom]:  # buffer size increase need
                self._starts[chrom].resize(pointer + size)
                self._ends[chrom].resize(pointer + size)
                self._reverses[chrom].resize(pointer + size)
                self._inserts[chrom].resize(pointer + size)
                self._sizes[chrom] = pointer + size

        # assign values
        self._starts[chrom][pointer:pointer + size] = starts
        self._ends[chrom][pointer:pointer + size] = ends
        self._reverses[chrom][pointer:pointer + size] = reverses
        self._inserts[chrom][pointer:pointer + size] = inserts
        self._pointers[chrom] = pointer + size - 1

        # add stats
        self._total += size
        self._total_length += np.sum(ends - starts)


    cpdef void add_track(self, TagTrack track):
        """
        Add a track to this track.
        """
        assert not self._closed

        # loop through all chroms
        chroms = track.get_chromosomes()
        for chrom in chroms:
            data = track.get_data(chrom)
            self.add_chromosome(chrom, data[0], data[1], data[2].astype(np.uint32), data[3])


    cpdef void close(self):
        """
        Close this track.
        Indicating that there would not be any tag added into this track.
        All buffer will be released.
        """
        cdef:
            str chrom
            uint32_t pointer
            uint32_t size

        assert not self._closed

        for chrom, pointer in self._pointers.items():
            # trim empty buffer size
            size = pointer + 1
            self._starts[chrom].resize(size)
            self._ends[chrom].resize(size)
            self._reverses[chrom].resize(size)
            self._inserts[chrom].resize(size)
            # max chrom position as last base here
            self._chroms[chrom] = np.max(self._ends[chrom])

        # release memory
        self._pointers = None
        self._sizes = None

        # mark as closed
        self._closed = True


    cpdef set_chromosome_sizes(self, dict chrom_sizes):
        """
        Set the chromosome size of this track.
        Note: Tags out of this chromosome list or chromosome regions is filter out.
        """
        assert self._closed, "Track not closed!"

        self._total = 0
        self._total_length = 0

        chroms = list(self._chroms.keys())
        # loop through chromosomes
        for chrom in chroms:
            if chrom in chrom_sizes:
                size = chrom_sizes[chrom]
                # filter out tags
                indices = (self._starts[chrom] < size) & (self._ends[chrom] <= size)
                self._starts[chrom] = self._starts[chrom][indices]
                self._ends[chrom] = self._ends[chrom][indices]
                self._reverses[chrom] = self._reverses[chrom][indices]
                self._inserts[chrom] = self._inserts[chrom][indices]
                # re-calculate total
                self._total += self._starts[chrom].shape[0]
                self._total_length += np.sum(self._ends[chrom] - self._starts[chrom])
                # set size
                self._chroms[chrom] = size
            else:
                # remove chromosomes
                self._chroms.pop(chrom)
                self._starts.pop(chrom)
                self._ends.pop(chrom)
                self._reverses.pop(chrom)
                self._inserts.pop(chrom)


    def get_chromosomes(self):
        """
        Chromosomes stored in this track.
        """
        assert self._closed

        return tuple(sorted(self._chroms.keys(), key=nasort_key))


    @property
    def total(self):
        """
        Tag counts in this track
        """
        assert self._closed

        return self._total


    @property
    def total_length(self):
        """
        Total tag length in this track
        """
        assert self._closed

        return self._total_length


    def get_chromosome_size(self, str chrom):
        """
        Return the length of the given chromosome
        """
        assert self._closed

        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")

        return self._chroms[chrom]


    def get_data(self, str chrom):
        """
        Return a 2D array containing all data of tags from the given chromosome.
        """
        assert self._closed

        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")

        return np.stack([
            self._starts[chrom], self._ends[chrom],
            self._reverses[chrom], self._inserts[chrom]
        ])


    cpdef TagTrack filter(self, bool reverse, uint32_t min_insert):
        """
        Filter tags by strand and insert size.
        Note: This func return a new instance
        """
        assert self._closed

        track = TagTrack(buffer_size=self.buffer_size)

        # loop through chromosomes
        chroms = self.get_chromosomes()
        for chrom in chroms:
            # filter indices
            indices_reverse = (self._reverses[chrom] == reverse)
            indices_insert = (self._inserts[chrom] > min_insert)
            indices = indices_reverse & indices_insert
            if not np.any(indices):  # not tags remain
                continue
            # extract data
            starts = self._starts[chrom][indices]
            ends = self._ends[chrom][indices]
            reverses = self._reverses[chrom][indices].astype(np.uint32)
            inserts = self._inserts[chrom][indices]
            # add data
            track.add_chromosome(chrom, starts, ends, reverses, inserts)

        track.close()

        # reset chromosome sizes
        chrom_sizes = {
            chrom: self.get_chromosome_size(chrom) for chrom in chroms
        }
        track.set_chromosome_sizes(chrom_sizes)

        return track


    cpdef void extend(self, uint32_t left, uint32_t right):
        """
        Extend tags stored in this track.
        """
        cdef:
            str chrom
            long i
            uint32_t size
            uint32_t start_limit
            uint32_t end_limit
            uint32_t[:] starts
            uint32_t[:] ends

        assert self._closed, "track not closed!"

        self._total_length = 0

        # loop through chromosomes
        for chrom in self.get_chromosomes():
            # do not extend beyond chromosome size
            size = self._chroms[chrom]
            start_limit = left
            end_limit = size - right

            starts = self._starts[chrom]
            ends = self._ends[chrom]

            # loop through tags
            for i in range(self._starts[chrom].size):
                # left extend
                if starts[i] <= start_limit:
                    starts[i] = 0
                else:
                    starts[i] -= left
                # right extend
                if ends[i] >= end_limit:
                    ends[i] = size
                else:
                    ends[i] += right
                # update total length
                self._total_length += ends[i] - starts[i]
            # end chromosome


    cpdef tuple make_windows(self, window_size=100_000):
        """
        Return the windows of chromosomes by the given windows size
        """
        return make_windows(self._chroms, window_size)


    cpdef uint32_t count(self, str chrom, uint32_t start, uint32_t end):
        assert self._closed, "track not closed!"
        # assert chrom in self._chroms, f"{chrom} not in track"
        assert start < end

        if chrom not in self._chroms:
            return 0
        return np.sum((self._starts[chrom] >= start) & (self._starts[chrom] < end))


    cpdef ndarray[uint32_t, ndim=1] count_windows(self, uint32_t window_size=100_000):
        """
        Count tag fall in windows make by the given window size.
        """
        cdef:
            tuple chrom_windows
            tuple chroms
            str chrom
            dict windows_counts
            uint32_t count
            uint32_t shift
            ndarray[uint32_t, ndim=1] windows
            ndarray[long, ndim=1] counts
            ndarray[uint32_t, ndim=1] total_counts

        # hom many windows?
        windows_counts = dict()
        chrom_windows = self.make_windows(window_size)
        for chrom, *_ in iter(chrom_windows):
            count = windows_counts.get(chrom, 0) + 1
            windows_counts[chrom] = count

        # init counts array
        total_counts = np.zeros(np.sum(list(windows_counts.values())), dtype=np.uint32)

        # loop through chromosomes
        shift = 0
        chroms = self.get_chromosomes()
        for chrom in iter(chroms):
            # tag fall into windows
            windows = self._starts[chrom] // window_size
            # count windows and assign to counts array
            windows, counts = np.unique(windows, return_counts=True)
            total_counts[windows + shift] += counts.astype(np.uint32)
            # shift to next chromosome
            count = windows_counts[chrom]
            shift += count

        return total_counts


    cpdef void to_bed(self, object writable, str name_prefix):
        """
        Write tags to a BED format file file.
        Columns: chrom, start, end, name, insert_size, strand
        """
        cdef:
            str chrom
            uint32_t[:] starts
            uint32_t[:] ends
            uint8_t[:] reverses
            uint32_t[:] inserts

            long order
            long i

            uint32_t start
            uint32_t end
            str strand
            uint32_t insert
            str name

        assert self._closed

        writable.write("#chrom\tstart\tend\tname\tinsert_size\tstrand\n")

        # loop through chromosomes
        for chrom in self.get_chromosomes():
            starts = self._starts[chrom]
            ends = self._ends[chrom]
            reverses = self._reverses[chrom]
            inserts = self._inserts[chrom]

            order =1
            # loop through tags
            for i in range(self._starts[chrom].size):
                start = starts[i]
                end = ends[i]
                strand = "-" if reverses[i] else "+"
                insert = inserts[i]
                name = f"{name_prefix}_{order}"

                writable.write(f"{chrom}\t{start}\t{end}\t{name}\t{insert}\t{strand}\n")

                order += 1
            # end chromosome
