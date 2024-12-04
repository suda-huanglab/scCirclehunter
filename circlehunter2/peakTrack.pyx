# cython: language_level=3

from collections import OrderedDict

import numpy as np

from numpy cimport uint8_t, uint32_t, float64_t, ndarray
cimport cython

from natsort import natsorted


cdef class PeakTrack:
    """
    Class representing Peaks.
    """

    def __init__(self, dict chrom_sizes, uint32_t buffer_size=256):
        self._buffer_size = buffer_size

        self._chroms = OrderedDict()
        chroms = natsorted(chrom_sizes.keys())
        for chrom in chroms:
            self._chroms[chrom] = chrom_sizes[chrom]

        self._pointers = dict()
        self._sizes = dict()
        self._tmp_starts = dict()
        self._tmp_ends = dict()
        self._tmp_max_logp = dict()
        self._tmp_max_fc = dict()
        self._tmp_mean_logp = dict()
        self._tmp_mean_fc = dict()
        self._counts = dict()

        self._closed = False


    cpdef void add_peak(self, str chrom, uint32_t start, uint32_t end,
                        float64_t max_logp, float64_t max_fc,
                        float64_t mean_logp, float64_t mean_fc):
        """
        Add a peak into this track.
        Note: peak can't overlap each other!
        """
        cdef long pointer

        # track not closed
        assert not self._closed, "track closed!"
        # chromosome size exists
        assert chrom in self._chroms, f"no such chromosome: {chrom}"
        # chromosome size match
        assert end <= self._chroms[chrom], f"end position {end} larger than chromosome size: {self._chroms[chrom]}"
        # start smaller than end
        assert start < end, f"start larger than end: {start}-{end}"

        if chrom not in self._pointers:  # new chrom
            pointer = 0
            self._tmp_starts[chrom] = np.empty(self._buffer_size, dtype=np.uint32)
            self._tmp_ends[chrom] = np.empty(self._buffer_size, dtype=np.uint32)
            self._tmp_max_logp[chrom] = np.empty(self._buffer_size, dtype=np.float64)
            self._tmp_max_fc[chrom] = np.empty(self._buffer_size, dtype=np.float64)
            self._tmp_mean_logp[chrom] = np.empty(self._buffer_size, dtype=np.float64)
            self._tmp_mean_fc[chrom] = np.empty(self._buffer_size, dtype=np.float64)
            self._sizes[chrom] = self._buffer_size
        else:
            pointer = self._pointers[chrom] + 1
            if pointer >= self._sizes[chrom]:  # buffer size increase needed
                self._tmp_starts[chrom].resize(pointer + self._buffer_size)
                self._tmp_ends[chrom].resize(pointer + self._buffer_size)
                self._tmp_max_logp[chrom].resize(pointer + self._buffer_size)
                self._tmp_max_fc[chrom].resize(pointer + self._buffer_size)
                self._tmp_mean_logp[chrom].resize(pointer + self._buffer_size)
                self._tmp_mean_fc[chrom].resize(pointer + self._buffer_size)
                self._sizes[chrom] += self._buffer_size

        # add peak
        self._tmp_starts[chrom][pointer] = start
        self._tmp_ends[chrom][pointer] = end
        self._tmp_max_logp[chrom][pointer] = max_logp
        self._tmp_max_fc[chrom][pointer] = max_fc
        self._tmp_mean_logp[chrom][pointer] = mean_logp
        self._tmp_mean_fc[chrom][pointer] = mean_fc
        self._pointers[chrom] = pointer


    cpdef close(self):
        """
        Close this track.
        Indicating that there would not be any peak added into this track.
        All buffer will be released.
        """
        cdef:
            long i
            str chrom
            uint32_t count
            uint32_t left
            uint32_t right
            uint32_t last
            uint32_t[:] positions
            float64_t[:] values

        assert not self._closed, "track closed!"

        # count data size
        count = 0
        for chrom, pointer in self._pointers.items():
            self._counts[chrom] = pointer + 1
            count += pointer + 1

        # init data memory
        self._starts = np.empty(count, dtype=np.uint32)
        self._ends = np.empty(count, dtype=np.uint32)
        self._max_logp = np.empty(count, dtype=np.float64)
        self._max_fc = np.empty(count, dtype=np.float64)
        self._mean_logp = np.empty(count, dtype=np.float64)
        self._mean_fc = np.empty(count, dtype=np.float64)

        chroms = self.get_chromosomes()
        last = 0
        left = 0
        # loop through chromosomes
        for i in range(len(chroms)):
            chrom = chroms[i]
            if chrom in self._counts:
                # assign values
                right = left + self._counts[chrom]
                count = self._counts[chrom]
                # kind of wire? but error if assign directly
                positions = self._tmp_starts[chrom][:count] + last
                self._starts[left:right] = positions
                positions = self._tmp_ends[chrom][:count] + last
                self._ends[left:right] = positions
                values = self._tmp_max_logp[chrom][:count]
                self._max_logp[left:right] = values
                values = self._tmp_max_fc[chrom][:count]
                self._max_fc[left:right] = values
                values = self._tmp_mean_logp[chrom][:count]
                self._mean_logp[left:right] = values
                values = self._tmp_mean_fc[chrom][:count]
                self._mean_fc[left:right] = values
                # shift
                left += self._counts[chrom]
            # end of chromosome
            last += self._chroms[chrom]

        # release buffers
        self._tmp_starts = None
        self._tmp_ends = None
        self._tmp_max_logp = None
        self._tmp_max_fc = None
        self._tmp_mean_logp = None
        self._tmp_mean_fc = None

        # sort by positions
        indices = np.lexsort((self._starts.base, self._ends.base))
        self._starts = self._starts.base[indices]
        self._ends = self._ends.base[indices]
        self._max_logp = self._max_logp.base[indices]
        self._max_fc = self._max_fc.base[indices]
        self._mean_logp = self._mean_logp.base[indices]
        self._mean_fc = self._mean_fc.base[indices]

        # check duplicates?
        groups, counts = np.unique(
            np.stack([self._starts, self._ends]), axis=1, return_counts=True
        )
        assert np.alltrue(counts == 1), "peak duplicated!"

        # mark as closed
        self._closed = True


    def get_chromosomes(self):
        """
        Chromosomes stored in this track.
        """
        return tuple(natsorted(self._chroms.keys()))


    def get_chromosome_size(self, str chrom):
        """
        Return the size of the given chromosome
        """
        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")
        return self._chroms[chrom]


    def get_chromosome_peaks(self, str chrom):
        """
        Return the peak counts of the given chromosome
        """
        assert self._closed

        if chrom not in self._counts:
            if chrom not in self._chroms:
                raise KeyError(f"No such chromosome: {chrom}")
            else:
                return 0

        return self._counts[chrom]


    def get_regions(self):
        """
        Return region of peaks stored in this track.
        """
        assert self._closed

        return np.stack([self._starts.base, self._ends.base])


    def get_significances(self):
        """
        Return significance of peaks stored in this track.
        """
        assert self._closed

        return np.stack([
            self._max_logp.base, self._max_fc.base,
            self._mean_logp.base, self._mean_fc.base
        ])


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef ndarray[uint32_t, ndim=1] count_overlaps(self, PeakTrack other):
        """
        Return how many overlap of each region in this track between the other.
        """
        cdef:
            long i_self
            long i_other
            long self_counts
            long region_counts
            uint32_t[:, :] regions
            uint32_t[:] counts

        assert self._closed

        # check chromosome sizes
        chroms1 = self.get_chromosomes()
        chroms2 = other.get_chromosomes()

        assert len(chroms1) == len(chroms2), "chromosomes not match"

        for chrom1, chrom2 in zip(chroms1, chroms2):
            assert chrom1 == chrom2

            size1 = self.get_chromosome_size(chrom1)
            size2 = other.get_chromosome_size(chrom2)

            assert size1 == size2
        # check pass

        # extract data
        regions = other.get_regions()

        self_counts = self._starts.size
        other_counts = regions.shape[1]

        counts = np.zeros_like(self._starts)

        i_self = 0
        i_other = 0

        # loop through regions
        while i_other < other_counts and i_self < self_counts:
            if regions[1, i_other] < self._starts[i_self]:  # region at left
                i_other += 1
            elif regions[0, i_other] > self._ends[i_self]:  # region at right
                i_self += 1
            else:  # overlap
                counts[i_self] += 1
                # move to next when region end first
                if regions[1, i_other] < self._ends[i_self]:
                    i_other += 1
                elif regions[1, i_other] > self._ends[i_self]:
                    i_self += 1
                else: # move to next together
                    i_other += 1
                    i_self += 1

        return counts.base


    cpdef float64_t coverage(self, str chrom, uint32_t start, uint32_t end):
        """
        Return coverage of a given region.
        """
        cdef:
            uint32_t length
            uint32_t covered

            long i

        assert self._closed

        assert end > start

        if chrom not in self._chroms:
            return 0.0

        assert end <= self._chroms[chrom], f"end position {end} larger than chromosome bound"

        length = end - start

        chroms = self.get_chromosomes()
        for i in range(len(chroms)):
            if chroms[i] == chrom:
                break
            else:
                start += self.get_chromosome_size(chroms[i])
                end += self.get_chromosome_size(chroms[i])

        overlaps = (self._starts.base < end) & (self._ends.base > start)
        starts = self._starts.base[overlaps]
        ends = self._ends.base[overlaps]

        starts = np.clip(starts, start, end)
        ends = np.clip(ends, start, end)
        covered = np.sum(ends - starts)

        return <float64_t> covered / <float64_t> length


    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep):
        """
        Filter peaks.
        Note: numpy bool array is uint8_t array!
        """
        assert self._closed

        assert self._starts.size == keep.size

        # filter peaks
        self._starts = self._starts.base[keep]
        self._ends = self._ends.base[keep]
        self._max_logp = self._max_logp.base[keep]
        self._max_fc = self._max_fc.base[keep]
        self._mean_logp = self._mean_logp.base[keep]
        self._mean_fc = self._mean_fc.base[keep]

        # re-count peaks fall into each chromosome
        chroms = self.get_chromosomes()
        rights = np.cumsum([
            self.get_chromosome_size(chrom) for chrom in chroms
        ])
        indices = np.searchsorted(self._starts, rights, side="right")
        counts = np.ediff1d(indices, to_begin=indices[0])
        self._counts = dict(zip(chroms, counts))


    cpdef void to_bed(self, object writable, str name_prefix):
        """
        Write peaks to a BED format file file.
        Columns: chrom, start, end, name, max_logp, strand, max_fc, mean_logp, mean_fc, peak
        """
        cdef:
            long i
            long count
            long last
            long shift
            long index
            long order
            str name

            uint32_t start
            uint32_t end
            float64_t max_logp
            float64_t max_fc
            float64_t mean_logp
            float64_t mean_fc

        writable.write("#chrom\tstart\tend\tname\tmax_logp\tstrand\tmax_fc\tmean_logp\tmean_fc\tpeak\n")

        last = 0
        shift = 0

        order = 1

        for chrom in self.get_chromosomes():
            count = self.get_chromosome_peaks(chrom)
            for i in range(count):
                index = last + i
                # extract data
                start = self._starts[index] - shift
                end = self._ends[index] - shift
                max_logp = self._max_logp[index]
                max_fc = self._max_fc[index]
                mean_logp = self._mean_logp[index]
                mean_fc = self._mean_fc[index]
                name = f"{name_prefix}_{order}"
                # write to file
                writable.write((
                    f"{chrom}\t{start}\t{end}\t{name}\t{max_logp:.6f}\t.\t"
                    f"{max_fc:.6f}\t{mean_logp:.6f}\t{mean_fc:.6f}\t{order}\n"
                ))
                # next peak
                order += 1
            # end one chromosome
            last += count
            shift += self.get_chromosome_size(chrom)
