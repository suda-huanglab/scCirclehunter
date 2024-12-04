# cython: language_level=3

from collections import OrderedDict

from numpy cimport uint32_t, float64_t, ndarray

from natsort import natsorted
import numpy as np


cdef class BedGraphTrack:
    """
    Class representing values associated with a genome region like bedGraph file.
    """

    def __init__(self, dict chrom_sizes, float64_t default=0.0):
        self._chroms = OrderedDict()
        chroms = natsorted(chrom_sizes.keys())
        for chrom in chroms:
            self._chroms[chrom] = chrom_sizes[chrom]
        self._tmp_positions = dict()
        self._tmp_values = dict()
        self._counts = dict()
        self._default = default
        self._closed = False


    cpdef add_chromosome(self,
                         str chrom,
                         ndarray[uint32_t, ndim=1] positions,
                         ndarray[float64_t, ndim=1] values):
        """
        Add a chromosome into this track.
        Positions is the continued end position of a region.
        Values should ordered as the position.
        Position can not be large than the give chromosome size!
        Exist chromosome will be replaced!
        """
        cdef:
            long i
            ndarray[long, ndim=1] indices

        # track not closed
        assert not self._closed, "track closed!"
        # chromosome size exists
        assert chrom in self._chroms, f"no such chromosome: {chrom}"
        # positions and values have same shape
        assert positions.shape[0] == values.shape[0], f"positions not the same shape with values"
        # positions is in continued order
        assert np.all(np.diff(positions) > 0), "positions is not continued"
        # positions match chromosome size
        assert positions[positions.shape[0] - 1] <= self._chroms[chrom], "positions not match the chromosome size"

        # last region of the chromosome is missing
        if positions[positions.shape[0] - 1] < self._chroms[chrom]:
            # extend the data array
            positions = np.resize(positions, positions.shape[0] + 1)
            values = np.resize(values, values.shape[0] + 1)
            # assign as default value
            positions[positions.shape[0] - 1] = self._chroms[chrom]
            values[positions.shape[0] - 1] = self._default

        # remove duplicates
        indices = np.where(np.diff(values) == 0)[0]
        positions = np.delete(positions, indices)
        values = np.delete(values, indices)

        # assign the values
        self._tmp_positions[chrom] = positions
        self._tmp_values[chrom] = values


    cpdef close(self):
        """
        Close this track.
        Indicating that there would not be any data added into this track.
        All buffer will be released.
        """
        cdef:
            long i
            uint32_t count
            uint32_t left
            uint32_t right
            uint32_t last
            uint32_t[:] positions
            float64_t[:] values

        assert not self._closed, "track closed!"

        chroms = self.get_chromosomes()

        # count chromosome regions
        count = 0
        for i in range(len(chroms)):
            if chroms[i] in self._tmp_positions:
                self._counts[chroms[i]] = self._tmp_positions[chroms[i]].shape[0]
                count += self._tmp_positions[chroms[i]].shape[0]
            else:
                self._counts[chroms[i]] = 1
                count += 1

        # init data array with region counts
        self._positions = np.zeros(count, dtype=np.uint32)
        self._values = np.zeros(count, dtype=np.float64)

        # assign values
        last = 0
        left = 0
        for i in range(len(chroms)):
            # extrac data from tmp
            if chroms[i] in self._tmp_positions:
                positions = self._tmp_positions[chroms[i]] + last
                values = self._tmp_values[chroms[i]]
            else:
                positions = np.repeat(last + self._chroms[chroms[i]], 1).astype(np.uint32)
                values = np.repeat(self._default, 1).astype(np.float64)

            # assign values
            right = left + self._counts[chroms[i]]
            self._positions[left:right] = positions
            self._values[left:right] = values

            last += self._chroms[chroms[i]]
            left += self._counts[chroms[i]]

        # release memory
        self._tmp_positions = None
        self._tmp_values = None

        # mark as closed
        self._closed = True


    def get_chromosomes(self):
        """
        Chromosomes stored in this track.
        """
        return tuple(self._chroms.keys())


    def get_chromosome_size(self, str chrom):
        """
        Return the size of the given chromosome
        """
        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")
        return self._chroms[chrom]


    def get_chromosome_regions(self, str chrom):
        """
        Return the regions count of the given chromosome
        """
        assert self._closed

        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")

        return self._counts[chrom]


    def get_positions(self):
        """
        Return the positions of this track
        """
        assert self._closed

        return self._positions.base


    def get_values(self):
        """
        Return the values of this track
        """
        assert self._closed

        return self._values.base


    cpdef BedGraphTrack new_values(self, ndarray[float64_t, ndim=1] values):
        """
        Return a new BedGraphTrack instance with values replaced
        """
        cdef:
            BedGraphTrack track
            str chrom
            uint32_t count
            uint32_t left
            uint32_t right
            uint32_t shift
            ndarray[uint32_t, ndim=1] positions
            ndarray[float64_t, ndim=1] new_values

        assert self._closed

        assert self._values.size == values.size, "values not the same shape"

        track = BedGraphTrack(chrom_sizes=dict(self._chroms), default=self._default)

        left = 0
        shift = 0
        # loop through chromosomes
        chroms = self.get_chromosomes()
        for chrom in chroms:
            # get regions
            count = self.get_chromosome_regions(chrom)
            right = left + count
            # extract regions positions and values
            positions = self._positions.base[left:right] - shift
            new_values = values[left:right]
            # assign values
            track.add_chromosome(chrom, positions, new_values)
            # shift positions
            shift += self.get_chromosome_size(chrom)
            left = right

        track.close()

        return track


    cpdef coverage(self, str chrom, uint32_t start, uint32_t end, uint32_t cutoff):
        cdef:
            long i
            long left
            long right

            uint32_t shift

        assert chrom in self._chroms, f"No such chromosome: {chrom}"
        assert end <= self._chroms[chrom], f"Positions {end} out of chromosome bound"
        assert start < end, f"region end {end} less than {start}"

        # find the chromosome and shift
        chroms = self.get_chromosomes()
        shift = 0
        for i in range(len(chroms)):
            if chroms[i] == chrom:
                break
            else:
                shift += self.get_chromosome_size(chroms[i])

        # shfit region
        start += shift
        end += shift

        # find regions that overlap
        left = 0
        right = 0

        i = 0
        while i < self._positions.size:
            if self._positions[i] > start:
                left = i
                break
            i += 1
        while i < self._positions.size:
            if self._positions[i] > end:
                right = i
                break
            i += 1

        # extract regions that overlap
        positions = np.empty(right - left + 1, dtype=np.uint32)
        values = np.empty(right - left + 1, dtype=np.float64)

        for i in range(left, right + 1):
            positions[i-left] = self._positions[i]
            values[i-left] = self._values[i]
        positions[right - left] = end

        # calculate region length
        lengths = np.ediff1d(positions, to_begin=positions[0] - start)

        # calculate depth mean
        mean = np.sum(values * lengths) / np.sum(lengths)

        # calculater high coverage fraction
        index = values >= cutoff
        coverage = np.sum(lengths[index]) / np.sum(lengths)

        return mean, coverage


    cpdef BedGraphTrack new_value(self, float64_t value):
        """
        Return a new BedGraphTrack with all region as the given single value.
        """
        track = BedGraphTrack(chrom_sizes=dict(self._chroms), default=value)
        track.close()

        return track


    cpdef float64_t mean(self):
        """
        Calculate the mean of the this track.
        May use as logp mu.
        """
        cdef:
            uint32_t last
            uint32_t current
            float64_t length
            float64_t total_length
            float64_t value
            float64_t total_value
            long i

        last = 0
        total_length = 0
        total_value = 0
        for i in range(len(self._positions)):
            length = <float64_t> (self._positions[i] - last)
            total_length += length
            value = self._values[i] * length
            total_value += value
            last = self._positions[i]
        return total_value / total_length
