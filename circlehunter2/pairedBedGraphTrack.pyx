# cython: language_level=3

from collections import OrderedDict

import numpy as np

from numpy cimport uint32_t, float64_t, ndarray
from cpython cimport bool
cimport cython

from circlehunter2.bedGraphTrack cimport BedGraphTrack
from circlehunter2.peakTrack cimport PeakTrack
from circlehunter2.utils cimport _logp as cal_logp


LOG_005 = -np.log(0.05)
LOG_010 = -np.log(0.10)


cdef class PairedBedGraphTrack:
    """
    Class representing two BedGraphTrack but has same positions.
    """

    cdef:
        object _chroms
        uint32_t[:] _positions
        uint32_t[:] _lengths
        float64_t[:] _values1
        float64_t[:] _values2
        dict _counts
        bool _compared
        float64_t[:] _logp
        float64_t[:] _fc


    def __init__(self, BedGraphTrack track1, BedGraphTrack track2=None):
        cdef:
            tuple chroms1
            str chrom1
            long size1
            tuple chroms2
            str chrom2
            long size2

        self._chroms = OrderedDict()
        self._counts = dict()
        self._compared = False

        if track2 is None:  # track2 is mean of track1
            track2 = track1.new_value(track1.mean())
            chroms1 = track1.get_chromosomes()
            for chrom1 in chroms1:
                self._chroms[chrom1] = track1.get_chromosome_size(chrom1)
        else:  # check chromosome sizes
            chroms1 = track1.get_chromosomes()
            chroms2 = track2.get_chromosomes()

            assert len(chroms1) == len(chroms2), "chromosomes not match"

            for chrom1, chrom2 in zip(chroms1, chroms2):

                assert chrom1 == chrom2

                size1 = track1.get_chromosome_size(chrom1)
                size2 = track2.get_chromosome_size(chrom2)

                assert size1 == size2

                # assign chromosome size
                self._chroms[chrom1] = size1

        self._pair_tracks(track1, track2)


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef void _pair_tracks(self, BedGraphTrack track1, BedGraphTrack track2):
        """
        Paired two BedGraphTrack, make them have the same positions.
        """
        cdef:
            uint32_t[:] positions1
            float64_t[:] values1

            uint32_t[:] positions2
            float64_t[:] values2

            uint32_t last_position

            long count
            long count1
            long count2

            long i1
            long i2

            long last_count
            long i
            long left
            long right

        # extract data
        positions1 = track1.get_positions()
        values1 = track1.get_values()
        positions2 = track2.get_positions()
        values2 = track2.get_values()

        # init var
        count = 0
        count1 = len(positions1)
        count2 = len(positions2)

        chroms = self.get_chromosomes()
        i = 0
        left = 0
        right = self.get_chromosome_size(chroms[i])
        last_count = 0

        # max possible size
        self._positions = np.empty(count1 + count2, dtype=np.uint32)
        self._lengths = np.empty(count1 + count2, dtype=np.uint32)
        self._values1 = np.empty(count1 + count2, dtype=np.float64)
        self._values2 = np.empty(count1 + count2, dtype=np.float64)

        last_position = 0

        i1 = 0
        i2 = 0

        # loop through positions
        while i1 < count1 and i2 < count2:
            if positions1[i1] < positions2[i2]:  # track one end
                self._positions[count] = positions1[i1]
                self._lengths[count] = positions1[i1] - last_position
                self._values1[count] = values1[i1]
                self._values2[count] = values2[i2]
                last_position = positions1[i1]
                count += 1
                i1 += 1
            elif positions1[i1] > positions2[i2]:  # track two end
                self._positions[count] = positions2[i2]
                self._lengths[count] = positions2[i2] - last_position
                self._values1[count] = values1[i1]
                self._values2[count] = values2[i2]
                last_position = positions2[i2]
                count += 1
                i2 += 1
            else:  # same end
                self._positions[count] = positions1[i1]
                self._lengths[count] = positions1[i1] - last_position
                self._values1[count] = values1[i1]
                self._values2[count] = values2[i2]
                count += 1
                # end of one chromosome, remember region counts
                if positions1[i1] == right:
                    # mark end of a chromosome
                    self._counts[chroms[i]] = count - last_count
                    left = right
                    last_count = count
                    # move to next chromosome
                    i += 1
                    if i < len(chroms):
                        right += self.get_chromosome_size(chroms[i])
                # end if
                last_position = positions1[i1]
                # next position
                i1 += 1
                i2 += 1

        # resize to save memory
        self._positions = np.resize(self._positions, count)
        self._lengths = np.resize(self._lengths, count)
        self._values1 = np.resize(self._values1, count)
        self._values2 = np.resize(self._values2, count)


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
        if chrom not in self._chroms:
            raise KeyError(f"No such chromosome: {chrom}")

        return self._counts[chrom]


    def get_positions(self):
        """
        Return the positions of this track
        """
        return self._positions.base


    def get_lengths(self):
        """
        Return the length of each region in this track
        """
        return self._lengths.base


    def get_values(self):
        """
        Return the values of this track
        """
        return np.vstack([self._values1, self._values2])

    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef void _compare(self):
        """
        Compare paired two BedGraphTrack, calculate logp and foldchange.
        """
        cdef:
            float64_t[:] logp_values
            float64_t[:] fc_values
            float64_t value1
            float64_t value2

        logp_values = np.empty(len(self._positions), dtype=np.float64)
        fc_values = np.empty(len(self._positions), dtype=np.float64)

        # loop through all regions, may be caching? or parallelism?
        for i in range(len(self._positions)):
            value1 = self._values1[i]
            value2 = self._values2[i]
            logp_values[i] = cal_logp(value1, value2)
            fc_values[i] = value1 / value2

        self._logp = logp_values
        self._fc = fc_values
        self._compared = True


    cpdef BedGraphTrack calculate_logp(self):
        """
        Calculate logp values between track1 and track2, return as a new BedGraphTrack.
        Something like -poisson.logpmf(x, y) in scipy.
        """
        cdef:
            BedGraphTrack track
            ndarray[float64_t, ndim=1] logp_values
            long i
            tuple chroms
            str chrom
            uint32_t count
            uint32_t left
            uint32_t right
            uint32_t shift
            ndarray[uint32_t, ndim = 1] positions
            ndarray[float64_t, ndim = 1] values

        if not self._compared:
            self._compare()

        logp_values = self._logp.base

        track = BedGraphTrack(chrom_sizes=dict(self._chroms))

        chroms = self.get_chromosomes()
        left = 0
        shift = 0

        # add logp to BedGraphTrack
        for i in range(len(chroms)):
            # get regions
            chrom = chroms[i]
            count = self.get_chromosome_regions(chrom)
            right = left + count
            # extract regions positions and values
            positions = self._positions.base[left:right] - shift
            values = logp_values[left:right]
            # assign values
            track.add_chromosome(chrom, positions, values)
            # shift positions
            shift += self.get_chromosome_size(chrom)
            left = right

        track.close()

        return track


    cpdef BedGraphTrack calculate_foldchange(self):
        """
        Calculate foldchange values between track1 and track2, return as a new BedGraphTrack.
        Something like x / y.
        """
        cdef:
            BedGraphTrack track
            ndarray[float64_t, ndim=1] fc_values
            long i
            tuple chroms
            str chrom
            uint32_t count
            uint32_t left
            uint32_t right
            uint32_t shift
            ndarray[uint32_t, ndim = 1] positions
            ndarray[float64_t, ndim = 1] values

        if not self._compared:
            self._compare()

        fc_values = self._fc.base

        track = BedGraphTrack(chrom_sizes=dict(self._chroms))

        chroms = self.get_chromosomes()
        left = 0
        shift = 0

        # add logp to BedGraphTrack
        for i in range(len(chroms)):
            # get regions
            chrom = chroms[i]
            count = self.get_chromosome_regions(chrom)
            right = left + count
            # extract regions positions and values
            positions = self._positions.base[left:right] - shift
            values = fc_values[left:right]
            # assign values
            track.add_chromosome(chrom, positions, values)
            # shift positions
            shift += self.get_chromosome_size(chrom)
            left = right

        track.close()

        return track


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef PeakTrack call_peaks(self, float64_t cutoff=LOG_005,  uint32_t max_gap=50, uint32_t min_length=200):
        """
        Call peaks with given cutoff.
        Peaks gap smaller than max_gap will be merged.
        Peak length smaller than min_length will be filter out.
        """
        cdef:
            PeakTrack track

            long i
            uint32_t left
            uint32_t size
            uint32_t right

            long[:] peak_lefts
            long[:] peak_sizes

            str chrom
            uint32_t start
            uint32_t end
            float64_t max_logp
            float64_t max_fc
            float64_t mean_logp
            float64_t mean_fc

            uint32_t length

            uint32_t last
            uint32_t last_size
            uint32_t shift

        if not self._compared:
            self._compare()

        # regions above cutoff
        indices = np.where(self._logp.base >= cutoff)[0]

        # extract peak regions
        ends = self._positions.base[indices]
        lengths = self._lengths.base[indices]
        starts = ends - lengths

        # gaps between peaks
        gaps = np.empty_like(lengths)
        gaps[1:lengths.size] = starts[1:lengths.size] - ends[0:lengths.size - 1]
        gaps[0] = 0

        # extend gaps between chromosomes
        chroms = self.get_chromosomes()
        rights = np.empty(len(chroms), dtype=np.uint32)
        for i in range(rights.size):
            rights[i] = self.get_chromosome_size(chroms[i])
        # end of each chromosome
        rights = np.cumsum(rights)
        # find gaps span two chromosomes, make them larger than given gap cutoff
        breakpoints = np.searchsorted(ends, rights, side="right")
        breakpoints = breakpoints[breakpoints < ends.size]
        gaps[breakpoints] += max_gap + 1

        # merge peaks with gap smaller than gap cutoff
        peak_groups = np.cumsum(gaps > max_gap)  # peaks with same group should be merged
        peak_groups, peak_lefts, peak_sizes = np.unique(
            peak_groups, return_index=True, return_counts=True
        )

        track = PeakTrack(chrom_sizes=dict(self._chroms))

        last = 0
        chrom = chroms[last]
        last_size = self.get_chromosome_size(chrom)
        shift = 0
        # loop through peaks
        for i in range(peak_groups.size):
            # merge peaks
            left = peak_lefts[i]
            size = peak_sizes[i]
            right = left + size

            start = starts[left]
            end = ends[right - 1]

            # filter by peak length
            length = end - start
            if length < min_length:
                continue

            # determine chromosome and shift
            while end > shift + last_size:
                shift += last_size
                last += 1
                chrom = chroms[last]
                last_size = self.get_chromosome_size(chrom)

            # shift position
            start -= shift
            end -= shift

            # extract data
            lengths = self._lengths.base[indices[left:right]]
            values1 = self._values1.base[indices[left:right]]
            values2 = self._values2.base[indices[left:right]]
            logp = self._logp.base[indices[left:right]]
            fc = self._fc.base[indices[left:right]]

            # max stats
            max_logp = np.max(logp)
            max_fc = np.max(fc)

            # mean stats
            mean_value1 = np.sum(values1 * lengths) / length
            mean_value2 = np.sum(values2 * lengths) / length
            mean_logp = cal_logp(mean_value1, mean_value2)
            mean_fc = mean_value1 / mean_value2

            # add peak
            track.add_peak(chrom, start, end, max_logp, max_fc, mean_logp, mean_fc)

        track.close()

        return track
