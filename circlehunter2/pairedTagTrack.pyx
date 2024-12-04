# cython: language_level=3

from natsort import natsorted
import numpy as np

from numpy cimport uint32_t, uint8_t, int64_t, ndarray
from libc.stdint cimport UINT32_MAX
cimport cython

from circlehunter2.bedTrack cimport BedTrack
from circlehunter2.pairedPeakTrack cimport PairedPeakTrack


cdef class PairedTagTrack:
    """
    Class representing PE reads alignment.
    """

    def __init__(self, dict chrom_sizes, uint32_t buffer_size=256):
        assert buffer_size > 0

        self._buffer_size = buffer_size
        self._chroms = dict()
        self._shifts = dict()
        shift = 0
        chroms = natsorted(chrom_sizes.keys())
        for chrom in chroms:
            self._chroms[chrom] = chrom_sizes[chrom]
            self._shifts[chrom] = shift
            shift += chrom_sizes[chrom]

        self._data = {
            "lefts1": np.empty(self._buffer_size, dtype=np.uint32),
            "reverses1": np.empty(self._buffer_size, dtype=np.uint8),
            "lefts2": np.empty(self._buffer_size, dtype=np.uint32),
            "reverses2": np.empty(self._buffer_size, dtype=np.uint8)
        }

        self._pointer = 0
        self._size = self._buffer_size

        self._closed = False


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


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void add_tag_pair(self,
                            str chrom1, uint32_t start1, bint reverse1,
                            str chrom2, uint32_t start2, bint reverse2):
        """
        Add one tag pair to this track.
        """
        cdef:
            long new_size
            uint32_t shift1
            uint32_t shift2

        assert not self._closed, "track closed!"

        assert chrom1 in self._chroms
        assert start1 < self._chroms[chrom1]

        assert chrom2 in self._chroms
        assert start2 < self._chroms[chrom2]

        if self._pointer >= self._size:  # buffer size increase needed
            new_size = self._pointer + self._buffer_size
            self._data["lefts1"] = np.resize(self._data["lefts1"], new_size)
            self._data["reverses1"] = np.resize(self._data["reverses1"], new_size)
            self._data["lefts2"] = np.resize(self._data["lefts2"], new_size)
            self._data["reverses2"] = np.resize(self._data["reverses2"], new_size)
            self._size = new_size

        # find shift
        shift1 = self._shifts[chrom1]
        shift2 = self._shifts[chrom2]

        # add tag
        self._data["lefts1"][self._pointer] = start1 + shift1
        self._data["reverses1"][self._pointer] = reverse1
        self._data["lefts2"][self._pointer] = start2 + shift2
        self._data["reverses2"][self._pointer] = reverse2

        # next pointer
        self._pointer += 1


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void add_track(self, PairedTagTrack track):
        """
        Add all tags into this track from given track.
        """
        assert not self._closed, "track closed!"

        # check chromosome sizes
        chroms1 = self.get_chromosomes()
        chroms2 = track.get_chromosomes()

        assert len(chroms1) == len(chroms2), "chromosomes not match"

        for chrom1, chrom2 in zip(chroms1, chroms2):
            assert chrom1 == chrom2
            size1 = self.get_chromosome_size(chrom1)
            size2 = track.get_chromosome_size(chrom2)
            assert size1 == size2

        data = track.get_data()

        size = data.size

        if self._pointer + data.size >= self._size:  # size increase needed
            # extend size
            new_size = self._pointer + data.size
            self._data["lefts1"] = np.resize(self._data["lefts1"], new_size)
            self._data["reverses1"] = np.resize(self._data["reverses1"], new_size)
            self._data["lefts2"] = np.resize(self._data["lefts2"], new_size)
            self._data["reverses2"] = np.resize(self._data["reverses2"], new_size)
            self._size = new_size

        insert_slice = slice(self._pointer, self._pointer + data.size)
        self._data["lefts1"][insert_slice] = data["lefts1"]
        self._data["reverses1"][insert_slice] = data["reverses1"]
        self._data["lefts2"][insert_slice] = data["lefts2"]
        self._data["reverses2"][insert_slice] = data["reverses2"]

        self._pointer = self._pointer + data.size


    cpdef void close(self):
        """
        Close this track.
        Indicating that there would not be any tag added into this track.
        All buffer will be released.
        """
        cdef:
            long new_size

        new_size = self._pointer
        # release memory
        self._data["lefts1"] = np.resize(self._data["lefts1"], new_size)
        self._data["reverses1"] = np.resize(self._data["reverses1"], new_size)
        self._data["lefts2"] = np.resize(self._data["lefts2"], new_size)
        self._data["reverses2"] = np.resize(self._data["reverses2"], new_size)
        self._size = new_size

        self._closed = True


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef set_chromosome_sizes(self, dict chrom_sizes):
        """
        Set the chromosome size of this track.
        Note: Tags out of this chromosome list or chromosome regions is filter out.
        """
        cdef:
            long old_chroms_len

            str chrom

            uint32_t[:] left_limits
            uint32_t[:] right_limits

            dict new_shifts
            uint32_t new_size

            uint32_t new_left_limit1
            uint32_t new_right_limit1
            uint32_t new_left_limit2
            uint32_t new_right_limit2

            uint32_t new_shift1
            uint32_t new_shift2

            uint32_t[:] new_lefts1
            uint8_t[:] new_reverse1
            uint32_t[:] new_lefts2
            uint8_t[:] new_reverse2

            long i
            long j
            long k

        assert self._closed, "Track not closed!"

        # calculate limits
        old_chroms = self.get_chromosomes()
        old_chroms_len = len(old_chroms)
        left_limits = np.empty(old_chroms_len, dtype=np.uint32)
        right_limits = np.empty(old_chroms_len, dtype=np.uint32)
        for i in range(old_chroms_len):
            left_limits[i] = self._shifts[old_chroms[i]]
            right_limits[i] = left_limits[i] + self._chroms[old_chroms[i]]

        # calculate new shift
        new_chroms = tuple(natsorted(chrom_sizes.keys()))
        new_shifts = dict()
        new_shift = 0
        for chrom in new_chroms:
            new_shifts[chrom] = new_shift
            new_shift += chrom_sizes[chrom]

        new_lefts1 = np.empty(self._size, dtype=np.uint32)
        new_reverse1 = np.empty(self._size, dtype=np.uint8)
        new_lefts2 = np.empty(self._size, dtype=np.uint32)
        new_reverse2 = np.empty(self._size, dtype=np.uint8)

        new_size = 0
        for i in range(self._size):
            # find old chrom1 and filter
            j = 0
            for j in range(old_chroms_len):
                if left_limits[j] <= self._data["lefts1"][i] < right_limits[j]:
                    break
            new_chrom1 = old_chroms[j]
            if new_chrom1 not in new_shifts:
                continue
            # find old chrom2 and filter
            k = 0
            for k in range(old_chroms_len):
                if left_limits[k] <= self._data["lefts2"][i] < right_limits[k]:
                    break
            new_chrom2 = old_chroms[k]
            if new_chrom2 not in new_shifts:
                continue

            # check if left1 out of chromosome bound
            new_left_limit1 = new_shifts[new_chrom1]
            new_right_limit1 = chrom_sizes[new_chrom1]
            new_right_limit1 += new_left_limit1
            if self._data["lefts1"][i] < new_left_limit1 or self._data["lefts1"][i] >= new_right_limit1:
                continue

            # check if left2 out of chromosome bound
            new_left_limit2 = new_shifts[new_chrom2]
            new_right_limit2 = chrom_sizes[new_chrom2]
            new_right_limit2 += new_left_limit2
            if self._data["lefts2"][i] < new_left_limit2 or self._data["lefts2"][i] >= new_right_limit2:
                continue

            # shift to new position
            new_shift1 = new_shifts[new_chrom1]
            new_shift2 = new_shifts[new_chrom2]

            new_lefts1[new_size] = self._data["lefts1"][i] - left_limits[j] + new_shift1
            new_reverse1[new_size] = self._data["reverses1"][i]
            new_lefts2[new_size] = self._data["lefts2"][i] - left_limits[k] + new_shift2
            new_reverse2[new_size] = self._data["reverses2"][i]

            new_size += 1

        # release unused memory
        self._data["lefts1"] = np.resize(new_lefts1, new_size)
        self._data["reverses1"] = np.resize(new_reverse1, new_size)
        self._data["lefts2"] = np.resize(new_lefts2, new_size)
        self._data["reverses2"] = np.resize(new_reverse2, new_size)
        self._size = new_size

        self._chroms = chrom_sizes
        self._shifts = new_shifts


    @property
    def total(self):
        """
        Return tag count of this track.
        """
        assert self._closed

        return self._size


    def get_data(self):
        """
        Return a 2D array containing all data of tags.
        """
        assert self._closed

        return np.rec.fromarrays([
            self._data["lefts1"], self._data["reverses1"],
            self._data["lefts2"], self._data["reverses2"]
        ], names=["lefts1", "reverses1", "lefts2", "reverses2"])


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef overlap(self, BedTrack track):
        """
        Check if tags overlap with the regions stored in the given bed track.
        """
        cdef:
            uint32_t[:] lefts
            uint32_t[:] rights
            int64_t[:] inserts
            long i

            uint8_t[:] overlap1
            uint8_t[:] overlap2

        assert self._closed

        # check chromosome sizes
        chroms1 = self.get_chromosomes()
        chroms2 = track.get_chromosomes()

        assert len(chroms1) == len(chroms2), "chromosomes not match"

        for chrom1, chrom2 in zip(chroms1, chroms2):
            assert chrom1 == chrom2
            size1 = self.get_chromosome_size(chrom1)
            size2 = track.get_chromosome_size(chrom2)
            assert size1 == size2

        data = track.get_data()

        # binary search left positions
        inserts = np.searchsorted(data["lefts"], self._data["lefts1"], side="right") - 1
        lefts = self._data["lefts1"]
        rights = data["rights"]
        overlap1 = np.empty(self._size, dtype=np.uint8)
        # loop through each tag if the tag fall in the region
        for i in range(self._size):
            if inserts[i] >= 0 and rights[inserts[i]] > lefts[i]:
                overlap1[i] = 1
            else:
                overlap1[i] = 0

        # binary search left positions
        inserts = np.searchsorted(data["lefts"], self._data["lefts2"], side="right") - 1
        lefts = self._data["lefts2"]
        rights = data["rights"]
        overlap2 = np.empty(self._size, dtype=np.uint8)
        # loop through each tag if the tag fall in the region
        for i in range(self._size):
            if inserts[i] >= 0 and rights[inserts[i]] > lefts[i]:
                overlap2[i] = 1
            else:
                overlap2[i] = 0

        return np.rec.fromarrays([overlap1.base, overlap2.base], names=["overlaps1", "overlaps2"])


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void filter(self, ndarray[uint8_t, ndim=1] keep):
        cdef:
            uint32_t new_size
            uint32_t[:] new_lefts1
            uint8_t[:] new_reverses1
            uint32_t[:] new_lefts2
            uint8_t[:] new_reverses2

            uint32_t[:] lefts1
            uint8_t[:] reverses1
            uint32_t[:] lefts2
            uint8_t[:] reverses2

            uint8_t[:] keep_view

            long i
            long j


        assert self._closed
        assert keep.size == self._size

        # new data
        new_size = np.sum(keep)
        new_lefts1 = np.empty(new_size, dtype=np.uint32)
        new_reverses1 = np.empty(new_size, dtype=np.uint8)
        new_lefts2 = np.empty(new_size, dtype=np.uint32)
        new_reverses2 = np.empty(new_size, dtype=np.uint8)

        # old_data
        lefts1 = self._data["lefts1"]
        reverses1 = self._data["reverses1"]
        lefts2 = self._data["lefts2"]
        reverses2 = self._data["reverses2"]
        keep_view = keep
        # move data to new_data
        i = j = 0
        for i in range(self._size):
            if keep_view[i] > 0:
                new_lefts1[j] = lefts1[i]
                new_reverses1[j] = reverses1[i]
                new_lefts2[j] = lefts2[i]
                new_reverses2[j] = reverses2[i]
                j += 1
        # assign data
        self._data["lefts1"] = new_lefts1.base
        self._data["reverses1"] = new_reverses1.base
        self._data["lefts2"] = new_lefts2.base
        self._data["reverses2"] = new_reverses2.base
        self._size = new_size


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef uint32_t[:] _merge_into_regions(self, uint32_t distance):
        """
        Merge all tags within distance into one region.
        Return an array indicate the region id each tag belong to.
        """
        cdef:
            uint32_t[:] right_limits

            uint32_t[:] lefts1
            uint8_t[:] reverses1
            uint32_t[:] lefts2
            uint8_t[:] reverses2

            uint32_t distance1
            uint32_t distance2

            uint32_t[:] regions
            uint32_t region

            uint32_t[:] queue
            long pointer

            long i
            long j
            long k

            long m
            long n

        # calculate limits
        chroms = self.get_chromosomes()
        chroms_len = len(chroms)
        right_limits = np.empty(chroms_len, dtype=np.uint32)
        for i in range(chroms_len):
            right_limits[i] = self._shifts[chroms[i]] + self._chroms[chroms[i]]

        lefts1 = self._data["lefts1"]
        reverses1 = self._data["reverses1"]
        lefts2 = self._data["lefts2"]
        reverses2 = self._data["reverses2"]

        regions = np.zeros(self._size, dtype=np.uint32)
        queue = np.empty(self._size, dtype=np.uint32)

        region = 1
        i = j = k = n = m = 0

        for i in range(self._size):  # loop through tags
            if regions[i] > 0:  # this tag has determined region, skip
                continue

            # a new region
            regions[i] = region

            # find all tags distance is less than given distance by DFS
            # these tags all belong to this region

            # push this tag into queue
            pointer = 0
            queue[pointer] = i

            while pointer >= 0:  # check all tags in the queue
                # pop from queue
                i = queue[pointer]
                pointer -= 1
                for j in range(i + 1, self._size):
                    if regions[j] > 0:  # this tag has determined region
                        continue
                    # tags orientation
                    if reverses1[i] != reverses1[j]:  # first tags have different orientation
                        continue
                    if reverses2[i] != reverses2[j]:  # next tags have different orientation
                        continue
                    # tags horizontal distance
                    if lefts1[i] < lefts1[j]:
                        distance1 = lefts1[j] - lefts1[i]
                    else:
                        distance1 = lefts1[i] - lefts1[j]
                    if distance1 > distance:  # tags horizontal distance is larger than distance
                        continue
                    # fist tags not in the same chromosome
                    for m in range(chroms_len):
                        if lefts1[i] < right_limits[m]:
                            break
                    for n in range(chroms_len):
                        if lefts1[j] < right_limits[n]:
                            break
                    if m != n:
                        continue
                    # tags vertical distance
                    if lefts2[i] < lefts2[j]:
                        distance2 = lefts2[j] - lefts2[i]
                    else:
                        distance2 = lefts2[i] - lefts2[j]
                    if distance2 > distance:  # tags vertical distance is larger than distance
                        continue
                    # next tags not in the same chromosome
                    for m in range(chroms_len):
                        if lefts2[i] < right_limits[m]:
                            break
                    for n in range(chroms_len):
                        if lefts2[j] < right_limits[n]:
                            break
                    if m != n:
                        continue
                    # tags pass all filter
                    # this tag belong to this region, too. push into the queue
                    pointer += 1
                    queue[pointer] = j
                    regions[j] = region
            # next region
            region += 1

        return regions


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef PairedPeakTrack call_peaks(self, uint32_t cutoff, uint32_t distance):
        """
        Call peaks.
        """
        cdef:
            uint32_t[:] regions
            uint32_t n

            uint32_t[:] left_limits
            uint32_t[:] right_limits

            uint32_t[:] lefts1
            uint8_t[:] reverses1
            uint32_t[:] lefts2
            uint8_t[:] reverses2

            uint32_t[:] peaks_left1
            uint32_t[:] peaks_right1
            uint8_t[:] peaks_reverse1
            uint32_t[:] peaks_left2
            uint32_t[:] peaks_right2
            uint8_t[:] peaks_reverse2
            uint32_t[:] peaks_tag_count

            uint32_t left
            uint32_t right

            str chrom1
            uint32_t start1
            uint32_t end1
            uint8_t reverse1
            str chrom2
            uint32_t start2
            uint32_t end2
            uint8_t reverse2
            uint32_t count

            long i
            long j

        assert self._closed
        
        # calculate limits
        chroms = self.get_chromosomes()
        chroms_len = len(chroms)
        left_limits = np.empty(chroms_len, dtype=np.uint32)
        right_limits = np.empty(chroms_len, dtype=np.uint32)
        for i in range(chroms_len):
            left_limits[i] = self._shifts[chroms[i]]
            right_limits[i] = left_limits[i] + self._chroms[chroms[i]]

        # merge tags into regions
        regions = self._merge_into_regions(distance)

        lefts1 = self._data["lefts1"]
        reverses1 = self._data["reverses1"]
        lefts2 = self._data["lefts2"]
        reverses2 = self._data["reverses2"]

        n = np.max(regions)
        peaks_left1 = np.full(n, UINT32_MAX, dtype=np.uint32)
        peaks_right1 = np.full(n, 0, dtype=np.uint32)
        peaks_reverse1 = np.empty(n, dtype=np.uint8)
        peaks_left2 = np.full(n, UINT32_MAX, dtype=np.uint32)
        peaks_right2 = np.full(n, 0, dtype=np.uint32)
        peaks_reverse2 = np.empty(n, dtype=np.uint8)
        peaks_tag_count = np.full(n, 0, dtype=np.uint32)

        i = j = 0
        for i in range(self._size):
            n = regions[i] - 1

            # first region
            left = lefts1[i]
            # find chromosome
            for j in range(chroms_len):
                if left < right_limits[j]:
                    break
            # extend a point to a region, but not out of chromosome bound
            if left >= right_limits[j] - distance:
                right = right_limits[j]
            else:
                right = left + distance + 1
            if left <= left_limits[j] + distance:
                left = left_limits[j]
            else:
                left -= distance
            # save to result
            if left < peaks_left1[n]:
                peaks_left1[n] = left
            if right > peaks_right1[n]:
                peaks_right1[n] = right
            peaks_reverse1[n] = reverses1[i]

            # next region
            left = lefts2[i]
            # find chromosome
            for j in range(chroms_len):
                if left < right_limits[j]:
                    break
            # extend a point to a region, but not out of chromosome bound
            if left >= right_limits[j] - distance:
                right = right_limits[j]
            else:
                right = left + distance + 1
            if left <= left_limits[j] + distance:
                left = left_limits[j]
            else:
                left -= distance
            # save to result
            if left < peaks_left2[n]:
                peaks_left2[n] = left
            if right > peaks_right2[n]:
                peaks_right2[n] = right
            peaks_reverse2[n] = reverses2[i]

            # count tags
            peaks_tag_count[n] += 1

        peaks = PairedPeakTrack(self._chroms, buffer_size=self._buffer_size)

        for i in range(peaks_tag_count.size):
            count = peaks_tag_count[i]
            # filter peaks
            if count < cutoff:
                continue

            # find chromosome and extract data for tag1
            for j in range(chroms_len):
                if peaks_left1[i] < right_limits[j]:
                    break
            chrom1 = chroms[j]
            start1 = peaks_left1[i] - left_limits[j]
            end1 = peaks_right1[i] - left_limits[j]
            reverse1 = peaks_reverse1[i]

            # find chromosome and extract data for tag2
            for j in range(chroms_len):
                if peaks_left2[i] < right_limits[j]:
                    break
            chrom2 = chroms[j]
            start2 = peaks_left2[i] - left_limits[j]
            end2 = peaks_right2[i] - left_limits[j]
            reverse2 = peaks_reverse2[i]

            peaks.add_peak_pair(chrom1, start1, end1, reverse1, chrom2, start2, end2, reverse2, count)

        peaks.close()

        return peaks


    @cython.nonecheck(False)
    @cython.initializedcheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cpdef void to_bedpe(self, object writable, uint32_t extend=50):
        """
        Write tags to a BEDPE format file.
        Columns: chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2
        """
        cdef:
            uint32_t[:] left_limits
            uint32_t[:] right_limits

            uint32_t[:] lefts1
            uint8_t[:] reverses1
            uint32_t[:] lefts2
            uint8_t[:] reverses2

            str chrom1
            uint32_t start1
            uint32_t end1
            str strand1

            str chrom2
            uint32_t start2
            uint32_t end2
            str strand2

            long order

            long i
            long j
            long k

        assert self._closed

        # find chromosome breakpoints
        chroms = self.get_chromosomes()
        chroms_len = len(chroms)
        left_limits = np.empty(chroms_len, dtype=np.uint32)
        right_limits = np.empty(chroms_len, dtype=np.uint32)
        for i in range(chroms_len):
            left_limits[i] = self._shifts[chroms[i]]
            right_limits[i] = left_limits[i] + self._chroms[chroms[i]]

        lefts1 = self._data["lefts1"]
        reverses1 = self._data["reverses1"]
        lefts2 = self._data["lefts2"]
        reverses2 = self._data["reverses2"]

        writable.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n")

        order = 1
        # loop through tags
        for i in range(self._size):
            # find chromosomes
            j = 0
            for j in range(chroms_len):
                if lefts1[i] < right_limits[j]:
                    break
            k = 0
            for k in range(chroms_len):
                if lefts2[i] < right_limits[k]:
                    break

            # extract data
            chrom1 = chroms[j]
            start1 = lefts1[i] - left_limits[j]
            end1 = start1 + extend
            strand1 = "-" if reverses1[i] else "+"

            chrom2 = chroms[k]
            start2 = lefts2[i] - left_limits[k]
            end2 = start2 + extend
            strand2 = "-" if reverses2[i] else "+"

            name = f"paired_tag_{order}"

            writable.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t.\t{strand1}\t{strand2}\n")

            order += 1
            # end tags
