# cython: language_level=3

import numpy as np
from libc.stdlib cimport atoi
from libc.math cimport log, lgamma
from numpy cimport uint32_t, float64_t, ndarray
cimport cython


cpdef tuple nasort_key(object obj):
    """
    Use this key func to sort chromosome names in nature order.
    Such as: chr1 -> chr2 -> chr10 -> chr11
    """
    cdef:
        int last
        int current
        int end
        bytes value

    if isinstance(obj, bytes):
        value = obj
    elif isinstance(obj, str):
        value = obj.encode("utf-8")
    else:
        return (obj,)


    keys = list()

    end = 0
    last = 0
    # loop through each char
    for end in range(len(value)):
        current = value[end]
        if 48 <= current <= 57:  # current char is digital
            continue
        else:  # not a digital char
            if end - last > 1:  # last chars is digital
                last = atoi(value[last:end + 1])  # digital chars merge as number
                keys.append(last)
            # add this char as integer by ASCII
            keys.append(current)
            last = end
    else:  # handle the last digital chars
        if end - last > 1:
            last = atoi(value[last:end + 1])
            keys.append(last)

    return tuple(keys)


cpdef tuple make_windows(dict chrom_sizes, long size):
        """
        Return the windows of chromosomes by the given windows size
        """
        cdef:
            list windows
            long chrom_size
            object chrom
            long last
            long current

        windows = list()

        chroms = sorted(chrom_sizes.keys(), key=nasort_key)
        # loop through chromosomes
        for chrom in iter(chroms):
            chrom_size = chrom_sizes[chrom]
            # trim the size to avoid judge in the loop
            chrom_size -= size
            last = 0
            current = 0
            while current < chrom_size:
                current = last + size
                windows.append((chrom, last, current))  # no need to judge
                last = current
            # add the trim out end
            windows.append((chrom, last, chrom_size + size))

        return tuple(windows)


cdef tuple _pileup(ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends):
    """
    Pileup the given tags and return the end positions and depth used by BedGraphTrack
    """
    cdef:
        long size_start
        long size_end
        uint32_t[:] start_positions
        uint32_t[:] end_positions
        long i_start
        long i_end
        uint32_t start
        uint32_t end
        long i
        float64_t  depth
        uint32_t last_start
        uint32_t last_end
        long pointer
        ndarray[uint32_t, ndim=1] result_position
        ndarray[float64_t, ndim=1] result_depth
        uint32_t[:] positions
        float64_t[:] depths

    assert starts.shape[0] == ends.shape[0], "starts and ends size do not match"

    start_positions = np.sort(starts)
    end_positions = np.sort(ends)

    result_position = np.empty(starts.shape[0] + ends.shape[0], dtype=np.uint32)
    result_depth = np.empty(starts.shape[0] + ends.shape[0], dtype=np.float64)
    positions = result_position
    depths = result_depth

    size_start = starts.shape[0]
    size_end = ends.shape[0]
    pointer = 0
    depth = 0
    i_start = 0
    i_end = 0
    last_start = 0
    last_end = 0

    # loop through all tags
    while i_start < size_start and i_end < size_end:
        start = start_positions[i_start]
        end = end_positions[i_end]
        if start < end:  # depth up
            if start != last_start:  # region changed
                # add last region
                positions[pointer] = start
                depths[pointer] = depth
                pointer += 1
                last_start = start
            # depth increase in next region
            depth += 1
            i_start += 1

        elif start > end:  # depth down
            if end != last_end:  # region changed
                # add last region
                positions[pointer] = end
                depths[pointer] = depth
                pointer += 1
                last_end = end
            # depth decrease in next region
            depth -= 1
            i_end += 1

        else:  # depth equal, skip
            i_start += 1
            i_end += 1

    # all start positions have processed, but some end positions last
    for i in range(i_end, size_end):
        end = end_positions[i]
        if end != last_end:
            positions[pointer] = end
            depths[pointer] = depth
            pointer += 1
            last_end = end
        depth -= 1

    # resize to save memory space
    result_position = np.resize(result_position, pointer)
    result_depth = np.resize(result_depth, pointer)

    return result_position, result_depth

# optimize performance by turn off these checks
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef tuple pileup(ndarray[uint32_t, ndim=1] starts, ndarray[uint32_t, ndim=1] ends):
    return _pileup(starts, ends)

cdef float64_t _logp(float64_t x, float64_t mu):
    return -(x * log(mu) - lgamma(x + 1) - mu)

cpdef float64_t logp(float64_t x, float64_t mu):
    return _logp(x, mu)
