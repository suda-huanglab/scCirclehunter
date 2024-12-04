# cython: language_level=3

from multiprocessing import Process, Queue
import re, queue

from tqdm import tqdm
import numpy as np

from numpy cimport uint8_t, uint16_t, uint32_t, float64_t
from pysam.libcalignmentfile cimport AlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile
from libc.stdint cimport UINT32_MAX
from cpython cimport bool
cimport cython

from circlehunter2.tagTrackCollection import TagTrackCollection
from circlehunter2.pairedTagTrackCollection import PairedTagTrackCollection
from circlehunter2.utils import make_windows as  make_windows_util

from circlehunter2.tagTrackCollection cimport TagTrackCollection
from circlehunter2.pairedTagTrackCollection cimport PairedTagTrackCollection
from circlehunter2.pairedTagTrack cimport PairedTagTrack


QUEUE_END = "QUEUE_END"
DIRECTION_INDEX = {False: -1, True: 0}


cpdef dict collect_chromosome_sizes(str filename):
    bam = AlignmentFile(filename, "rb")

    # check bam file and index
    assert bam.is_bam, f"Not a bam file: {bam.filename}"
    assert bam.check_index(), f"Index not exist for: {bam.filename}"

    # collect chromosome size
    chrom_sizes = dict()
    for chrom, chrom_size in zip(bam.references, bam.lengths):
        chrom_sizes[chrom] = chrom_size

    return chrom_sizes


cpdef tuple make_windows(str filename, long window_size):
        """
        Return the windows of chromosomes by the given windows size
        """
        cdef:
            list windows
            dict chrom_sizes
            long chrom_size
            str chrom
            long last
            long current

        chrom_sizes = collect_chromosome_sizes(filename)

        return make_windows_util(chrom_sizes, window_size)


cdef bint pass_filter(AlignedSegment read,
                      uint16_t flag_include=1, uint16_t flag_exclude=1548,
                      uint8_t mapq=10, float64_t mismatch_ratio=0.02,
                      str barcode_tag="CB"):
    """
    Judge whether the read pass filters
    """
    if not read.flag & flag_include:
        return 0
    if read.flag & flag_exclude:
        return 0
    if read.mapq < mapq:
        return 0
    if barcode_tag is not None and not read.has_tag(barcode_tag):
        return 0
    if read.has_tag('MD'):  # judge mismatch ratio by MD tag
        mismatch = len(re.findall(r'[ATCG]', read.get_tag('MD')))
        match = sum(
            map(int, (re.findall(r'(\d+)M', read.cigarstring)))
        )
        if mismatch > match * mismatch_ratio:
            return 0
    return 1


cpdef tuple build_tag_track_single_process(
        str filename, uint32_t min_insert_size=5000,
        uint16_t flag_include=1, uint16_t flag_exclude=1548,
        uint8_t mapq=10, float64_t mismatch_ratio=0.02, str barcode_tag="CB",
        uint32_t buffer_size=256, bool no_progressbar=True):
    """
    Read BAM file and build as a TagTrackCollection
    """
    cdef:
        AlignmentFile bam
        AlignedSegment read

        TagTrackCollection tag_collection
        PairedTagTrackCollection paired_tag_collection

        str barcode

        str chrom1
        uint32_t start1
        uint32_t end1
        bint reverse1

        str chrom2
        uint32_t start2
        uint32_t end2
        bint reverse2

        uint32_t insert

    bam = AlignmentFile(filename, "rb")

    # check bam file and index
    assert bam.is_bam, f"Not a bam file: {bam.filename}"
    assert bam.check_index(), f"Index not exist for: {bam.filename}"

    chrom_sizes = collect_chromosome_sizes(filename)

    tag_collection = TagTrackCollection(buffer_size=buffer_size)
    paired_tag_collection = PairedTagTrackCollection(
        chrom_sizes=chrom_sizes, buffer_size=buffer_size
    )

    pbar = tqdm(
        disable=no_progressbar,
        desc="reading reads from bam file", unit=" reads"
    )
    # loop through reads
    for read in bam:
        if pass_filter(read, flag_include, flag_exclude, mapq, mismatch_ratio, barcode_tag):
            if barcode_tag is not None:
                barcode = read.get_tag(barcode_tag)
            else:
                barcode = "NO_BARCODE"

            chrom1 = read.reference_name
            start1 = read.reference_start
            end1 = read.reference_end
            reverse1 = read.is_reverse

            insert = abs(read.template_length)

            # add tag
            tag_collection.add_tag(barcode, chrom1, start1, end1, reverse1, insert)

            # discordant reads
            if read.mate_is_mapped:
                chrom2 = read.next_reference_name
                start2 = read.next_reference_start
                reverse2 = read.mate_is_reverse
                if read.reference_name == read.next_reference_name:  # same chromosome
                    # may be negative
                    if reverse1 == reverse2:  # same direction
                        paired_tag_collection.add_tag_pair(
                            barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                        )
                    elif insert >= min_insert_size:  # large insert size
                        paired_tag_collection.add_tag_pair(
                            barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                        )
                else:  # different chromosome
                    paired_tag_collection.add_tag_pair(
                        barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                    )
            # end discordant
        pbar.update(1)
    pbar.close()

    tag_collection.close()
    paired_tag_collection.close()

    return tag_collection, paired_tag_collection


cpdef PairedTagTrack build_paired_tag_track_single_process(
        str filename, uint32_t min_insert_size=5000,
        uint16_t flag_include=1, uint16_t flag_exclude=1548,
        uint8_t mapq=10, float64_t mismatch_ratio=0.02, str barcode_tag="CB",
        uint32_t buffer_size=256, bool no_progressbar=True):
    """
    Read BAM file and build as a PairedTagTrack.
    Only used for test!!!
    """
    cdef:
        AlignmentFile bam
        PairedTagTrack track

        str chrom1
        uint32_t start1
        uint32_t end1
        bint reverse1

        str chrom2
        uint32_t start2
        uint32_t end2
        bint reverse2

        uint32_t insert

    bam = AlignmentFile(filename, "rb")

    # check bam file and index
    assert bam.is_bam, f"Not a bam file: {bam.filename}"
    assert bam.check_index(), f"Index not exist for: {bam.filename}"

    chrom_sizes = collect_chromosome_sizes(filename)

    track = PairedTagTrack(chrom_sizes=chrom_sizes, buffer_size=buffer_size)

    pbar = tqdm(
        disable=no_progressbar,
        desc="reading reads from bam file", unit=" reads"
    )
    # loop through reads
    for read in bam:
        if pass_filter(read, flag_include, flag_exclude, mapq, mismatch_ratio, barcode_tag):
            # barcode = read.get_tag(barcode_tag)
            chrom1 = read.reference_name
            start1 = read.reference_start
            reverse1 = read.is_reverse
            chrom2 = read.next_reference_name
            start2 = read.next_reference_start
            reverse2 = read.mate_is_reverse
            # may be negative, no need
            if read.reference_name == read.next_reference_name:
                # may be negative
                insert = abs(read.template_length)
            else:  # different reference
                insert = UINT32_MAX

            if insert >= min_insert_size:
                # add tag
                track.add_tag_pair(chrom1, start1, reverse1, chrom2, start2, reverse2)
        pbar.update(1)
    pbar.close()

    track.close()

    return track


cdef tuple _fetch_tags(
        dict chrom_sizes, AlignmentFile bam,
        str chrom, uint32_t start, uint32_t end, uint32_t min_insert_size=5000,
        uint16_t flag_include=1, uint16_t flag_exclude=1548,
        uint8_t mapq=10, float64_t mismatch_ratio=0.02, str barcode_tag="CB",
        uint32_t buffer_size=256):
    """
    Fetch read from BAM file and build TagTrackCollection.
    Used at multi process func only!
    """
    cdef:
        AlignedSegment read

        TagTrackCollection tag_collection
        PairedTagTrackCollection paired_tag_collection

        str barcode

        str chrom1
        uint32_t start1
        uint32_t end1
        bint reverse1

        str chrom2
        uint32_t start2
        uint32_t end2
        bint reverse2

        uint32_t insert

    tag_collection = TagTrackCollection(buffer_size=buffer_size)
    paired_tag_collection = PairedTagTrackCollection(
        chrom_sizes=chrom_sizes, buffer_size=buffer_size
    )

    for read in bam.fetch(chrom, start, end):
        if read.reference_start < start:
            continue
        if not pass_filter(read, flag_include, flag_exclude, mapq, mismatch_ratio, barcode_tag):
            continue

        if barcode_tag is not None:
            barcode = read.get_tag(barcode_tag)
        else:
            barcode = "NO_BARCODE"

        chrom1 = read.reference_name
        start1 = read.reference_start
        end1 = read.reference_end
        reverse1 = read.is_reverse

        insert = abs(read.template_length)

        # add tag
        tag_collection.add_tag(barcode, chrom1, start1, end1, reverse1, insert)

        # discordant reads?
        if read.mate_is_mapped:
            chrom2 = read.next_reference_name
            start2 = read.next_reference_start
            reverse2 = read.mate_is_reverse
            if read.reference_name == read.next_reference_name:  # same chromosome
                # may be negative
                if reverse1 == reverse2:  # same direction
                    paired_tag_collection.add_tag_pair(
                        barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                    )
                elif insert >= min_insert_size:  # large insert size
                    paired_tag_collection.add_tag_pair(
                        barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                    )
            else:  # different chromosome
                paired_tag_collection.add_tag_pair(
                    barcode, chrom1, start1, reverse1, chrom2, start2, reverse2
                )
        # end discordant

    tag_collection.close()
    paired_tag_collection.close()

    return tag_collection, paired_tag_collection

def _fetch_tags_process(filename, regions, results, min_insert_size=5000,
                        flag_include=1, flag_exclude=1548,
                        mapq=10, mismatch_ratio=0.02,
                        barcode_tag="CB",
                        buffer_size=256):
    """
    Multi processes target func to fetch tags as TagTrackCollection.
    Used at multi processed func only!
    """
    chrom_sizes = collect_chromosome_sizes(filename)

    bam = AlignmentFile(filename, 'rb')

    while True:
        # consume the input regions
        try:
            region = regions.get_nowait()
        except queue.Empty:  # queue not ready, keep waiting
            continue
        if region != QUEUE_END:  # got a region
            chrom, start, end = region
            # fetch tags
            collections = _fetch_tags(
                chrom_sizes, bam,
                chrom, start, end, min_insert_size,
                flag_include, flag_exclude,
                mapq, mismatch_ratio,
                barcode_tag,
                buffer_size
            )
            # produce result
            results.put(collections)
        else:  # queue end
            # put the queue end mark into queue to stop other processes
            regions.put(QUEUE_END)
            break
            # process end


def build_tag_track_multi_process(str filename,
                                  uint32_t min_insert_size=5000,
                                  uint16_t flag_include=1, uint16_t flag_exclude=1548,
                                  uint8_t mapq=10, float64_t mismatch_ratio=0.02,
                                  str barcode_tag="CB",
                                  uint32_t buffer_size=256,
                                  processes=2, window_size=64_000_000,
                                  no_progressbar=True):
    """
    Read BAM file and build as a TagTrackCollection using multi processes.
    Each process will process a given window size of the genome.
    """
    assert processes > 1, "processes less than 2, use single cell method!"

    chrom_sizes = collect_chromosome_sizes(filename)

    tag_collection = TagTrackCollection(buffer_size=buffer_size)
    paired_tag_collection = PairedTagTrackCollection(
        chrom_sizes=chrom_sizes, buffer_size=buffer_size
    )

    # create queue
    queue_regions, queue_results = Queue(), Queue()

    # make windows and put into the queue as jobs
    windows = make_windows(filename, window_size)
    for chrom, start, end in iter(windows):
        queue_regions.put((chrom, start, end))
    queue_regions.put(QUEUE_END)  # end the queue

    pbar = tqdm(
        total=len(windows), disable=no_progressbar,
        desc="reading chunks from bam file", unit=" chunks"
    )
    # start processes
    for _ in range(processes):
        Process(
            target=_fetch_tags_process,
            args=(
                filename, queue_regions, queue_results, min_insert_size,
                flag_include, flag_exclude,
                mapq, mismatch_ratio,
                barcode_tag,
                buffer_size
            )
        ).start()

    # processing results
    finished = 0
    while True:
        # check how many regions has finished
        size = queue_results.qsize()
        if size == 0:  # no result now, skip
            continue

        # collect results
        for _ in range(size):
            collections = queue_results.get()
            if collections[0].total > 0:
                tag_collection.add_collection(collections[0])
            if collections[1].total > 0:
                paired_tag_collection.add_collection(collections[1])
        pbar.update(size)

        # check if all regions has been processed
        finished += size
        if finished == len(windows):
            break  # finished
    pbar.close()

    tag_collection.close()
    paired_tag_collection.close()

    return tag_collection, paired_tag_collection


@cython.nonecheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef tuple fetch_breakpoints(AlignmentFile bam,
                              str chrom1, uint32_t start1, uint32_t end1, bool reverse1,
                              str chrom2, uint32_t start2, uint32_t end2, bool reverse2,
                              uint16_t flag_include=1, uint16_t flag_exclude=1548,
                              uint8_t mapq=10, float64_t mismatch_ratio=0.02,
                              str barcode_tag="CB", uint32_t buffer_size=256
):
    """
    Fetch breakpoint positions for the given paired peak.
    """
    cdef:
        AlignedSegment read

        uint32_t[:] discordant_pos
        uint32_t discordant_size
        uint32_t discordant_pointer

        uint32_t[:] clipped_pos
        uint32_t clipped_size
        uint32_t clipped_pointer


    discordant_pos = np.empty(buffer_size, dtype=np.uint32)
    discordant_size = buffer_size
    discordant_pointer = 0

    clipped_pos = np.empty(buffer_size, dtype=np.uint32)
    clipped_size = buffer_size
    clipped_pointer = 0

    for read in bam.fetch(chrom1, start1, end1):
        if read.reference_start < start1:
            continue
        if not pass_filter(read, flag_include, flag_exclude, mapq, mismatch_ratio, barcode_tag):
            continue
        if read.reference_end > end1:  # unmapped read do no have reference_end
            continue
        if read.is_reverse != reverse1:  # may be a clipped reads
            if (
                # should propre pair
                read.is_proper_pair
                # reads end with a clipped cigar
                and read.cigartuples[DIRECTION_INDEX[reverse1]][0] in {4, 5}
            ):  # this is a clipped reads
                if clipped_pointer >= clipped_size:  # buffer increase needed
                    clipped_size += buffer_size
                    clipped_pos = np.resize(clipped_pos, clipped_size)

                if read.is_reverse:
                    clipped_pos[clipped_pointer] = read.reference_end
                else:
                    clipped_pos[clipped_pointer] = read.reference_start
                clipped_pointer += 1

            # this can't be a discordant reads
            continue
        if read.mate_is_reverse != reverse2:
            continue
        if read.next_reference_name != chrom2:
            continue
        if read.next_reference_start < start2 or read.next_reference_start > end2:
            continue

        # this is a discordant read
        if discordant_pointer >= discordant_size:  # buffer increase needed
            discordant_size += buffer_size
            discordant_pos = np.resize(discordant_pos, discordant_size)
        # add discordant position
        if reverse1:
            discordant_pos[discordant_pointer] = read.reference_start
        else:
            discordant_pos[discordant_pointer] = read.reference_end
        discordant_pointer += 1

        # check if this is a clipped reads, too
        if read.cigartuples[DIRECTION_INDEX[read.is_reverse]][0] in {4, 5}:
            if clipped_pointer >= clipped_size:  # buffer increase needed
                clipped_size += buffer_size
                clipped_pos = np.resize(clipped_pos, clipped_size)

            if read.is_reverse:
                clipped_pos[clipped_pointer] = read.reference_start
            else:
                clipped_pos[clipped_pointer] = read.reference_end
            clipped_pointer += 1

    discordant_pos = np.resize(discordant_pos, discordant_pointer)
    clipped_pos = np.resize(clipped_pos, clipped_pointer)
    return discordant_pos.base, clipped_pos.base


def count_cross_reads(AlignmentFile bam,
                         str chrom, uint32_t pos, bool reverse,
                         uint32_t error=5, uint32_t extend=500,
                         uint16_t flag_include=1, uint16_t flag_exclude=1548,
                         uint8_t mapq=10, float64_t mismatch_ratio=0.02,
                         str barcode_tag="CB"):
    """
    Count reads cross a breakpoint or not.
    """
    fetch_start = pos - error if reverse else pos - extend
    fetch_end = pos + extend if reverse else pos + error

    cross = set()
    disjoint = set()

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        # count only proper paired reads
        if not read.is_proper_pair:
            continue

        if not pass_filter(read, flag_include, flag_exclude, mapq, mismatch_ratio, barcode_tag):
            continue

        # extract paired reads region
        if read.is_reverse:
            end = read.reference_end
            start = end + read.template_length
        else:
            start = read.reference_start
            end = start + read.template_length

        if start <= pos - error and end >= pos + error:  # cross
            cross.add(read.query_name)
        elif reverse and end >= pos - error:  # disjoint reverse
            disjoint.add(read.query_name)
        elif not reverse and end <= pos + error:  # disjoint forward
            disjoint.add(read.query_name)

    return len(cross), len(disjoint)
