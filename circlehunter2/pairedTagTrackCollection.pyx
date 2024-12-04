# cython: language_level=3

from natsort import natsorted

from numpy cimport uint32_t

from circlehunter2.pairedTagTrack cimport PairedTagTrack

cdef class PairedTagTrackCollection:
    """
    Class representing all paired tags from multi samples or cells.
    """
    def __init__(self, dict chrom_sizes, uint32_t buffer_size=256):
        self._buffer_size = buffer_size

        self._chroms = dict()
        self._shifts = dict()
        shift = 0
        chroms = natsorted(chrom_sizes.keys())
        for chrom in chroms:
            self._chroms[chrom] = chrom_sizes[chrom]
            self._shifts[chrom] = shift
            shift += chrom_sizes[chrom]

        self._tracks = dict()

        self._total = 0
        self._closed = False


    cpdef void add_tag_pair(self, str barcode,
                            str chrom1, uint32_t start1, bint reverse1,
                            str chrom2, uint32_t start2, bint reverse2):
        """
        Add tag pair to this collection
        """
        assert not self._closed

        if barcode not in self._tracks:
            self._tracks[barcode] = PairedTagTrack(
                chrom_sizes=self._chroms, buffer_size=self._buffer_size
            )
        self._tracks[barcode].add_tag_pair(
            chrom1, start1, reverse1, chrom2, start2, reverse2
        )


    cpdef void add_track(self, str barcode, PairedTagTrack track):
        """
        Add track to this collection
        """
        assert not self._closed

        if barcode not in self._tracks:  # new barcode
            self._tracks[barcode] = PairedTagTrack(
                chrom_sizes=self._chroms, buffer_size=self._buffer_size
            )
        # add track
        self._tracks[barcode].add_track(track)


    cpdef void add_collection(self, PairedTagTrackCollection collection):
        """
        Add collection to this collection
        """
        assert not self._closed

        barcodes = collection.get_barcodes()
        # loop through barcodes
        for barcode in barcodes:
            track = collection.get_track(barcode)
            self.add_track(barcode, track)


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
            # stats
            self._total += track.total

        # mark as closed
        self._closed = True


    cpdef set_chromosome_sizes(self, dict chrom_sizes):
        """
        Set the chromosome size of this collection.
        Note: Tags out of this chromosome list or chromosome regions is filter out.
        """
        assert self._closed, "Track not closed!"

        self._total = 0

        for barcode in self.get_barcodes():
            self._tracks[barcode].set_chromosome_sizes(chrom_sizes)
            self._total += self._tracks[barcode].total

        # set size of this collection
        chroms = natsorted(self._chroms.keys())
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

        for barcode in self.get_barcodes():
            if barcode not in input_barcodes:
                track = self._tracks.pop(barcode)
                self._total -= track.total
                track = None


    @property
    def total(self):
        """
        Tag pair counts in this collection
        """
        assert self._closed, "track not closed!"

        return self._total


    def get_barcodes(self):
        """
        Chromosomes stored in this collection.
        """
        assert self._closed

        return tuple(natsorted(self._tracks.keys()))


    def get_track(self, str barcode):
        """
        Return the track associated with given barcode
        """
        assert self._closed, "track not closed!"

        if barcode not in self._tracks:
            raise KeyError(f"No such track associated with barcode: {barcode}")

        return self._tracks[barcode]


    cpdef PairedTagTrack pool_barcodes(self, list barcodes):
        """
        Pool tracks within the given barcodes and return a new PairedTagTrack instance.
        """
        assert self._closed, "track not closed!"

        # barcodes intersect
        pool = set(barcodes) & set(self.get_barcodes())

        track = PairedTagTrack(self._chroms, buffer_size=self._buffer_size)

        for barcode in pool:
            track.add_track(self._tracks[barcode])

        track.close()

        return track
