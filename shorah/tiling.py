from abc import ABC, abstractmethod
from typing import List, Tuple

class TilingStrategy(ABC):
    """An abstract class that defines the tiling strategy i.e. how a reference
    gnome is split into windows for later analysis.

    Only one methods needs to specified by concrete classes.
    """

    @abstractmethod
    def get_window_tilings(self) -> List[Tuple[int, int]]:
        """
        Returns:
            A list of tuples that indicate the starting point and length of
            each window. It is assumed that this list is **sorted** by starting
            points ascending. First integer is the **1-based** starting position,
            second the window length. For example (with fixed window length 201)::

                [(1, 201), (68, 201), (135, 201), (202, 201)]

        """
        pass

    @abstractmethod
    def get_reference_name() -> str:
        """
        Returns:
            Name of the reference genome.
        """
        pass

    @abstractmethod
    def get_region_end() -> str:
        pass


class EquispacedTilingStrategy(TilingStrategy):
    """Implements a tiling strategy that puts windows of fixed length at equally
    spaced distances.

    Default instantiation with `EquispacedTilingStrategy(region)` will result in
    a window of length 201. In that case, each window will have an overlap of
    2/3 with the previous window.

    Attributes:
        reference_name: Name of the reference genome.
        start: Start of region in reference genome (1-based like samtools).
        end: End of region (inclusive).
        window_length: Constant number of bases considered at once per loop.
        incr: Increment between each window.
        exact_conformance_overlap_at_boundary: The old implementation in C++
            assumes that the `incr` is always one third of the `window_length`.
            This leads to issues on the boundary. Set to `True` for testing vs.
            the old implementation.
        use_full_reference_as_region: Assume that the region string includes the
            reference genome in its full length.
    """

    def __init__(self, region, window_length=201, incr=201//3,
        exact_conformance_overlap_at_boundary=False,
        use_full_reference_as_region=False) -> None:

        if window_length%incr != 0 or incr <= 0 or window_length <= 0:
            raise ValueError("window_length has to be divisible by incr")
        if region == "":
            raise ValueError("empty region string is not allowed, use e.g. chr1:50-150")
        if exact_conformance_overlap_at_boundary == True and use_full_reference_as_region == True:
            raise ValueError("this combination of arguments is not allowed")

        self.window_length = window_length
        self.incr = incr
        self.exact_conformance_overlap_at_boundary = exact_conformance_overlap_at_boundary
        self.use_full_reference_as_region = use_full_reference_as_region

        reference_name, start, end = self.__parse_region(region)
        self.reference_name = reference_name
        self.start = start
        self.end = end

    def __parse_region(self, region: str) -> Tuple[str, int, int]:
        """
        Args:
            region: A genomic sequence compatibile with samtools.
                Example: "chr1:10000-20000"
        Returns:
            start: Start of region in reference genome (1-based like samtools).
            end: End of region (inclusive).
        """
        tmp = region.split(":")
        reference_name = tmp[0]
        tmp = tmp[1].split("-")
        start = int(tmp[0])
        end = int(tmp[1])
        return reference_name, start, end

    def get_window_tilings(self) -> List[Tuple[int, int]]:
        """Implements :meth:`~shorah.tiling.TilingStrategy.get_window_tilings`.
        """
        if self.use_full_reference_as_region == True:
            window_positions = list(range(
                1,
                self.end - 1,
                self.incr
            ))
            #while window_positions[-1] + self.window_length >= self.end:
            #    del window_positions[-1] # FIXME uncommented to create one single window

        else:
            window_positions = list(range(
                self.start - self.incr * 3 if self.exact_conformance_overlap_at_boundary
                    else self.start - self.window_length, # this is 1-based
                self.end + 1 - (self.window_length//self.incr - 3) * self.incr if self.exact_conformance_overlap_at_boundary
                    else self.end + 1, # TODO why +1
                self.incr
            ))

        # add one more window at the end
        if self.exact_conformance_overlap_at_boundary == True:
            window_positions.append(window_positions[-1] + self.incr)

        return [(i, self.window_length) for i in window_positions]

    def get_reference_name(self):
        return self.reference_name

    def get_region_end(self):
        return self.end


class PrimerTilingStrategy(TilingStrategy):
    """Implements a tiling strategy that it is based on the primer scheme used
    for sequencing.

    Attributes:
        amplicons: A data structure containing the coordinates of each primer in
            the scheme relative to the reference genome.

            See more information here:
            https://artic.readthedocs.io/en/latest/primer-schemes/

            The datastructure is based upon the `foobar.insert.bed` file
            described in the documentation. See examples here:
            https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019/V3

            https://primalscheme.com
            https://github.com/aresti/primalscheme
    """

    def __init__(self, insert_bed_path: str) -> None:
        self.amplicons: List[Tuple[int, int]] = []
        with open(insert_bed_path) as f:
            last_ref_seq = None
            for line in f:
                L = line.strip().split()
                # 0-based, exclusive
                self.amplicons.append((int(L[1]), int(L[2])))

                if last_ref_seq != L[0] and last_ref_seq != None:
                    raise InputError("Insert files with more than on reference \
                        sequence are not supported.")
                else:
                    last_ref_seq = L[0]
            self.reference_name = last_ref_seq

    def get_window_tilings(self) -> List[Tuple[int, int]]:
        rv = []
        #start = self.amplicons[0][0]
        #end = self.amplicons[-1][1]
        for amplicon in self.amplicons:
            rv.append( (amplicon[0], amplicon[1] - amplicon[0]) )
            # add two more windows +/-50bp
            #if amplicon[0]-50 > start:
            #    rv.append( (amplicon[0]-50, amplicon[1] - amplicon[0]) )
            #if amplicon[1]+50 < end:
            #    rv.append( (amplicon[0]+50, amplicon[1] - amplicon[0]) )
            #print(rv)

            # TODO: for amplicon mode add more windows ???

        return rv

    def get_reference_name(self) -> str:
        return self.reference_name

    def get_region_end(self):
        return self.amplicons[-1][1]
