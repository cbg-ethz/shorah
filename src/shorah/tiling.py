from abc import ABC, abstractmethod
from typing import List, Tuple

class TilingStrategy(ABC):
    """An abstract class that defines the tiling strategy i.e. how a reference
    gnome is split into windows for later analysis. 

    Only one methods needs to specified by concrete classes.
    """

    @abstractmethod
    def get_window_tilings(self, start: int, end: int) -> List[Tuple[int, int]]:
        """
        Args:
            start: Start of region in reference genome (1 based like samtools).  
            end: End of region (inclusive).

        Returns:
            A list of tuples that indicate the starting point and length of
            each window. It assumed that this list is **sorted** by starting 
            points ascending. First integer is the starting position, second the 
            window length. For example (with fixed window length 201)::

                [(0, 201), (67, 201), (134, 201), (201, 201)]

        """
        pass

class EquispacedTilingStrategy(TilingStrategy):
    """Implements a tiling strategy that puts windows of fixed length at equally
    spaced distances. 

    Default instantiation with `EquispacedTilingStrategy()` will result in   
    a window of length 201. In that case, each window will have an overlap of 
    2/3 with the previous window.  

    Attributes:
        window_length: Constant number of bases considered at once per loop.
        incr: Increment between each window.
        exact_conformance_overlap_at_boundary: The old implementation in C++
            assumes that the `incr` is always one third of the `window_length`.
            This leads to issues on the boundary. Set to `Ture` for testing vs. 
            the old implementation. 
    """

    def __init__(self, window_length=201, incr=201//3, 
        exact_conformance_overlap_at_boundary=False) -> None:

        if window_length%incr != 0 or incr <= 0 or window_length <= 0:
            raise ValueError("window_length has to be divisible by incr")
        
        self.window_length = window_length
        self.incr = incr
        self.exact_conformance_overlap_at_boundary = exact_conformance_overlap_at_boundary

    def get_window_tilings(self, start: int, end: int) -> List[Tuple[int, int]]:
        """Implements :meth:`~shorah.tiling.TilingStrategy.get_window_tilings`.
        """
        
        window_positions = list(range(
            start - self.incr * 3 if self.exact_conformance_overlap_at_boundary else start - self.window_length,
            end + 1 - (self.window_length//self.incr - 3) * self.incr if self.exact_conformance_overlap_at_boundary else end + 1, 
            self.incr 
        ))

        # add one more window at the end
        if self.exact_conformance_overlap_at_boundary == True:
            window_positions.append(window_positions[-1] + self.incr)
        
        return [(i, self.window_length) for i in window_positions]

class PrimerTilingStrategy(TilingStrategy):
    def __init__(self, primer_path: str) -> None:
        pass