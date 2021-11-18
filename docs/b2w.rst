b2w
===

.. danger::
   The `b2w.py` differs from the deprecated `b2w.cpp` in the following points:
   
      * Can no longer run with a samtools compatible region
      * etc # TODO

Window tiling strategies
************************

Equispaced
----------
The window moves along the region in equally spaced steps.  

Based on primer scheme
-------------------
Information on the postion of amplicons during the sequencing processes is exploited to place windows
in a fashion so that partial overlapping is avoided. 

User-defined
------------
Users may implement their own tiling strategy by deriving their own class from
from the `TilingStrategy` abstract class.  

API
*** 

.. autofunction:: shorah.b2w.build_windows

.. autoclass:: shorah.tiling.TilingStrategy
   :members:

.. autoclass:: shorah.tiling.EquispacedTilingStrategy
   :members:

.. autoclass:: shorah.tiling.PrimerTilingStrategy
   :members: