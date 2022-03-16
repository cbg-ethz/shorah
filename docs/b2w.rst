Window building (b2w)
=====================

.. caution::
   The new ``b2w.py`` differs from the deprecated ``b2w.cpp`` in the following points:

      * ``b2w.cpp`` should be called directly in Python through pybind11 through ``libshorah.b2w(...)`` an appropriate arguments. 
      * ``reads.fas`` is a file not conformant to FASTA that is generated during a run of each version. The second and third column of ``reads.fas`` cannot be exactly reproduced (off-by-one bug).

.. note::
   The pysam library uses 0-based indexing which comes more natural to Python. This project, however, uses *1-based indexing* 
   in all (ephemeral) output files if not specified otherwise.

Window tiling strategies
************************

Equispaced
----------
The window moves along the region in equally spaced steps. This is the *default* behavior.   

Based on primer scheme
----------------------
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