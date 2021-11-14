b2w
===

.. danger::
   The `b2w.py` differs from the deprecated `b2w.cpp` in the following points:
   
      * Can no longer run with a samtools compatible region
      * etc

Window tiling strategies
************************

Equispaced
----------
The window moves along the region in equally spaced steps.  

User-defined
------------
The user provides an additional file (i.e. an "insert file") that specifies where the amplification was 
performed during sequencing.

Adaptive
--------
Window starting positions are calculated endogenously.

API
***

.. caution::
   ``reads.fas`` does not comply with the FASTA format.    

.. autofunction:: shorah.b2w.build_windows

.. autoclass:: shorah.tiling.TilingStrategy
   :members:

.. autoclass:: shorah.tiling.EquispacedTilingStrategy
   :members: