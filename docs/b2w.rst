b2w
===

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

.. autofunction:: shorah.b2w.b2w