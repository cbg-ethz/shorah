CLI
===

.. caution::
    * The CLI option ``-of vcf`` will yield a ``NotImplementedError``. Only ``csv`` output is supported at the moment.
    * The ``csv`` output format changed. The old format was to restricitive because it forced the user to use three windows.  
    * ``shorah amplicon`` mode has been deprecated and is no longer available.
    * ``shorah shotgun`` might be misleading name because not only shotgun sequencing data is supported.

Output files
************

* ``raw_snv.tsv``
* ``raw_snv_<sigma>.tsv``
* ``raw_snv_<sigma>_final.csv``

Replace ``<sigma>`` with the appropriate value passed in through the CLI. These files will be in the folder ``snv`` after
a successful of ``shorah shotgun``.

Interface
*********

.. argparse::
   :module: shorah.cli
   :func: all_parsers
   :prog: shorah
   :path: shotgun
