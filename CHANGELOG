New in version 0.5.1
---------------------

+ bug fixes in bams2msa.py

+ default behaviour of s2f.py is now not to propagate gaps

New in version 0.5
---------------------

+ SAM/BAM format added! Users can now produce MSA files
  from their alignments in SAM format

+ README and INSTALL have been updated, and more
  information is available on the web, at the ShoRAH
  documentation page
  https://wiki-bsse.ethz.ch/display/ShoRAH/Documentation

New in version 0.4
---------------------

+ dpm_sampler now computes Hamming distance in
  time proportional to the distance rather than
  to the length of the sequences.

+ dpm_sampler now clubs identical reads into
  objects. This method achieves a significant speed-up
  with Illumina datasets. Up to 65000 Illumina reads were
  analysed in a single window with good results.

+ freqEst sorts the output according to the frequency


New in version 0.3
---------------------

+ dpm_sampler is now C++! Thanks to the structures
  defined in C++ libraries (map and multimap), we
  are able to run requiring much less memory

+ Parallelization! dec.py now calls diri_sampler using
  a pool of independent workers, exploiting all the
  available compuational power, as well as s2f.py does

+ The sampling has been improved, now we can have a more
  reliable estimate of the quality of our local
  haplotype reconstruction

+ The output of diri_sampler is now better organized
  (if -k is given, all intermediate files are saved in
  subdirectories of the current)

+ The alignment program runs now in linear time (with
  respect to the number of reads), and deals with indels
  in a more clever way

+ All python programs now use the logging module to
  write logs of their operations

+ plot_sampling.py and plot_stat.py can be used to produce
  graph showing the behaviour of the Gibbs sampling


New in version 0.2
---------------------

+ New method to assign the reads after the sampling

+ The alignment is now provided by a separate
  program (step2far.py), so that the user can input
  his/her own alignment and install EMBOSS only if
  really necessary

+ Fixed numerical bugs when dealing with very high or
  very low probabilities
