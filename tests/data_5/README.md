### Test to check fil.cpp implementation accounting for long deletions

Test files `SNV.txt` and `SNVs_0.010000.txt` are obtained by running `shorah shutgun`, e.g:

```
shorah shotgun -a 0.1 -w 42 -x 100000 -p 0.9 -c 0 -r REF:42-272 -R 42 -b test_aln.cram -f ref.fasta
```

The test script `test_long_deletions.py` uses `pysam` and `NumPy`, which can be installed using pip or conda.  
