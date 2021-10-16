import filecmp
import os
import glob
from shorah import b2w

def test_cmp_raw():
    b2w.main()
    p = os.path.dirname(__file__)

    spec_files = [os.path.basename(x) for x in glob.glob(os.path.join(p, 'data/*.fas'))]
    spec_files.extend(['coverage.txt'])

    match, mismatch, errors = filecmp.cmpfiles(
        os.path.join(p, 'data'),
        p,
        spec_files,
        shallow=False
    ) 
    print(match)
    print(mismatch)
    print(errors)

    assert len(mismatch) == len(errors) == 0
