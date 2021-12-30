import subprocess
import os
import helper_long_deletions
# TODO run through Python directly

dir = "./data_5"

def test_long_deletions():
    original = subprocess.run(
        "./shotgun_prepare.sh", shell=True, check=True, cwd=dir
    )
    assert original.returncode == 0

    p = os.path.dirname(__file__)
    os.chdir(os.path.join(p, dir))
    print(os.getcwd())
    # Input data
    bamfile = "test_aln.cram"
    snvsfile = "snv/SNV.txt"
    outfile = "snv/SNVs_0.010000.txt"

    helper_long_deletions.main(bamfile=bamfile, snvsfile=snvsfile, outfile=outfile)

    os.chdir(p)
