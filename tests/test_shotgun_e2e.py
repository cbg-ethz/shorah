import subprocess
import filecmp

def test_e2e_shorah():
    original = subprocess.run(
        "./shotgun_test.sh", shell=True, check=True, cwd="./data_1"
    )
    assert original.returncode == 0

    assert filecmp.cmp(
        "./data_1/test.csv",
        "./data_1/snv/SNVs_0.010000_final.csv",
        shallow=False
    )
