import libshorah

def test_exec_dpm_sampler():
    retcode = libshorah.exec_dpm_sampler(
        "data_4/w-HXB2-2335-2535.reads.fas", 
        1305,
        0.1,
        261,
        K_cluster_start=20, 
        R_seed=42
    )

    assert retcode == 0