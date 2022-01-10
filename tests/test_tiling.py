from shorah import tiling
import pytest

def test_equispaced():
    strategy = tiling.EquispacedTilingStrategy("HBX2:201-268", 201, 67)
    actual = strategy.get_window_tilings()

    expected = [(0, 201), (67, 201), (134, 201), (201, 201), (268, 201)]

    assert len(actual) == len(expected)
    assert all([a == b for a, b in zip(actual, expected)])

def test_equispaced_with_HBX2():
    end = 3713
    strategy = tiling.EquispacedTilingStrategy(f"HBX2:2469-{end}", 201, 67)
    actual = strategy.get_window_tilings()
    print(actual)

    assert actual[0][0] == 2268
    assert actual[-1][0] == 3675 
    assert actual[-1][0] + actual[-1][1] >= end


def test_equispaced_wrong_incr():
    with pytest.raises(ValueError):
        tiling.EquispacedTilingStrategy("HBX:100-200", 201, 68)


def test_primer_init():
    strategy = tiling.PrimerTilingStrategy("./data_1/scheme.insert.bed")
    window_tilings = strategy.get_window_tilings()
    first = window_tilings[0]
    last = window_tilings[-1]

    assert first[0] == 34
    assert first[1] == 373

    assert last[0] == 9338 - 1
    assert last[1] == 357

    assert len(window_tilings()) == 93