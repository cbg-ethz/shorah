from shorah import tiling
import pytest

def test_equispaced():
    strategy = tiling.EquispacedTilingStrategy(201, 67)
    actual = strategy.get_window_tilings(201, 268)

    expected = [(0, 201), (67, 201), (134, 201), (201, 201), (268, 201)]

    assert len(actual) == len(expected)
    assert all([a == b for a, b in zip(actual, expected)])

def test_equispaced_with_HBX2():
    strategy = tiling.EquispacedTilingStrategy(201, 67)
    end = 3713
    actual = strategy.get_window_tilings(2469, end)
    print(actual)

    assert actual[0][0] == 2268
    assert actual[-1][0] == 3675 
    assert actual[-1][0] + actual[-1][1] >= end


def test_equispaced_wrong_incr():
    with pytest.raises(ValueError):
        tiling.EquispacedTilingStrategy(201, 68)