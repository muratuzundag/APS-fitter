import numpy as np
import pytest

from compact_spacing import Splitting, correct_to_m0, fit_sequence, gls_line


def test_known_weighted_line():
    k = np.array([1.0, 2.0, 5.0, 9.0, 11.0, 12.0, 14.0, 16.0])
    p = np.array([440.200, 457.583, 513.200, 585.728, 618.124, 638.109, 679.172, 714.105])
    e = np.array([1.5, 1.0, 2.4, 0.5, 1.0, 0.1, 0.5, 2.6])
    r = fit_sequence("l1", k, p, e, n_samples=100, seed=2)
    assert r.spacing_gls_s == pytest.approx(18.156890313753294)
    assert r.spacing_gls_error_s == pytest.approx(0.06838209249044948)
    assert r.dof == 6
    assert r.p16_s < r.spacing_s < r.p84_s


def test_frequency_correction():
    p = np.array([440.2, 457.583, 513.2])
    m = np.array([-1.0, 0.0, 1.0])
    split = Splitting(2.0, 0.1)
    obs = 1e6 / (1e6 / p + m * split.value_uhz)
    corrected, _, derivative = correct_to_m0(obs, np.full(3, 0.1), m, split)
    assert np.allclose(corrected, p)
    assert derivative[1] == 0


def test_reproducibility_and_order_offset():
    k = np.array([1, 2, 5, 9])
    p = np.array([440.2, 457.583, 513.2, 585.728])
    kwargs = {"period_s": p, "period_error_s": np.ones(4), "n_samples": 100, "seed": 21}
    a = fit_sequence("a", k=k, **kwargs)
    b = fit_sequence("a", k=k, **kwargs)
    c = fit_sequence("a", k=k + 10, **kwargs)
    assert np.array_equal(a.samples_s, b.samples_s)
    assert a.spacing_gls_s == pytest.approx(c.spacing_gls_s)


@pytest.mark.parametrize("errors", [[1, 0, 1], [1, -1, 1], [1, np.nan, 1]])
def test_invalid_errors(errors):
    with pytest.raises(ValueError):
        fit_sequence("bad", [1, 2, 3], [440, 460, 480], errors, n_samples=100)


def test_nonzero_m_without_splitting():
    with pytest.raises(ValueError):
        fit_sequence("bad", [1, 2, 3], [440, 460, 480], [1, 1, 1],
                     m=[-1, 0, 1], n_samples=100)


def test_degenerate_orders():
    with pytest.raises(ValueError):
        fit_sequence("bad", [1, 1, 1], [440, 460, 480], [1, 1, 1], n_samples=100)


def test_invalid_covariance():
    with pytest.raises(ValueError):
        gls_line([1, 2, 3], [440, 460, 480], np.zeros((3, 3)))


def test_invalid_splitting():
    with pytest.raises(ValueError):
        Splitting(-1)
