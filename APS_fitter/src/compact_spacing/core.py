"""Weighted/GLS fitting and Monte Carlo propagation for mean period spacings."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Splitting:
    """First-order rotational splitting in frequency units (microhertz)."""
    value_uhz: float
    error_uhz: float = 0.0


@dataclass
class FitResult:
    sequence: str
    n_modes: int
    intercept_s: float
    spacing_gls_s: float
    spacing_gls_error_s: float
    spacing_s: float
    error_minus_s: float
    error_plus_s: float
    p16_s: float
    p84_s: float
    chi2: float
    dof: int
    reduced_chi2: float
    corrected_period_s: np.ndarray
    corrected_error_s: np.ndarray
    residual_s: np.ndarray
    samples_s: np.ndarray


def validate_inputs(k, period_s, period_error_s, m=None):
    k = np.asarray(k, dtype=float)
    p = np.asarray(period_s, dtype=float)
    e = np.asarray(period_error_s, dtype=float)
    if not (k.ndim == p.ndim == e.ndim == 1 and len(k) == len(p) == len(e)):
        raise ValueError("k, period_s, and period_error_s must be equal-length 1-D arrays")
    if len(k) < 3:
        raise ValueError("At least three modes are required")
    if not (np.all(np.isfinite(k)) and np.all(np.isfinite(p)) and np.all(np.isfinite(e))):
        raise ValueError("Inputs must be finite")
    if np.unique(k).size < 2:
        raise ValueError("At least two distinct radial orders are required")
    if np.any(p <= 0) or np.any(e <= 0):
        raise ValueError("Periods and period errors must be positive")
    if m is None:
        m_arr = np.zeros(len(k), dtype=float)
    else:
        m_arr = np.asarray(m, dtype=float)
        if m_arr.shape != k.shape or not np.all(np.isfinite(m_arr)):
            raise ValueError("m must be finite and have the same shape as k")
    return k, p, e, m_arr


def correct_to_m0(period_s, period_error_s, m, splitting=None):
    """Correct observed periods to m=0 using nu_obs = nu_0 + m*delta_nu.

    Returns corrected periods, propagated independent period errors, and the
    derivative d(P_corr)/d(delta_nu) used for shared splitting covariance.
    """
    p = np.asarray(period_s, dtype=float)
    e = np.asarray(period_error_s, dtype=float)
    m = np.asarray(m, dtype=float)
    if splitting is None:
        return p.copy(), e.copy(), np.zeros_like(p)
    if splitting.value_uhz < 0 or splitting.error_uhz < 0:
        raise ValueError("Splitting and its error must be non-negative")
    nu_obs_uhz = 1.0e6 / p
    nu0_uhz = nu_obs_uhz - m * splitting.value_uhz
    if np.any(nu0_uhz <= 0):
        raise ValueError("Rotational correction produced a non-positive m=0 frequency")
    pc = 1.0e6 / nu0_uhz
    dp_dpobs = (pc / p) ** 2
    ec = np.abs(dp_dpobs) * e
    dp_dsplit = m * pc**2 / 1.0e6
    return pc, ec, dp_dsplit


def gls_line(k, period_s, covariance):
    """Fit period = intercept + spacing*k with a full known covariance."""
    k = np.asarray(k, dtype=float)
    y = np.asarray(period_s, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    x = np.column_stack((np.ones_like(k), k))
    try:
        cinv_x = np.linalg.solve(cov, x)
        cinv_y = np.linalg.solve(cov, y)
        normal = x.T @ cinv_x
        beta_cov = np.linalg.inv(normal)
    except np.linalg.LinAlgError as exc:
        raise ValueError("The covariance or regression matrix is singular") from exc
    beta = beta_cov @ (x.T @ cinv_y)
    residual = y - x @ beta
    chi2 = float(residual @ np.linalg.solve(cov, residual))
    return beta, beta_cov, residual, chi2


def fit_sequence(
    sequence,
    k,
    period_s,
    period_error_s,
    m=None,
    splitting=None,
    intrinsic_scatter_s=0.0,
    n_samples=20_000,
    seed=1729,
    resample_modes=False,
):
    """Fit one sequence and propagate errors with parametric Monte Carlo.

    Measurement errors are treated as independent Gaussian 1-sigma errors.
    A splitting uncertainty is drawn once per realization, preserving the
    correlation shared by all modes in a sequence. Optional mode resampling
    estimates sensitivity to the finite set of identified modes.
    """
    k, p, e, m = validate_inputs(k, period_s, period_error_s, m)
    if splitting is None and np.any(m != 0):
        raise ValueError("Non-zero m values require a rotational splitting")
    if intrinsic_scatter_s < 0:
        raise ValueError("intrinsic_scatter_s must be non-negative")
    if n_samples < 100:
        raise ValueError("n_samples must be at least 100")
    pc, ec, dp_dsplit = correct_to_m0(p, e, m, splitting)
    variance = ec**2 + intrinsic_scatter_s**2
    covariance = np.diag(variance)
    if splitting is not None and splitting.error_uhz > 0:
        covariance += np.outer(dp_dsplit, dp_dsplit) * splitting.error_uhz**2

    beta, beta_cov, residual, chi2 = gls_line(k, pc, covariance)
    dof = len(k) - 2
    rng = np.random.default_rng(seed)
    samples = np.empty(n_samples, dtype=float)
    accepted = 0
    attempts = 0
    max_attempts = max(10 * n_samples, n_samples + 1000)

    while accepted < n_samples and attempts < max_attempts:
        attempts += 1
        if resample_modes:
            idx = rng.integers(0, len(k), len(k))
            if np.unique(k[idx]).size < 2:
                continue
        else:
            idx = np.arange(len(k))

        p_draw = rng.normal(p, e)
        if intrinsic_scatter_s > 0:
            p_draw += rng.normal(0.0, intrinsic_scatter_s, len(p))
        if np.any(p_draw <= 0):
            continue

        if splitting is None:
            split_draw = None
        else:
            split_value = (
                rng.normal(splitting.value_uhz, splitting.error_uhz)
                if splitting.error_uhz > 0 else splitting.value_uhz
            )
            if split_value < 0:
                continue
            split_draw = Splitting(split_value, 0.0)

        pcd, ecd, _ = correct_to_m0(p_draw, e, m, split_draw)
        var_draw = ecd**2 + intrinsic_scatter_s**2
        try:
            b, _, _, _ = gls_line(k[idx], pcd[idx], np.diag(var_draw[idx]))
        except ValueError:
            continue
        samples[accepted] = b[1]
        accepted += 1

    if accepted < n_samples:
        raise RuntimeError("Could not generate enough valid Monte Carlo samples")

    p16, median, p84 = np.percentile(samples, [16.0, 50.0, 84.0])
    reduced = chi2 / dof if dof > 0 else np.nan
    return FitResult(
        sequence=str(sequence), n_modes=len(k), intercept_s=float(beta[0]),
        spacing_gls_s=float(beta[1]), spacing_gls_error_s=float(np.sqrt(beta_cov[1, 1])),
        spacing_s=float(median), error_minus_s=float(median - p16),
        error_plus_s=float(p84 - median), p16_s=float(p16), p84_s=float(p84),
        chi2=chi2, dof=dof, reduced_chi2=float(reduced),
        corrected_period_s=pc, corrected_error_s=ec, residual_s=residual,
        samples_s=samples,
    )
