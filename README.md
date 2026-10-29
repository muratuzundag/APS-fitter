# APS-fitter

![APS-fitter logo](assets/APS_fitter.png)

**Mean period-spacing fits for pulsating compact stars, with explicit uncertainty propagation.**

APS-fitter fits identified mode sequences using periods, one-standard-deviation period errors, and assigned relative radial orders. It produces a mean spacing with asymmetric Monte Carlo errors, reference-fit statistics, residual plots, and saved slope samples.

Suitable applications include approximately uniform g-mode period sequences in white dwarfs and hot subdwarfs. The software fits an **assigned sequence**: it does not automatically determine radial orders, identify angular degrees, search for an unknown spacing, or resolve ambiguous multiplets.

**Status:** Version `0.1.0`, initial release. Read the statistical limitations before interpreting the formal errors as a complete scientific uncertainty.

## Installation

Requires Python 3.10 or later, NumPy 1.24 or later, and Matplotlib 3.7 or later. NumPy and Matplotlib are installed automatically; SciPy is not required.

Extract the source bundle and open a terminal in the folder containing `pyproject.toml`:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
compact-spacing --help
```

On Windows PowerShell, replace the activation command with:

```powershell
.venv\Scripts\Activate.ps1
```

If the executable is not on your `PATH`, use:

```bash
python -m compact_spacing --help
```

There is no requirement to install a wheel separately when installing this source directory. A PyPI release is not assumed by these instructions.

## Quick start

From the repository root:

```bash
compact-spacing examples/example_modes.csv --output-dir results --samples 5000 --seed 1729
```

This fits `l1` and `l2` independently, with no rotational correction or extra scatter. It saves one three-panel PDF per sequence, a results table, and the Monte Carlo slope samples. To save PNG figures instead:

```bash
compact-spacing examples/example_modes.csv --output-dir results_png --format png
```

Use `--samples 50000` for a larger run. More samples reduce numerical Monte Carlo noise, not the observational uncertainty or uncertainty in the mode identification. Plots are saved without displaying a window unless `--show` is supplied.

## Input CSV

Use a plain comma-separated file with these exact headers:

| Column | Required | Meaning |
|---|---|---|
| `sequence` | Yes | Independently fitted sequence label, e.g. `l1` or `starA_l2` |
| `k` | Yes | Assigned relative radial order, treated as exact |
| `period_s` | Yes | Positive observed period, in seconds |
| `period_err_s` | Yes | Positive, finite 1-sigma period error, in seconds |
| `m` | No | Identified azimuthal order; omitted or blank means zero |

For example, the first sequence in the bundled file is:

```csv
sequence,k,period_s,period_err_s,m
l1,1,440.200,1.5,0
l1,2,457.583,1.0,0
l1,5,513.200,2.4,0
l1,9,585.728,0.5,0
l1,11,618.124,1.0,0
l1,12,638.109,0.1,0
l1,14,679.172,0.5,0
l1,16,714.105,2.6,0
```

The full [example file](examples/example_modes.csv) also contains a second sequence. Its zero `m` values are a no-correction demonstration, not an independent determination of central-component identifications.

- Supply at least three rows and at least two distinct radial orders for each sequence.
- Preserve radial-order gaps; missing modes should not be erased by renumbering detections consecutively.
- A common offset to all radial orders changes the intercept, not the slope.
- Use physically meaningful integer orders even though the reader accepts general numeric values.
- Use separate labels for different stars and different independently fitted sequences.
- A label such as `l1` does not cause the program to determine or validate angular degree.
- Normally use one independent period per radial order and angular degree. Multiplet components are not automatically combined into a centroid.
- Avoid comments, metadata rows, missing required fields, NaN, and infinities.
- Use simple labels; labels that map to the same output filename are rejected.

If converting a frequency and its small Gaussian error from microhertz to a period and period error, use:

```math
P = \frac{10^6}{\nu}, \qquad
\sigma_P \simeq \frac{10^6\sigma_\nu}{\nu^2}.
```

Here the frequency quantities are numerical values in microhertz and periods are in seconds. This is first-order error propagation; large or asymmetric frequency errors need additional care.

## Command-line options

| Option | Default | Meaning |
|---|---|---|
| `INPUT` | Required | Input CSV path |
| `--output-dir DIR` | `period_spacing_results` | Results directory |
| `--samples N` | `20000` | Accepted Monte Carlo draws; minimum 100 |
| `--seed N` | `1729` | Base random seed |
| `--splitting SEQ=VALUE[,ERROR]` | None | Adjacent-component frequency splitting and optional 1-sigma error in microhertz |
| `--intrinsic-scatter SEQ=VALUE` | Zero | Fixed extra scatter parameter in seconds |
| `--resample-modes` | Off | Also resample rows with replacement |
| `--format {pdf,png}` | `pdf` | Figure format |
| `--show` | Off | Also display figures |
| `--version` | — | Print the software version |
| `--help` | — | Print usage information |

Repeat splitting or scatter options for different sequences. Labels must match the CSV exactly. Repeated options for the same sequence are rejected. Existing output files with the same names are overwritten, so use distinct output directories when comparing runs.

## Rotational multiplets

### Correction and units

For identified azimuthal components, the adopted convention is:

```math
\nu_{k\ell m}=\nu_{k\ell0}+m\,\delta\nu_\ell,
\qquad
P_{k\ell0}=\frac{10^6}{10^6/P_{\mathrm{obs}}-m\,\delta\nu_\ell}.
```

Periods are in seconds; the splitting is in **microhertz**, not seconds or cycles per day. Positive `m` means a component at higher frequency under this convention. Translate assignments from a different convention before entering them.

The splitting is the separation for **one step in `m`**, not necessarily the separation between any two observed components and not a rotation frequency to which a Ledoux correction remains to be applied. The correction implements an adopted first-order symmetric splitting description; it is not a model of higher-order rotational effects.

After entering identified `m` values, an example command is:

```bash
compact-spacing modes.csv --output-dir corrected_results --splitting l1=2.0,0.1 --splitting l2=3.0,0.1
```

These splitting numbers illustrate syntax; replace them with measurements appropriate to your target. Supplying `l1=2.0` treats the value as fixed. A nonzero `m` without splitting raises an error; all-zero `m` values produce no correction and trigger a warning if a splitting was supplied.

One splitting value is drawn **per sequence per realization**, preserving its shared effect on all modes in that sequence. The splitting and period measurements are otherwise assumed independent. Covariance arising when they were estimated from the same modes is not automatically included.

Negative splitting draws and invalid periods or corrected frequencies are rejected. Use Gaussian splitting errors only when small compared with the positive splitting; large uncertainties can make rejection alter the effective distribution.

### Unknown identifications

A blank `m` means zero numerically, **not** an unknown identification that the program marginalizes over. The code does not randomly assign `m` and does not validate the allowed `m` range against an angular degree.

For ambiguous assignments, compare explicit hypotheses and describe their effects separately. Formal marginalization requires a model of their probabilities; this is not supplied here.

## Statistical method

### Reference line

Each assigned sequence is modeled as:

```math
P_i=P_0+k_i\,\Delta\Pi+\epsilon_i.
```

The reference slope estimates the sequence's mean spacing. The intercept depends on the adopted zero point of relative radial order and is not automatically an absolute asymptotic phase parameter.

The reference fit uses generalized least squares (GLS):

```math
\widehat{\boldsymbol\beta}=(X^T C^{-1}X)^{-1}X^T C^{-1}\boldsymbol P,
\qquad
\operatorname{Cov}(\widehat{\boldsymbol\beta})=(X^T C^{-1}X)^{-1}.
```

The design matrix contains an intercept and radial order; the coefficients are intercept and slope. Without splitting uncertainty the covariance is diagonal, giving inverse-variance weighted least squares. With splitting uncertainty the reference covariance includes the shared contribution calculated with first-order derivatives.

The supplied errors are treated as absolute errors; the covariance is not automatically multiplied by reduced chi-square. The convention is analogous to absolute-error weighting described in the [SciPy fitting documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html), although APS-fitter implements the linear algebra directly with NumPy.

### Monte Carlo errors

By default, each realization:

1. Perturbs each measured period with its Gaussian 1-sigma observational error.
2. Adds the optional independent Gaussian scatter term.
3. Draws the optional shared splitting and transforms to inferred central periods.
4. Fits the realized sequence conditionally using diagonal weights.
5. Stores the slope.

`spacing_s` is the median of the slope distribution; the two reported errors reach its 16th and 84th percentiles. These are propagated percentile intervals, **not automatically Bayesian credible intervals** or guaranteed-coverage confidence intervals. The draws perturb the measured periods rather than simulating a new data set about the fitted line.

The reference GLS fit includes the shared splitting covariance, but the conditional fit in each realization uses diagonal weights after drawing a splitting. These are different estimators when splitting uncertainty is important; do not describe the code as full-covariance GLS in every realization. A joint model or matched-estimator propagation is preferable when that distinction matters.

### Extra scatter and resampling

A fixed extra scatter parameter can be requested with:

```bash
compact-spacing modes.csv --intrinsic-scatter l1=0.5 --intrinsic-scatter l2=0.5
```

The software neither estimates this parameter nor propagates uncertainty in its value. Without splitting correction it is added in quadrature to measurement errors and contributes an independent Gaussian perturbation. With correction, the perturbation is added before transformation while the same numerical variance is added to weights after transformation; those conventions are only approximately equivalent for corrections close to identity. Review or extend the model before combining substantial corrections with this option.

Mode trapping may create structured residuals, not independent Gaussian scatter. Do not select an extra scatter value merely to obtain preferred error bars.

`--resample-modes` additionally resamples rows with replacement. It can assess sensitivity to the mode set, but combining resampling and measurement perturbations can count overlapping uncertainty contributions; calibrated interval coverage is not established. Small samples and uneven radial-order coverage can cause unstable results. Realizations with fewer than two distinct orders are rejected. See [Andrae (2010)](https://arxiv.org/abs/1009.2755) for background on error propagation and resampling.

## Output files

For `l1` and `l2`, the output directory contains:

```text
results/
├── period_spacing_results.csv
├── l1_spacing.pdf
├── l1_samples.npy
├── l2_spacing.pdf
└── l2_samples.npy
```

The table records:

| Field | Interpretation |
|---|---|
| `sequence`, `n_modes` | Label and input row count |
| `spacing_s` | Median sampled slope, seconds |
| `error_minus_s`, `error_plus_s` | Distances to the 16th and 84th percentiles, seconds |
| `p16_s`, `p84_s` | The percentile endpoints, seconds |
| `spacing_gls_s`, `spacing_gls_error_s` | Reference GLS slope and formal standard error, seconds |
| `intercept_s` | Reference intercept, seconds |
| `chi2`, `dof`, `reduced_chi2` | Reference chi-square, number of rows minus two, and their ratio |

Each figure shows corrected periods and the reference line, residuals against radial order, and the sampled slope distribution. The vertical error bars show **independent measurement errors only**: they omit the shared splitting covariance and optional extra scatter. No fitted-line confidence band is drawn.

Each `.npy` file contains sampled slopes in seconds, not an MCMC chain or intercept samples. Load one with `np.load("results/l1_samples.npy")`.

### Interpretation of the example

For the supplied periods without rotational correction or extra scatter, reference weighted-fit slopes are approximately 18.1569 s for `l1` and 10.8692 s for `l2`. Their reduced chi-square values are approximately 18.15 and 8.64. The Monte Carlo medians are close to these values but vary slightly with the seed and number of draws.

The high reduced chi-square values demonstrate that the specified error budget and uniform-spacing line do not explain all the scatter. They do not uniquely identify its cause. Formal measurement-only errors should not be reported as total astrophysical uncertainties. Do not divide the slope-distribution width by the square root of the number of draws.

## Python API

```python
import numpy as np

from compact_spacing import fit_sequence
from compact_spacing.plotting import plot_result

k = np.array([1, 2, 5, 9, 11, 12, 14, 16])
period_s = np.array([
    440.200, 457.583, 513.200, 585.728,
    618.124, 638.109, 679.172, 714.105,
])
period_error_s = np.array([1.5, 1.0, 2.4, 0.5, 1.0, 0.1, 0.5, 2.6])

result = fit_sequence(
    sequence="l1", k=k, period_s=period_s,
    period_error_s=period_error_s, n_samples=20_000, seed=1729,
)
print(
    f"DeltaPi = {result.spacing_s:.4f} "
    f"+{result.error_plus_s:.4f} -{result.error_minus_s:.4f} s"
)
plot_result(k, result, "l1_spacing.pdf")
```

**Naming difference:** CSV uses `period_err_s`; Python uses `period_error_s`.

To correct identified multiplets, additionally import `Splitting` from `compact_spacing` and pass `m=...` and `splitting=Splitting(value_uhz=..., error_uhz=...)`. The `m` array must match the periods.

The low-level `gls_line(k, period_s, covariance)` accepts a covariance matrix, but the high-level fitter and CLI do not accept arbitrary measurement covariance. Using the low-level function does not automatically propagate that covariance through the sampler.

## Limitations and troubleshooting

Radial orders, sequence membership, and identifications are fixed. Period errors are modeled as independent Gaussian errors before correction. The code does not fit trapping signatures, curvature, structural models, outlier mixtures, or a joint degree-dependent asymptotic relation across sequences. Splitting is fixed in form within each sequence; cross-sequence covariance is not modeled.

| Problem | Action |
|---|---|
| Executable not found | Activate the environment or run `python -m compact_spacing` |
| Input file missing | Check the working directory or use an absolute path |
| Missing columns or invalid fields | Check exact headers, units, finite positive errors, and row formatting |
| Unknown sequence option | Match labels and capitalization exactly |
| Nonpositive corrected frequency | Check splitting units, magnitude, and `m` convention |
| Sampling failure | Review excessive errors, invalid corrections, and degenerate samples |
| Large reduced chi-square | Inspect the residuals and identification/error model rather than automatically inflating errors |

On a headless Linux server:

```bash
MPLBACKEND=Agg compact-spacing examples/example_modes.csv --output-dir results
```

Avoid `--show` if no graphical display is available.

## Reproducibility

Retain the input CSV, full command, software version or commit, environment, and seed:

```bash
python --version
compact-spacing --version
python -m pip freeze > requirements-used.txt
```

The first sequence uses the base seed and subsequent sequences use incremented seeds in their order of first appearance in the CSV. Reordering sequences can therefore change their random streams. A fixed seed is not a guarantee of bit-for-bit reproducibility across NumPy versions or platforms; see the [NumPy compatibility policy](https://numpy.org/doc/2.2/reference/random/compatibility.html).

In publications, specify radial-order and multiplet identification procedures, scatter and splitting assumptions, whether resampling was enabled, and the meaning of the reported percentile intervals.

## Development, citation, and license

Install development tools and run checks:

```bash
python -m pip install -e ".[dev]"
python -m pytest
python -m ruff check .
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for contributions. Tests verify the included behavior, not universal statistical adequacy for arbitrary data.

Cite the software version used; [CITATION.cff](CITATION.cff) contains citation metadata. Use an archived-release DOI if one is later assigned, and do not imply a DOI already exists.

Copyright (c) 2026 Murat Uzundag. Software is distributed under the [MIT License](LICENSE). The generated logo and its provenance are described in [assets/README.md](assets/README.md); no endorsement by an organization associated with a reference illustration is implied.

For initial repository upload instructions, see [UPLOAD_TO_GITHUB.md](UPLOAD_TO_GITHUB.md).
