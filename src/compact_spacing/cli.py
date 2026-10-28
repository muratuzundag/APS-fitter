"""Command-line interface."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

from .core import Splitting, fit_sequence
from .io import read_modes, write_results
from .plotting import plot_result


def parse_mapping(values, kind):
    mapping = {}
    for text in values or []:
        if "=" not in text:
            raise argparse.ArgumentTypeError(f"{kind} must use SEQUENCE=VALUE syntax")
        name, value = text.split("=", 1)
        name = name.strip()
        if kind == "splitting":
            pieces = [x.strip() for x in value.split(",")]
            if len(pieces) not in (1, 2):
                raise argparse.ArgumentTypeError("splitting must be SEQUENCE=VALUE[,ERROR]")
            mapping[name] = Splitting(float(pieces[0]), float(pieces[1]) if len(pieces) == 2 else 0.0)
        else:
            mapping[name] = float(value)
    return mapping


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "sequence"


def build_parser():
    p = argparse.ArgumentParser(
        prog="compact-spacing",
        description="Fit mean compact-star period spacings with propagated uncertainties.")
    p.add_argument("input", help="CSV with sequence,k,period_s,period_err_s[,m]")
    p.add_argument("--output-dir", default="period_spacing_results")
    p.add_argument("--samples", type=int, default=20_000)
    p.add_argument("--seed", type=int, default=1729)
    p.add_argument("--splitting", action="append", default=[], metavar="SEQ=UHZ[,ERR]",
                   help="Known rotational frequency splitting and 1-sigma error; repeat per sequence")
    p.add_argument("--intrinsic-scatter", action="append", default=[], metavar="SEQ=SECONDS",
                   help="Optional independent extra scatter; repeat per sequence")
    p.add_argument("--resample-modes", action="store_true",
                   help="Also bootstrap rows; use for mode-set sensitivity, not measurement error alone")
    p.add_argument("--format", choices=("pdf", "png"), default="pdf")
    p.add_argument("--show", action="store_true")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        split_map = parse_mapping(args.splitting, "splitting")
        scatter_map = parse_mapping(args.intrinsic_scatter, "scatter")
        groups = read_modes(args.input)
    except (ValueError, argparse.ArgumentTypeError) as exc:
        raise SystemExit(f"Input error: {exc}") from exc

    unknown = (set(split_map) | set(scatter_map)) - set(groups)
    if unknown:
        raise SystemExit("Unknown sequence(s) in options: " + ", ".join(sorted(unknown)))

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for offset, (name, data) in enumerate(groups.items()):
        if name in split_map and np.all(data["m"] == 0):
            print(f"Warning: splitting supplied for {name}, but all m values are zero.")
        result = fit_sequence(
            name, **data, splitting=split_map.get(name),
            intrinsic_scatter_s=scatter_map.get(name, 0.0),
            n_samples=args.samples, seed=args.seed + offset,
            resample_modes=args.resample_modes)
        results.append(result)
        plot_result(data["k"], result, out / f"{safe_name(name)}_spacing.{args.format}",
                    show=args.show)
        np.save(out / f"{safe_name(name)}_samples.npy", result.samples_s)
        print(f"{name}: DeltaPi = {result.spacing_s:.4f} "
              f"+{result.error_plus_s:.4f} -{result.error_minus_s:.4f} s; "
              f"chi2/dof = {result.chi2:.2f}/{result.dof}")
    write_results(out / "period_spacing_results.csv", results)
    return 0
