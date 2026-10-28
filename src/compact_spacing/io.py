"""CSV input and output helpers."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

REQUIRED = {"sequence", "k", "period_s", "period_err_s"}


def read_modes(path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header")
        missing = REQUIRED - set(reader.fieldnames)
        if missing:
            raise ValueError("Missing CSV columns: " + ", ".join(sorted(missing)))
        rows = list(reader)
    if not rows:
        raise ValueError("Input CSV contains no modes")

    groups = {}
    for line, row in enumerate(rows, start=2):
        name = row["sequence"].strip()
        if not name:
            raise ValueError(f"Empty sequence at CSV line {line}")
        try:
            item = {
                "k": float(row["k"]),
                "period_s": float(row["period_s"]),
                "period_error_s": float(row["period_err_s"]),
                "m": float(row.get("m", "0") or 0),
            }
        except ValueError as exc:
            raise ValueError(f"Non-numeric value at CSV line {line}") from exc
        groups.setdefault(name, []).append(item)

    output = {}
    for name, items in groups.items():
        items.sort(key=lambda x: x["k"])
        output[name] = {key: np.array([r[key] for r in items], dtype=float)
                        for key in ("k", "period_s", "period_error_s", "m")}
    return output


def write_results(path, results):
    fields = ["sequence", "n_modes", "spacing_s", "error_minus_s", "error_plus_s",
              "p16_s", "p84_s", "spacing_gls_s", "spacing_gls_error_s",
              "intercept_s", "chi2", "dof", "reduced_chi2"]
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow({field: getattr(r, field) for field in fields})
