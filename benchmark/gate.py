#!/usr/bin/env python3
"""Benchmark regression gate for fortnum.

Reads a benchmark run emitted by ``bench_main --json`` and compares each
benchmark against a checked-in baseline. A primal regression beyond the
tolerance factor fails the gate (nonzero exit). Derivative-product overheads
(jvp/vjp/grad/hvp relative to primal) are reported but, by default, do not
fail the gate: runner timing noise on those ratios is not yet characterized.

Run JSON schema (one object per benchmark under ``benchmarks``):

    {
      "benchmarks": [
        {
          "name": "version_string_len",
          "reps": 50000000,
          "ns_per_call": 3.39,
          "backend": "primal",          # analytic|implicit|trace|generated|primal
          "deriv_ns_per_call": null,    # number once a derivative kernel exists
          "jvp_primal": null,           # deriv/primal overhead ratios; null=no kernel
          "vjp_primal": null,
          "grad_primal": null,
          "hvp_primal": null
        }
      ]
    }

The baseline (baseline.json) uses the same schema. The gate matches by
``name``; benchmarks present in the run but missing from the baseline are
reported as new and ignored (not a regression). Benchmarks present in the
baseline but absent from the run fail the gate (a benchmark vanished).

Usage:
    bench_main --json | gate.py --baseline baseline.json
    gate.py --baseline baseline.json --run run.json
    gate.py --baseline baseline.json --run run.json --factor 2.0

Refreshing the baseline: see README.md, "Regression gate".
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Primal regression factor: a run slower than baseline * factor fails.
DEFAULT_FACTOR = 2.0

# Derivative-overhead fields tracked per benchmark. These stay informational
# (non-blocking) until runner noise on the ratios is characterized; flip with
# --gate-derivative to make a derivative regression fail too.
DERIV_FIELDS = ("jvp_primal", "vjp_primal", "grad_primal", "hvp_primal")


def load_run(text: str) -> dict:
    data = json.loads(text)
    if "benchmarks" not in data or not isinstance(data["benchmarks"], list):
        raise ValueError("benchmark JSON must contain a 'benchmarks' list")
    by_name: dict[str, dict] = {}
    for entry in data["benchmarks"]:
        name = entry.get("name")
        if not name:
            raise ValueError("benchmark entry missing 'name'")
        by_name[name] = entry
    return by_name


def fnum(value):
    """Return value as float, or None for missing/null entries."""
    if value is None:
        return None
    return float(value)


def gate(baseline: dict, run: dict, factor: float, gate_derivative: bool):
    """Compare run against baseline. Returns (failed, lines)."""
    lines: list[str] = []
    failed = False

    for name, base in sorted(baseline.items()):
        if name not in run:
            lines.append(f"FAIL  {name}: present in baseline, missing from run")
            failed = True
            continue

        cur = run[name]
        base_ns = fnum(base.get("ns_per_call"))
        cur_ns = fnum(cur.get("ns_per_call"))

        if base_ns is None or base_ns <= 0.0:
            lines.append(f"WARN  {name}: baseline ns_per_call missing or <= 0; skipping primal check")
        elif cur_ns is None:
            lines.append(f"FAIL  {name}: run ns_per_call missing")
            failed = True
        else:
            ratio = cur_ns / base_ns
            limit = base_ns * factor
            status = "FAIL" if cur_ns > limit else "ok"
            if status == "FAIL":
                failed = True
            lines.append(
                f"{status:5s} {name}: primal {cur_ns:.3f} ns vs baseline {base_ns:.3f} ns "
                f"(x{ratio:.2f}, limit x{factor:.2f})"
            )

        # Derivative-product overheads: informational unless --gate-derivative.
        for field in DERIV_FIELDS:
            base_r = fnum(base.get(field))
            cur_r = fnum(cur.get(field))
            if base_r is None and cur_r is None:
                continue
            if base_r is None or base_r <= 0.0:
                lines.append(f"INFO  {name}: {field} run={cur_r} (no baseline)")
                continue
            if cur_r is None:
                lines.append(f"INFO  {name}: {field} missing in run (baseline {base_r:.3f})")
                continue
            dratio = cur_r / base_r
            over = cur_r > base_r * factor
            tag = "FAIL" if (over and gate_derivative) else "INFO"
            if tag == "FAIL":
                failed = True
            lines.append(
                f"{tag:5s} {name}: {field} {cur_r:.3f} vs baseline {base_r:.3f} "
                f"(x{dratio:.2f}, limit x{factor:.2f})"
                + ("" if gate_derivative else " [non-blocking]")
            )

    new = sorted(set(run) - set(baseline))
    for name in new:
        lines.append(f"INFO  {name}: new benchmark, not in baseline (ignored)")

    return failed, lines


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="fortnum benchmark regression gate")
    parser.add_argument(
        "--baseline",
        type=Path,
        default=Path(__file__).resolve().parent / "baseline.json",
        help="baseline JSON (default: benchmark/baseline.json)",
    )
    parser.add_argument(
        "--run",
        type=Path,
        default=None,
        help="run JSON (default: read stdin)",
    )
    parser.add_argument(
        "--factor",
        type=float,
        default=DEFAULT_FACTOR,
        help=f"regression factor; slower than baseline*factor fails (default {DEFAULT_FACTOR})",
    )
    parser.add_argument(
        "--gate-derivative",
        action="store_true",
        help="make derivative-overhead regressions blocking (default: informational)",
    )
    args = parser.parse_args(argv)

    baseline = load_run(args.baseline.read_text())
    run_text = args.run.read_text() if args.run else sys.stdin.read()
    run = load_run(run_text)

    failed, lines = gate(baseline, run, args.factor, args.gate_derivative)
    for line in lines:
        print(line)

    if failed:
        print("benchmark gate: FAIL")
        return 1
    print("benchmark gate: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
