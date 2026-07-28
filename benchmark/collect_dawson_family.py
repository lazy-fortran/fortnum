#!/usr/bin/env python3
"""Collect a portable Dawson fused-versus-separate selection record."""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path


PRODUCTS = ("jvp", "vjp")
CANDIDATES = ("fused", "separate")
PERF_EVENTS = ("cycles", "instructions", "cache-references", "cache-misses")


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=True)


def cpu_model() -> str:
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def compiler_version() -> str:
    compiler = os.environ.get("FC", "gfortran")
    result = run([compiler, "--version"])
    return result.stdout.splitlines()[0]


def source_revision(repo: Path) -> str:
    return run(["git", "-C", str(repo), "rev-parse", "HEAD"]).stdout.strip()


def measure(executable: Path, candidate: str, product: str) -> dict[str, object]:
    result = run([str(executable), candidate, product])
    samples = [float(line) for line in result.stdout.splitlines() if line.strip()]
    if len(samples) != 31:
        raise RuntimeError(f"expected 31 samples, received {len(samples)}")
    median = statistics.median(samples)
    mad = statistics.median(abs(value - median) for value in samples)

    rss_result = run([str(executable), candidate, product, "--peak-rss"])
    rss_lines = [line for line in rss_result.stdout.splitlines() if line.strip()]
    if len(rss_lines) != 1:
        raise RuntimeError("peak-RSS mode did not emit exactly one value")

    return {
        "product": product,
        "candidate": candidate,
        "median_ns_per_workload": median,
        "mad_ns_per_workload": mad,
        "peak_rss_bytes": int(rss_lines[0]),
        "samples_ns_per_workload": samples,
    }


def perf_counters(executable: Path, candidate: str, product: str) -> dict[str, object]:
    perf = shutil.which("perf")
    if perf is None:
        return {"available": False, "reason": "perf executable not found"}

    with tempfile.NamedTemporaryFile() as output:
        command = [
            perf,
            "stat",
            "-x",
            ",",
            "-e",
            ",".join(PERF_EVENTS),
            "-o",
            output.name,
            str(executable),
            candidate,
            product,
        ]
        result = run(command, check=False)
        if result.returncode != 0:
            reason = result.stderr.strip().splitlines()
            return {
                "available": False,
                "reason": reason[-1] if reason else "perf stat failed",
            }
        output.seek(0)
        counters: dict[str, int] = {}
        for raw_line in output.read().decode().splitlines():
            fields = raw_line.split(",")
            if len(fields) < 3 or fields[0] in ("<not counted>", "<not supported>"):
                continue
            event = fields[2].split(":", 1)[0]
            if event in PERF_EVENTS:
                counters[event] = int(float(fields[0]))
        if len(counters) != len(PERF_EVENTS):
            return {
                "available": False,
                "reason": "one or more requested counters were unavailable",
                **counters,
            }
        return {"available": True, **counters}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    executable = args.executable.resolve()
    repo = Path(__file__).resolve().parent.parent
    run([str(executable), "fused", "jvp", "--validate"])

    measurements = [
        measure(executable, candidate, product)
        for product in PRODUCTS
        for candidate in CANDIDATES
    ]
    counters = {
        f"{product}_{candidate}": perf_counters(executable, candidate, product)
        for product in PRODUCTS
        for candidate in CANDIDATES
    }
    selections: dict[str, object] = {}
    for product in PRODUCTS:
        product_results = [
            item for item in measurements if item["product"] == product
        ]
        winner = min(product_results, key=lambda item: item["median_ns_per_workload"])
        loser = max(product_results, key=lambda item: item["median_ns_per_workload"])
        selections[product] = {
            "selected": winner["candidate"],
            "speedup_over_other": (
                loser["median_ns_per_workload"] / winner["median_ns_per_workload"]
            ),
            "selection_basis": "minimum measured median wall-clock time",
        }

    record = {
        "schema_version": 1,
        "source_revision": source_revision(repo),
        "operation": "Dawson outer generated value/JVP/VJP family",
        "workload": "one value and one contracted product per scalar call",
        "hardware": {
            "cpu": cpu_model(),
            "machine": platform.machine(),
            "runner": os.environ.get("RUNNER_NAME", "local"),
        },
        "toolchain": {
            "compiler": compiler_version(),
            "compiler_flags": "-O3",
        },
        "sampling": {
            "warmup_rounds": 3,
            "samples": 31,
            "calls_per_sample": 10_000_000,
        },
        "validation": {
            "status": "passed",
            "jvp_oracle": "central difference of sin(dawson(x)) + dawson(x)^2",
            "vjp_oracle": "independent scalar adjoint identity",
        },
        "candidates": measurements,
        "perf_stat": counters,
        "selection": selections,
    }
    args.output.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
