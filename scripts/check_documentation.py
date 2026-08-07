#!/usr/bin/env python3
"""Validate maintained documentation against repository-owned facts."""

from __future__ import annotations

import csv
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path


LEGACY_TERMS = (
    "analytic_rule",
    "implicit_rule",
    "trace_rule",
    "transparent_ad",
    "primal_only",
)
STALE_PATTERNS = (
    r"\bno derivative code ships\b",
    r"\bderivatives? (?:do not|does not) ship\b",
    r"\breserved for (?:issue )?#\d+\b",
    r"\bissue #\d+\b",
    r"\bstatus: accepted \(issue\b",
    r"\bcurrent requirement that every procedure have exactly one\b",
)
API_PREFIXES = (
    "bessel_",
    "bspline_",
    "dawson",
    "det2",
    "det3",
    "fft_",
    "fixed_point_",
    "fortnum_",
    "gamma_",
    "gauss_",
    "grid_search",
    "hyperg_",
    "integrate",
    "inv2",
    "inv3",
    "jacobian_",
    "lagrange_",
    "layout_",
    "levin_",
    "linear_solve_",
    "lu_",
    "multiroot_",
    "ode_",
    "oracle_",
    "pack_block",
    "rng_",
    "root_",
    "status_",
    "unpack_block",
)
MODULE_EXCLUSIONS = {
    "fortnum_active_vector",
    "fortnum_ad_interfaces",
    "fortnum_build_selection",
    "fortnum_derivative_registry",
    "fortnum_kinds",
    "fortnum_status",
    "fortnum_version",
}


def maintained_markdown(root: Path) -> list[Path]:
    paths = [
        root / "README.md",
        root / "CONTRIBUTING.md",
        root / "ROADMAP.md",
        root / "benchmark/README.md",
    ]
    paths.extend(sorted((root / "docs").rglob("*.md")))
    paths.extend(sorted((root / ".github").rglob("*.md")))
    return paths


def github_slug(heading: str) -> str:
    slug = heading.strip().lower()
    slug = re.sub(r"[^\w\-\s]", "", slug)
    return re.sub(r"[\s]+", "-", slug)


def anchors(path: Path) -> set[str]:
    result: set[str] = set()
    counts: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^#{1,6}\s+(.+?)\s*#*\s*$", line)
        if match is None:
            continue
        base = github_slug(match.group(1))
        count = counts.get(base, 0)
        counts[base] = count + 1
        result.add(base if count == 0 else f"{base}-{count}")
    return result


def check_links(root: Path, paths: list[Path]) -> list[str]:
    errors: list[str] = []
    link_pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    for path in paths:
        if not path.is_file():
            errors.append(f"missing maintained document: {path.relative_to(root)}")
            continue
        text = path.read_text(encoding="utf-8")
        for raw_target in link_pattern.findall(text):
            target = raw_target.strip().strip("<>")
            if re.match(r"^(?:https?|mailto):", target):
                continue
            file_part, separator, fragment = target.partition("#")
            destination = path if not file_part else (path.parent / file_part)
            if not destination.exists():
                errors.append(
                    f"{path.relative_to(root)}: missing link target {file_part!r}"
                )
                continue
            if separator and fragment and destination.suffix.lower() == ".md":
                if fragment not in anchors(destination):
                    errors.append(
                        f"{path.relative_to(root)}: missing anchor "
                        f"{destination.relative_to(root)}#{fragment}"
                    )
    return errors


def check_docs_map(root: Path) -> list[str]:
    map_path = root / "docs/README.md"
    text = map_path.read_text(encoding="utf-8")
    linked = set()
    for target in re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", text):
        file_part = target.partition("#")[0]
        if not file_part or re.match(r"^(?:https?|mailto):", file_part):
            continue
        linked.add((map_path.parent / file_part).resolve())
    expected = {
        path.resolve()
        for path in (root / "docs").rglob("*.md")
        if path != map_path
    }
    missing = expected - linked
    return [
        f"docs/README.md does not index {path.relative_to(root)}"
        for path in sorted(missing)
    ]


def source_text(root: Path) -> str:
    parts = []
    for path in sorted((root / "src").rglob("*.f90")):
        parts.append(path.read_text(encoding="utf-8"))
    header = root / "include/fortnum.h"
    if header.is_file():
        parts.append(header.read_text(encoding="utf-8"))
    return "\n".join(parts).lower()


def check_api(root: Path) -> list[str]:
    errors: list[str] = []
    api_path = root / "docs/api.md"
    api = api_path.read_text(encoding="utf-8").lower()
    sources = source_text(root)

    for path in sorted((root / "src").rglob("*.f90")):
        relative = path.relative_to(root).as_posix()
        if any(part in relative for part in ("src/generated/", "src/bindings/", "src/testing/")):
            continue
        match = re.search(
            r"^\s*module\s+([a-z][a-z0-9_]*)\b",
            path.read_text(encoding="utf-8"),
            re.IGNORECASE | re.MULTILINE,
        )
        if match is None:
            continue
        module = match.group(1).lower()
        if module in MODULE_EXCLUSIONS:
            continue
        if re.search(rf"\b{re.escape(module)}\b", api) is None:
            errors.append(f"docs/api.md does not mention module {module}")

    for span in re.findall(r"`([a-z][a-z0-9_]*)`", api):
        if not span.startswith(API_PREFIXES):
            continue
        if re.search(rf"\b{re.escape(span)}\b", sources) is None:
            errors.append(f"docs/api.md names absent source symbol {span}")
    return errors


def check_terminology(paths: list[Path], root: Path) -> list[str]:
    errors: list[str] = []
    prose_paths = [path for path in paths if path.name != "ROADMAP.md"]
    sha_pattern = re.compile(r"\b[0-9a-f]{40}\b")
    for path in prose_paths:
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        for term in LEGACY_TERMS:
            if term in lowered:
                errors.append(f"{path.relative_to(root)}: legacy term {term}")
        for pattern in STALE_PATTERNS:
            if re.search(pattern, lowered):
                errors.append(
                    f"{path.relative_to(root)}: stale claim matches {pattern!r}"
                )
        if sha_pattern.search(lowered):
            errors.append(
                f"{path.relative_to(root)}: hard-coded revision; link its owner"
            )
    return errors


def check_generated_revision(root: Path) -> list[str]:
    errors: list[str] = []
    lock_path = root / "tools/codegen/fortsym.lock"
    revision = lock_path.read_text(encoding="utf-8").strip()
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        return ["tools/codegen/fortsym.lock is not one full hexadecimal revision"]
    rk54_lock_path = root / "tools/codegen/fortsym-rk54.lock"
    rk54_revision = rk54_lock_path.read_text(encoding="utf-8").strip()
    if re.fullmatch(r"[0-9a-f]{40}", rk54_revision) is None:
        return [
            "tools/codegen/fortsym-rk54.lock is not one full hexadecimal revision"
        ]
    # Only fortsym's own output carries a fortsym revision. src/generated also
    # holds kernels fortad produced, and the inventory is what says which is
    # which, so the classification there is the authority rather than the
    # directory a file happens to sit in.
    fortad = set()
    inventory = root / "docs/design/derivative_kernel_inventory.csv"
    if inventory.exists():
        with inventory.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if row["classification"].strip() == "fortad-generated":
                    fortad.add(row["path"].strip())
    for path in sorted((root / "src/generated").glob("*.f90")):
        if str(path.relative_to(root)) in fortad:
            continue
        text = path.read_text(encoding="utf-8")
        path_revision = (
            rk54_revision if "Generator: gen_rk54_device" in text else revision
        )
        expected = f"Generator revision: fortsym@{path_revision}"
        if expected not in text:
            lock_name = (
                "fortsym-rk54.lock"
                if "Generator: gen_rk54_device" in text
                else "fortsym.lock"
            )
            errors.append(
                f"{path.relative_to(root)} does not match tools/codegen/{lock_name}"
            )
    return errors


def check_report(root: Path) -> list[str]:
    errors: list[str] = []
    csv_path = root / "benchmark/report/data/mechanism_tournaments.csv"
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return ["mechanism tournament CSV is empty"]

    mechanisms = ("analytical", "autodiff", "hybrid", "diagnostic")
    selected = {
        mechanism: sum(row["selected_mechanism"] == mechanism for row in rows)
        for mechanism in mechanisms
    }
    fastest = {
        mechanism: sum(row["fastest_mechanism"] == mechanism for row in rows)
        for mechanism in mechanisms
    }
    ratios = []
    for row in rows:
        first = float(row["fastest_ns"])
        second = float(row["runner_up_ns"])
        if first <= 0.0 or second < first:
            errors.append(f"invalid tournament timing in {row['workload']}")
            continue
        ratios.append(second / first)
        record = root / "benchmark/reference" / row["source_record"]
        if not record.is_file():
            errors.append(f"missing tournament source record {row['source_record']}")
    if errors:
        return errors

    report = (root / "docs/design/differentiation_report.md").read_text(
        encoding="utf-8"
    )
    required = [
        f"contains {len(rows)} derivative",
        f"| `analytical` | {selected['analytical']} | "
        f"{100.0 * selected['analytical'] / len(rows):.1f}% | "
        f"{fastest['analytical']} |",
        f"| `autodiff` | {selected['autodiff']} | "
        f"{100.0 * selected['autodiff'] / len(rows):.1f}% | "
        f"{fastest['autodiff']} |",
        f"| `hybrid` | {selected['hybrid']} | "
        f"{100.0 * selected['hybrid'] / len(rows):.1f}% | "
        f"{fastest['hybrid']} |",
        f"| diagnostic | {selected['diagnostic']} | "
        f"{100.0 * selected['diagnostic'] / len(rows):.1f}% | "
        f"{fastest['diagnostic']} |",
        f"ranges from {min(ratios):.3f} to {max(ratios):,.3f}",
        f"Its median is {statistics.median(ratios):.3f}",
        f"geometric mean is {math.exp(statistics.fmean(map(math.log, ratios))):.3f}",
        f"{sum(value <= 1.2 for value in ratios)} workloads have a",
    ]
    for phrase in required:
        if phrase not in report:
            errors.append(f"differentiation report is missing current value: {phrase}")
    return errors


def check_no_figures(root: Path) -> list[str]:
    errors = []
    for base in (root / "docs", root / "benchmark/report"):
        for path in base.rglob("*.png"):
            if "build" not in path.relative_to(base).parts:
                errors.append(f"generated figure committed at {path.relative_to(root)}")
    return errors


def check_inventory(root: Path) -> list[str]:
    checker = root / "scripts/check_derivative_kernel_inventory.py"
    result = subprocess.run(
        [sys.executable, str(checker), "--root", str(root)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode == 0:
        return []
    detail = result.stderr.strip() or result.stdout.strip()
    return [f"derivative ownership inventory drift:\n{detail}"]


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    paths = maintained_markdown(root)
    errors = []
    errors.extend(check_links(root, paths))
    errors.extend(check_docs_map(root))
    errors.extend(check_api(root))
    errors.extend(check_terminology(paths, root))
    errors.extend(check_generated_revision(root))
    errors.extend(check_report(root))
    errors.extend(check_no_figures(root))
    errors.extend(check_inventory(root))
    if errors:
        print("documentation drift check failed:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
