#!/usr/bin/env python3
"""Check that every production derivative kernel has an ownership class."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


CLASSES = {
    "fortsym-generated",
    "fortad-generated",
    "hand-written algorithmic",
    "hand-written stable recurrence",
    "implicit solve",
    "frozen trace",
    "generated hybrid boundary",
}
MECHANISMS = {"autodiff", "analytical", "hybrid"}
PRODUCT_MARKERS = ("jvp", "vjp", "adjoint", "sensitivity")
VALUE_ONLY_GENERATED_SUFFIX = "_value_kernel.f90"
WRAPPER_GENERATOR = "tools/codegen/app/gen_enzyme_scalar_wrappers.f90"


def procedure_definitions(path: Path) -> set[str]:
    """Return procedure definitions, excluding declarations in interfaces."""
    definitions: set[str] = set()
    interface_depth = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        lowered = line.lower()
        if re.match(r"^(abstract\s+)?interface(?:\s|$)", lowered):
            interface_depth += 1
            continue
        if re.match(r"^end\s+interface(?:\s|$)", lowered):
            interface_depth = max(0, interface_depth - 1)
            continue
        if interface_depth or not line:
            continue
        tokens = re.findall(r"[A-Za-z][A-Za-z0-9_]*", line)
        lowered_tokens = [token.lower() for token in tokens]
        for keyword in ("subroutine", "function"):
            if keyword not in lowered_tokens:
                continue
            index = lowered_tokens.index(keyword)
            if index + 1 < len(tokens):
                definitions.add(tokens[index + 1].lower())
            break
    return definitions


def discover_source_kernels(root: Path) -> set[tuple[str, str]]:
    kernels: set[tuple[str, str]] = set()
    for path in sorted((root / "src").rglob("*.f90")):
        relative = path.relative_to(root).as_posix()
        definitions = procedure_definitions(path)
        for symbol in definitions:
            if any(marker in symbol for marker in PRODUCT_MARKERS):
                kernels.add((relative, symbol))
        if relative.startswith("src/generated/"):
            if relative.endswith(VALUE_ONLY_GENERATED_SUFFIX):
                continue
            for symbol in definitions:
                if symbol.startswith("fortnum_"):
                    kernels.add((relative, symbol))
    return kernels


def discover_wrapper_artifacts(root: Path) -> set[str]:
    path = root / WRAPPER_GENERATOR
    text = path.read_text(encoding="utf-8")
    artifacts = {
        Path(match).name
        for match in re.findall(r'"([^"]*fortnum_enzyme_[A-Za-z0-9_]+\.f90)"', text)
    }
    scalar_loop = re.search(
        r"do\s+active_inputs\s*=\s*(\d+)\s*,\s*(\d+)", text, re.IGNORECASE
    )
    scalar_template = "fortnum_enzyme_scalar_\"//suffix//\".f90"
    if scalar_loop is None or scalar_template not in text:
        raise ValueError("cannot identify generated scalar-wrapper range")
    first, last = (int(value) for value in scalar_loop.groups())
    artifacts.update(
        f"fortnum_enzyme_scalar_p{count}.f90"
        for count in range(first, last + 1)
    )
    return artifacts


def load_inventory(path: Path) -> tuple[set[tuple[str, str]], set[str], list[str]]:
    source_entries: set[tuple[str, str]] = set()
    wrapper_artifacts: set[str] = set()
    errors: list[str] = []
    seen: set[tuple[str, str]] = set()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected = {"path", "symbols", "classification", "mechanisms", "reason"}
        if set(reader.fieldnames or ()) != expected:
            return set(), set(), [f"invalid columns: expected {sorted(expected)}"]
        for line_number, row in enumerate(reader, start=2):
            source_path = row["path"].strip()
            classification = row["classification"].strip()
            mechanisms = {
                item.strip() for item in row["mechanisms"].split(";") if item.strip()
            }
            symbols = {
                item.strip().lower()
                for item in row["symbols"].split(";")
                if item.strip()
            }
            if classification not in CLASSES:
                errors.append(
                    f"line {line_number}: invalid classification {classification!r}"
                )
            if not mechanisms or not mechanisms <= MECHANISMS:
                errors.append(
                    f"line {line_number}: invalid mechanisms {sorted(mechanisms)}"
                )
            if not row["reason"].strip():
                errors.append(f"line {line_number}: reason is empty")
            for symbol in symbols:
                key = (source_path, symbol)
                if key in seen:
                    errors.append(f"line {line_number}: duplicate {source_path}:{symbol}")
                seen.add(key)
            if classification == "generated hybrid boundary":
                if source_path != WRAPPER_GENERATOR:
                    errors.append(
                        f"line {line_number}: hybrid boundary has unexpected generator"
                    )
                wrapper_artifacts.update(symbols)
            else:
                source_entries.update((source_path, symbol) for symbol in symbols)
    return source_entries, wrapper_artifacts, errors


def format_entries(entries: set[tuple[str, str]]) -> str:
    return "\n".join(f"  {path}:{symbol}" for path, symbol in sorted(entries))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path("docs/design/derivative_kernel_inventory.csv"),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    inventory_path = args.inventory
    if not inventory_path.is_absolute():
        inventory_path = root / inventory_path

    inventory, wrapper_inventory, errors = load_inventory(inventory_path)
    discovered = discover_source_kernels(root)
    missing = discovered - inventory
    stale = inventory - discovered
    if missing:
        errors.append("unclassified derivative kernels:\n" + format_entries(missing))
    if stale:
        errors.append("inventory entries without a kernel:\n" + format_entries(stale))

    discovered_wrappers = discover_wrapper_artifacts(root)
    missing_wrappers = discovered_wrappers - wrapper_inventory
    stale_wrappers = wrapper_inventory - discovered_wrappers
    if missing_wrappers:
        errors.append(
            "unclassified generated hybrid boundaries:\n  "
            + "\n  ".join(sorted(missing_wrappers))
        )
    if stale_wrappers:
        errors.append(
            "inventory hybrid boundaries not emitted:\n  "
            + "\n  ".join(sorted(stale_wrappers))
        )

    if errors:
        print("derivative-kernel inventory check failed:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
