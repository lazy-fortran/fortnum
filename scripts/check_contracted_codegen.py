#!/usr/bin/env python3
"""Enforce direct contracted generation for committed JVP kernels."""

from __future__ import annotations

import re
import sys
from pathlib import Path


CONTRACTED_APIS = re.compile(r"\b(?:jvp|directional_derivative)\s*\(", re.I)
JACOBIAN_CALL = re.compile(r"\bjacobian\s*\(", re.I)
GENERATOR_BANNER = re.compile(r"^!\s*Generator:\s*(\S+)\s*$", re.M)


def code_without_comments(path: Path) -> str:
    return "\n".join(
        line.split("!", 1)[0] for line in path.read_text(encoding="utf-8").splitlines()
    )


def check(root: Path) -> list[str]:
    errors: list[str] = []
    app = root / "tools/codegen/app"
    generators = {path.stem: path for path in app.glob("*.f90")}

    for name, path in sorted(generators.items()):
        if JACOBIAN_CALL.search(code_without_comments(path)):
            errors.append(f"{path.relative_to(root)} materializes a Jacobian")

    for kernel in sorted((root / "src/generated").glob("*_jvp_kernel.f90")):
        text = kernel.read_text(encoding="utf-8")
        match = GENERATOR_BANNER.search(text)
        if match is None:
            errors.append(f"{kernel.relative_to(root)} has no generator banner")
            continue
        generator = generators.get(match.group(1))
        if generator is None:
            errors.append(
                f"{kernel.relative_to(root)} names absent generator {match.group(1)}"
            )
            continue
        source = generator.read_text(encoding="utf-8")
        generator_code = code_without_comments(generator)
        if not CONTRACTED_APIS.search(generator_code):
            if "contracted-jvp-direct" not in source.lower():
                errors.append(
                    f"{generator.relative_to(root)} lacks a contracted JVP "
                    f"construction for {kernel.name}"
                )
    return errors


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    errors = check(root)
    if errors:
        print("contracted-codegen check failed:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
