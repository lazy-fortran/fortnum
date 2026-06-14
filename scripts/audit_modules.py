#!/usr/bin/env python3
"""Audit fortnum module public surfaces against docs/api.md.

Checks:
  - every module uses 'private' default (no default-public leakage)
  - each explicitly public name is documented in docs/api.md
  - reports undocumented public names as warnings
  - exits non-zero if any violation is found

Usage:
    python3 scripts/audit_modules.py [--src SRC_DIR] [--api API_MD]
"""

import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Fortran parser helpers
# ---------------------------------------------------------------------------

def _strip_fortran_comments(line: str) -> str:
    """Remove inline Fortran comment (! ...) from a line."""
    in_str = False
    char = None
    for i, c in enumerate(line):
        if in_str:
            if c == char:
                in_str = False
        else:
            if c in ('"', "'"):
                in_str = True
                char = c
            elif c == "!":
                return line[:i]
    return line


def _continuation_lines(path: Path) -> list:
    """Return logical lines: free-form Fortran continuation (&) joined."""
    raw = path.read_text(errors="replace").splitlines()
    logical = []
    buf = ""
    for raw_line in raw:
        stripped = _strip_fortran_comments(raw_line).rstrip()
        if stripped.endswith("&"):
            buf += stripped[:-1] + " "
        else:
            buf += stripped
            logical.append(buf)
            buf = ""
    if buf:
        logical.append(buf)
    return logical


def parse_module(path: Path) -> list:
    """Parse a single .f90 file; return list of info dicts, one per module."""
    lines = _continuation_lines(path)
    modules = []
    current = None

    re_module_start = re.compile(
        r"^\s*module\s+(\w+)\s*$", re.IGNORECASE
    )
    re_module_end = re.compile(
        r"^\s*end\s+module(?:\s+\w+)?\s*$", re.IGNORECASE
    )
    re_private_stmt = re.compile(
        r"^\s*private\s*$", re.IGNORECASE
    )
    re_public_list = re.compile(
        r"^\s*public\s*::\s*(.+)$", re.IGNORECASE
    )
    re_inline_public = re.compile(
        r"^\s*(?:\w[\w\s,()=*:]+),\s*public\s*(?:::|,)", re.IGNORECASE
    )
    re_type_public = re.compile(
        r"^\s*type\s*,\s*public\s*::\s*(\w+)", re.IGNORECASE
    )
    re_abstract_iface = re.compile(
        r"^\s*abstract\s+interface\s*$", re.IGNORECASE
    )
    re_end_iface = re.compile(
        r"^\s*end\s+interface", re.IGNORECASE
    )

    in_iface = 0

    for line in lines:
        if current is None:
            m = re_module_start.match(line)
            if m:
                name = m.group(1).lower()
                current = {
                    "name": name,
                    "file": path,
                    "has_private_default": False,
                    "public_names": [],
                }
                in_iface = 0
            continue

        if re_module_end.match(line):
            modules.append(current)
            current = None
            in_iface = 0
            continue

        if re_abstract_iface.match(line):
            in_iface += 1
            continue
        if re_end_iface.match(line):
            if in_iface > 0:
                in_iface -= 1
            continue

        if re_private_stmt.match(line):
            current["has_private_default"] = True
            continue

        if in_iface:
            continue

        m = re_public_list.match(line)
        if m:
            raw = m.group(1)
            names = [n.strip().lower() for n in raw.split(",") if n.strip()]
            current["public_names"].extend(names)
            continue

        m = re_type_public.match(line)
        if m:
            current["public_names"].append(m.group(1).lower())
            continue

        if re_inline_public.match(line):
            colon_pos = line.find("::", line.lower().find("public"))
            if colon_pos >= 0:
                rest = line[colon_pos + 2:]
                name_part = rest.split("=")[0].strip().lower()
                if re.match(r"^\w+$", name_part):
                    current["public_names"].append(name_part)

    return modules


# ---------------------------------------------------------------------------
# docs/api.md parser
# ---------------------------------------------------------------------------

def parse_api_md(path: Path) -> set:
    """Extract documented names from api.md.

    Pulls names from:
      - inline backtick: `name`
      - Markdown heading: ### `name(...)`
      - Fortran code-fence signature lines:
          subroutine name(
          function name(
          type :: name
          integer, parameter :: NAME
    """
    text = path.read_text(errors="replace")
    documented = set()

    # 1. All backtick-quoted identifiers
    for m in re.finditer(r"`(\w+)`", text):
        documented.add(m.group(1).lower())

    # 2. Fortran code block signatures
    in_code = False
    for line in text.splitlines():
        if line.strip().startswith("```"):
            in_code = not in_code
            continue
        if not in_code:
            # Heading: ### `name(...)` or #### `name(...)`
            m = re.match(r"^#{1,6}\s+`(\w+)", line)
            if m:
                documented.add(m.group(1).lower())
            continue
        # Inside a Fortran code block
        # subroutine / function / pure subroutine / elemental function ...
        m = re.match(
            r"^\s*(?:pure\s+|elemental\s+)?(?:subroutine|function)\s+(\w+)\s*\(",
            line, re.IGNORECASE
        )
        if m:
            documented.add(m.group(1).lower())
            continue
        # abstract interface name
        m = re.match(r"^\s*(?:function|subroutine)\s+(\w+)\s*\(", line, re.IGNORECASE)
        if m:
            documented.add(m.group(1).lower())
            continue
        # type :: name_t
        m = re.match(r"^\s*type\s*::\s*(\w+)", line, re.IGNORECASE)
        if m:
            documented.add(m.group(1).lower())
            continue
        # integer, parameter :: NAME = value
        m = re.match(r"^\s*integer,\s*parameter\s*::\s*(\w+)", line, re.IGNORECASE)
        if m:
            documented.add(m.group(1).lower())
            continue

    return documented


# ---------------------------------------------------------------------------
# Module classification
# ---------------------------------------------------------------------------

# Intentionally internal helper modules: used only by sibling fortnum modules.
PRIVATE_HELPER_MODULES = {
    "fortnum_ode_cash_karp",
    "fortnum_ode_events",
    "fortnum_ode_wrapper",
    "fortnum_special_bessel",
    "fortnum_special_dawson",
    "fortnum_special_gamma",
}

# Documented public API modules.
PUBLIC_API_MODULES = {
    "fortnum_kinds",
    "fortnum_status",
    "fortnum_special",
    "fortnum_special_complex_bessel",
    "fortnum_special_hypergeometric_1f1",
    "fortnum_special_erf_cbind",
    "fortnum_fft",
    "fortnum_quadrature",
    "fortnum_integrate_gk",
    "fortnum_integrate",
    "fortnum_levin",
    "fortnum_ode",
    "fortnum_ode_dop853",
    "fortnum_roots",
    "fortnum_multiroot",
    "fortnum_rng",
    "fortnum_interp",
    "fortnum_polynomial",
    "fortnum_bspline",
    "fortnum_oracle",
}

# C ABI binding modules: every public name is a bind(c) wrapper documented in
# docs/migration_libneo.md (the C-consumer surface), not in the Fortran api.md.
# They keep the bare 'private' default; the audit does not require api.md
# entries for the C-callable symbol names.
CBINDING_MODULES = {
    "fortnum_capi",
    "fortnum_capi_bspline",
}

# Internal utility modules; not expected in api.md.
INTERNAL_MODULES = {
    "fortnum_version",
}


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def audit(src_dir: Path, api_md: Path) -> int:
    """Run the audit; return exit code (0 = clean, 1 = violations)."""
    documented = parse_api_md(api_md)

    all_modules = []
    for f90 in sorted(src_dir.rglob("*.f90")):
        all_modules.extend(parse_module(f90))

    violations = 0
    warnings = 0

    print(f"fortnum module audit  (src={src_dir},  api={api_md})")
    print("=" * 72)

    for mod in all_modules:
        name = mod["name"]
        rel = mod["file"].relative_to(src_dir.parent)

        issues = []
        warns = []

        # 1. Every module must declare 'private' as the default accessibility.
        if not mod["has_private_default"]:
            issues.append(
                f"VIOLATION: module '{name}' missing bare 'private' statement"
                f" — is default-public and may leak helpers"
            )

        is_helper = name in PRIVATE_HELPER_MODULES
        is_internal = name in INTERNAL_MODULES
        is_public_api = name in PUBLIC_API_MODULES
        is_cbinding = name in CBINDING_MODULES

        # 2. Public-API modules: every public name must be in api.md.
        if is_public_api:
            for pname in mod["public_names"]:
                if pname not in documented:
                    warns.append(
                        f"WARN: public name '{pname}' not documented in api.md"
                    )

        # 3. Internal modules: flag exports missing from api.md for reviewers.
        if is_internal:
            for pname in mod["public_names"]:
                if pname not in documented:
                    warns.append(
                        f"WARN: internal module '{name}' exports '{pname}'"
                        f" (not in api.md; intentionally internal)"
                    )

        # Combine and print.
        if issues or warns:
            print(f"\n{rel}  [{name}]")
            for msg in issues:
                print(f"  {msg}")
                violations += 1
            for msg in warns:
                print(f"  {msg}")
                warnings += 1
        else:
            category = (
                "helper"   if is_helper
                else "internal" if is_internal
                else "c-binding" if is_cbinding
                else "public-api"
            )
            pub = ", ".join(mod["public_names"]) if mod["public_names"] else "(none)"
            print(
                f"  OK  [{category:10s}]  {name}  public={pub}"
            )

    print()
    print("=" * 72)
    print(f"Modules checked : {len(all_modules)}")
    print(f"Violations      : {violations}")
    print(f"Warnings        : {warnings}")

    if violations:
        print("\nAudit FAILED: violations must be fixed.")
        return 1
    if warnings:
        print("\nAudit PASSED with warnings (undocumented public names).")
        return 0
    print("\nAudit PASSED — all public surfaces match api.md.")
    return 0


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--src",
        type=Path,
        default=repo_root / "src",
        help="path to src/ directory (default: repo_root/src)",
    )
    parser.add_argument(
        "--api",
        type=Path,
        default=repo_root / "docs" / "api.md",
        help="path to docs/api.md (default: repo_root/docs/api.md)",
    )
    args = parser.parse_args()
    sys.exit(audit(args.src, args.api))


if __name__ == "__main__":
    main()
