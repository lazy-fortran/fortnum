#!/usr/bin/env python3
"""Count floating-point instructions in a compiled artifact, from disassembly.

This is the ``N_machine`` level of the four-count cost chain (``N_sym``,
``N_emit``, ``N_machine``, measured time).  ``N_emit`` is what the generator
wrote; ``N_machine`` is what the compiler produced.  The gap is the
compiler's doing and moves in both directions:

- **down**, when FMA fusion turns two floating-point operations into one
  instruction;
- **up**, when register pressure forces spills, adding loads and stores that
  never appeared in the source.

The measurement is taken from the compiled artifact, never from source text:

- **CPU**:  ``objdump -d`` on the object file
- **CUDA**: ``cuobjdump -sass`` / ``nvdisasm`` for SASS (not PTX; ptxas
  rewrites PTX substantially)

Reported separately and never conflated:

- **instruction count** and **FLOP count**.  One FMA is one instruction and
  two FLOPs.
- **spill loads and stores**, counted apart from ordinary memory traffic
- **register count**, where the tool reports it (``nvdisasm`` does; objdump
  and ``cuobjdump -sass`` do not)

The binary utilities are not present everywhere.  When a tool is absent the
corresponding level is simply not measured: the manifest field is absent,
never ``"skipped"`` and never zero.  Absence is the honest encoding of "no
evidence".

``N_emit`` is parsed from the generated Fortran source when ``--source`` is
given, so the acceptance comparison can be produced by one command.  This is
the count of floating-point operations the generator wrote, and it is kept in
a separate ``n_emit`` block; ``n_machine`` never contains source-derived
numbers.

Usage::

    python3 benchmark/collect_machine_count.py \
        --object /path/to/kernel.o \
        --source src/generated/fortnum_det2_jvp_kernel.f90 \
        --target cpu \
        --output benchmark/reference/<host>_det2_machine_count.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------------
# Floating-point instruction classification.
#
# x86-64 (objdump) uses the SSE/AVX mnemonics; every FP instruction below is
# scalar or packed add/sub/mul/div/sqrt plus the FMA family.  One FMA is one
# instruction and two FLOPs; every other FP instruction is one FLOP.
# ---------------------------------------------------------------------------

X86_FP_PREFIXES = (
    # scalar SSE
    "addsd", "subsd", "mulsd", "divsd", "sqrtsd",
    "addss", "subss", "mulss", "divss", "sqrtss",
    # packed SSE
    "addpd", "subpd", "mulpd", "divpd", "sqrtpd",
    "addps", "subps", "mulps", "divps", "sqrtps",
    # AVX (VEX-encoded) scalar and packed
    "vaddsd", "vsubsd", "vmulsd", "vdivsd", "vsqrtsd",
    "vaddss", "vsubss", "vmulss", "vdivss", "vsqrtss",
    "vaddpd", "vsubpd", "vmulpd", "vdivpd", "vsqrtpd",
    "vaddps", "vsubps", "vmulps", "vdivps", "vsqrtps",
)

X86_FMA_PREFIXES = (
    "vfmadd", "vfnmadd", "vfmsub", "vfnmsub",
    "fmadd", "fnmadd", "fmsub", "fnmsub",
)

# Spill loads/stores on x86-64 go through the stack frame.  They are FP value
# moves that reference the stack pointer (or the frame pointer when a frame
# is materialized).  The integer `mov 0x28(%rsp),%rax` argument-pointer
# reloads are not FP spills and are not counted.
X86_FP_MOVE_PREFIXES = (
    "movsd", "vmovsd", "movss", "vmovss",
    "movapd", "vmovapd", "movaps", "vmovaps",
    "movupd", "vmovupd", "movups", "vmovups",
    "movq", "vmovq", "movdqa", "vmovdqa", "movdqu", "vmovdqu",
)

# SASS (NVIDIA) floating-point instructions.  FFMA is the fused multiply-add:
# one instruction, two FLOPs.  The rest are one FLOP each.
SASS_FP_OPS = {
    "FADD": 1, "FSUB": 1, "FMUL": 1, "FMNMX": 1,
    "FDIV": 1, "FSQRT": 1, "FRND": 1,
    "FADD32I": 1, "FMUL32I": 1,
    "FFMA": 2, "FFMA32I": 2,
    "FMNMX32I": 1,
}

# SASS local-memory loads and stores are the spill traffic: the compiler
# places register spills in local (stack) memory.  LDL/STL are therefore
# spill loads/stores, counted apart from ordinary global memory traffic.
SASS_SPILL_LOADS = ("LDL", "LDLU")
SASS_SPILL_STORES = ("STL", "STLU")

TRANSCENDENTAL_CALLS = (
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
    "exp", "log", "log10", "sqrt", "abs",
)

INSTRUCTION_RE = re.compile(
    r"^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2,}\s+)+([A-Za-z0-9.]+)", re.M
)


# ---------------------------------------------------------------------------
# N_emit: the floating-point operations the generator wrote, parsed from the
# emitted Fortran source.  This is "what the generator wrote"; it is kept
# strictly separate from the disassembly-derived N_machine block.
# ---------------------------------------------------------------------------

def count_emit_operations(source_text: str) -> dict[str, int]:
    """Count FP operations in emitted Fortran assignment expressions.

    Each ``+``, ``-``, ``*``, ``/`` is one instruction and one FLOP;
    ``a**k`` expands to ``k - 1`` multiplications; each transcendental call
    is one instruction (counted in a separate ``transcendental`` bucket and
    not as a FLOP, because the generator cannot know the cost the library
    call has at machine level).
    """
    body = source_text.split("subroutine", 1)[-1]
    # Strip directives, declarations, and comments so only assignment RHS
    # arithmetic is counted.
    lines = []
    for raw in body.splitlines():
        line = raw.split("!")[0].strip()
        if not line:
            continue
        if line.startswith("$"):
            continue
        if re.match(r"^(use|implicit|real|integer|intent|type|contains|end|" +
                    r"subroutine|function|public|private|module|interface)", line):
            continue
        lines.append(line)

    text = " ".join(lines)
    # Remove line-continuation ampersands and Fortran kind suffixes.
    text = re.sub(r"&\s*", " ", text)
    text = re.sub(r"_\w+", "", text)

    # A leading unary minus on an assignment right-hand side is negation,
    # not a binary subtraction; strip exactly the leading `-` after each `=`.
    text = re.sub(r"=\s*-", "=", text)

    adds = text.count("+")
    subs = text.count("-")
    # Each `**` contains two `*` characters, so remove them from the single
    # multiplication count and then add the exponentiation cost separately.
    muls = text.count("*") - 2 * text.count("**")
    divs = text.count("/")
    transcendental = 0
    for name in TRANSCENDENTAL_CALLS:
        transcendental += len(re.findall(rf"\b{name}\s*\(", text))

    # exponentiation: a**k is k-1 multiplications
    for match in re.finditer(r"\*\*(\d+)", text):
        muls += int(match.group(1)) - 1

    return {
        "add": max(0, adds),
        "sub": max(0, subs),
        "mul": max(0, muls),
        "div": max(0, divs),
        "transcendental": transcendental,
    }


def emit_block(breakdown: dict[str, int]) -> dict[str, int]:
    instructions = sum(breakdown.values())
    flops = (
        breakdown["add"] + breakdown["sub"] + breakdown["mul"] +
        breakdown["div"]
    )
    return {
        "instruction_count": instructions,
        "flop_count": flops,
        "breakdown": breakdown,
    }


# ---------------------------------------------------------------------------
# N_machine: disassembly-derived counts.
# ---------------------------------------------------------------------------

def count_cpu_object(objdump_path: Path, object_path: Path) -> dict[str, int]:
    """Count FP instructions in an x86-64 object with ``objdump -d``."""
    result = subprocess.run(
        [str(objdump_path), "-d", str(object_path)],
        check=True, text=True, capture_output=True,
    )
    disassembly = result.stdout

    fp_instructions = 0
    flops = 0
    fma_count = 0
    spill_loads = 0
    spill_stores = 0

    for match in INSTRUCTION_RE.finditer(disassembly):
        opcode = match.group(1)
        if any(opcode.startswith(p) for p in X86_FMA_PREFIXES):
            fma_count += 1
            fp_instructions += 1
            flops += 2
            continue
        if any(opcode.startswith(p) for p in X86_FP_PREFIXES):
            fp_instructions += 1
            flops += 1
            continue
        if any(opcode.startswith(p) for p in X86_FP_MOVE_PREFIXES):
            line = disassembly[match.start():disassembly.find("\n", match.start())]
            if re.search(r"%(?:rsp|rbp)", line):
                # objdump uses AT&T syntax: `opcode src, dst`.  The final
                # comma-separated operand is the destination, so a stack
                # destination is a spill store and a register destination fed
                # from the stack is a spill load.
                destination = line.split(",")[-1]
                if re.search(r"%(?:rsp|rbp)", destination):
                    spill_stores += 1
                else:
                    spill_loads += 1

    return {
        "instruction_count": fp_instructions,
        "flop_count": flops,
        "fma_count": fma_count,
        "spill_loads": spill_loads,
        "spill_stores": spill_stores,
    }


def count_sass_tool(tool_path: Path, object_path: Path) -> str:
    """Run the SASS disassembler and return its text output."""
    name = tool_path.name
    if "cuobjdump" in name:
        return subprocess.run(
            [str(tool_path), "-sass", str(object_path)],
            check=True, text=True, capture_output=True,
        ).stdout
    # nvdisasm
    return subprocess.run(
        [str(tool_path), "-c", str(object_path)],
        check=True, text=True, capture_output=True,
    ).stdout


def count_sass(text: str, tool_name: str) -> dict[str, int]:
    """Count FP instructions and spills in SASS output."""
    fp_instructions = 0
    flops = 0
    fma_count = 0
    spill_loads = 0
    spill_stores = 0

    for line in text.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        # SASS lines look like "        /*0000*/ FFMA R0, R1, R2, R3 ;"
        opcode = None
        for token in tokens:
            if token.startswith("/*"):
                continue
            candidate = token.rstrip(";")
            if candidate.isupper() and candidate in SASS_FP_OPS:
                opcode = candidate
                break
            if candidate in ("LDL", "LDLU", "STL", "STLU"):
                opcode = candidate
                break
        if opcode is None:
            continue
        if opcode in SASS_FP_OPS:
            fp_instructions += 1
            flops += SASS_FP_OPS[opcode]
            if opcode.startswith("FFMA"):
                fma_count += 1
        elif opcode in SASS_SPILL_LOADS:
            spill_loads += 1
        elif opcode in SASS_SPILL_STORES:
            spill_stores += 1

    result = {
        "instruction_count": fp_instructions,
        "flop_count": flops,
        "fma_count": fma_count,
        "spill_loads": spill_loads,
        "spill_stores": spill_stores,
    }

    # Register count: nvdisasm reports the register usage; cuobjdump -sass
    # does not.  Only include it when the tool reports it.
    if "nvdisasm" in tool_name:
        registers = _parse_sass_register_count(text)
        if registers is not None:
            result["register_count"] = registers
    return result


def _parse_sass_register_count(text: str) -> int | None:
    """nvdisasm reports register usage as e.g. ``REG:24`` or a header line."""
    match = re.search(r"REG:\s*(\d+)", text)
    if match:
        return int(match.group(1))
    match = re.search(r"\.reg\s*\.b32\s+%r<(\d+)>", text)
    if match:
        return int(match.group(1))
    return None


# ---------------------------------------------------------------------------
# Top-level measurement and manifest construction.
# ---------------------------------------------------------------------------

def source_revision(repo: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True, text=True, capture_output=True,
    ).stdout.strip()


def compiler_version() -> str:
    compiler = os.environ.get("FC", "gfortran")
    result = subprocess.run(
        [compiler, "--version"], check=False, text=True, capture_output=True
    )
    if result.returncode != 0:
        return compiler
    return result.stdout.splitlines()[0]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count floating-point instructions in a compiled artifact"
    )
    parser.add_argument("--object", type=Path, required=True,
                        help="compiled object or binary to disassemble")
    parser.add_argument("--source", type=Path, default=None,
                        help="generated Fortran source for the N_emit count")
    parser.add_argument("--target", choices=("cpu", "cuda", "auto"),
                        default="auto",
                        help="disassembly target (default: auto-detect)")
    parser.add_argument("--output", type=Path, default=None,
                        help="JSON output path (default: print to stdout)")
    args = parser.parse_args()

    object_path = args.object.resolve()
    repo = Path(__file__).resolve().parent.parent

    target = args.target
    if target == "auto":
        target = "cuda" if object_path.suffix in (".cubin", ".cuda") else "cpu"

    # Pick the tool.  Absence of the tool means the level is not measured:
    # the n_machine block is omitted from the manifest, never "skipped" and
    # never zeroed.
    n_machine: dict[str, object] | None = None
    toolchain: dict[str, str] = {}
    if target == "cpu":
        objdump = shutil.which("objdump")
        if objdump is not None:
            n_machine = count_cpu_object(Path(objdump), object_path)
            toolchain["disassembler"] = subprocess.run(
                [objdump, "--version"], check=True, text=True,
                capture_output=True,
            ).stdout.splitlines()[0]
            toolchain["disassembly"] = "objdump -d on the object file"
    elif target == "cuda":
        tool = None
        for candidate in ("cuobjdump", "nvdisasm"):
            found = shutil.which(candidate)
            if found is not None:
                tool = Path(found)
                break
        if tool is not None:
            text = count_sass_tool(tool, object_path)
            n_machine = count_sass(text, tool.name)
            toolchain["disassembler"] = tool.name
            toolchain["disassembly"] = (
                "SASS from cuobjdump/nvdisasm; PTX is not machine code and "
                "ptxas rewrites it"
            )

    # N_emit: "what the generator wrote", parsed from the emitted source.
    n_emit: dict[str, object] | None = None
    if args.source is not None:
        source_text = args.source.read_text(encoding="utf-8")
        breakdown = count_emit_operations(source_text)
        n_emit = emit_block(breakdown)

    record: dict[str, object] = {
        "schema_version": 1,
        "source_revision": source_revision(repo),
        "operation": "machine-code floating-point instruction count "
                     "(N_machine) from disassembly",
        "target": target,
        "object": str(object_path),
        "toolchain": {
            "compiler": compiler_version(),
            **toolchain,
        },
    }
    if n_emit is not None:
        record["n_emit"] = n_emit
    if n_machine is not None:
        record["n_machine"] = n_machine
        if n_emit is not None:
            record["fusion"] = {
                "emitted_instructions": n_emit["instruction_count"],
                "machine_instructions": n_machine["instruction_count"],
                "instruction_reduction": (
                    n_emit["instruction_count"] -
                    n_machine["instruction_count"]
                ),
                "emitted_flops": n_emit["flop_count"],
                "machine_flops": n_machine["flop_count"],
                "flop_count_preserved": (
                    n_emit["flop_count"] == n_machine["flop_count"]
                ),
            }

    text = json.dumps(record, indent=2) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
