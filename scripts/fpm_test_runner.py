#!/usr/bin/env python3
"""Run an fpm test with the same per-test arguments used by fo."""

import os
from pathlib import Path
import sys
import tomllib


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: fpm_test_runner.py EXECUTABLE [ARG ...]")
    executable = str(Path(sys.argv[1]).resolve())
    root = Path(__file__).resolve().parents[1]
    with (root / "fpm.toml").open("rb") as manifest:
        test_args = tomllib.load(manifest)["extra"]["fo"]["test-args"]
    args = test_args.get(Path(executable).stem, [])
    os.chdir(root)
    os.execv(executable, [executable, *args, *sys.argv[2:]])


if __name__ == "__main__":
    main()
