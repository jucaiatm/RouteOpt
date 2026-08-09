#!/usr/bin/env python3
"""Build and run the exact-pricing regression without a Gurobi installation."""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys


def run(command: list[str], cwd: pathlib.Path) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    app_dir = pathlib.Path(__file__).resolve().parents[1]
    source_dir = app_dir / "src"
    test_dir = app_dir / "tests"
    stub_dir = test_dir / "gurobi_stub"
    build_dir = app_dir / "build" / "stub-tests"
    build_dir.mkdir(parents=True, exist_ok=True)

    compiler = os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("No C++20 compiler was found; set CXX explicitly")

    executable = build_dir / ("test_exact_pricing.exe" if os.name == "nt" else "test_exact_pricing")
    command = [
        compiler,
        "-std=c++20",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        f"-I{source_dir}",
        f"-I{stub_dir}",
        str(source_dir / "vurpm.cpp"),
        str(stub_dir / "gurobi_stub.cpp"),
        str(test_dir / "test_exact_pricing.cpp"),
        "-o",
        str(executable),
    ]
    run(command, app_dir)
    run([str(executable)], app_dir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f"VURP-M pricing regression failed: {error}", file=sys.stderr)
        raise SystemExit(1)
