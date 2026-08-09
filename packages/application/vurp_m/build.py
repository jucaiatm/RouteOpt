#!/usr/bin/env python3
"""Configure and build the RouteOpt VURP-M BPC application."""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


def run(command: list[str], cwd: pathlib.Path) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-type", default="Release")
    parser.add_argument("--jobs", type=int, default=0)
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()

    app_dir = pathlib.Path(__file__).resolve().parent
    build_dir = app_dir / "build"
    if args.clean and build_dir.exists():
        import shutil
        shutil.rmtree(build_dir)

    run([
        "cmake", "-S", str(app_dir), "-B", str(build_dir),
        f"-DCMAKE_BUILD_TYPE={args.build_type}",
    ], app_dir)
    build_command = ["cmake", "--build", str(build_dir)]
    if args.jobs > 0:
        build_command.extend(["-j", str(args.jobs)])
    run(build_command, app_dir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        print(f"Build failed with exit code {error.returncode}", file=sys.stderr)
        raise SystemExit(error.returncode)
