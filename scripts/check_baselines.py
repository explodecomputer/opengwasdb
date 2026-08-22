#!/usr/bin/env python3
"""Fail when a change adds ruff or mypy findings.

`CONTRIBUTING.md` records pre-existing counts rather than pretending the tree
is clean, and asks that a change add none. This enforces that against
`.baselines.json`.

Counts, not a diff. A diff is the better check -- it catches a change that
fixes one finding while introducing another -- but it is fragile in CI because
line numbers move. This catches the regression that matters (the count going
up) and is worth having; the diff recipe in CONTRIBUTING.md remains the right
local check before pushing.

A count *below* baseline is reported, not failed, with the number to commit.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

BASELINES = Path(".baselines.json")
TOOLS = {
    "ruff": ["ruff", "check", ".", "--output-format=concise"],
    "mypy": ["mypy", "opengwasdb"],
}
_COUNT = re.compile(r"Found (\d+) error")


def count(tool: str) -> int:
    out = subprocess.run(TOOLS[tool], capture_output=True, text=True).stdout
    m = _COUNT.search(out)
    if m:
        return int(m.group(1))
    # Both tools report success without a "Found N errors" line.
    if "All checks passed" in out or "Success: no issues found" in out:
        return 0
    print(f"{tool}: could not parse output:\n{out[-2000:]}", file=sys.stderr)
    raise SystemExit(2)


def main() -> int:
    baselines = json.loads(BASELINES.read_text())
    failed = False
    for tool in TOOLS:
        expected, actual = baselines[tool], count(tool)
        if actual > expected:
            print(
                f"{tool}: {actual} findings, baseline {expected} "
                f"(+{actual - expected}). Fix them, or if they are genuinely "
                f"pre-existing update {BASELINES} and say why in the commit message.",
                file=sys.stderr,
            )
            failed = True
        elif actual < expected:
            print(
                f"{tool}: {actual} findings, below the baseline of {expected} — "
                f"nice. Set \"{tool}\": {actual} in {BASELINES} to lock it in."
            )
        else:
            print(f"{tool}: {actual} findings, at baseline.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
