#!/usr/bin/env python3
"""Fail when a change to `opengwasdb/` leaves CHANGELOG.md's Unreleased section untouched.

`CONTRIBUTING.md` says every user-visible change is recorded under
`Unreleased` as it lands, so that cutting a version is a matter of dating a
section rather than reconstructing two months of history from the log. That is
a convention until something checks it.

Checking merely that `CHANGELOG.md` was *modified* would pass on a typo fix in
an old entry, so this compares the extracted `Unreleased` section before and
after: the section itself has to differ.

Runs locally as well as in CI:

    pixi run -e dev python scripts/check_changelog.py --base-ref origin/dev
"""

from __future__ import annotations

import argparse
import subprocess
import sys

CHANGELOG = "CHANGELOG.md"
WATCHED = ("opengwasdb/",)
SKIP_LABEL = "no-changelog"


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], capture_output=True, text=True, check=True
    ).stdout


def changed_files(base_ref: str) -> list[str]:
    merge_base = _git("merge-base", base_ref, "HEAD").strip()
    return [p for p in _git("diff", "--name-only", merge_base, "HEAD").splitlines() if p]


def unreleased_section(text: str) -> str:
    """The `## [Unreleased]` block, up to the next `## [` heading.

    Returns "" when the heading is absent, which makes a CHANGELOG that has
    lost its Unreleased section compare unequal to one that has it -- the
    right answer, since that is itself a problem.
    """
    lines = text.splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("## [Unreleased]"))
    except StopIteration:
        return ""
    rest = lines[start + 1 :]
    end = next((i for i, ln in enumerate(rest) if ln.startswith("## [")), len(rest))
    return "\n".join(rest[:end]).strip()


def file_at(ref: str, path: str) -> str:
    try:
        return _git("show", f"{ref}:{path}")
    except subprocess.CalledProcessError:
        return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-ref", required=True, help="branch this change targets, e.g. origin/dev")
    args = ap.parse_args()

    changed = changed_files(args.base_ref)
    watched = [p for p in changed if p.startswith(WATCHED)]
    if not watched:
        print(f"No changes under {'/'.join(WATCHED)} — changelog entry not required.")
        return 0

    merge_base = _git("merge-base", args.base_ref, "HEAD").strip()
    before = unreleased_section(file_at(merge_base, CHANGELOG))
    after = unreleased_section(open(CHANGELOG, encoding="utf-8").read())

    if before != after:
        print(f"CHANGELOG.md Unreleased section updated alongside {len(watched)} source file(s).")
        return 0

    print(
        f"{CHANGELOG}: the Unreleased section is unchanged, but this branch changes "
        f"{len(watched)} file(s) under {'/'.join(WATCHED)}:\n"
        + "".join(f"    {p}\n" for p in watched[:10])
        + (f"    ... and {len(watched) - 10} more\n" if len(watched) > 10 else "")
        + "\nAdd an entry under '## [Unreleased]' describing what a user would notice.\n"
        f"If this change genuinely has no user-visible effect, apply the "
        f"'{SKIP_LABEL}' label to the pull request.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
