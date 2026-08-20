"""The Reference-Completion ancestry-match filter (ADR 0028, issue #108).

`derive_impute_analysis_ids` decides which Analyses an LD panel may impute.
Its no-ancestry-anywhere fallback -- impute everything -- is deliberate
backward compatibility, but it is also the one path where a genuinely
mixed-ancestry store gets imputed against a single panel's LD/EAF with no
signal to the operator, so it must announce itself.
"""

from __future__ import annotations

import logging

from opengwasdb.completion.ancestry_filter import derive_impute_analysis_ids


def _rows(*pairs: tuple[str, str]) -> list[dict[str, str]]:
    return [{"analysis_id": aid, "assigned_ancestry": anc} for aid, anc in pairs]


def test_matching_subset_is_selected_when_ancestry_is_known():
    got = derive_impute_analysis_ids(
        _rows(("a1", "EUR"), ("a2", "EAS"), ("a3", "EUR")), "EUR"
    )
    assert got == {"a1", "a3"}


def test_no_ancestry_anywhere_still_imputes_everything():
    # Backward compatibility (pre-ADR-0028 releases): None means "impute all".
    assert derive_impute_analysis_ids(_rows(("a1", ""), ("a2", "")), "EUR") is None


def test_no_ancestry_anywhere_warns_that_one_panel_is_applied_to_all(caplog):
    """The fallback must be loud: silently imputing every Analysis against one
    panel is indistinguishable, in the logs, from correctly imputing a
    genuinely single-ancestry store (issue #108)."""
    with caplog.at_level(logging.WARNING, logger="opengwasdb.completion.ancestry_filter"):
        derive_impute_analysis_ids(_rows(("a1", ""), ("a2", ""), ("a3", "")), "EUR")

    warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert warnings, "no-ancestry fallback imputed every Analysis without warning"
    message = warnings[0].getMessage()
    assert "3" in message  # how many Analyses are affected
    assert "EUR" in message  # which panel they are all being imputed against


def test_known_ancestry_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING, logger="opengwasdb.completion.ancestry_filter"):
        derive_impute_analysis_ids(_rows(("a1", "EUR"), ("a2", "EAS")), "EUR")
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


def test_ancestry_known_but_none_match_returns_empty_not_none():
    # Distinct from the fallback: ancestry *is* known, nothing matches, so
    # nothing should be imputed -- returning None here would impute everything
    # against a panel every Analysis is known not to match.
    assert derive_impute_analysis_ids(_rows(("a1", "EAS"), ("a2", "AFR")), "EUR") == set()
