"""GWAS-VCF SourceReader (issue #19).

Wraps the existing bcftools-based `opengwasdb.build.vcf_source` logic behind
the `SourceReader` interface with no behaviour change: `stream_associations`
is a thin re-shaping of `stream_vcf_associations`'s existing tuples.
`extract_at_sites` is new surface -- a combined AF+SE-at-sites lookup did not
exist as a single function before this reader -- built from the same
bcftools-query-plus-`orient_to_canonical` pattern
`opengwasdb.ancestry.extract.extract_af_at_sites` already established for AF
alone, reusing its regions-file, palindrome-filtering, and AF-parsing helpers
rather than duplicating them.
"""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Iterable, Iterator, KeysView
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.ancestry.extract import _parse_af, is_palindromic, write_regions_file
from opengwasdb.build.vcf_source import _require_bcftools, stream_vcf_associations
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics
from opengwasdb.variants.normalise import VariantNormalisationError, orient_to_canonical

GWAS_VCF_CAPABILITY = "opengwasdb.gwas-vcf"


@dataclass(frozen=True)
class GwasVcfReader:
    """SourceReader for one GWAS-VCF file."""

    path: str | Path

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        for chrom, pos, ref, alt, z, se, scale in stream_vcf_associations(self.path):
            yield ReaderAssociation(
                chromosome=chrom,
                position=pos,
                ref=ref,
                alt=alt,
                z=z,
                se=se,
                stored_effect_scale=scale,
            )

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        wanted = (
            alids if isinstance(alids, (set, frozenset, dict, KeysView)) else set(alids)
        )
        if not wanted:
            return {}
        bcftools = _require_bcftools()
        with tempfile.TemporaryDirectory() as tmp_dir:
            regions_path = write_regions_file(wanted, Path(tmp_dir) / "regions.tsv")
            cmd = [
                bcftools,
                "query",
                "-R",
                str(regions_path),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t[%AF]\t[%SE]\n",
                str(self.path),
            ]
            out: dict[str, SiteMetrics] = {}
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
            try:
                for line in proc.stdout:  # type: ignore[union-attr]
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) != 6:
                        continue
                    chrom, pos, ref, alt, af_str, se_str = parts
                    if "," in alt or is_palindromic(ref, alt):
                        continue
                    af = _parse_af(af_str)
                    se = _parse_positive_float(se_str)
                    if af is None or se is None:
                        continue
                    try:
                        orientation = orient_to_canonical(chrom, pos, alt, ref)
                    except VariantNormalisationError:
                        continue
                    alid = orientation.variant.alid
                    if alid not in wanted:
                        continue
                    out[alid] = SiteMetrics(
                        af=(1.0 - af) if orientation.flipped else af,
                        se=se,
                    )
            finally:
                proc.stdout.close()  # type: ignore[union-attr]
                proc.wait()
        return out


def _parse_positive_float(value: str) -> float | None:
    if value in {".", ""}:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if parsed > 0.0 else None
