"""OpenGWASDB command line interface."""

import json
from pathlib import Path

import numpy as np
import typer

from opengwasdb.build.observed import build_dense_observed_from_sources
from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.layouts.ragged.build_besd import build_ragged_from_besd
from opengwasdb.query import query_store
from opengwasdb.store import open_store
from opengwasdb.validation import validate_store

app = typer.Typer(no_args_is_help=True)


@app.command()
def info(store_path: Path) -> None:
    """Print basic manifest information for a local Store Release."""

    store = open_store(store_path)
    manifest = store.manifest
    typer.echo(f"store_id: {manifest.store_id}")
    typer.echo(f"release_id: {manifest.release_id}")
    typer.echo(f"format_version: {manifest.format_version}")
    typer.echo(f"primary_layout: {manifest.primary_layout.value}")
    typer.echo(f"association_coverage: {manifest.association_coverage.value}")
    typer.echo(f"completion_state: {manifest.completion_state.value}")
    typer.echo(f"reference_assembly: {manifest.reference_assembly}")


@app.command("validate")
def validate_command(store_path: Path) -> None:
    """Validate a local Store Release."""

    result = validate_store(store_path)
    if result.ok:
        typer.echo("valid")
        return
    for error in result.errors:
        typer.echo(f"error: {error}", err=True)
    raise typer.Exit(1)


@app.command("build-dense")
def build_dense_command(
    source_path: Path,
    output_path: Path,
    store_id: str = typer.Option(...),
    release_id: str = typer.Option(...),
    reference_assembly: str = typer.Option("GRCh37"),
    overwrite: bool = typer.Option(False),
) -> None:
    """Build a Dense Observed-Only store from a tiny TSV/CSV source."""

    result = build_dense_observed_from_sources(
        [source_path],
        output_path,
        store_id=store_id,
        release_id=release_id,
        reference_assembly=reference_assembly,
        overwrite=overwrite,
    )
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
            },
            sort_keys=True,
        )
    )


@app.command("build-dense-vcf")
def build_dense_vcf_command(
    manifest_path: Path,
    output_path: Path,
    store_id: str = typer.Option(...),
    release_id: str = typer.Option(...),
    overwrite: bool = typer.Option(False),
) -> None:
    """Build a Dense Observed-Only store from a manifest of GWAS-VCF files.

    MANIFEST_PATH is a TSV with columns: trait_id, file_path, trait_name, n.
    VCF files must be in GRCh37/hg19 coordinates; liftover to hg38 is applied inline.
    """

    result = build_dense_from_vcf_manifest(
        manifest_path,
        output_path,
        store_id=store_id,
        release_id=release_id,
        overwrite=overwrite,
    )
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
            },
            sort_keys=True,
        )
    )


@app.command("build-ragged-besd")
def build_ragged_besd_command(
    besd_prefix: Path,
    output_path: Path,
    store_id: str = typer.Option(...),
    release_id: str = typer.Option(...),
    tissue: str = typer.Option(None),
    overwrite: bool = typer.Option(False),
) -> None:
    """Build a Ragged Observed-Only store from BESD files.

    BESD_PREFIX is the path without extension (.esi, .epi, .besd are appended).
    The BESD dataset must already be in GRCh38; no liftover is applied.
    """
    result = build_ragged_from_besd(
        besd_prefix,
        output_path,
        store_id=store_id,
        release_id=release_id,
        tissue=tissue or None,
        overwrite=overwrite,
    )
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
                "n_associations": result.n_associations,
            },
            sort_keys=True,
        )
    )


@app.command("query-phewas")
def query_phewas_command(store_path: Path, identifier: str) -> None:
    """Extract one variant across all analyses (PheWAS)."""

    _emit_results(query_store(store_path).phewas(identifier))


@app.command("query-range")
def query_range_command(
    store_path: Path,
    chromosome: str,
    start: int,
    end: int,
) -> None:
    """Query a genomic range."""

    _emit_results(query_store(store_path).range(chromosome, start, end))


@app.command("query-analysis")
def query_analysis_command(store_path: Path, analysis_id: str) -> None:
    """Extract all finite associations for one analysis."""

    _emit_results(query_store(store_path).analysis(analysis_id))


@app.command("query-lookup")
def query_lookup_command(store_path: Path, identifiers: str, analysis_ids: str) -> None:
    """Query comma-separated variants against comma-separated analyses."""

    _emit_results(
        query_store(store_path).lookup(
            [item for item in identifiers.split(",") if item],
            [item for item in analysis_ids.split(",") if item],
        )
    )


@app.command("query-top-hits")
def query_top_hits_command(
    store_path: Path,
    threshold: float = typer.Option(5e-8),
    limit: int | None = typer.Option(None),
) -> None:
    """Return ranked top-hit associations."""

    _emit_results(query_store(store_path).top_hits(threshold=threshold, limit=limit))


def _emit_results(result: dict[str, np.ndarray]) -> None:
    rows = [
        {
            "variant_index": int(vi),
            "analysis_index": int(ai),
            "z": float(z),
            "se": float(se),
        }
        for vi, ai, z, se in zip(
            result["variant_index"],
            result["analysis_index"],
            result["z"],
            result["se"],
            strict=True,
        )
    ]
    typer.echo(json.dumps(rows, sort_keys=True))
