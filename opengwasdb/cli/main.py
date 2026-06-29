"""OpenGWASDB command line interface."""

import json
from pathlib import Path

import typer

from opengwasdb.build.observed import build_dense_observed_from_sources
from opengwasdb.query import AssociationResult, query_store
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


@app.command("query-variant")
def query_variant_command(store_path: Path, identifier: str) -> None:
    """Query one variant by canonical ALID or alias across analyses."""

    _emit_results(query_store(store_path).variant(identifier))


@app.command("query-range")
def query_range_command(
    store_path: Path,
    chromosome: str,
    start: int,
    end: int,
    analysis_id: str | None = typer.Option(None),
) -> None:
    """Query a genomic range."""

    _emit_results(query_store(store_path).range(chromosome, start, end, analysis_id=analysis_id))


@app.command("query-analysis")
def query_analysis_command(store_path: Path, analysis_id: str) -> None:
    """Extract all finite associations for one analysis."""

    _emit_results(query_store(store_path).analysis(analysis_id))


@app.command("query-phewas")
def query_phewas_command(store_path: Path, identifier: str) -> None:
    """Extract one variant across all analyses."""

    _emit_results(query_store(store_path).phewas(identifier))


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


def _emit_results(results: list[AssociationResult]) -> None:
    typer.echo(json.dumps([result.to_dict() for result in results], sort_keys=True))
