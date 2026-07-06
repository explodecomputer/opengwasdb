"""OpenGWASDB command line interface."""

import json
import logging
from pathlib import Path

import numpy as np
import typer

from opengwasdb.build.observed import build_dense_observed_from_sources
from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.layouts.dense.complete import (
    complete_dense_store,
    resume_dense_completion,
)
from opengwasdb.layouts.dense.constants import DEFAULT_CHUNK_SHAPE
from opengwasdb.layouts.dense.top_hits import build_top_hit_indexes
from opengwasdb.layouts.ragged.build_besd import build_ragged_from_besd
from opengwasdb.layouts.ragged.complete import complete_ragged_store
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.query import query_store
from opengwasdb.store import open_store
from opengwasdb.validation import validate_store

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

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
    n_workers: int = typer.Option(1, help="Fork-based process pool size for Pass 1 and Pass 2"),
    chunk_variants: int = typer.Option(
        DEFAULT_CHUNK_SHAPE[0], help="Zarr chunk size along the variant axis"
    ),
    chunk_analyses: int = typer.Option(
        DEFAULT_CHUNK_SHAPE[1], help="Zarr chunk size along the analysis (trait) axis"
    ),
) -> None:
    """Build a Dense Observed-Only store from a manifest of GWAS-VCF files.

    MANIFEST_PATH is a TSV with columns: trait_id, file_path, trait_name, n.
    VCF files must be in GRCh37/hg19 coordinates; liftover to hg38 is applied inline.

    The zarr chunk shape (default 1000x1000) can be tuned with --chunk-variants /
    --chunk-analyses; a narrower analysis chunk speeds up per-analysis (bulk) reads
    at the cost of larger per-variant (phewas) reads.
    """

    result = build_dense_from_vcf_manifest(
        manifest_path,
        output_path,
        store_id=store_id,
        release_id=release_id,
        overwrite=overwrite,
        n_workers=n_workers,
        chunk_shape=(chunk_variants, chunk_analyses),
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
    source_build: str = typer.Option("hg38"),
    overwrite: bool = typer.Option(False),
) -> None:
    """Build a Ragged Observed-Only store from BESD files.

    BESD_PREFIX is the path without extension (.esi, .epi, .besd are appended).
    Use --source-build hg19 to liftover coordinates to hg38 inline.
    """
    result = build_ragged_from_besd(
        besd_prefix,
        output_path,
        store_id=store_id,
        release_id=release_id,
        tissue=tissue or None,
        source_build=source_build,
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


@app.command("complete-ragged")
def complete_ragged_command(
    source_path: Path,
    dest_path: Path,
    ld_panel: Path = typer.Option(..., help="Root of LD panel (ld_dir/{ancestry}/{chr}/...)"),
    ancestry: str = typer.Option("EUR"),
    cis_window_bp: int = typer.Option(1_000_000),
    min_cor: float = typer.Option(0.7),
    release_id: str = typer.Option(None),
    overwrite: bool = typer.Option(False),
) -> None:
    """Produce a Reference-Completed ragged store from an observed-only store."""
    import time
    t0 = time.time()
    result = complete_ragged_store(
        source_path,
        dest_path,
        ld_panel,
        ancestry=ancestry,
        cis_window_bp=cis_window_bp,
        min_cor=min_cor,
        release_id=release_id or None,
        overwrite=overwrite,
    )
    elapsed = time.time() - t0
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
                "n_associations": result.n_associations,
                "n_imputed": result.n_imputed,
                "n_missing": result.n_missing,
                "elapsed_s": round(elapsed, 1),
            },
            sort_keys=True,
        )
    )


@app.command("complete-dense")
def complete_dense_command(
    source_path: Path,
    dest_path: Path,
    ld_panel: Path = typer.Option(..., help="Root of LD panel (ld_dir/{ancestry}/{chr}/...)"),
    ancestry: str = typer.Option("EUR"),
    min_cor: float = typer.Option(0.7),
    thresh: float = typer.Option(0.9),
    release_id: str = typer.Option(None),
    n_workers: int = typer.Option(1, help="LD-block process-pool size"),
    overwrite: bool = typer.Option(False),
) -> None:
    """Produce a Reference-Completed Dense store from a Full Coverage observed-only store."""
    import time
    t0 = time.time()
    result = complete_dense_store(
        source_path,
        dest_path,
        ld_panel,
        ancestry=ancestry,
        min_cor=min_cor,
        thresh=thresh,
        release_id=release_id or None,
        n_workers=n_workers,
        overwrite=overwrite,
    )
    elapsed = time.time() - t0
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
                "n_imputed": result.n_imputed,
                "n_missing_off_panel": result.n_missing_off_panel,
                "n_missing_imputation_failed": result.n_missing_imputation_failed,
                "elapsed_s": round(elapsed, 1),
            },
            sort_keys=True,
        )
    )


@app.command("complete-dense-resume")
def complete_dense_resume_command(
    checkpoint_dir: Path,
    n_workers: int = typer.Option(1, help="LD-block process-pool size"),
) -> None:
    """Resume an interrupted complete-dense run from its checkpoint directory.

    Takes only the checkpoint directory path; all other build parameters are
    loaded from the build_params.json written by the original run.
    """
    import time
    t0 = time.time()
    result = resume_dense_completion(checkpoint_dir, n_workers=n_workers)
    elapsed = time.time() - t0
    typer.echo(
        json.dumps(
            {
                "output_path": str(result.output_path),
                "n_variants": result.n_variants,
                "n_analyses": result.n_analyses,
                "n_imputed": result.n_imputed,
                "n_missing_off_panel": result.n_missing_off_panel,
                "n_missing_imputation_failed": result.n_missing_imputation_failed,
                "elapsed_s": round(elapsed, 1),
            },
            sort_keys=True,
        )
    )


@app.command("build-ragged-top-hits")
def build_ragged_top_hits_command(store_path: Path) -> None:
    """Build (or rebuild) the top-hit index for a Ragged store."""
    build_ragged_top_hit_indexes(store_path)
    typer.echo("done")


@app.command("build-dense-top-hits")
def build_dense_top_hits_command(store_path: Path) -> None:
    """Build (or rebuild) the top-hit index for a Dense store."""
    build_top_hit_indexes(store_path)
    typer.echo("done")


@app.command("query-phewas")
def query_phewas_command(store_path: Path, identifier: str) -> None:
    """Extract one variant across all analyses (PheWAS)."""

    _emit_results(query_store(store_path).phewas(identifier))


@app.command("query-range-phewas")
def query_range_phewas_command(
    store_path: Path,
    chromosome: str,
    start: int,
    end: int,
) -> None:
    """Regional PheWAS: all variants in a genomic range across all analyses."""

    _emit_results(query_store(store_path).range_phewas(chromosome, start, end))


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
