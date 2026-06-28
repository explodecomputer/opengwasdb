"""OpenGWASDB command line interface."""

from pathlib import Path

import typer

from opengwasdb.store import open_store

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

