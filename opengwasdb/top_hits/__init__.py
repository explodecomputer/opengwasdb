"""The Top-Hit Index format, shared by every Primary Storage Layout.

Owns the on-disk schema -- ``analysis_offsets``, ``variant_index``,
``analysis_index``, ``abs_z``, ``z``, ``se``, ``p_value``, optionally
``imputed``, one Zarr group per threshold tier under ``data.zarr/top_hits/``
-- end to end: ``writer.write()`` builds it from pre-collected candidate
cells, ``reader.TopHitReader``/``reader.counts()`` read it back, and
``validation.validate_group_structure()`` checks it. Dense, Ragged, and
Hybrid all harvest candidate cells from their own physical storage (a band
scan of the dense matrix; a scan of the Ragged CSR) and hand them to
``writer.write()`` -- that harvesting stays layout-specific, in
``opengwasdb.layouts.dense.top_hits``/``opengwasdb.layouts.ragged.top_hits``,
since it is the one part of the pipeline that genuinely differs by layout.
"""
