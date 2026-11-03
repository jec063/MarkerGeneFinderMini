# Planned v3 public contract — frozen in v2.9

This is the supported contract to retain when v3 removes `src.*`, `_core`'s
legacy facade, and the flat CLI. v2.9 retains these compatibility interfaces.
No new statistical methods or default changes are part of this transition.

The reviewed [machine-readable snapshot](public-api-v3.json) freezes root
exports, callable parameter names/kinds/defaults, service entry points,
configuration defaults, capability defaults, result column order, schema
versions, and CLI command argument destinations. A regression test compares
the implementation to this snapshot; changes require intentional review.

## Python

- `analyze(expression, config=None, *, registry=...) -> MarkerResult` and
  `analyze_file(input_path, output_path, config=None, *, registry=...,
  hash_input_file=False) -> MarkerResult` are the primary analysis APIs.
- `AnalysisConfig` is immutable. JSON lists become tuples for cluster/gene
  selectors. `to_dict`, `to_json`, `from_dict`, and `from_json` round-trip the
  schema-2 document. Defaults are frozen in the snapshot, including Native,
  top 5, pseudocount 0.1, one worker, and 1,000-gene chunks.
- `MarkerResult` contains `markers` (DataFrame), `config`, and `metadata`.
  `to_csv(path)` writes without an index; `to_bundle(path=None)` returns ZIP
  bytes and optionally writes them. Metadata includes app version, row count,
  and versioned provenance: requested/resolved engine, engine version,
  effective config, expression source, and scientific dependency versions.
  Single-input hashing remains opt-in. Runtime durations are not deterministic.
- Public services remain available in `analysis`, `inspection`, `batch`,
  `benchmark`, `expression_summary`, `filter_presets`, `profiles`, `settings`,
  `report_bundle`, and `progress`. The migration guide lists legacy replacements.
  Stable constants include the existing service `__all__` exports,
  `FILTER_DEFAULTS`, `FILTER_PRESETS`, `PROFILE_DEFAULTS`, `PROFILE_STATE_KEYS`,
  settings schema constants, and `ProgressCallback`.
- `MarkerFinderError` is the typed API error base; `ConfigurationError`,
  `InputValidationError`, and `SchemaError` also subclass `ValueError`.
  Lower-level service/readers retain their documented `ValueError` behavior;
  filesystem failures can raise `OSError`. Private modules/helpers are not API.

Built-in result columns and order are frozen in the snapshot. Markers keep
cluster selection order and deterministic ranking/tie handling. Native and
Scanpy must each match their own regression fixtures, not one another. Floating
regression tolerances are `rtol=1e-6`, `atol=1e-10`; counts, labels, and ordering
must match exactly. Input expression-source selection never normalizes data.

## Engines

An engine exposes `name` and implements `analyze_dataframe(expression, config)`
and `analyze_file(input_path, output_path, config)`, returning `EngineOutput`.
Register instances with `EngineRegistry`; discovery is opt-in through
`load_engine_plugins` and the `marker_gene_finder.engines` entry-point group.
Names start with a lowercase letter and contain lowercase letters, digits,
underscores, or hyphens. Discovery does not mutate the default registry.

`EngineCapabilities` declares `input_formats` (csv/h5ad/dataframe),
`supports_layers`, `supports_raw`, and `supports_workers`. Unsupported requests
fail before engine invocation. To preserve v2 plugin behavior, omitted
capabilities retain the permissive defaults in the snapshot; explicit
declarations are recommended. Optional string `version` is recorded in
provenance; unknown versions remain null.

EngineOutput contains a DataFrame and finite JSON-compatible diagnostics.
Columns must be unique and include non-null `cluster` and `gene`. Recognized
numeric columns must be finite; probabilities/fractions are in [0,1], rank is a
positive integer, and counts are non-negative integers. Plugins may supply
additional columns. The typed API owns authoritative app-version/row-count
metadata, validates output, and publishes returned marker rows atomically.
The plugin receives a staging output path; it must not assume the user's final
filename or create sidecar artifacts beside that staging path.

## CLI and file behavior

`marker-gene-finder` and `python -m marker_gene_finder` share the same CLI:
`analyze`, `inspect`, `plan`, `batch`, `batch-plan`, `verify`, `benchmark`,
`compare-benchmarks`, and `list-engines`. Analysis settings use `--config`.
Exit 0 means success; exit 1 denotes a failed workflow result (for example,
unexecutable plan, batch failure, invalid verification, or benchmark regression);
exit 2 denotes argument/config/input errors. Commands write their reports
before exit 1. The command-specific flags are frozen in the snapshot.

Settings schema 1 remains readable via migration to settings schema 2; config
schema 2, manifest schema 2, and benchmark schema 1 remain supported. Profiles
and result bundles retain their readers. A major application version does not
automatically change data schemas. Unsupported artifacts receive explicit
schema errors; the migration guide describes recovery.

Writes publish complete individual files by atomic replacement where supported
by the filesystem. This is not a multi-file transaction or power-loss durability
guarantee. Progress callbacks keep `(completed, total, cluster)`, one-based
counts, cluster order, caller-thread execution, and propagated callback errors.
