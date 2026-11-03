# Migrating to the supported API and CLI

v2.9 keeps legacy interfaces available while preparing their removal in v3.
Change imports and commands now; statistical defaults do not change. The
[planned v3 contract](public-api-v3.md) defines the supported surface.

## Python imports

| Legacy module | Supported replacement |
| --- | --- |
| `src.marker_finder` | Services listed below; prefer `marker_gene_finder.analyze` / `analyze_file` for new analyses |
| `src.analysis_profiles` | `marker_gene_finder.profiles` |
| `src.expression_summary` | `marker_gene_finder.expression_summary` |
| `src.filter_presets` | `marker_gene_finder.filter_presets` |
| `src.report_bundle` | `marker_gene_finder.report_bundle` |
| `src.settings_import` | `marker_gene_finder.settings` |
| `src.version` | `marker_gene_finder.version` (or root `APP_VERSION`) |

The six helper modules retain the same public names and signatures. Import
`PROFILE_DEFAULTS`, `PROFILE_STATE_KEYS`, profile functions, presets, settings
readers, and bundle/summary helpers from their replacement module.

`src.marker_finder` is split by responsibility:

| Names | Replacement module |
| --- | --- |
| `run`, `load_expression`, `find_markers`, `find_markers_scanpy`, `RESULT_COLUMNS`, `RANK_BY_CHOICES`, `DEFAULT_GENE_CHUNK_SIZE`, `DEFAULT_MAX_DENSE_MEMORY_MB` | `marker_gene_finder.analysis` |
| `inspect_input`, `build_analysis_plan`, `estimate_dense_expression_bytes`, `PLAN_SCHEMA_VERSION`, `CSV_INSPECTION_CHUNK_ROWS` | `marker_gene_finder.inspection` |
| `run_batch`, `build_batch_plan`, `verify_batch_manifest`, `BATCH_MANIFEST_SCHEMA_VERSION`, `BATCH_PLAN_SCHEMA_VERSION`, `BATCH_VERIFICATION_SCHEMA_VERSION` | `marker_gene_finder.batch` |
| `benchmark_analysis`, `compare_benchmarks`, `BENCHMARK_SCHEMA_VERSION`, `BENCHMARK_COMPARISON_SCHEMA_VERSION` | `marker_gene_finder.benchmark` |
| `ProgressCallback` | `marker_gene_finder.progress` |
| `APP_VERSION` | `marker_gene_finder.version` |
| `main`, `build_parser` | `marker_gene_finder.cli`; arguments now use subcommands, as below |

Underscore-prefixed helpers and incidental imported dependencies are internal,
not supported public imports. Patch or inject the documented engine interface
in integrations; legacy `_core` patch forwarding is transitional compatibility.
`src` itself is an empty package, and `src._compat` is internal warning machinery.

Existing DataFrame-returning calls can change only their import:

```python
from marker_gene_finder.analysis import run
markers = run("expression.csv", "markers.csv", top_n=10)
```

For configuration, plugins, and reproducibility metadata, use the typed API:

```python
from marker_gene_finder import AnalysisConfig, analyze_file
config = AnalysisConfig(top_n=10, target_clusters=("A",))
result = analyze_file("expression.csv", "markers.csv", config)
result.to_bundle("analysis.zip")
config.to_json("analysis-config.json")
```

`analyze` accepts a DataFrame and returns `MarkerResult`; `analyze_file` writes
CSV and returns the same result type. The lower-level `run` retains its existing
DataFrame result and diagnostics/progress parameters. `load_expression` keeps
its historical `cell_column="cell_id"` default; `run` and `AnalysisConfig` use
`"cell"`. Pass the column explicitly when moving between these interfaces.

## Commands

Install with `python -m pip install '.[web,test]'`. Use either
`marker-gene-finder` or `python -m marker_gene_finder`. Both invoke the same
subcommand parser. Replace `python src/marker_finder.py`,
`python -m src.marker_finder`, and the historical flat installed command.
Calling `python -m marker_gene_finder._core` is also transitional, not a v3 API.

| Legacy mode | Replacement |
| --- | --- |
| `--input INPUT --output OUTPUT` | `marker-gene-finder analyze INPUT OUTPUT --config config.json` |
| `--input INPUT --inspect-json REPORT` | `marker-gene-finder inspect INPUT --output REPORT --config config.json` |
| `--input INPUT --plan-json REPORT` | `marker-gene-finder plan INPUT --output REPORT --config config.json` |
| `--input-dir INPUTS --output-dir OUTPUTS` | `marker-gene-finder batch INPUTS OUTPUTS --config config.json` |
| Batch plus `--batch-plan-json REPORT` | `marker-gene-finder batch-plan INPUTS OUTPUTS --output REPORT --config config.json` |
| `--verify-batch-manifest MANIFEST --verification-json REPORT` | `marker-gene-finder verify MANIFEST --output REPORT` |
| `--input INPUT --benchmark-json REPORT --benchmark-repeats N` | `marker-gene-finder benchmark INPUT --output REPORT --repeats N --config config.json` |
| `--compare-benchmarks BASE CANDIDATE --comparison-json REPORT` | `marker-gene-finder compare-benchmarks BASE CANDIDATE --output REPORT` |

Configuration is optional; omission selects the unchanged defaults. For example,
save `{"top_n": 10, "engine": "native"}` as `config.json`, then run the first
command above. Full documents may include `schema_name` equal to
`marker-gene-finder.analysis-config` and `schema_version` equal to `2`.

`--diagnostics-json` becomes analyze's `--metadata`, which includes provenance
alongside runtime diagnostics. `--quiet` needs no replacement: subcommands do
not print per-cluster progress. The Python `analysis.run` callback remains
available when progress is needed. Paths for input, config, output, metadata,
and bundle must be distinct in the analyze command.

## Saved files

| Artifact | Reader and migration behavior |
| --- | --- |
| Analysis config schema 2 | `AnalysisConfig.from_json` / `from_dict`; omitted schema metadata is accepted; unknown fields or schemas raise `SchemaError` |
| Flat settings/profile schema 1 | `migrate_settings_document` converts to schema 2 in memory; settings/profile upload readers also migrate automatically |
| Named settings/profile schema 2 | `load_analysis_settings_upload` / `load_analysis_profile_upload`; JSON or ZIP with `analysis_settings.json` |
| Existing result bundles | Settings readers retain legacy ZIP support; new `MarkerResult.to_bundle` includes replayable `analysis_config.json` and provenance |
| Batch manifest schema 2 | `verify_batch_manifest` and CLI `verify`; `batch --resume` checks fingerprints and hashes before skipping |
| Benchmark schema 1 | `compare_benchmarks` and CLI `compare-benchmarks` |

Application version and data schema version are separate. Do not change a saved
schema number to match the application version. Unsupported settings schemas
receive an error naming supported versions. Verification rejects unsupported
batch schemas; rerun the batch into a new output directory to generate a current
manifest instead of editing hashes by hand. A changed config fingerprint may
cause a safe rerun during migration, even when scientific settings are equivalent.

Profiles contain portable UI settings, not dataset identity or a complete
`AnalysisConfig`. Apply them through the app or profile readers. To replay a new
result bundle through the CLI, extract its `analysis_config.json` and supply the
original dataset separately. Do not pass an `analysis_settings.json` envelope as
`--config`. Upload readers return UI state, not an AnalysisConfig document.

## Executable examples and verification

Run `python examples/migration.py /tmp/mgf-migration` from an installed checkout.
The [example source](../examples/migration.py) creates its own data and exercises
CSV/API equivalence, H5AD X/layers/raw in both engines, a reference plugin,
configuration round trips, profiles, bundles, all CLI workflows, batch resume,
verification, and benchmark comparison. Tests execute this same file.

Run the primary suite with legacy imports blocked:

```sh
MGF_DISABLE_LEGACY=1 python -m pytest -q --ignore=tests/compatibility
python -m pytest -q tests/compatibility
```

CI builds and installs wheels on Python 3.10, 3.11, and 3.12, then copies the
example, Streamlit app, and wheel-check script outside the repository. It checks
the imported package path, blocks legacy imports, executes the example and app,
and smoke-tests the installed command. This v2.9 preparation does not replace
the v3 release-candidate gate: repeat the examples, scientific regressions,
installation matrix, and artifact migration checks after the actual removals.

## Complete flat-option mapping

Use snake_case fields in the JSON passed to `--config`. Repeated target/gene/prefix
flags become JSON arrays. The same config is accepted by analyze, inspect, plan,
batch, batch-plan, and benchmark (each consumes its relevant settings).

| Legacy option | Replacement |
| --- | --- |
| `-h`, `--help` | `marker-gene-finder --help` or command `--help` |
| `--version` | `marker-gene-finder --version` |
| `--input` | Input positional argument on analyze/inspect/plan/benchmark |
| `--input-dir` | Input directory positional argument on batch/batch-plan |
| `--verify-batch-manifest` | `verify MANIFEST` |
| `--compare-benchmarks` | `compare-benchmarks BASE CANDIDATE` |
| `--output` | Output positional argument on analyze |
| `--output-dir` | Output directory positional argument on batch/batch-plan |
| `--verification-json` | verify `--output REPORT` |
| `--comparison-json` | compare-benchmarks `--output REPORT` |
| `--engine` | Config `engine` |
| `--cluster-column` | Config `cluster_column` |
| `--cell-column` | Config `cell_column` |
| `--max-dense-memory-mb` | Config `max_dense_memory_mb`; translate legacy 0 (disabled) to JSON null |
| `--gene-chunk-size` | Config `gene_chunk_size` |
| `--workers` | Config `workers` |
| `--diagnostics-json` | analyze `--metadata REPORT` |
| `--inspect-json` | `inspect INPUT --output REPORT` |
| `--plan-json` | `plan INPUT --output REPORT` |
| `--benchmark-json` | `benchmark INPUT --output REPORT` |
| `--benchmark-repeats` | benchmark `--repeats N` |
| `--max-regression-percent` | compare-benchmarks `--max-regression-percent N` |
| `--resume` | batch/batch-plan `--resume` |
| `--retry-failed` | batch/batch-plan `--retry-failed` |
| `--fail-fast` | batch `--fail-fast` |
| `--include-pattern` | batch/batch-plan repeated `--include-pattern GLOB` |
| `--exclude-pattern` | batch/batch-plan repeated `--exclude-pattern GLOB` |
| `--batch-plan-json` | `batch-plan INPUTS OUTPUTS --output REPORT` |
| `--quiet` | No flag needed; subcommands do not emit per-cluster progress |
| `--target-cluster` | Config `target_clusters` array of string cluster labels |
| `--min-cluster-size` | Config `min_cluster_size` |
| `--exclude-gene` | Config `excluded_genes` array |
| `--exclude-prefix` | Config `excluded_prefixes` array |
| `--top-n` | Config `top_n` |
| `--pseudocount` | Config `pseudocount` |
| `--min-pct` | Config `min_pct` |
| `--min-cells-expressing` | Config `min_cells_expressing` |
| `--min-cluster-mean` | Config `min_cluster_mean` |
| `--min-mean-difference` | Config `min_mean_difference` |
| `--max-other-mean` | Config `max_other_mean` |
| `--max-cells-expressing-out` | Config `max_cells_expressing_out` |
| `--min-log2fc` | Config `min_log2fc` |
| `--max-pct-out` | Config `max_pct_out` |
| `--rank-by` | Config `rank_by` |
| `--min-pct-difference` | Config `min_pct_difference` |
| `--max-p-adj` | Config `max_p_adj` |
| `--layer` | Config `layer` |
| `--use-raw` | Config `use_raw` |
