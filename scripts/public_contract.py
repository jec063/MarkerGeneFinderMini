"""Describe the planned v3 surface; reviewed snapshot lives in docs/."""
import argparse
from dataclasses import asdict
import importlib
import inspect

import marker_gene_finder as mgf
from marker_gene_finder.analysis import RESULT_COLUMNS
from marker_gene_finder.cli import build_parser


SERVICES = {
    "analysis": ["run", "load_expression", "find_markers", "find_markers_scanpy"],
    "inspection": ["inspect_input", "build_analysis_plan", "estimate_dense_expression_bytes"],
    "batch": ["run_batch", "build_batch_plan", "verify_batch_manifest"],
    "benchmark": ["benchmark_analysis", "compare_benchmarks"],
    "expression_summary": ["summarize_expression", "summarize_clusters"],
    "filter_presets": ["get_filter_preset"],
    "profiles": ["reset_profile_state", "build_analysis_profile", "compare_profile_state", "is_newer_app_version"],
    "settings": ["build_settings_document", "migrate_settings_document", "extract_analysis_settings_upload", "load_analysis_settings_upload", "load_analysis_profile_upload", "load_analysis_settings"],
    "report_bundle": ["build_analysis_bundle"],
    "cli": ["main", "build_parser"],
    "progress": ["console_progress"],
}


def parameters(function):
    result = []
    for p in inspect.signature(function).parameters.values():
        entry = {"name": p.name, "kind": p.kind.name}
        if p.default is not inspect.Parameter.empty:
            value = p.default
            entry["default"] = (value if value is None or isinstance(value, (str, int, float, bool, list, tuple, dict))
                                else type(value).__module__ + "." + type(value).__qualname__)
        result.append(entry)
    return result


def contract():
    callables = {}
    for name in mgf.__all__:
        value = getattr(mgf, name)
        if inspect.isfunction(value):
            callables[name] = parameters(value)
    for cls in (mgf.AnalysisConfig, mgf.MarkerResult, mgf.EngineCapabilities,
                mgf.EngineOutput, mgf.EngineInfo, mgf.EngineRegistry):
        callables[cls.__name__] = parameters(cls)
    for cls, methods in [
        (mgf.AnalysisConfig, ["to_dict", "to_json", "from_dict", "from_json", "to_analysis_kwargs"]),
        (mgf.MarkerResult, ["to_csv", "to_bundle"]),
        (mgf.EngineRegistry, ["get", "register", "names", "copy"]),
        (mgf.AnalysisEngine, ["analyze_dataframe", "analyze_file"]),
    ]:
        for name in methods:
            callables[cls.__name__ + "." + name] = parameters(getattr(cls, name))
    for module, names in SERVICES.items():
        imported = importlib.import_module("marker_gene_finder." + module)
        for name in names:
            callables[module + "." + name] = parameters(getattr(imported, name))
    parser = build_parser()
    commands = next(a.choices for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    service_exports = {}
    schemas = {}
    for module in SERVICES:
        imported = importlib.import_module("marker_gene_finder." + module)
        if hasattr(imported, "__all__"):
            service_exports[module] = sorted(imported.__all__)
        for name, value in vars(imported).items():
            if name.endswith("_SCHEMA_VERSION"):
                schemas[module + "." + name] = value
    return {
        "exports": sorted(mgf.__all__), "callables": callables,
        "service_exports": service_exports, "schemas": schemas,
        "analysis_defaults": asdict(mgf.AnalysisConfig()),
        "engine_capability_defaults": asdict(mgf.EngineCapabilities()),
        "result_columns": RESULT_COLUMNS,
        "commands": {name: [a.dest for a in command._actions] for name, command in commands.items()},
        "config_schema": mgf.CONFIG_SCHEMA_VERSION,
        "settings_schema": mgf.SETTINGS_SCHEMA_VERSION,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(contract(), indent=2))
