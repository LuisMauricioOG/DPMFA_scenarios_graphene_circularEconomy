# -*- coding: utf-8 -*-
"""
Prospective DPMFA Local Sensitivity Analysis (direct differential sensitivity analysis ) for S1 (Waste-as-Usual) and S3 (Circular Economy-oriented)
@author: Luis Mauricio Ortiz-Galvez (with help of ChatGPT 5.5)

Based on Bornhöft, N.A., Nowack, B. & Hilty, L.M. (2021).
Representation, Propagation, and Interpretation of Uncertain Knowledge in Dynamic Probabilistic Material Flow Models.
Environ Model Assess 26, 709–721. https://doi.org/10.1007/s10666-021-09775-5



Purpose
-------
1) build mean-based local-sensitivity tables:
   - mean cumulative sink values
   - relative change of means
   - paired run-level Q5-Q95 relative-change interval
   - normalized sensitivity ratio where input_change_pct is defined

2) create:
   - Excel sensitivity table
   - Sensitivity heatmap matrices
   - tornado ranking plots for ±20% scalar TC variants
   - annual sensitivity tables for selected sinks in 2030 and 2050
   - cumulative sink bubble maps for 2014-2050
   - S1 BattWDism_dynamic annual output-response curves
   - S1 BattWDism_dynamic annual normalized S(t), using TC values extracted from model files
   - cross-scenario S1_BattWDism_dynamic vs S3_base comparison,
     using S1_BattWDism_dynamic as the reference.

Important methodological note
-----------------------------
The normalized local sensitivity coefficient is calculated only when the input
change is a scalar percentage change:

    S = (%Delta output) / (%Delta input)

Allocation variants and dynamic trajectory variants are still evaluated as
output-side changes, but not as scalar normalized sensitivity ratios unless
an effective input change is explicitly available.
For BattWDism_dynamic, this script extracts the yearly TC trajectory from the
S1 model files and calculates an optional annual S(t).
"""

from __future__ import annotations

import ast
import math
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Publication-style plotting defaults. Arial is attempted first; common fallbacks
# are used if Arial is not installed on the system running the script.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

PROJECT_ROOT = Path(
    r"C:\Users\ortl\Documents\Empa\Articles\pDPMFA_AdMa_EU_CE"
)

RUNS_SOURCE_ROOT = (
    PROJECT_ROOT
    / "ProspectiveMFA_scenarios"
    / "Extended_DPMFA_scenarios_graphene"
    / "SensitivityAnalysis_runs"
)

RUN_DATE = datetime.now().strftime("%Y%m%d")
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

RESULTS_ROOT = PROJECT_ROOT / f"{RUN_DATE}_SensitivityAnalysis_results"

# Some existing runner scripts may still write raw outputs to an old fixed folder.
# The master script checks all candidate raw-output folders after each run.
RAW_RUNNER_OUTPUT_DIRS = [
    RESULTS_ROOT / "_raw_runner_outputs",
    PROJECT_ROOT / "20260622_SensitivityAnalysis_results" / "_raw_runner_outputs",
]

SCENARIOS_TO_RUN = ["S1", "S3"]

# Set False if you already have organized result folders under RESULTS_ROOT.
RUN_RUNNERS_FIRST = True

# If True, the script renames an existing destination folder as *_backup_<timestamp>
# before moving new results into S1/variant or S3/variant.
BACKUP_EXISTING_VARIANT_FOLDERS = True

PYTHON_EXE = sys.executable

START_YEAR = 2014
END_YEAR = 2050
YEARS = list(range(START_YEAR, END_YEAR + 1))

CUMULATIVE_PERIODS = [
    {"label": "cumulative_full_boundary_2014_2050", "from_year": 2014, "to_year": 2050},
    {"label": "cumulative_policy_period_2025_2050", "from_year": 2025, "to_year": 2050},
]
SELECTED_YEARS = [2030, 2040, 2050]

# Practical relevance threshold used in output-side tables and plots.
PRACTICAL_EFFECT_THRESHOLD_PCT = 20.0

# Normalized sensitivity outputs are usually clearest for final sinks.
SINK_OUTPUTS = {
    "EoW_alternatives": [
        "Results_SinkDynamic_EoW_alternatives*.csv",
        "Results_SinkDynamic_End-of-Waste*.csv",
    ],
    "Reuse": ["Results_SinkDynamic_Reuse*.csv"],
    "Landfill": ["Results_SinkDynamic_Landfill*.csv"],
    "Elimination": ["Results_SinkDynamic_Elimination*.csv"],
    "Export": ["Results_SinkDynamic_Export*.csv"],
    "Unknown": ["Results_SinkDynamic_Unknown*.csv"],
    "Roadside": ["Results_SinkDynamic_Roadside*.csv"],
    "STsoil": [
        "Results_SinkDynamic_STsoil*.csv",
        "Results_SinkDynamic_Sludge treated soil*.csv",
    ],
    "NUsoil": [
        "Results_SinkDynamic_NUsoil*.csv",
        "Results_SinkDynamic_Natural*Urban soil*.csv",
    ],
    "Surfacewater": ["Results_SinkDynamic_Surfacewater*.csv"],
    "Subsurface": ["Results_SinkDynamic_Subsurface*.csv"],
}

SINK_LABELS = {
    "EoW_alternatives": "End-of-Waste alternatives",
    "Reuse": "Reuse",
    "Landfill": "Landfill",
    "Elimination": "Elimination",
    "Export": "Export",
    "Unknown": "Unknown",
    "Roadside": "Roadside",
    "STsoil": "Sludge-treated soil",
    "NUsoil": "Natural & urban soil",
    "Surfacewater": "Surface water",
    "Subsurface": "Subsurface",
}

MAIN_SINKS_FOR_HEATMAP = [
    "EoW_alternatives",
    "Reuse",
    "Landfill",
    "Elimination",
    "Export",
    "Unknown",
    "Roadside",
    "STsoil",
    "NUsoil",
    "Surfacewater",
    "Subsurface",
]

MAIN_SINKS_FOR_TORNADO = [
    "EoW_alternatives",
    "Landfill",
    "Elimination",
    "Unknown",
]

# Annual sensitivity table and simple bubble maps.
# These focus on the same sink compartments requested for article/SI reporting.
ANNUAL_TABLE_YEARS = [2030, 2050]
ANNUAL_TABLE_SINKS = [
    "Landfill",
    "EoW_alternatives",
    "Unknown",
    "Elimination",
]

KEY_SINKS_DYNAMIC_CURVE = [
    "EoW_alternatives",
    "Landfill",
    "Elimination",
    "Unknown",
]

# Use only variants whose scalar input change is exactly ±20% in the tornado ranking.
TORNADO_ONLY_ABS_20_PERCENT = True

SAVE_DPI = 600


# =============================================================================
# 2. INPUT-CHANGE METADATA
# =============================================================================

def infer_input_change_pct(scenario: str, variant: str) -> Optional[float]:
    """
    Input-change convention:
        input_change_pct = ((variant input - baseline input) / baseline input) * 100

    Use None when the variant is structural, categorical, or time-dependent
    without one scalar perturbation.
    """
    v = variant.lower()

    if v in {"base", "baseline", "reference", "original"}:
        return None

    # Allocation/fate changes.
    allocation_variants = {
        "battw2lf",
        "battrec_unk2lf",
        "battrec_unk2lf",
        "battw_recres2lf",
        "battw_recres2lf",
    }
    if v in allocation_variants or "2lf" in v:
        return None

    # Dynamic trajectory case.
    if v == "battwdism_dynamic":
        return None


    # Static ±20% tire TC variants.
    if v.endswith("_high") and ("tiresenrec" in v or "tiresrec" in v):
        return 20.0
    if v.endswith("_low") and ("tiresenrec" in v or "tiresrec" in v):
        return -20.0

    # S3 direct recycling deterministic variants.
    # Baseline TC = 1.0. Low = 0.40, high = 0.95.
    if scenario == "S3" and v == "battdirrec_low":
        return -60.0
    if scenario == "S3" and v == "battdirrec_high":
        return -5.0

    return None


# =============================================================================
# 3. GENERAL HELPERS
# =============================================================================

def ensure_dirs() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    for folder in RAW_RUNNER_OUTPUT_DIRS:
        folder.mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / "_sensitivity_comparison_tables").mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / "_normalized_sensitivity_outputs").mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / "_cross_scenario_S1dynamic_vs_S3base").mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / "_annual_sensitivity_tables").mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / "_cumulative_sink_bubble_maps").mkdir(parents=True, exist_ok=True)


def period_index(year: int) -> int:
    idx = year - START_YEAR
    if idx < 0 or year > END_YEAR:
        raise ValueError(f"Year {year} is outside {START_YEAR}-{END_YEAR}.")
    return idx


def clean_filename(text: str) -> str:
    text = str(text)
    for ch in ['\\', '/', ':', '*', '?', '"', '<', '>', '|']:
        text = text.replace(ch, "_")
    return text.replace(" ", "_")


def format_seconds(seconds: float) -> str:
    seconds = int(seconds)
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def discover_runner_jobs() -> List[Dict[str, Any]]:
    """
    Discover runner scripts from:
        SensitivityAnalysis_runs/S1_SA
        SensitivityAnalysis_runs/S3_SA

    Expected filename pattern:
        DPMFA_graphene_TBC_runner_S1_base.py
        DPMFA_graphene_TBC_runner_S1_BattWDism_dynamic.py
        DPMFA_graphene_TBC_runner_S3_TiresRec_low.py
    """
    jobs: List[Dict[str, Any]] = []

    for scenario in SCENARIOS_TO_RUN:
        scenario_folder = RUNS_SOURCE_ROOT / f"{scenario}_SA"
        if not scenario_folder.exists():
            print(f"Warning: scenario folder not found: {scenario_folder}")
            continue

        pattern = f"DPMFA_graphene_TBC_runner_{scenario}_*.py"
        scripts = sorted(scenario_folder.glob(pattern))

        for script in scripts:
            if "backup" in script.name.lower():
                continue

            prefix = f"DPMFA_graphene_TBC_runner_{scenario}_"
            variant = script.stem.replace(prefix, "")

            # Exclude accidental utility files that match weakly.
            if not variant or variant.lower().startswith("_"):
                continue

            jobs.append({
                "scenario": scenario,
                "variant": variant,
                "script": script,
                "input_change_pct": infer_input_change_pct(scenario, variant),
            })

    # Stable useful order: base first, then the rest alphabetically.
    jobs = sorted(jobs, key=lambda j: (j["scenario"], 0 if j["variant"].lower() == "base" else 1, j["variant"]))
    return jobs


def snapshot_raw_output_dirs() -> set[Path]:
    before: set[Path] = set()
    for raw_dir in RAW_RUNNER_OUTPUT_DIRS:
        if raw_dir.exists():
            before.update([p.resolve() for p in raw_dir.iterdir() if p.is_dir()])
    return before


def find_newest_new_raw_folder(before: set[Path]) -> Path:
    after: List[Path] = []
    for raw_dir in RAW_RUNNER_OUTPUT_DIRS:
        if raw_dir.exists():
            after.extend([p for p in raw_dir.iterdir() if p.is_dir() and p.resolve() not in before])

    if not after:
        # Fallback: if no new folder was detected, use newest folder across candidates.
        all_candidates: List[Path] = []
        for raw_dir in RAW_RUNNER_OUTPUT_DIRS:
            if raw_dir.exists():
                all_candidates.extend([p for p in raw_dir.iterdir() if p.is_dir()])
        if not all_candidates:
            raise FileNotFoundError(
                "No runner output folder found in any raw-output candidate folder."
            )
        return max(all_candidates, key=lambda p: p.stat().st_mtime)

    return max(after, key=lambda p: p.stat().st_mtime)


def run_runner_jobs(jobs: List[Dict[str, Any]]) -> None:
    if not jobs:
        raise ValueError("No runner jobs discovered.")

    ensure_dirs()
    print(f"Results root: {RESULTS_ROOT}")
    print(f"Runner source root: {RUNS_SOURCE_ROOT}")
    print(f"Discovered {len(jobs)} runner jobs.")

    total_jobs = len(jobs)
    batch_start = time.perf_counter()
    completed_durations: List[float] = []

    for i, job in enumerate(jobs, start=1):
        scenario = job["scenario"]
        variant = job["variant"]
        script = Path(job["script"])

        if not script.exists():
            raise FileNotFoundError(f"Runner script not found: {script}")

        elapsed = time.perf_counter() - batch_start
        if completed_durations:
            avg = sum(completed_durations) / len(completed_durations)
            remaining = avg * (total_jobs - i + 1)
            remaining_text = format_seconds(remaining)
        else:
            remaining_text = "available after first completed job"

        print("\n" + "=" * 80)
        print(f"Running job {i}/{total_jobs}: {scenario} / {variant}")
        print(f"Elapsed: {format_seconds(elapsed)}")
        print(f"Estimated remaining: {remaining_text}")
        print(f"Script: {script}")
        print("=" * 80)

        before = snapshot_raw_output_dirs()
        start = time.perf_counter()

        completed = subprocess.run(
            [PYTHON_EXE, str(script)],
            cwd=str(script.parent),
            text=True,
            capture_output=True,
        )

        if completed.stdout:
            print("\n--- Runner STDOUT ---")
            print(completed.stdout)
        if completed.stderr:
            print("\n--- Runner STDERR ---")
            print(completed.stderr)

        if completed.returncode != 0:
            raise RuntimeError(
                f"Runner failed for {scenario} / {variant}\n"
                f"Script: {script}\n"
                f"Return code: {completed.returncode}"
            )

        raw_result_folder = find_newest_new_raw_folder(before)
        destination = RESULTS_ROOT / scenario / variant
        destination.parent.mkdir(parents=True, exist_ok=True)

        if destination.exists():
            if BACKUP_EXISTING_VARIANT_FOLDERS:
                backup = destination.with_name(destination.name + f"_backup_{TIMESTAMP}")
                print(f"Existing destination found. Moving old folder to: {backup}")
                shutil.move(str(destination), str(backup))
            else:
                print(f"Removing existing destination: {destination}")
                shutil.rmtree(destination)

        print("Moving result folder:")
        print(f"  from: {raw_result_folder}")
        print(f"  to:   {destination}")
        shutil.move(str(raw_result_folder), str(destination))

        duration = time.perf_counter() - start
        completed_durations.append(duration)
        print(f"Finished job {i}/{total_jobs}: {scenario} / {variant}")
        print(f"Job duration: {format_seconds(duration)}")

    print("\nAll runner jobs completed and organized.")


# =============================================================================
# 4. READING AND COMPARING RUN-BY-YEAR MATRICES
# =============================================================================

def find_latest_matching_file(folder: Path, patterns: List[str]) -> Optional[Path]:
    matches: List[Path] = []
    for pattern in patterns:
        matches.extend(folder.rglob(pattern))
    matches = [p for p in matches if p.is_file()]
    if not matches:
        return None
    return max(matches, key=lambda p: p.stat().st_mtime)


def read_run_by_year_matrix(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    first_col = str(df.columns[0])
    if first_col.startswith("Unnamed") or first_col.lower() in ["index", "run"]:
        df = df.drop(columns=df.columns[0])

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")

    if df.empty:
        raise ValueError(f"No numeric data found in {csv_path}")

    return df


def extract_indicator_values(
    matrix: pd.DataFrame,
    indicator: str,
    year: Optional[int] = None,
    cumulative_from_year: Optional[int] = None,
    cumulative_to_year: Optional[int] = None,
) -> np.ndarray:
    values = matrix.to_numpy(dtype=float)

    if indicator == "cumulative":
        if cumulative_from_year is None or cumulative_to_year is None:
            raise ValueError("cumulative_from_year and cumulative_to_year are required.")
        start_idx = period_index(cumulative_from_year)
        end_idx = period_index(cumulative_to_year)
        return np.nansum(values[:, start_idx:end_idx + 1], axis=1)

    if indicator == "year":
        if year is None:
            raise ValueError("year must be provided for indicator='year'.")
        idx = period_index(year)
        return values[:, idx]

    raise ValueError(f"Unknown indicator: {indicator}")


def safe_nanmean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.nanmean(values)) if values.size else np.nan


def safe_nanmedian(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.nanmedian(values)) if values.size else np.nan


def safe_nanpercentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.nanpercentile(values, q)) if values.size else np.nan


def compare_variant_to_base(
    base_values: np.ndarray,
    variant_values: np.ndarray,
    input_change_pct: Optional[float] = None,
) -> Dict[str, Any]:
    base_values = np.asarray(base_values, dtype=float)
    variant_values = np.asarray(variant_values, dtype=float)

    n_paired = min(base_values.size, variant_values.size)
    base_paired = base_values[:n_paired]
    variant_paired = variant_values[:n_paired]
    finite_pair = np.isfinite(base_paired) & np.isfinite(variant_paired)

    base_finite = base_paired[finite_pair]
    variant_finite = variant_paired[finite_pair]
    absolute_change = variant_finite - base_finite

    mean_base = safe_nanmean(base_finite)
    mean_variant = safe_nanmean(variant_finite)
    absolute_change_of_means = mean_variant - mean_base if np.isfinite(mean_base) and np.isfinite(mean_variant) else np.nan

    if np.isfinite(mean_base) and mean_base != 0:
        relative_change_of_means_pct = (absolute_change_of_means / mean_base) * 100.0
    else:
        relative_change_of_means_pct = np.nan

    valid_relative = finite_pair & (base_paired != 0)
    if valid_relative.any():
        relative_change = ((variant_paired[valid_relative] - base_paired[valid_relative]) / base_paired[valid_relative]) * 100.0
    else:
        relative_change = np.array([], dtype=float)

    mean_relative_change_pct = safe_nanmean(relative_change)
    median_relative_change_pct = safe_nanmedian(relative_change)
    q05_relative_change_pct = safe_nanpercentile(relative_change, 5)
    q95_relative_change_pct = safe_nanpercentile(relative_change, 95)

    if np.isfinite(q05_relative_change_pct) and np.isfinite(q95_relative_change_pct):
        q05_q95_crosses_zero = bool(q05_relative_change_pct <= 0 <= q95_relative_change_pct)
    else:
        q05_q95_crosses_zero = True

    if input_change_pct is not None and np.isfinite(input_change_pct) and input_change_pct != 0:
        sensitivity_ratio_mean = relative_change_of_means_pct / input_change_pct
        sensitivity_ratio_median = median_relative_change_pct / input_change_pct
        sensitivity_ratio_available = True
    else:
        sensitivity_ratio_mean = np.nan
        sensitivity_ratio_median = np.nan
        sensitivity_ratio_available = False

    robust_direction = not q05_q95_crosses_zero
    passes_threshold_mean = (
        bool(abs(relative_change_of_means_pct) >= PRACTICAL_EFFECT_THRESHOLD_PCT)
        if np.isfinite(relative_change_of_means_pct)
        else False
    )

    return {
        "n_paired_runs": int(n_paired),
        "n_finite_pairs": int(finite_pair.sum()),
        "n_runs_used_for_relative_change": int(relative_change.size),
        "n_zero_baseline_excluded_from_relative_change": int(finite_pair.sum() - relative_change.size),

        "mean_base": mean_base,
        "mean_variant": mean_variant,
        "absolute_change_of_means": absolute_change_of_means,
        "relative_change_of_means_pct": relative_change_of_means_pct,
        "abs_relative_change_of_means_pct": abs(relative_change_of_means_pct) if np.isfinite(relative_change_of_means_pct) else np.nan,

        "mean_relative_change_pct": mean_relative_change_pct,
        "median_relative_change_pct": median_relative_change_pct,
        "q05_relative_change_pct": q05_relative_change_pct,
        "q95_relative_change_pct": q95_relative_change_pct,
        "q05_q95_crosses_zero": q05_q95_crosses_zero,
        "robust_direction": robust_direction,
        "passes_practical_threshold_mean": passes_threshold_mean,
        "robust_and_practically_relevant_mean": bool(robust_direction and passes_threshold_mean),

        "input_change_pct": input_change_pct,
        "sensitivity_ratio_available": sensitivity_ratio_available,
        "sensitivity_ratio_mean": sensitivity_ratio_mean,
        "sensitivity_ratio_median": sensitivity_ratio_median,
    }


def discover_organized_results() -> Dict[str, Dict[str, Path]]:
    scenarios: Dict[str, Dict[str, Path]] = {}

    for scenario in SCENARIOS_TO_RUN:
        scenario_dir = RESULTS_ROOT / scenario
        if not scenario_dir.exists():
            print(f"Warning: organized scenario folder not found: {scenario_dir}")
            continue

        variants = {}
        for p in sorted(scenario_dir.iterdir()):
            if p.is_dir() and not p.name.startswith("_") and "backup" not in p.name.lower():
                variants[p.name] = p

        if variants:
            scenarios[scenario] = variants

    return scenarios


def find_baseline_variant(variant_map: Dict[str, Path]) -> Tuple[str, Path]:
    for name, path in variant_map.items():
        if name.lower() in {"base", "baseline", "reference", "original"}:
            return name, path
    raise ValueError(f"Could not find baseline folder among: {list(variant_map.keys())}")


def build_sensitivity_table(jobs: List[Dict[str, Any]]) -> pd.DataFrame:
    input_change_by_variant = {
        (job["scenario"], job["variant"]): job.get("input_change_pct")
        for job in jobs
    }

    scenarios = discover_organized_results()
    rows: List[Dict[str, Any]] = []

    for scenario, variants in scenarios.items():
        base_name, base_folder = find_baseline_variant(variants)

        for variant, variant_folder in variants.items():
            if variant == base_name:
                continue

            input_change_pct = input_change_by_variant.get(
                (scenario, variant),
                infer_input_change_pct(scenario, variant),
            )

            for output_name, patterns in SINK_OUTPUTS.items():
                base_file = find_latest_matching_file(base_folder, patterns)
                variant_file = find_latest_matching_file(variant_folder, patterns)

                if base_file is None or variant_file is None:
                    print(
                        f"Skipping missing output: {scenario}/{variant}/{output_name} "
                        f"(base={base_file}, variant={variant_file})"
                    )
                    continue

                base_matrix = read_run_by_year_matrix(base_file)
                variant_matrix = read_run_by_year_matrix(variant_file)

                # Cumulative periods.
                for period in CUMULATIVE_PERIODS:
                    base_values = extract_indicator_values(
                        base_matrix,
                        "cumulative",
                        cumulative_from_year=period["from_year"],
                        cumulative_to_year=period["to_year"],
                    )
                    variant_values = extract_indicator_values(
                        variant_matrix,
                        "cumulative",
                        cumulative_from_year=period["from_year"],
                        cumulative_to_year=period["to_year"],
                    )

                    metrics = compare_variant_to_base(
                        base_values,
                        variant_values,
                        input_change_pct=input_change_pct,
                    )

                    rows.append({
                        "scenario": scenario,
                        "baseline": base_name,
                        "variant": variant,
                        "output": output_name,
                        "output_label": SINK_LABELS.get(output_name, output_name),
                        "indicator": period["label"],
                        "indicator_type": "cumulative",
                        "from_year": period["from_year"],
                        "to_year": period["to_year"],
                        "year": np.nan,
                        "base_file": str(base_file),
                        "variant_file": str(variant_file),
                        **metrics,
                    })

                # Annual selected years.
                for year in SELECTED_YEARS:
                    base_values = extract_indicator_values(base_matrix, "year", year=year)
                    variant_values = extract_indicator_values(variant_matrix, "year", year=year)

                    metrics = compare_variant_to_base(
                        base_values,
                        variant_values,
                        input_change_pct=input_change_pct,
                    )

                    rows.append({
                        "scenario": scenario,
                        "baseline": base_name,
                        "variant": variant,
                        "output": output_name,
                        "output_label": SINK_LABELS.get(output_name, output_name),
                        "indicator": f"year_{year}",
                        "indicator_type": "year",
                        "from_year": np.nan,
                        "to_year": np.nan,
                        "year": year,
                        "base_file": str(base_file),
                        "variant_file": str(variant_file),
                        **metrics,
                    })

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("Sensitivity table is empty. Check that runner outputs exist and file patterns match.")

    df["rank_by_abs_mean_effect"] = (
        df.groupby(["scenario", "output", "indicator"])["abs_relative_change_of_means_pct"]
        .rank(method="dense", ascending=False)
    )

    df["rank_by_abs_sensitivity_ratio"] = (
        df.groupby(["scenario", "output", "indicator"])["sensitivity_ratio_mean"]
        .transform(lambda s: s.abs().rank(method="dense", ascending=False))
    )

    return df


def save_sensitivity_table(df: pd.DataFrame) -> Path:
    output_dir = RESULTS_ROOT / "_sensitivity_comparison_tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    path = output_dir / f"DPMFA_sensitivity_comparison_ALL_{TIMESTAMP}.xlsx"

    compact_cols = [
        "scenario",
        "baseline",
        "variant",
        "output",
        "output_label",
        "indicator",
        "mean_base",
        "mean_variant",
        "absolute_change_of_means",
        "relative_change_of_means_pct",
        "q05_relative_change_pct",
        "q95_relative_change_pct",
        "q05_q95_crosses_zero",
        "robust_and_practically_relevant_mean",
        "input_change_pct",
        "sensitivity_ratio_available",
        "sensitivity_ratio_mean",
        "rank_by_abs_mean_effect",
        "rank_by_abs_sensitivity_ratio",
    ]

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="all_results", index=False)
        df[compact_cols].to_excel(writer, sheet_name="summary_compact", index=False)

    print(f"Saved sensitivity table: {path}")
    return path


# =============================================================================
# 5. BORNHOFT-STYLE HEATMAP AND TORNADO PLOTS
# =============================================================================

def make_normalized_sensitivity_heatmap_excel(df: pd.DataFrame) -> Path:
    output_dir = RESULTS_ROOT / "_normalized_sensitivity_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    path = output_dir / f"Normalized_TC_sensitivity_heatmaps_{TIMESTAMP}.xlsx"

    subset = df[
        (df["indicator"] == "cumulative_full_boundary_2014_2050")
        & (df["sensitivity_ratio_available"] == True)
        & (df["output"].isin(MAIN_SINKS_FOR_HEATMAP))
    ].copy()

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        subset.to_excel(writer, sheet_name="source_data", index=False)

        for scenario in sorted(subset["scenario"].unique()):
            ssub = subset[subset["scenario"] == scenario].copy()
            if ssub.empty:
                continue

            matrix = ssub.pivot_table(
                index="variant",
                columns="output_label",
                values="sensitivity_ratio_mean",
                aggfunc="mean",
            )
            matrix.to_excel(writer, sheet_name=f"{scenario}_S_matrix")

            matrix_abs = matrix.abs()
            matrix_abs.to_excel(writer, sheet_name=f"{scenario}_abs_S_matrix")

    print(f"Saved normalized sensitivity heatmap workbook: {path}")
    return path


def plot_heatmaps(df: pd.DataFrame) -> None:
    output_dir = RESULTS_ROOT / "_normalized_sensitivity_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    subset = df[
        (df["indicator"] == "cumulative_full_boundary_2014_2050")
        & (df["sensitivity_ratio_available"] == True)
        & (df["output"].isin(MAIN_SINKS_FOR_HEATMAP))
    ].copy()

    for scenario in sorted(subset["scenario"].unique()):
        ssub = subset[subset["scenario"] == scenario].copy()
        if ssub.empty:
            continue

        matrix = ssub.pivot_table(
            index="variant",
            columns="output_label",
            values="sensitivity_ratio_mean",
            aggfunc="mean",
        )

        # Stable column order.
        ordered_cols = [SINK_LABELS[o] for o in MAIN_SINKS_FOR_HEATMAP if SINK_LABELS[o] in matrix.columns]
        matrix = matrix[ordered_cols]

        fig_w = max(8, 0.55 * len(matrix.columns) + 3.5)
        fig_h = max(3.5, 0.45 * len(matrix.index) + 2.2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        vals = matrix.to_numpy(dtype=float)
        vmax = np.nanmax(np.abs(vals)) if np.isfinite(vals).any() else 1.0
        vmax = max(vmax, 1.0)

        im = ax.imshow(vals, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")

        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(matrix.index)))
        ax.set_yticklabels(matrix.index)
        ax.set_xlabel("Sink compartment")
        ax.set_ylabel("Scalar TC variant")
        ax.set_title(f"{scenario}: Normalized local sensitivity coefficient")

        for i in range(vals.shape[0]):
            for j in range(vals.shape[1]):
                val = vals[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("S = (%Δ mean cumulative sink flow) / (%Δ TC)")

        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"normalized_TC_sensitivity_heatmap_{scenario}.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)


def plot_tornado_rankings(df: pd.DataFrame) -> None:
    output_dir = RESULTS_ROOT / "_normalized_sensitivity_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    subset = df[
        (df["indicator"] == "cumulative_full_boundary_2014_2050")
        & (df["sensitivity_ratio_available"] == True)
        & (df["output"].isin(MAIN_SINKS_FOR_TORNADO))
    ].copy()

    if TORNADO_ONLY_ABS_20_PERCENT:
        subset = subset[np.isclose(subset["input_change_pct"].abs(), 20.0)].copy()

    for scenario in sorted(subset["scenario"].unique()):
        ssub = subset[subset["scenario"] == scenario].copy()
        if ssub.empty:
            continue

        # Rank variants by the maximum absolute S across selected sinks.
        rows = []
        for variant, group in ssub.groupby("variant"):
            idx = group["sensitivity_ratio_mean"].abs().idxmax()
            dominant = group.loc[idx]
            rows.append({
                "variant": variant,
                "max_abs_S": float(abs(dominant["sensitivity_ratio_mean"])),
                "dominant_signed_S": float(dominant["sensitivity_ratio_mean"]),
                "dominant_sink": dominant["output_label"],
                "input_change_pct": dominant["input_change_pct"],
            })

        rank_df = pd.DataFrame(rows).sort_values("max_abs_S", ascending=True)
        if rank_df.empty:
            continue

        fig_h = max(3.2, 0.45 * len(rank_df) + 1.8)
        fig, ax = plt.subplots(figsize=(8.0, fig_h))

        y = np.arange(len(rank_df))
        ax.barh(y, rank_df["dominant_signed_S"])
        ax.axvline(0, color="black", linewidth=0.8)

        labels = [
            f"{row.variant} ({row.dominant_sink})"
            for row in rank_df.itertuples(index=False)
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Dominant normalized sensitivity coefficient, S")
        ax.set_title(
            f"{scenario}: tornado ranking of scalar TC variants"
            + (" (±20% only)" if TORNADO_ONLY_ABS_20_PERCENT else "")
        )

        for yi, val in zip(y, rank_df["dominant_signed_S"]):
            ax.text(
                val + (0.02 if val >= 0 else -0.02),
                yi,
                f"{val:.2f}",
                va="center",
                ha="left" if val >= 0 else "right",
                fontsize=8,
            )

        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"scalar_TC_tornado_{scenario}.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)

        rank_path = output_dir / f"scalar_TC_tornado_ranking_{scenario}.csv"
        rank_df.sort_values("max_abs_S", ascending=False).to_csv(rank_path, index=False)
        print(f"Saved: {rank_path}")



# =============================================================================
# 6. ANNUAL SENSITIVITY TABLES AND CUMULATIVE BUBBLE MAPS
# =============================================================================

def _column_label_for_sink_year(sink: str, year: int) -> str:
    return f"{SINK_LABELS.get(sink, sink)} {year}"


def _ordered_scalar_variants(df: pd.DataFrame, scenario: str) -> List[str]:
    sub = df[(df["scenario"] == scenario) & (df["sensitivity_ratio_available"] == True)].copy()
    if sub.empty:
        return []

    # Keep a deterministic and readable order: direct recycling first when present,
    # then tire energy recovery and tire recycling variants.
    preferred = [
        "BattDirRec_high",
        "BattDirRec_low",
        "TiresEnRec_high",
        "TiresEnRec_low",
        "TiresRec_high",
        "TiresRec_low",
    ]
    present = list(sub["variant"].dropna().astype(str).unique())
    ordered = [v for v in preferred if v in present]
    ordered += sorted([v for v in present if v not in ordered])
    return ordered


def make_annual_sensitivity_tables(df: pd.DataFrame) -> Path:
    """
    Create annual sensitivity tables for selected sinks and years.

    This table is separate from the cumulative bubble map. It focuses only on:
        Landfill, End-of-Waste alternatives, Unknown, and Elimination
    for:
        2030 and 2050 annual mass flows.

    The workbook contains, for each scenario:
      1) mean relative change (%) for all non-baseline variants;
      2) absolute change (t/y) for all non-baseline variants;
      3) normalized S coefficients for scalar TC variants only;
      4) baseline annual mean mass flows (t/y).
    """
    output_dir = RESULTS_ROOT / "_annual_sensitivity_tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"Annual_sink_sensitivity_tables_2030_2050_{TIMESTAMP}.xlsx"

    indicators = [f"year_{year}" for year in ANNUAL_TABLE_YEARS]
    subset = df[
        (df["indicator"].isin(indicators))
        & (df["output"].isin(ANNUAL_TABLE_SINKS))
    ].copy()

    if subset.empty:
        raise ValueError("No annual data found for ANNUAL_TABLE_YEARS and ANNUAL_TABLE_SINKS.")

    def columns_for_table():
        cols = []
        for year in ANNUAL_TABLE_YEARS:
            for sink in ANNUAL_TABLE_SINKS:
                cols.append((sink, year, _column_label_for_sink_year(sink, year)))
        return cols

    def build_base_row(ssub: pd.DataFrame) -> Dict[str, Any]:
        row = {"Variant": "Basic scenario (t/y)", "Input change (%)": np.nan}
        for sink, year, col_label in columns_for_table():
            vals = ssub[(ssub["output"] == sink) & (ssub["indicator"] == f"year_{year}")]["mean_base"].dropna()
            row[col_label] = float(vals.iloc[0]) if len(vals) else np.nan
        return row

    def ordered_all_variants(ssub: pd.DataFrame) -> List[str]:
        preferred = [
            "BattWDism_dynamic",
            "BattW2Lf",
            "BattRec_Unk2Lf",
            "BattW_RecRes2Lf",
            "BattDirRec_high",
            "BattDirRec_low",
            "TiresEnRec_high",
            "TiresEnRec_low",
            "TiresRec_high",
            "TiresRec_low",
        ]
        present = list(ssub["variant"].dropna().astype(str).unique())
        ordered = [v for v in preferred if v in present]
        ordered += sorted([v for v in present if v not in ordered])
        return ordered

    def build_wide_table(
        scenario: str,
        value_col: str,
        row_label: str,
        scalar_only: bool = False,
        include_base_row: bool = True,
    ) -> pd.DataFrame:
        ssub = subset[subset["scenario"] == scenario].copy()
        if scalar_only:
            row_variants = _ordered_scalar_variants(ssub, scenario)
        else:
            row_variants = ordered_all_variants(ssub)

        rows = []
        if include_base_row:
            rows.append(build_base_row(ssub))

        for variant in row_variants:
            vsub = ssub[ssub["variant"] == variant].copy()
            if vsub.empty:
                continue

            row = {"Variant": variant}
            input_vals = vsub["input_change_pct"].dropna().unique()
            row["Input change (%)"] = float(input_vals[0]) if len(input_vals) else np.nan

            for sink, year, col_label in columns_for_table():
                vals = vsub[(vsub["output"] == sink) & (vsub["indicator"] == f"year_{year}")][value_col].dropna()
                row[col_label] = float(vals.iloc[0]) if len(vals) else np.nan
            rows.append(row)

        out = pd.DataFrame(rows)
        out.insert(0, "Scenario", scenario)
        out.insert(2, "Metric", row_label)
        return out

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        subset.to_excel(writer, sheet_name="source_long_all_variants", index=False)

        for scenario in sorted(subset["scenario"].dropna().unique()):
            base_table = build_wide_table(
                scenario,
                value_col="mean_base",
                row_label="Baseline annual mean mass flow (t/y)",
                scalar_only=False,
                include_base_row=True,
            ).iloc[[0]].copy()

            rel_all = build_wide_table(
                scenario,
                value_col="relative_change_of_means_pct",
                row_label="Mean relative change vs baseline (%)",
                scalar_only=False,
                include_base_row=True,
            )

            abs_all = build_wide_table(
                scenario,
                value_col="absolute_change_of_means",
                row_label="Absolute change vs baseline (t/y)",
                scalar_only=False,
                include_base_row=True,
            )

            s_scalar = build_wide_table(
                scenario,
                value_col="sensitivity_ratio_mean",
                row_label="S = (%Delta annual mean sink flow) / (%Delta TC)",
                scalar_only=True,
                include_base_row=True,
            )

            base_table.to_excel(writer, sheet_name=f"{scenario}_base_annual_t", index=False)
            rel_all.to_excel(writer, sheet_name=f"{scenario}_rel_pct_all", index=False)
            abs_all.to_excel(writer, sheet_name=f"{scenario}_abs_t_all", index=False)
            s_scalar.to_excel(writer, sheet_name=f"{scenario}_S_scalar_TC", index=False)

    print(f"Saved annual sensitivity tables: {path}")
    return path


def _colorblind_diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "blue_gray_orange",
        ["#3B6FB6", "#F7F7F7", "#D55E00"],
        N=256,
    )


def _bubble_sizes(values: np.ndarray, min_size: float = 35, max_size: float = 950) -> np.ndarray:
    vals = np.abs(np.asarray(values, dtype=float))
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return np.full_like(vals, min_size)
    vmax = np.nanpercentile(finite, 95)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = np.nanmax(finite)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    scaled = np.minimum(vals, vmax) / vmax
    return min_size + (max_size - min_size) * np.sqrt(scaled)


def _ordered_variants_for_bubble(df: pd.DataFrame, scenario: str) -> List[str]:
    preferred = [
        "BattWDism_dynamic",
        "BattW2Lf",
        "BattRec_Unk2Lf",
        "BattW_RecRes2Lf",
        "BattDirRec_high",
        "BattDirRec_low",
        "TiresEnRec_high",
        "TiresEnRec_low",
        "TiresRec_high",
        "TiresRec_low",
    ]
    present = list(df[df["scenario"] == scenario]["variant"].dropna().astype(str).unique())
    ordered = [v for v in preferred if v in present]
    ordered += sorted([v for v in present if v not in ordered])
    return ordered


def plot_cumulative_sink_bubble_map(df: pd.DataFrame) -> None:
    """
    Create the main bubble map for cumulative sink changes until 2050.

    This is separate from the annual table. It uses all selected sink compartments
    and all non-baseline sensitivity variants available in the results.

    Encodings:
        color = signed mean relative change in cumulative sink flow vs baseline (%)
        size  = absolute mean relative change in cumulative sink flow vs baseline (%)

    No filled/open robustness encoding is used in this figure. Q5-Q95 values remain
    available in the Excel comparison table for uncertainty checking.
    """
    from matplotlib.colors import SymLogNorm, TwoSlopeNorm

    output_dir = RESULTS_ROOT / "_cumulative_sink_bubble_maps"
    output_dir.mkdir(parents=True, exist_ok=True)

    indicator = "cumulative_full_boundary_2014_2050"
    plot_df = df[
        (df["indicator"] == indicator)
        & (df["output"].isin(MAIN_SINKS_FOR_HEATMAP))
    ].copy()

    if plot_df.empty:
        print("No cumulative full-boundary data found for the cumulative bubble map.")
        return

    plot_df["effect"] = pd.to_numeric(plot_df["relative_change_of_means_pct"], errors="coerce")
    plot_df = plot_df[np.isfinite(plot_df["effect"])].copy()
    if plot_df.empty:
        print("No finite mean relative changes found for the cumulative bubble map.")
        return

    scenarios = [s for s in SCENARIOS_TO_RUN if s in set(plot_df["scenario"].astype(str))]
    if not scenarios:
        scenarios = sorted(plot_df["scenario"].dropna().astype(str).unique())

    sinks = [s for s in MAIN_SINKS_FOR_HEATMAP if s in set(plot_df["output"].astype(str))]
    cmap = _colorblind_diverging_cmap()

    max_abs = np.nanmax(np.abs(plot_df["effect"].to_numpy(dtype=float)))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    if max_abs > 100:
        norm = SymLogNorm(linthresh=10, vmin=-max_abs, vmax=max_abs, base=10)
    else:
        norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    fig_height = max(7.0, 3.6 * len(scenarios))
    fig_width = 10.5
    fig, axes = plt.subplots(nrows=len(scenarios), ncols=1, figsize=(fig_width, fig_height), sharey=True, constrained_layout=False)
    if len(scenarios) == 1:
        axes = [axes]

    last_scatter = None
    panel_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for ax, scenario, letter in zip(axes, scenarios, panel_letters):
        sub = plot_df[plot_df["scenario"].astype(str) == scenario].copy()
        variants = _ordered_variants_for_bubble(sub, scenario)
        if not variants:
            ax.set_visible(False)
            continue

        sub = sub[sub["variant"].isin(variants)].copy()
        sub["x"] = sub["variant"].map({v: i for i, v in enumerate(variants)})
        sub["y"] = sub["output"].map({s: i for i, s in enumerate(sinks)})
        sub = sub[np.isfinite(sub["x"]) & np.isfinite(sub["y"])].copy()

        values = sub["effect"].to_numpy(dtype=float)
        last_scatter = ax.scatter(
            sub["x"], sub["y"], s=_bubble_sizes(values), c=values,
            cmap=cmap, norm=norm, edgecolors="0.25", linewidths=0.8, alpha=0.95,
        )

        for _, row in sub.iterrows():
            val = float(row["effect"])
            if abs(val) >= PRACTICAL_EFFECT_THRESHOLD_PCT:
                ax.annotate(
                    f"{val:.0f}%", xy=(row["x"], row["y"]), xytext=(4, 4), textcoords="offset points",
                    fontsize=7, ha="left", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                )

        ax.set_xticks(range(len(variants)))
        ax.set_xticklabels(variants, rotation=45, ha="right")
        ax.set_yticks(range(len(sinks)))
        ax.set_yticklabels([SINK_LABELS.get(s, s) for s in sinks])
        ax.invert_yaxis()
        ax.set_xlabel("Sensitivity variant")
        ax.set_ylabel("Sink compartment")
        ax.grid(True, linestyle=":", linewidth=0.6, color="0.80")
        ax.set_axisbelow(True)
        ax.text(1.02, 0.98, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")
        ax.set_title(f"{scenario} — cumulative 2014-2050")

    if last_scatter is not None:
        fig.subplots_adjust(left=0.24, right=0.83, bottom=0.08, top=0.94, hspace=0.55)
        cax = fig.add_axes([0.86, 0.18, 0.025, 0.68])
        cbar = fig.colorbar(last_scatter, cax=cax)
        cbar.set_label("Mean relative change vs baseline (%)")

    stem = "cumulative_sink_bubble_map_2014_2050_all_variants"
    for ext in ["png", "pdf", "svg"]:
        out = output_dir / f"{stem}.{ext}"
        fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)

    for scenario in scenarios:
        sub = plot_df[plot_df["scenario"].astype(str) == scenario].copy()
        variants = _ordered_variants_for_bubble(sub, scenario)
        if not variants:
            continue
        sub = sub[sub["variant"].isin(variants)].copy()
        sub["x"] = sub["variant"].map({v: i for i, v in enumerate(variants)})
        sub["y"] = sub["output"].map({s: i for i, s in enumerate(sinks)})
        sub = sub[np.isfinite(sub["x"]) & np.isfinite(sub["y"])].copy()
        values = sub["effect"].to_numpy(dtype=float)

        fig, ax = plt.subplots(figsize=(max(7.5, 0.75 * len(variants) + 2.8), max(5.5, 0.45 * len(sinks) + 2.0)))
        sc = ax.scatter(
            sub["x"], sub["y"], s=_bubble_sizes(values), c=values,
            cmap=cmap, norm=norm, edgecolors="0.25", linewidths=0.8, alpha=0.95,
        )
        for _, row in sub.iterrows():
            val = float(row["effect"])
            if abs(val) >= PRACTICAL_EFFECT_THRESHOLD_PCT:
                ax.annotate(
                    f"{val:.0f}%", xy=(row["x"], row["y"]), xytext=(4, 4), textcoords="offset points",
                    fontsize=7, ha="left", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                )
        ax.set_xticks(range(len(variants)))
        ax.set_xticklabels(variants, rotation=45, ha="right")
        ax.set_yticks(range(len(sinks)))
        ax.set_yticklabels([SINK_LABELS.get(s, s) for s in sinks])
        ax.invert_yaxis()
        ax.set_xlabel("Sensitivity variant")
        ax.set_ylabel("Sink compartment")
        ax.set_title(f"{scenario} — cumulative 2014-2050")
        ax.grid(True, linestyle=":", linewidth=0.6, color="0.80")
        ax.set_axisbelow(True)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Mean relative change vs baseline (%)")
        fig.tight_layout()
        stem = f"cumulative_sink_bubble_map_{scenario}_2014_2050"
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"{stem}.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)


# =============================================================================
# 7. DYNAMIC BATTERY COLLECTION/DISMANTLING TRAJECTORY
# =============================================================================

def safe_eval_numeric(expr: str) -> Optional[float]:
    """
    Safely evaluate simple numeric expressions such as:
        0.4
        0.4*0.5
        1
    """
    expr = expr.strip()
    try:
        node = ast.parse(expr, mode="eval")
    except SyntaxError:
        return None

    allowed = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Num,
        ast.Constant,
        ast.Mult,
        ast.Div,
        ast.Add,
        ast.Sub,
        ast.Pow,
        ast.USub,
        ast.UAdd,
    )
    if not all(isinstance(n, allowed) for n in ast.walk(node)):
        return None

    try:
        return float(eval(compile(node, "<expr>", "eval"), {"__builtins__": {}}, {}))
    except Exception:
        return None


def extract_bw_s_tc_values(model_file: Path) -> Dict[int, float]:
    """
    Extract BW_S_y0 ... BW_S_y36 central TC values from a model file.

    Supports:
        BW_S_y0 = cp.TransferConstant(0.4)
        BW_S_y1 = cp.TransferDistribution(tf.TriangTruncDet,
                   [0.4*0.5, 0.4, 0.4*1.5, 1, 0, 1])

    The central value is taken as:
        - TransferConstant argument, or
        - second value inside the distribution list, because this is the mode/central value.
    """
    if not model_file.exists():
        raise FileNotFoundError(f"Model file not found: {model_file}")

    values: Dict[int, float] = {}
    text = model_file.read_text(encoding="utf-8", errors="ignore").splitlines()

    for line in text:
        if "BW_S_y" not in line or "=" not in line:
            continue

        m = re.search(r"BW_S_y(\d+)\s*=", line)
        if not m:
            continue

        idx = int(m.group(1))
        year = START_YEAR + idx

        # Case 1: TransferConstant(...)
        m_const = re.search(r"TransferConstant\s*\(\s*([^)]+)\)", line)
        if m_const:
            val = safe_eval_numeric(m_const.group(1))
            if val is not None:
                values[year] = val
                continue

        # Case 2: list inside brackets [...]
        m_list = re.search(r"\[([^\]]+)\]", line)
        if m_list:
            parts = [p.strip() for p in m_list.group(1).split(",")]
            if len(parts) >= 2:
                # Central/mode value is the second parameter.
                val = safe_eval_numeric(parts[1])
                if val is not None:
                    values[year] = val
                    continue

    if not values:
        raise ValueError(f"No BW_S_y values found in model file: {model_file}")

    missing = [year for year in YEARS if year not in values]
    if missing:
        raise ValueError(
            f"Missing BW_S_y values for years {missing} in {model_file}. "
            f"Found years: {sorted(values)}"
        )

    return values


def discover_s1_model_files() -> Tuple[Path, Path]:
    s1_source = RUNS_SOURCE_ROOT / "S1_SA"
    base = s1_source / "DPMFA_graphene_TBC_model_S1_base.py"
    dynamic = s1_source / "DPMFA_graphene_TBC_model_S1_BattWDism_dynamic.py"
    return base, dynamic


def make_dynamic_battery_curves() -> None:
    output_dir = RESULTS_ROOT / "_normalized_sensitivity_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    base_model, dynamic_model = discover_s1_model_files()

    try:
        tc_base = extract_bw_s_tc_values(base_model)
        tc_dynamic = extract_bw_s_tc_values(dynamic_model)
    except Exception as exc:
        print(f"Could not extract dynamic battery TC values: {exc}")
        return

    tc_rows = []
    for year in YEARS:
        b = tc_base[year]
        d = tc_dynamic[year]
        input_change_pct = ((d - b) / b * 100.0) if b != 0 else np.nan
        tc_rows.append({
            "year": year,
            "tc_base": b,
            "tc_dynamic_variant": d,
            "input_change_pct": input_change_pct,
        })

    tc_df = pd.DataFrame(tc_rows)

    tc_path = output_dir / "S1_BattWDism_dynamic_extracted_TC_trajectory.xlsx"
    tc_df.to_excel(tc_path, index=False)
    print(f"Saved extracted TC trajectory: {tc_path}")

    # TC trajectory figure.
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(tc_df["year"], tc_df["tc_base"], marker="o", linewidth=1.7, label="S1 base")
    ax.plot(tc_df["year"], tc_df["tc_dynamic_variant"], marker="o", linewidth=1.7, label="S1 BattWDism_dynamic")
    ax.set_xlabel("Year")
    ax.set_ylabel("Battery collection/dismantling TC")
    ax.set_title("S1 battery collection/dismantling transfer-coefficient trajectory")
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle=":", linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        out = output_dir / f"S1_BattWDism_dynamic_TC_trajectory.{ext}"
        fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)

    # Annual output response and normalized annual S(t).
    base_folder = RESULTS_ROOT / "S1" / "base"
    dynamic_folder = RESULTS_ROOT / "S1" / "BattWDism_dynamic"

    if not base_folder.exists() or not dynamic_folder.exists():
        print("S1 base or BattWDism_dynamic result folder not found; skipping annual dynamic output curves.")
        return

    response_rows = []

    for sink, patterns in SINK_OUTPUTS.items():
        base_file = find_latest_matching_file(base_folder, patterns)
        dyn_file = find_latest_matching_file(dynamic_folder, patterns)

        if base_file is None or dyn_file is None:
            continue

        base_matrix = read_run_by_year_matrix(base_file)
        dyn_matrix = read_run_by_year_matrix(dyn_file)

        base_mean_by_year = np.nanmean(base_matrix.to_numpy(dtype=float), axis=0)
        dyn_mean_by_year = np.nanmean(dyn_matrix.to_numpy(dtype=float), axis=0)

        for i, year in enumerate(YEARS):
            y_base = base_mean_by_year[i]
            y_dyn = dyn_mean_by_year[i]
            output_change_pct = ((y_dyn - y_base) / y_base * 100.0) if y_base != 0 else np.nan

            input_change_pct = tc_df.loc[tc_df["year"] == year, "input_change_pct"].iloc[0]
            normalized_s_t = output_change_pct / input_change_pct if np.isfinite(output_change_pct) and np.isfinite(input_change_pct) and input_change_pct != 0 else np.nan

            response_rows.append({
                "year": year,
                "sink": sink,
                "sink_label": SINK_LABELS.get(sink, sink),
                "mean_base": y_base,
                "mean_dynamic_variant": y_dyn,
                "absolute_change": y_dyn - y_base,
                "relative_change_pct": output_change_pct,
                "input_change_pct": input_change_pct,
                "normalized_S_t": normalized_s_t,
            })

    response_df = pd.DataFrame(response_rows)
    response_path = output_dir / "S1_BattWDism_dynamic_annual_output_response_and_S_t.xlsx"
    response_df.to_excel(response_path, index=False)
    print(f"Saved dynamic annual response table: {response_path}")

    # Plot output-side annual response for key sinks.
    plot_df = response_df[response_df["sink"].isin(KEY_SINKS_DYNAMIC_CURVE)].copy()

    if not plot_df.empty:
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        for sink in KEY_SINKS_DYNAMIC_CURVE:
            ssub = plot_df[plot_df["sink"] == sink]
            if ssub.empty:
                continue
            ax.plot(ssub["year"], ssub["relative_change_pct"], marker="o", linewidth=1.5, label=SINK_LABELS.get(sink, sink))

        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Year")
        ax.set_ylabel("Annual mean output change vs S1 base (%)")
        ax.set_title("S1 BattWDism_dynamic: annual output response")
        ax.grid(True, linestyle=":", linewidth=0.6)
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"S1_BattWDism_dynamic_annual_output_response.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        for sink in KEY_SINKS_DYNAMIC_CURVE:
            ssub = plot_df[plot_df["sink"] == sink]
            if ssub.empty:
                continue
            ax.plot(ssub["year"], ssub["normalized_S_t"], marker="o", linewidth=1.5, label=SINK_LABELS.get(sink, sink))

        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Year")
        ax.set_ylabel("Normalized annual sensitivity coefficient S(t)")
        ax.set_title("S1 BattWDism_dynamic: normalized annual sensitivity")
        ax.grid(True, linestyle=":", linewidth=0.6)
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"S1_BattWDism_dynamic_normalized_S_t.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)


# =============================================================================
# 8. CROSS-SCENARIO COMPARISON: S1 DYNAMIC VS S3 BASE
# =============================================================================

def make_cross_scenario_comparison() -> None:
    """
    Compare S3_base against S1_BattWDism_dynamic.

    Reference scenario:
        S1_BattWDism_dynamic

    Relative difference:
        (S3_base - S1_BattWDism_dynamic) / S1_BattWDism_dynamic * 100

    This is not a local sensitivity coefficient. It is a diagnostic comparison
    of mean sink outputs under two cases that share an idealized battery
    collection/recycling trajectory but differ in broader scenario assumptions.
    """
    output_dir = RESULTS_ROOT / "_cross_scenario_S1dynamic_vs_S3base"
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_folder = RESULTS_ROOT / "S1" / "BattWDism_dynamic"
    comp_folder = RESULTS_ROOT / "S3" / "base"

    if not ref_folder.exists() or not comp_folder.exists():
        print("Cross-scenario folders not found; skipping cross-scenario comparison.")
        print(f"Reference folder: {ref_folder}")
        print(f"Comparison folder: {comp_folder}")
        return

    cumulative_rows = []
    annual_rows = []

    for sink, patterns in SINK_OUTPUTS.items():
        ref_file = find_latest_matching_file(ref_folder, patterns)
        comp_file = find_latest_matching_file(comp_folder, patterns)

        if ref_file is None or comp_file is None:
            continue

        ref_matrix = read_run_by_year_matrix(ref_file)
        comp_matrix = read_run_by_year_matrix(comp_file)

        # Full cumulative period, reference = S1 dynamic.
        ref_cum = extract_indicator_values(
            ref_matrix, "cumulative", cumulative_from_year=2014, cumulative_to_year=2050
        )
        comp_cum = extract_indicator_values(
            comp_matrix, "cumulative", cumulative_from_year=2014, cumulative_to_year=2050
        )

        ref_mean = safe_nanmean(ref_cum)
        comp_mean = safe_nanmean(comp_cum)
        abs_diff = comp_mean - ref_mean
        rel_diff_pct = (abs_diff / ref_mean * 100.0) if ref_mean != 0 else np.nan

        cumulative_rows.append({
            "sink": sink,
            "sink_label": SINK_LABELS.get(sink, sink),
            "reference_case": "S1_BattWDism_dynamic",
            "comparison_case": "S3_base",
            "reference_mean_cumulative_t": ref_mean,
            "comparison_mean_cumulative_t": comp_mean,
            "absolute_difference_t": abs_diff,
            "relative_difference_pct_vs_S1_BattWDism_dynamic": rel_diff_pct,
        })

        ref_mean_by_year = np.nanmean(ref_matrix.to_numpy(dtype=float), axis=0)
        comp_mean_by_year = np.nanmean(comp_matrix.to_numpy(dtype=float), axis=0)

        for i, year in enumerate(YEARS):
            ref_y = ref_mean_by_year[i]
            comp_y = comp_mean_by_year[i]
            abs_y = comp_y - ref_y
            rel_y = (abs_y / ref_y * 100.0) if ref_y != 0 else np.nan

            annual_rows.append({
                "year": year,
                "sink": sink,
                "sink_label": SINK_LABELS.get(sink, sink),
                "reference_case": "S1_BattWDism_dynamic",
                "comparison_case": "S3_base",
                "reference_mean_t": ref_y,
                "comparison_mean_t": comp_y,
                "absolute_difference_t": abs_y,
                "relative_difference_pct_vs_S1_BattWDism_dynamic": rel_y,
            })

    cumulative_df = pd.DataFrame(cumulative_rows)
    annual_df = pd.DataFrame(annual_rows)

    xlsx_path = output_dir / f"cross_scenario_S1_BattWDism_dynamic_vs_S3_base_{TIMESTAMP}.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        cumulative_df.to_excel(writer, sheet_name="cumulative_2014_2050", index=False)
        annual_df.to_excel(writer, sheet_name="annual_mean", index=False)
    print(f"Saved cross-scenario workbook: {xlsx_path}")

    if cumulative_df.empty:
        return

    # Point-range plot: reference and comparison cumulative mean masses.
    plot_df = cumulative_df.sort_values("reference_mean_cumulative_t", ascending=True)
    y = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(8.0, max(4.8, 0.38 * len(plot_df) + 1.8)))
    ax.hlines(
        y,
        plot_df["reference_mean_cumulative_t"],
        plot_df["comparison_mean_cumulative_t"],
        color="0.65",
        linewidth=1.2,
    )
    ax.scatter(plot_df["reference_mean_cumulative_t"], y, marker="o", label="S1 BattWDism_dynamic")
    ax.scatter(plot_df["comparison_mean_cumulative_t"], y, marker="s", label="S3 base")

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["sink_label"])
    ax.set_xlabel("Mean cumulative sink mass, 2014–2050 (t)")
    ax.set_ylabel("Sink compartment")
    ax.set_title("Cross-scenario comparison: S3 base vs S1 BattWDism_dynamic")
    ax.grid(True, axis="x", linestyle=":", linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        out = output_dir / f"cross_scenario_S1_BattWDism_dynamic_vs_S3_base_point_range.{ext}"
        fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)

    # Annual relative-difference curves for key sinks.
    plot_annual = annual_df[annual_df["sink"].isin(KEY_SINKS_DYNAMIC_CURVE)].copy()
    if not plot_annual.empty:
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        for sink in KEY_SINKS_DYNAMIC_CURVE:
            ssub = plot_annual[plot_annual["sink"] == sink]
            if ssub.empty:
                continue
            ax.plot(
                ssub["year"],
                ssub["relative_difference_pct_vs_S1_BattWDism_dynamic"],
                marker="o",
                linewidth=1.5,
                label=SINK_LABELS.get(sink, sink),
            )

        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Year")
        ax.set_ylabel("S3 base relative difference vs S1 BattWDism_dynamic (%)")
        ax.set_title("Annual mean sink differences, reference = S1 BattWDism_dynamic")
        ax.grid(True, linestyle=":", linewidth=0.6)
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"cross_scenario_annual_relative_difference_reference_S1dynamic.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            print(f"Saved: {out}")
        plt.close(fig)


# =============================================================================
# 9. MAIN
# =============================================================================

def main() -> None:
    ensure_dirs()

    print("\n" + "=" * 80)
    print("MASTER DPMFA LOCAL SENSITIVITY WORKFLOW")
    print("=" * 80)
    print(f"Project root:        {PROJECT_ROOT}")
    print(f"Runner source root:  {RUNS_SOURCE_ROOT}")
    print(f"New results root:    {RESULTS_ROOT}")
    print(f"Scenarios:           {SCENARIOS_TO_RUN}")
    print("=" * 80)

    jobs = discover_runner_jobs()
    if not jobs:
        raise ValueError(f"No runner jobs found in {RUNS_SOURCE_ROOT}")

    print("\nDiscovered jobs:")
    for job in jobs:
        print(
            f"  - {job['scenario']} / {job['variant']} "
            f"(input_change_pct={job['input_change_pct']})"
        )

    if RUN_RUNNERS_FIRST:
        run_runner_jobs(jobs)
    else:
        print("\nRUN_RUNNERS_FIRST=False; using existing organized result folders.")

    print("\nBuilding sensitivity-comparison table...")
    sensitivity_df = build_sensitivity_table(jobs)
    save_sensitivity_table(sensitivity_df)

    print("\nCreating annual sensitivity tables for 2030 and 2050...")
    make_annual_sensitivity_tables(sensitivity_df)

    print("\nCreating cumulative sink bubble map for 2014-2050...")
    plot_cumulative_sink_bubble_map(sensitivity_df)

    print("\nCreating normalized sensitivity heatmap matrices and plots...")
    make_normalized_sensitivity_heatmap_excel(sensitivity_df)
    plot_heatmaps(sensitivity_df)

    print("\nCreating tornado plots for scalar TC variants...")
    plot_tornado_rankings(sensitivity_df)

    print("\nCreating S1 BattWDism_dynamic annual curves...")
    make_dynamic_battery_curves()

    print("\nCreating cross-scenario comparison...")
    make_cross_scenario_comparison()

    print("\nDone.")
    print(f"All outputs saved under: {RESULTS_ROOT}")


if __name__ == "__main__":
    main()
