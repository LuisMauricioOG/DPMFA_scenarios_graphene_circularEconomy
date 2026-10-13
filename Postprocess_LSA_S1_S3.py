# -*- coding: utf-8 -*-
"""
Post-process local sensitivity-analysis results.

This script DOES NOT rerun the DPMFA simulations. It reads the latest
DPMFA_sensitivity_comparison_ALL_*.xlsx workbook created by the LSA workflow and creates:

1) one cumulative bubble map for 2014-2050 using all sink compartments and all
   non-baseline sensitivity variants, excluding any variants listed in EXCLUDE_VARIANTS;

2) a supplementary annual sensitivity table for annual mass flows in 2030 and 2050
   for Landfill, End-of-Waste alternatives, Unknown, and Elimination;

3) the exact long-format data used for the figure and tables.

Figure encoding:
    color = signed mean relative change vs corresponding scenario baseline (%)
            using a colorblind-friendly blue-neutral-orange diverging palette
    size  = absolute mean relative change vs corresponding scenario baseline (%)

The Q5-Q95 robustness columns remain in the source table, but the figure uses
only bubble size and color.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm, TwoSlopeNorm


# =============================================================================
# USER SETTINGS
# =============================================================================

PROJECT_ROOT = Path(
    r"C:\Users\ortl\Documents\Empa\Articles\pDPMFA_AdMa_EU_CE"
)

# Leave as None to automatically use the latest workbook found.
# You can also paste a specific path, for example:
# SENSITIVITY_WORKBOOK = Path(r"C:\...\DPMFA_sensitivity_comparison_ALL_20260807_165000.xlsx")
SENSITIVITY_WORKBOOK: Optional[Path] = None

# Prefer today's dated result folder, but fall back to any *_SensitivityAnalysis_results folder.
RUN_DATE = datetime.now().strftime("%Y%m%d")
PREFERRED_RESULTS_ROOT = PROJECT_ROOT / f"{RUN_DATE}_SensitivityAnalysis_results"

# Exclude variants that should not appear anywhere after rerunning/cleaning.
EXCLUDE_VARIANTS = {
    "BattDirRec_Trap",
}

SCENARIOS_TO_PLOT = ["S1", "S3"]

CUMULATIVE_INDICATOR = "cumulative_full_boundary_2014_2050"

ALL_SINKS_ORDER = [
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

ALL_SINK_LABELS = {
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

ANNUAL_TABLE_YEARS = [2030, 2050]
ANNUAL_TABLE_SINKS = [
    "Landfill",
    "EoW_alternatives",
    "Unknown",
    "Elimination",
]

PRACTICAL_LABEL_THRESHOLD_PCT = 20.0
SAVE_DPI = 600

# Colorblind-friendly palette.
# Diverging values:
#   negative change -> blue
#   near zero       -> light neutral
#   positive change -> orange
# This avoids red-green contrasts.
COLOR_NEGATIVE = "#3B6FB6"
COLOR_NEUTRAL = "#F7F7F7"
COLOR_POSITIVE = "#D55E00"
EDGE_COLOR = "#333333"
GRID_COLOR = "#D0D0D0"

# Arial first; matplotlib will fall back if Arial is unavailable on the machine.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def find_latest_sensitivity_workbook() -> Path:
    if SENSITIVITY_WORKBOOK is not None:
        path = Path(SENSITIVITY_WORKBOOK)
        if not path.exists():
            raise FileNotFoundError(
                f"SENSITIVITY_WORKBOOK does not exist: {path}"
            )
        return path

    candidate_patterns = []

    # 1) Today/current result folder first.
    candidate_patterns.append(
        PREFERRED_RESULTS_ROOT
        / "_sensitivity_comparison_tables"
        / "DPMFA_sensitivity_comparison_ALL_*.xlsx"
    )

    # 2) Any dated sensitivity-analysis result folder.
    candidate_patterns.append(
        PROJECT_ROOT
        / "*_SensitivityAnalysis_results"
        / "_sensitivity_comparison_tables"
        / "DPMFA_sensitivity_comparison_ALL_*.xlsx"
    )

    candidates: List[Path] = []

    for pattern in candidate_patterns:
        candidates.extend(pattern.parent.glob(pattern.name))

    candidates = [
        p for p in candidates
        if p.exists() and p.is_file()
    ]

    if not candidates:
        raise FileNotFoundError(
            "No DPMFA_sensitivity_comparison_ALL_*.xlsx workbook found."
        )

    return max(
            candidates,
            key=lambda p: p.stat().st_mtime
        )

def infer_results_root_from_workbook(workbook: Path) -> Path:
    # Expected:
    # RESULTS_ROOT / "_sensitivity_comparison_tables" / workbook.xlsx
    if workbook.parent.name == "_sensitivity_comparison_tables":
        return workbook.parent.parent
    return PREFERRED_RESULTS_ROOT


# =============================================================================
# DATA CLEANING
# =============================================================================

def normalize_sink(value: Any) -> str:
    s = str(value).strip()
    low = s.lower().replace("_", " ").replace("-", " ")

    mapping = {
        "eow alternatives": "EoW_alternatives",
        "end of waste alternatives": "EoW_alternatives",
        "end of waste alternative": "EoW_alternatives",
        "reuse": "Reuse",
        "landfill": "Landfill",
        "elimination": "Elimination",
        "export": "Export",
        "unknown": "Unknown",
        "roadside": "Roadside",
        "stsoil": "STsoil",
        "st soil": "STsoil",
        "sludge treated soil": "STsoil",
        "nu soil": "NUsoil",
        "nusoil": "NUsoil",
        "natural urban soil": "NUsoil",
        "natural & urban soil": "NUsoil",
        "surfacewater": "Surfacewater",
        "surface water": "Surfacewater",
        "subsurface": "Subsurface",
    }

    return mapping.get(low, s)


def load_sensitivity_table(workbook: Path) -> pd.DataFrame:
    # Prefer all_results because it keeps every column.
    xls = pd.ExcelFile(workbook)
    sheet = "all_results" if "all_results" in xls.sheet_names else "summary_compact"
    df = pd.read_excel(workbook, sheet_name=sheet)

    required = [
        "scenario",
        "variant",
        "output",
        "indicator",
        "mean_base",
        "mean_variant",
        "absolute_change_of_means",
        "relative_change_of_means_pct",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"The workbook is missing required columns: {missing}. "
            f"Workbook: {workbook}, sheet: {sheet}"
        )

    df["scenario"] = df["scenario"].astype(str)
    df["variant"] = df["variant"].astype(str)
    df["output"] = df["output"].map(normalize_sink)
    df["output_label"] = df["output"].map(lambda x: ALL_SINK_LABELS.get(x, x))

    df = df[~df["variant"].isin(EXCLUDE_VARIANTS)].copy()

    numeric_cols = [
        "mean_base",
        "mean_variant",
        "absolute_change_of_means",
        "relative_change_of_means_pct",
        "input_change_pct",
        "sensitivity_ratio_mean",
        "q05_relative_change_pct",
        "q95_relative_change_pct",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


# =============================================================================
# ANNUAL TABLES
# =============================================================================

def wide_columns() -> List[tuple[str, int, str]]:
    cols = []
    for year in ANNUAL_TABLE_YEARS:
        for sink in ANNUAL_TABLE_SINKS:
            label = f"{ALL_SINK_LABELS.get(sink, sink)} {year}"
            cols.append((sink, year, label))
    return cols


def ordered_variants(df: pd.DataFrame, scenario: str, scalar_only: bool = False) -> List[str]:
    sub = df[df["scenario"] == scenario].copy()
    if scalar_only and "sensitivity_ratio_mean" in sub.columns:
        sub = sub[np.isfinite(sub["sensitivity_ratio_mean"])].copy()

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

    present = list(sub["variant"].dropna().astype(str).unique())
    out = [v for v in preferred if v in present]
    out += sorted([v for v in present if v not in out])
    return out


def make_wide_table(
    df: pd.DataFrame,
    scenario: str,
    value_col: str,
    metric_label: str,
    scalar_only: bool = False,
    include_baseline_row: bool = True,
) -> pd.DataFrame:
    ssub = df[df["scenario"] == scenario].copy()
    rows: List[Dict[str, Any]] = []

    if include_baseline_row:
        row = {
            "Scenario": scenario,
            "Metric": "Baseline annual mean mass flow (t/y)",
            "Variant": "Basic scenario",
            "Input change (%)": np.nan,
        }
        for sink, year, label in wide_columns():
            vals = ssub[
                (ssub["output"] == sink) &
                (ssub["indicator"] == f"year_{year}")
            ]["mean_base"].dropna()
            row[label] = float(vals.iloc[0]) if len(vals) else np.nan
        rows.append(row)

    for variant in ordered_variants(ssub, scenario, scalar_only=scalar_only):
        vsub = ssub[ssub["variant"] == variant].copy()
        if vsub.empty:
            continue

        row = {
            "Scenario": scenario,
            "Metric": metric_label,
            "Variant": variant,
            "Input change (%)": np.nan,
        }
        if "input_change_pct" in vsub.columns:
            vals = vsub["input_change_pct"].dropna().unique()
            if len(vals):
                row["Input change (%)"] = float(vals[0])

        for sink, year, label in wide_columns():
            vals = vsub[
                (vsub["output"] == sink) &
                (vsub["indicator"] == f"year_{year}")
            ][value_col].dropna()
            row[label] = float(vals.iloc[0]) if len(vals) else np.nan
        rows.append(row)

    return pd.DataFrame(rows)


def save_annual_tables(df: pd.DataFrame, output_root: Path) -> Path:
    output_dir = output_root / "_annual_sink_sensitivity_tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    annual_df = df[
        (df["indicator"].isin([f"year_{y}" for y in ANNUAL_TABLE_YEARS])) &
        (df["output"].isin(ANNUAL_TABLE_SINKS))
    ].copy()

    if annual_df.empty:
        raise ValueError(
            "No annual rows found for the requested years/sinks. "
            "Expected indicators like year_2030 and year_2050."
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = output_dir / f"annual_sink_sensitivity_2030_2050_{timestamp}.xlsx"

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        annual_df.to_excel(writer, sheet_name="source_long", index=False)

        for scenario in SCENARIOS_TO_PLOT:
            if scenario not in set(annual_df["scenario"].astype(str)):
                continue

            base = make_wide_table(
                annual_df, scenario,
                value_col="mean_base",
                metric_label="Baseline annual mean mass flow (t/y)",
                scalar_only=False,
                include_baseline_row=True,
            ).iloc[[0]]

            rel = make_wide_table(
                annual_df, scenario,
                value_col="relative_change_of_means_pct",
                metric_label="Mean relative change vs baseline (%)",
                scalar_only=False,
                include_baseline_row=True,
            )

            abs_change = make_wide_table(
                annual_df, scenario,
                value_col="absolute_change_of_means",
                metric_label="Absolute change vs baseline (t/y)",
                scalar_only=False,
                include_baseline_row=True,
            )

            if "sensitivity_ratio_mean" in annual_df.columns:
                scalar_s = make_wide_table(
                    annual_df, scenario,
                    value_col="sensitivity_ratio_mean",
                    metric_label="S = (%Δ annual mean sink flow) / (%Δ TC)",
                    scalar_only=True,
                    include_baseline_row=True,
                )
            else:
                scalar_s = pd.DataFrame()

            base.to_excel(writer, sheet_name=f"{scenario}_baseline_t_y", index=False)
            rel.to_excel(writer, sheet_name=f"{scenario}_rel_change_pct", index=False)
            abs_change.to_excel(writer, sheet_name=f"{scenario}_abs_change_t_y", index=False)
            if not scalar_s.empty:
                scalar_s.to_excel(writer, sheet_name=f"{scenario}_S_scalar_TC", index=False)

    return path


# =============================================================================
# CUMULATIVE BUBBLE MAP
# =============================================================================

def colorblind_diverging_cmap():
    """Colorblind-friendly diverging map for signed relative changes."""
    return LinearSegmentedColormap.from_list(
        "blue_neutral_orange_colorblind",
        [COLOR_NEGATIVE, COLOR_NEUTRAL, COLOR_POSITIVE],
        N=256,
    )


def bubble_sizes(values: np.ndarray, min_size: float = 35, max_size: float = 1050) -> np.ndarray:
    vals = np.abs(np.asarray(values, dtype=float))
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return np.full(vals.shape, min_size)

    # Clip at 95th percentile so a single extreme value does not hide all other bubbles.
    vmax = np.nanpercentile(finite, 95)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = np.nanmax(finite)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    scaled = np.minimum(vals, vmax) / vmax
    return min_size + (max_size - min_size) * np.sqrt(scaled)


def make_effect_norm(values: np.ndarray):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

    max_abs = float(np.nanmax(np.abs(finite)))
    if max_abs <= 0:
        max_abs = 1.0

    # Use symmetric log when extreme values exist, e.g. +9714%.
    if max_abs > 300:
        return SymLogNorm(linthresh=20, linscale=1.0, vmin=-max_abs, vmax=max_abs, base=10)
    return TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)


def prepare_cumulative_bubble_data(df: pd.DataFrame) -> pd.DataFrame:
    plot_df = df[
        (df["indicator"] == CUMULATIVE_INDICATOR) &
        (df["output"].isin(ALL_SINKS_ORDER))
    ].copy()

    if plot_df.empty:
        available = sorted(df["indicator"].dropna().astype(str).unique())
        raise ValueError(
            f"No rows found for indicator '{CUMULATIVE_INDICATOR}'. "
            f"Available indicators: {available}"
        )

    plot_df["effect"] = pd.to_numeric(plot_df["relative_change_of_means_pct"], errors="coerce")
    plot_df = plot_df[np.isfinite(plot_df["effect"])].copy()

    if plot_df.empty:
        raise ValueError("Cumulative bubble-map data exists but all effect values are non-finite.")

    return plot_df


def plot_cumulative_bubble_map(df: pd.DataFrame, output_root: Path) -> List[Path]:
    output_dir = output_root / "_cumulative_sink_bubble_maps"
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_df = prepare_cumulative_bubble_data(df)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    data_path = output_dir / f"cumulative_sink_bubble_map_data_2014_2050_{timestamp}.xlsx"
    plot_df.to_excel(data_path, index=False)

    scenarios = [s for s in SCENARIOS_TO_PLOT if s in set(plot_df["scenario"].astype(str))]
    if not scenarios:
        scenarios = sorted(plot_df["scenario"].dropna().astype(str).unique())

    sinks = [s for s in ALL_SINKS_ORDER if s in set(plot_df["output"].astype(str))]
    sink_to_y = {sink: i for i, sink in enumerate(sinks)}

    cmap = colorblind_diverging_cmap()
    norm = make_effect_norm(plot_df["effect"].to_numpy(dtype=float))

    saved_paths: List[Path] = []

    fig_height = max(7.5, 3.8 * len(scenarios))
    fig_width = 10.8
    fig, axes = plt.subplots(
        nrows=len(scenarios),
        ncols=1,
        figsize=(fig_width, fig_height),
        sharey=True,
        constrained_layout=False,
    )
    if len(scenarios) == 1:
        axes = [axes]

    panel_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    last_scatter = None

    for ax, scenario, letter in zip(axes, scenarios, panel_letters):
        sub = plot_df[plot_df["scenario"] == scenario].copy()
        variants = ordered_variants(sub, scenario, scalar_only=False)
        if not variants:
            ax.set_visible(False)
            continue

        variant_to_x = {variant: i for i, variant in enumerate(variants)}
        sub = sub[sub["variant"].isin(variants)].copy()
        sub["x"] = sub["variant"].map(variant_to_x)
        sub["y"] = sub["output"].map(sink_to_y)
        sub = sub[np.isfinite(sub["x"]) & np.isfinite(sub["y"])].copy()

        values = sub["effect"].to_numpy(dtype=float)
        last_scatter = ax.scatter(
            sub["x"],
            sub["y"],
            s=bubble_sizes(values),
            c=values,
            cmap=cmap,
            norm=norm,
            edgecolors=EDGE_COLOR,
            linewidths=0.8,
            alpha=0.95,
        )

        # Keep labels only for effects large enough to matter visually.
        for _, row in sub.iterrows():
            val = float(row["effect"])
            if abs(val) >= PRACTICAL_LABEL_THRESHOLD_PCT:
                ax.annotate(
                    f"{val:.0f}%",
                    xy=(row["x"], row["y"]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                    ha="left",
                    va="bottom",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                )

        ax.set_xticks(range(len(variants)))
        ax.set_xticklabels(variants, rotation=45, ha="right")
        ax.set_yticks(range(len(sinks)))
        ax.set_yticklabels([ALL_SINK_LABELS.get(s, s) for s in sinks])
        ax.invert_yaxis()
        ax.set_xlabel("Sensitivity variant")
        ax.set_ylabel("Sink compartment")
        ax.grid(True, linestyle=":", linewidth=0.6, color=GRID_COLOR)
        ax.set_axisbelow(True)
        ax.text(
            1.02,
            0.98,
            letter,
            transform=ax.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
            ha="left",
        )

    if last_scatter is not None:
        fig.subplots_adjust(left=0.24, right=0.83, bottom=0.08, top=0.97, hspace=0.55)
        cax = fig.add_axes([0.86, 0.18, 0.025, 0.68])
        cbar = fig.colorbar(last_scatter, cax=cax)
        cbar.set_label("Mean relative change vs baseline (%)")

    stem = f"cumulative_sink_bubble_map_all_sinks_2014_2050_{timestamp}"
    for ext in ["png", "pdf", "svg"]:
        out = output_dir / f"{stem}.{ext}"
        fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
        saved_paths.append(out)
    plt.close(fig)

    # Also save one figure per scenario for editing.
    for scenario in scenarios:
        sub = plot_df[plot_df["scenario"] == scenario].copy()
        variants = ordered_variants(sub, scenario, scalar_only=False)
        if not variants:
            continue

        variant_to_x = {variant: i for i, variant in enumerate(variants)}
        sub = sub[sub["variant"].isin(variants)].copy()
        sub["x"] = sub["variant"].map(variant_to_x)
        sub["y"] = sub["output"].map(sink_to_y)
        sub = sub[np.isfinite(sub["x"]) & np.isfinite(sub["y"])].copy()

        values = sub["effect"].to_numpy(dtype=float)

        fig, ax = plt.subplots(
            figsize=(max(7.8, 0.75 * len(variants) + 2.8), max(5.8, 0.45 * len(sinks) + 2.0))
        )
        sc = ax.scatter(
            sub["x"],
            sub["y"],
            s=bubble_sizes(values),
            c=values,
            cmap=cmap,
            norm=norm,
            edgecolors=EDGE_COLOR,
            linewidths=0.8,
            alpha=0.95,
        )

        for _, row in sub.iterrows():
            val = float(row["effect"])
            if abs(val) >= PRACTICAL_LABEL_THRESHOLD_PCT:
                ax.annotate(
                    f"{val:.0f}%",
                    xy=(row["x"], row["y"]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                    ha="left",
                    va="bottom",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                )

        ax.set_xticks(range(len(variants)))
        ax.set_xticklabels(variants, rotation=45, ha="right")
        ax.set_yticks(range(len(sinks)))
        ax.set_yticklabels([ALL_SINK_LABELS.get(s, s) for s in sinks])
        ax.invert_yaxis()
        ax.set_xlabel("Sensitivity variant")
        ax.set_ylabel("Sink compartment")
        ax.grid(True, linestyle=":", linewidth=0.6, color=GRID_COLOR)
        ax.set_axisbelow(True)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Mean relative change vs baseline (%)")
        fig.tight_layout()

        stem = f"cumulative_sink_bubble_map_{scenario}_all_sinks_2014_2050_{timestamp}"
        for ext in ["png", "pdf", "svg"]:
            out = output_dir / f"{stem}.{ext}"
            fig.savefig(out, dpi=SAVE_DPI if ext == "png" else None, bbox_inches="tight")
            saved_paths.append(out)
        plt.close(fig)

    saved_paths.append(data_path)
    return saved_paths


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    workbook = find_latest_sensitivity_workbook()
    output_root = infer_results_root_from_workbook(workbook)

    print("=" * 80)
    print("POST-PROCESSING LOCAL SENSITIVITY RESULTS")
    print("=" * 80)
    print(f"Input workbook: {workbook}")
    print(f"Output root:    {output_root}")

    df = load_sensitivity_table(workbook)
    print(f"Rows loaded after exclusions: {len(df)}")
    print(f"Scenarios found: {sorted(df['scenario'].dropna().astype(str).unique())}")
    print(f"Indicators found: {sorted(df['indicator'].dropna().astype(str).unique())}")

    annual_path = save_annual_tables(df, output_root)
    bubble_paths = plot_cumulative_bubble_map(df, output_root)

    manifest_path = output_root / "postprocess_LSA_outputs_manifest.txt"
    with open(manifest_path, "w", encoding="utf-8") as f:
        f.write("Post-processing local sensitivity-analysis outputs\n")
        f.write(f"Input workbook: {workbook}\n")
        f.write(f"Output root: {output_root}\n\n")
        f.write("Annual table:\n")
        f.write(f"{annual_path}\n\n")
        f.write("Bubble-map outputs:\n")
        for path in bubble_paths:
            f.write(f"{path}\n")

    print("\nSaved annual table:")
    print(f"  {annual_path}")
    print("\nSaved cumulative bubble-map outputs:")
    for path in bubble_paths:
        print(f"  {path}")
    print("\nManifest:")
    print(f"  {manifest_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
