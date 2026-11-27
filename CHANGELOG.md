# CHANGELOG – local sensitivity analysis (S1, S3)

Changes made 2026-11-24. Nothing was re-run when these edits were made.

## Model fixes (SensitivityAnalysis_runs)
- S1 `Hyd_loss_100/high/low`: `impurities -> EoW_alternatives` transfer active again.
- Hyd_loss runner labels fixed.
- `Hyd_loss_high/low` ranges now ±50 %, matching the base.
- Open (not changed): S1 `Hyd_loss_high/low` still enable the unused `battW_RecMat_GBM` compartment (harmless); `TiresRec_high` cementFill mean 0.040 vs computed 0.0388 (rounding).

## Portability of the repository copy (2026-11-27)
- All hard-coded absolute (machine-specific) paths were removed from the repository copy of the scripts (models, runners, LSA, post-processing).
- Inputs: models find `for_py_graphene_TBC.csv` and the shared modules by looking in the script's folder and up to two levels above it.
- Outputs: runners, the LSA script and the post-processing script write to `<repository>/results_runner/` (git-ignored). Set the environment variable `DPMFA_OUTPUT_ROOT` to change this.
- Top-level runners now write the raw per-flow `.txt` files inside the run's own results folder (`results_runner/<timestamp>_<year>/Flows_txt/RawData_S*/`) instead of the fixed folder `results_runner/2025_10_16_Results/Flows_txt/`, and the `'\\'` Windows separators were replaced by `os.path.join`.
- Not re-run: the edits were checked for syntax and for path resolution only (the LSA discovers all 24 variant runners), not by re-running the simulations.

## Transfer coefficient changed per variant
| Variant | TC changed (base -> variant) |
|---|---|
| `TiresEnRec` (S1, S3) | `tiresW_EnergyRec -> energyRecovery`: 0.09 -> 0.108 / 0.072 (±20 %) |
| `TiresRec` (S1, S3) | Mechanical-recycling share of `tiresW_Recycling`: 0.796 -> 0.955 / 0.637 (±20 %) |
| `Hyd_loss_high / low` (S1, S3) | `battW_HydroLeaching -> _losses`: 0.09 -> 0.108 / 0.072 (±20 %) |
| `Hyd_loss_50` / `Hyd_loss_100` | Same TC set to 0.50 / 1.0 |
| `BattW2Lf` (S1) | `battW_Unidentified`: Landfill 0.29 / Unknown 0.71 -> Landfill 1 |
| `BattRec_Unk2Lf` (S1), `BattW_RecRes2Lf` (S3) | `battW_Recycling_Residue`: Unknown -> Landfill |
| `BattWDism_dynamic` (S1) | `battW` yearly TC y7–y36: flat 0.5 -> S3 trajectory |
| `BattDirRec_high / low` (S3) | `battW_DirectRecycling -> RecMat_GBM`: 1 -> 0.95 (truncated triangular) / 0.40 |

## LocalSensitivityAnalysis_S1_S3.py
- `infer_input_change_pct(scenario, variant, profile)`:
  - `original`: unchanged values (`BattDirRec_high` = -5 %).
  - `corrected`: `BattDirRec_high` = -17.79 % (analytic mean of the triangular [0.475, 0.95, 1.425] truncated at 1 is 0.822). `BattDirRec_low` stays -60 %.
  - `Hyd_loss_high/low` = ±20 %; `Hyd_loss_50` = +455.6 % (scalar S); `Hyd_loss_100` = no S (relative change only).
- Six views (`VIEWS`), computed from the stored table, no re-simulation:

  | View | Profile | Hyd_loss variants | Outputs |
  |---|---|---|---|
  | A | original | none | tornado, heatmaps, bubble maps, annual tables |
  | B | corrected | none | same |
  | C | original | ±20 % | same |
  | D | corrected | ±20 % | same |
  | E | original | ±20 %, 50, 100 | same |
  | F | corrected | ±20 %, 50, 100 | same |

  Outputs go to `<RESULTS_ROOT>/_views/view_<X>_<profile>_<Hyd>/`, filenames tagged `_view<X>`. Legacy output folders now hold view A.
- Tornado rule: only variants with `|Δinput| <= 20 %` (the ±20 % variants or anything within that range, e.g. `BattDirRec_high` at -5 % / -17.8 %, labelled with its real Δinput). `BattDirRec_low` (-60 %) and `Hyd_loss_50` (+456 %) are not in the tornado. The tornado ranks by max |S|. (Switches `TORNADO_EXPLICIT_ALLOW` / `TORNADO_EXTENDED_RANGE` exist but are empty.) Every view (A–F) gets the full set of outputs; E/F differ from C/D only by the extra Hyd_loss_50/100 rows in bubble maps, heatmaps (Hyd_loss_50 only; _100 has no S) and annual tables.
- Workbook: `all_results` and `summary_compact` unchanged in structure; new sheets `views_definition`, `view_A` … `view_F`.
- Module docstring rewritten: finite-difference secant on Monte Carlo means (ratio of means), "inspired by" Bornhöft et al., %Δinput = relative change of the perturbed TC (effective mean, siblings rebalanced by TC normalization), per-variant input-change rules, seed/pairing caveat. Annual-table label is now "%Δ perturbed TC"; the heatmap gets a one-line footnote defining S.
- Heatmap colorbar label is now just an upright "S". Figure font is Arial (fallbacks Liberation Sans, DejaVu Sans).

## Postprocess_LSA_S1_S3.py
- Imports the view logic from `LocalSensitivityAnalysis_S1_S3.py`; default outputs = view A; per-view outputs in `_views_postprocess/`.
- `Hyd_loss_*` added to `ordered_variants`.

## Expected result changes
- `BattDirRec_high` (Unknown sink): S goes from about -8.1 (original) to about -2.3 (corrected), close to `BattDirRec_low` (-2.28). "By far the largest sensitivity" (manuscript lines 308–309) no longer holds under the corrected profile.

## Manuscript notes (PDFs not edited)
- Method is a finite-difference normalized sensitivity on Monte Carlo mean outputs, inspired by Bornhöft et al. (2021), not identical.
- Table S5: BattDirRec shown as "TC = 0.95" / trapezoidal; code uses triangular (truncated triangular for high).
- Naming: `BattRec_Unk2Lf` (code) vs `BattW_RecRes2Lf` (S1, Table S5).
- Line 317: dynamic-variant baseline is about 0.5, not 1.0.
- Add Hyd_loss rows to Table S5; include them in the ±20 % sentence at line 206.
