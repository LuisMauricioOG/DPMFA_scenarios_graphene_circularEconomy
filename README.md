# DPMFA_scenarios_graphene_circularEconomy
These simulated scenarios are developed as part of the following publication: 
"Advanced materials in a circular economy: Insights from a prospective material flow scenario analysis of graphene in Europe"

Authors: Luis Mauricio Ortiz-Galvez, Lora Dameska, and Bernd Nowack

Journal: Resources, Conservation, and Recycling

These three prospective scenarios generate probability distributions associated with the flows of graphene within the lifecycle of three product categories: batteries, tires, and construction, and potential releases to the environment, at the European level (EU27, Norway, the United Kingdom, and Switzerland) from 2014 to 2050 with yearly intervals. 

The tool is written using the DPMFA Python package (github.com/empa-tsl/dpmfa), consisting of:
* Components.py: the different classes of compartments and flows are defined;
* Model.py: the required methods to link the compartments together with the flows are defined;
* Simulator.py: includes the methods required to run the simulation and extract results;
* The main scripts with each scenario modeled are: DPMFA_graphene_TBC_model_S1.py and DPMFA_graphene_TBC_runner_S1.py

Each DPMFA_graphene_TBC_model_S.py takes as input a CSV file containing the probability distributions of the produced amounts of graphene in each year in the corresponding product categories: "for_py_graphene_TBC.csv"

Additionally, it includes a Local Sensitivity Analysis (LSA) to test specific assumptions and evaluate their impact on the results in the Waste-as-Usual and Circular Economy-oriented scenarios (compared to the corresponding base scenarios). The analysis uses direct differential sensitivity analysis to assess the robustness of DPMFA outputs with respect to variations in individual model parameters using mean parameter values. The variant scenarios are inside the folder SensitivityAnalysis_runs, while the LSA and post-processing results scripts are within the same folder as the main scenarios used in this study.

## Getting started

**Requirements:** Python 3.10 or newer with `numpy`, `pandas`, `scipy` and `matplotlib` (the LSA workbook export also needs `openpyxl`), plus the DPMFA classes included here (`components.py`, `Model.py`, `Simulator.py`).

**Folder layout** (no path editing is needed; every script finds its files relative to its own location):
```
<repository root>
├── for_py_graphene_TBC.csv                 input data (graphene production distributions)
├── components.py, Model.py, Simulator.py, TriangularBalanced.py, TruncatingFunctions.py, Cumulative_sum_script.py
├── DPMFA_graphene_TBC_model_S1/S2/S3.py    scenario models
├── DPMFA_graphene_TBC_runner_S1/S2/S3.py   run a scenario and write results
├── LocalSensitivityAnalysis_S1_S3.py       master script: runs all variants, builds sensitivity tables and figures
├── Postprocess_LSA_S1_S3.py                re-plots the LSA results (run after the LSA)
├── SensitivityAnalysis_runs/S1_SA, S3_SA   one model + one runner per sensitivity variant
└── results_runner/                         created on first run; all outputs go here (ignored by git)
```

**Running**
1. Single scenario: `python DPMFA_graphene_TBC_runner_S1.py` (or S2 / S3). Results are written to `results_runner/<timestamp>_<year>/`.
2. Sensitivity analysis: `python LocalSensitivityAnalysis_S1_S3.py` runs every variant in `SensitivityAnalysis_runs/` and writes tables and figures to `results_runner/<date>_SensitivityAnalysis_results/`. Then run `python Postprocess_LSA_S1_S3.py`.
3. To write results somewhere else, set the environment variable `DPMFA_OUTPUT_ROOT` before running (for example `set DPMFA_OUTPUT_ROOT=D:\results` on Windows). The runners, the LSA script and the post-processing script all use the same rule, so they always agree.

The LSA runs 24 model variants with Monte Carlo simulations and can take a long time. Set `RUN_RUNNERS_FIRST = False` in the user-settings block of `LocalSensitivityAnalysis_S1_S3.py` if the variant results already exist.

See `CHANGELOG.md` for the latest changes to the sensitivity variants.

For further details on each scenario, the user is referred to the publication above.
