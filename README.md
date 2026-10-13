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

For further details on each scenario, the user is referred to the publication above.
