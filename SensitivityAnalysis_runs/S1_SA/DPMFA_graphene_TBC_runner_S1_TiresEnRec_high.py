# -*- coding: utf-8 -*-
#######################################
# @author: Luis Mauricio Ortiz-Galvez
# Last modification on October 16, 2025
 
startYear = 2014 #from 2014 to 2050
Tperiods = 37     # total time range considered, also for np.arange, important for development of  Accumulation, release, Concentration
                    # if Tperiods is changed, the import settings of the production data in the model must also be changed (if Tperiods get higher)
RUNs = 1000 # e.g., 1000 RUNs means 1000 simulations
Speriod = 36 # specific period for concentrations and Mass-flows
########################################



# In[2]:
# from datetime import datetime
import time
import os
import DPMFA_graphene_TBC_model_S1_TiresEnRec_high as ScI_Model # here space is not allowed in model name
import Simulator as sc
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nr
from scipy.stats import gaussian_kde
import pandas as pd
import matplotlib.ticker as ticker
import csv
import Cumulative_sum_script as cs


# In[3]:
## prepare time print for folders
t = time.strftime("%Y-%m-%d-%H-%M_%S")
tYear = str(startYear + Speriod)

# Sensitivity-analysis output root
RESULTS_ROOT = r"C:\\Users\\ortl\\Documents\\Empa\\Articles\\pDPMFA_AdMa_EU_CE\\20260622_SensitivityAnalysis_results\\_raw_runner_outputs"

SCENARIO_LABEL = "S1"
VARIANT_LABEL = "TiresEnRec_high"

dirsname = os.path.join(
    RESULTS_ROOT,
    f"{SCENARIO_LABEL}_{VARIANT_LABEL}_{t}_{tYear}"
)

print("dirsname: " + dirsname)

if not os.path.exists(dirsname):
    os.makedirs(dirsname)

def withOpen_under_directory(filename, mode, delimiter_setting):
    filepath = os.path.join(dirsname, filename + "_" + t + ".csv")
    with open(filepath, mode, newline="") as RM:
        a = csv.writer(RM, delimiter=delimiter_setting)
        data = np.asarray(Flows)
        a.writerows(data)


def toCSV_under_directory(dataname, filename, index_value, seperator):
    filepath = os.path.join(dirsname, filename + "_" + t + ".csv")
    dataname.to_csv(filepath, index=index_value, sep=seperator)


def clean_filename(text):
    """
    Make compartment names safe for Windows filenames.
    """
    text = str(text)
    for char in ['\\', '/', ':', '*', '?', '"', '<', '>', '|', '&']:
        text = text.replace(char, '_')
    text = text.replace(' ', '_')
    return text


def export_cumulative_sink_results(
    sink_matrices,
    scenario_label,
    output_dir,
    start_year,
    t,
    cumulative_start_year=None
):
    """
    Export cumulative mass to sink compartments over time.

    Parameters
    ----------
    sink_matrices : dict
        Dictionary where keys are sink names and values are matrices with shape:
        RUNs x Tperiods.

        Example:
        {
            "Landfill": Landfill_sink,
            "Surfacewater": Surfacewater_sink
        }

    scenario_label : str
        Example: "S1", "S3", "S1_recycling_high".

    output_dir : str
        Same output folder used by the runner: dirsname.

    start_year : int
        Model start year, e.g. 2014.

    t : str
        Timestamp already used in your runner.

    cumulative_start_year : int or None
        If None, cumulative sums start from the model start year.
        If 2025, cumulative sums start from 2025 onward.
    """

    years = list(range(start_year, start_year + Tperiods))

    if cumulative_start_year is None:
        start_idx = 0
        cumulative_period_label = f"{start_year}_{years[-1]}"
    else:
        start_idx = cumulative_start_year - start_year
        cumulative_period_label = f"{cumulative_start_year}_{years[-1]}"

    if start_idx < 0 or start_idx >= Tperiods:
        raise ValueError("cumulative_start_year is outside the model time range.")

    summary_rows = []

    # These dictionaries will become wide tables:
    # rows = sink names, columns = years.
    wide_tables = {
        "mean": {},
        "sd": {},
        "min": {},
        "q05": {},
        "q15": {},
        "q25": {},
        "q50": {},
        "q75": {},
        "q85": {},
        "q95": {},
        "max": {},
    }

    for sink_name, annual_matrix in sink_matrices.items():

        annual_matrix = np.asarray(annual_matrix, dtype=float)

        if annual_matrix.ndim != 2:
            print(f"Skipping {sink_name}: matrix is not 2D.")
            continue

        # Create a matrix with same shape, but cumulative only from selected start year.
        selected_annual = annual_matrix[:, start_idx:]

        cumulative_selected = np.cumsum(selected_annual, axis=1)

        # Put the cumulative result back into full model-time columns.
        cumulative_full = np.full_like(annual_matrix, np.nan, dtype=float)
        cumulative_full[:, start_idx:] = cumulative_selected

        # Export raw cumulative run-by-year matrix.
        cumulative_df = pd.DataFrame(cumulative_full, columns=years)

        raw_filename = (
            f"Results_SinkCumulativeDynamic_"
            f"{clean_filename(sink_name)}_{scenario_label}_{cumulative_period_label}_{t}.csv"
        )

        cumulative_df.to_csv(
            os.path.join(output_dir, raw_filename),
            index=True
        )

        # Calculate yearly statistics across MC runs.
        stats_by_year = {
            "mean": [],
            "sd": [],
            "min": [],
            "q05": [],
            "q15": [],
            "q25": [],
            "q50": [],
            "q75": [],
            "q85": [],
            "q95": [],
            "max": [],
        }

        for period, year in enumerate(years):
            values = cumulative_full[:, period]
            values = values[~np.isnan(values)]

            if values.size == 0:
                stats = {
                    "mean": np.nan,
                    "sd": np.nan,
                    "min": np.nan,
                    "q05": np.nan,
                    "q15": np.nan,
                    "q25": np.nan,
                    "q50": np.nan,
                    "q75": np.nan,
                    "q85": np.nan,
                    "q95": np.nan,
                    "max": np.nan,
                }
            else:
                stats = {
                    "mean": np.mean(values),
                    "sd": np.std(values),
                    "min": np.min(values),
                    "q05": np.percentile(values, 5),
                    "q15": np.percentile(values, 15),
                    "q25": np.percentile(values, 25),
                    "q50": np.percentile(values, 50),
                    "q75": np.percentile(values, 75),
                    "q85": np.percentile(values, 85),
                    "q95": np.percentile(values, 95),
                    "max": np.max(values),
                }

            for stat_name, stat_value in stats.items():
                stats_by_year[stat_name].append(stat_value)

            summary_rows.append({
                "scenario": scenario_label,
                "sink": sink_name,
                "year": year,
                "cumulative_from_year": cumulative_start_year if cumulative_start_year is not None else start_year,
                "mean": stats["mean"],
                "sd": stats["sd"],
                "min": stats["min"],
                "q05": stats["q05"],
                "q15": stats["q15"],
                "q25": stats["q25"],
                "q50": stats["q50"],
                "q75": stats["q75"],
                "q85": stats["q85"],
                "q95": stats["q95"],
                "max": stats["max"],
            })

        for stat_name in wide_tables.keys():
            wide_tables[stat_name][sink_name] = stats_by_year[stat_name]

    # Long-format summary table: one row per sink-year.
    summary_df = pd.DataFrame(summary_rows)

    summary_filename = (
        f"RESULTS_CUMULATIVE_To_Sinks_summary_"
        f"{scenario_label}_{cumulative_period_label}_{t}.csv"
    )

    summary_df.to_csv(
        os.path.join(output_dir, summary_filename),
        index=False
    )

    # Wide-format tables: rows = sinks, columns = years.
    for stat_name, stat_dict in wide_tables.items():
        wide_df = pd.DataFrame(stat_dict, index=years).T
        wide_df.columns = years

        wide_filename = (
            f"RESULTS_CUMULATIVE_To_Sinks_{stat_name}_"
            f"{scenario_label}_{cumulative_period_label}_{t}.csv"
        )

        wide_df.to_csv(
            os.path.join(output_dir, wide_filename),
            index=True
        )

    print("Cumulative sink results exported.")
    print(os.path.join(output_dir, summary_filename))

def export_cumulative_sinks_master_summary(
    sink_matrices,
    scenario_label,
    output_dir,
    start_year,
    cumulative_from_year,
    cumulative_to_year,
    t,
    save_run_level_vectors=True
):
    """
    Export one master summary table with cumulative mass to all sink compartments.

    This replaces the separate post-analysis script that read sink CSV files
    and created ALL_SINKS_cumulative_summary_MASTER.csv.

    The calculation is done directly from the annual sink matrices:
        rows = Monte Carlo runs
        columns = years / model periods

    For each sink:
        cumulative_mass_per_run = sum from cumulative_from_year to cumulative_to_year
    """

    start_idx = cumulative_from_year - start_year
    end_idx = cumulative_to_year - start_year

    if start_idx < 0 or end_idx >= Tperiods:
        raise ValueError(
            f"Cumulative period {cumulative_from_year}-{cumulative_to_year} "
            f"is outside the model range {start_year}-{start_year + Tperiods - 1}."
        )

    summary_rows = []
    run_level_dict = {}

    for sink_name, annual_matrix in sink_matrices.items():

        annual_matrix = np.asarray(annual_matrix, dtype=float)

        if annual_matrix.ndim != 2:
            print(f"Skipping {sink_name}: matrix is not 2D.")
            continue

        # Sum annual sink inflows per Monte Carlo run over the selected period.
        cumulative_per_run = np.nansum(
            annual_matrix[:, start_idx:end_idx + 1],
            axis=1
        )

        clean_sink_name = clean_filename(sink_name)
        run_level_dict[clean_sink_name] = pd.Series(cumulative_per_run)

        summary_rows.append({
            "scenario": scenario_label,
            "sink": sink_name,
            "n_runs": cumulative_per_run.shape[0],
            "columns_used": end_idx - start_idx + 1,
            "years_used": f"{cumulative_from_year}-{cumulative_to_year}",
            "mean": np.mean(cumulative_per_run),
            "median": np.median(cumulative_per_run),
            "std": np.std(cumulative_per_run, ddof=1),
            "min": np.min(cumulative_per_run),
            "max": np.max(cumulative_per_run),
            "p05": np.percentile(cumulative_per_run, 5),
            "p15": np.percentile(cumulative_per_run, 15),
            "p25": np.percentile(cumulative_per_run, 25),
            "p75": np.percentile(cumulative_per_run, 75),
            "p85": np.percentile(cumulative_per_run, 85),
            "p95": np.percentile(cumulative_per_run, 95),
        })

    master_df = pd.DataFrame(summary_rows)

    master_csv_path = os.path.join(
        output_dir,
        f"ALL_SINKS_cumulative_summary_MASTER_{scenario_label}_{cumulative_from_year}_{cumulative_to_year}_{t}.csv"
    )

    master_xlsx_path = os.path.join(
        output_dir,
        f"ALL_SINKS_cumulative_summary_MASTER_{scenario_label}_{cumulative_from_year}_{cumulative_to_year}_{t}.xlsx"
    )

    master_df.to_csv(master_csv_path, index=False)

    # Optional Excel version
    try:
        master_df.to_excel(master_xlsx_path, index=False)
    except Exception as exc:
        print(f"Could not save Excel master summary: {exc}")

    if save_run_level_vectors:
        run_level_df = pd.DataFrame(run_level_dict)

        run_level_path = os.path.join(
            output_dir,
            f"ALL_SINKS_cumulative_run_level_{scenario_label}_{cumulative_from_year}_{cumulative_to_year}_{t}.csv"
        )

        run_level_df.to_csv(run_level_path, index=False)

    print("Master cumulative sink summary exported:")
    print(master_csv_path)

## Exporting likelihood distributions
## This saves the empirical probability distribution pool used by RandomChoiceTransfer
def export_likelihood_distributions(model_module, scenario_label, output_dir):
    """
    Export all likelihood-distribution arrays created in the model.

    It automatically finds variables ending in '_dist', for example:
        cementKilns_dist
        Pyrolisis_dist
        tiresW_MechanicalRecycling_dist
        battW_SizeReduction_dist

    It saves:
        1. the raw sampled distribution arrays
        2. summary statistics for each distribution
    """

    dist_names = [
        name for name in dir(model_module)
        if name.endswith("_dist")
    ]

    if len(dist_names) == 0:
        print("No variables ending in '_dist' were found in the model.")
        return

    samples_dict = {}
    summary_rows = []

    for name in dist_names:
        values = getattr(model_module, name)

        try:
            values = np.asarray(values, dtype=float).ravel()
        except Exception:
            print(f"Skipping {name}: could not convert to numeric array.")
            continue

        if values.size == 0:
            print(f"Skipping {name}: empty array.")
            continue

        samples_dict[name] = pd.Series(values)

        summary_rows.append({
            "scenario": scenario_label,
            "distribution_name": name,
            "n_values": values.size,
            "min": np.min(values),
            "q05": np.percentile(values, 5),
            "q25": np.percentile(values, 25),
            "median": np.percentile(values, 50),
            "mean": np.mean(values),
            "q75": np.percentile(values, 75),
            "q95": np.percentile(values, 95),
            "max": np.max(values),
            "sd": np.std(values)
        })

    samples_df = pd.DataFrame(samples_dict)
    summary_df = pd.DataFrame(summary_rows)

    samples_path = os.path.join(
        output_dir,
        f"LIKELIHOOD_distribution_samples_{scenario_label}_{t}.csv"
    )

    summary_path = os.path.join(
        output_dir,
        f"LIKELIHOOD_distribution_summary_{scenario_label}_{t}.csv"
    )

    samples_df.to_csv(samples_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    print("Likelihood distributions exported:")
    print(samples_path)
    print(summary_path)


# In[4]:
print('Results here in console for period: ' + str(Speriod))
plt.close()

model = ScI_Model.model
model.checkModelValidity()

export_likelihood_distributions(
    model_module=ScI_Model,
    scenario_label=f"{SCENARIO_LABEL}_{VARIANT_LABEL}",
    output_dir=dirsname
)

simulator = sc.Simulator(RUNs, Tperiods, 2250, True, True)
simulator.setModel(model)
simulator.runSimulation()

# In[6]:

sinks = simulator.getSinks()
stocks = simulator.getStocks()
# In[7]:
#print "logged Inflows called"

Comp_with_loggedInflows = simulator.getLoggedInflows() # compartment with loggedInflows

# In[8]:
#print "loggedInflows" + str(loggedInflows)
Comp_with_loggedOutflows = simulator.getLoggedOutflows() # compartment with LoggedOutflows

RoundV = 6 #rounding value
# In[9]:
## this part is used to display all flows!

print('logged Outflows: ')
for Comp in Comp_with_loggedOutflows: # loggedOutflows is the compartment list of compartmensts with loggedoutflows
    print('Flows from ' + Comp.name + ':' )
    for Target_name, value in Comp.outflowRecord.items(): # in this case name is the key, value is the matrix(data), in this case .items is needed
        print(Comp.name + ' --> ' + str(Target_name))
       # Saving raw data
        print(value[:,Speriod])
        raw = value[:, Speriod] 
        raw_dir = os.path.join(dirsname, "Flows_txt", f"RawData_{SCENARIO_LABEL}_{VARIANT_LABEL}")
        os.makedirs(raw_dir, exist_ok=True)

        raw_filename = clean_filename(f"{Comp.name}_TO_{str(Target_name)}") + ".txt"

        np.savetxt(
            os.path.join(raw_dir, raw_filename),
            raw,
            delimiter=","
        )         

        print('Mean--> ' + str(round(np.mean(value[:,Speriod]),RoundV)))
        print('SD--> ' + str(round(np.std(value[:,Speriod]),RoundV)))
        print('Q5--> ' + str(round(np.percentile(value[:,Speriod],5),RoundV)))
        print('Q15--> ' + str(round(np.percentile(value[:,Speriod],15),RoundV)))
        print('Q25--> ' + str(round(np.percentile(value[:,Speriod],25),RoundV)))
        print('Q50--> ' + str(round(np.percentile(value[:,Speriod],50),RoundV)))
        print('Q75--> ' + str(round(np.percentile(value[:,Speriod],75),RoundV)))
        print('Q85--> ' + str(round(np.percentile(value[:,Speriod],85),RoundV)))
        print('Q95--> ' + str(round(np.percentile(value[:,Speriod],95),RoundV)))
        print(' ')
    print(' ')
print(' ')


# In[10]:
############### yearly inflow  ######################
inflowMatrix = Comp_with_loggedInflows['Total Inflow']
print('flow record: ' + str('Total Inflow'))
print('mean: ' + str(round(np.mean(inflowMatrix[:,Speriod]),RoundV)))
print('SD: ' + str(round(np.std(inflowMatrix[:,Speriod]),RoundV)))
#print('mode:'+str(round(peakValue,4)))
print('Q5: ' + str(round(np.percentile(inflowMatrix[:,Speriod],5),RoundV)))
print('Q15: ' + str(round(np.percentile(inflowMatrix[:,Speriod],15),RoundV)))
print('Q25: ' + str(round(np.percentile(inflowMatrix[:,Speriod],25),RoundV)))
print('Q50: ' + str(round(np.percentile(inflowMatrix[:,Speriod],50),RoundV)))
print('Q75: ' + str(round(np.percentile(inflowMatrix[:,Speriod],75),RoundV)))
print('Q85: ' + str(round(np.percentile(inflowMatrix[:,Speriod],85),RoundV)))
print('Q95: ' + str(round(np.percentile(inflowMatrix[:,Speriod],95),RoundV)))


# In[11]:

############### ALL DISAGGREGATED FLOWS ######################
import csv

# export all flows to an excel - means
Flows = []
for Comp in Comp_with_loggedOutflows: # loggedOutflows is the compartment list of compartmensts with loggedoutflows
    for Target_name, value in Comp.outflowRecord.items(): # in this case name is the key, value is the matrix(data), in this case .items is needed
        x = []
        for i in np.arange(0,Tperiods):
            x.append((np.mean(value[:,i])))
            z = [Comp.name, 'to',Target_name]
            y= x + z
        Flows.append(y)

    withOpen_under_directory('RESULTS_AllFlows_S1_mean_2014_2050', 'w', ' ')

# export all flows to an excel - standard deviations
Flows = []
for Comp in Comp_with_loggedOutflows: # loggedOutflows is the compartment list of compartmensts with loggedoutflows
    for Target_name, value in Comp.outflowRecord.items(): # in this case name is the key, value is the matrix(data), in this case .items is needed
        x = []
        for i in np.arange(0,Tperiods):
            x.append((np.std(value[:,i])))
            z = [Comp.name, 'to',Target_name]
            y= x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_std_2014_2050', 'w', ' ')


# In[11]:

# VA: export all flows to an excel - Q5
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 5))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q5_2014_2050', 'w', ' ')

# In[13]:
# VA: export all flows to an excel - Q15
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 15))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q15_2014_2050', 'w', ' ')    
    
    
# In[14]:
# VA: export all flows to an excel - Q25
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 25))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q25_2014_2050', 'w', ' ')

# In[15]:
# VA: export all flows to an excel - Q75
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 75))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q75_2014_2050', 'w', ' ')

# In[16]:
# VA: export all flows to an excel - Q85
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 85))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q85_2014_2050', 'w', ' ')
    
# In[17]:
# VA: export all flows to an excel - Q95
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 95))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_Q95_2014_2050', 'w', ' ')

# In[18]:
# VA: export all flows to an excel - min
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(min(value[:,i]))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_min_2014_2050', 'w', ' ')

# In[19]:
# VA: export all flows to an excel - max
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(max(value[:,i]))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_max_2014_2050', 'w', ' ')

# In[20]:
# VA: export all flows to an excel - medians
Flows = []
for Comp in Comp_with_loggedOutflows:
    for Target_name, value in Comp.outflowRecord.items():
        x = []
        for i in np.arange(0, Tperiods):
            x.append(np.percentile(value[:,i], 50))
            z = [Comp.name, 'to', Target_name]
            y = x + z
        Flows.append(y)
    withOpen_under_directory('RESULTS_AllFlows_S1_medians_2014_2050', 'w', ' ')


# In[21]:

##### for yearly stock in consumption; PMC = production, manufacture, and consumption stage (as one compartment)

print('PMC flow to technical and environmental compartments Inflows')
compartmentlist = simulator.getCompartments()

PMCrellist = [comp for comp in compartmentlist if ('PMCrel' in comp.categories)]
PMCrel2TE = simulator.getCombinedOutflows(PMCrellist) # PMC2TE is a dictionary

print('From Use and Waste categories to tech. and env. compartments')

for key,value in PMCrel2TE.items():
    Tinflow_Mean = np.mean(value[:,Speriod])
    Tinflow_SD = np.std(value[:,Speriod])
    Tinflow_Q50 = np.percentile(value[:,Speriod], 50)
    Tinflow_Q5 = np.percentile(value[:,Speriod], 5)
    Tinflow_Q15 = np.percentile(value[:,Speriod], 15)
    Tinflow_Q25 = np.percentile(value[:,Speriod], 25)
    Tinflow_Q75 = np.percentile(value[:,Speriod], 75)
    Tinflow_Q85 = np.percentile(value[:,Speriod], 85)
    Tinflow_Q95 = np.percentile(value[:,Speriod], 95)
    Tinflow_min = min(value[:,Speriod])
    Tinflow_max = max(value[:,Speriod])
    print(str(key) + ': ')
    print('Mean --> ' + str(round(Tinflow_Mean, 4)))
    print('SD --> ' + str(round(Tinflow_SD, 4)))
    print('Q50 --> ' + str(round(Tinflow_Q50, 4)))
    print('Q5 --> ' + str(round(Tinflow_Q5, 4)))
    print('Q15 --> ' + str(round(Tinflow_Q15, 4)))
    print('Q25 --> ' + str(round(Tinflow_Q25, 4)))
    print('Q75 --> ' + str(round(Tinflow_Q75, 4)))
    print('Q85 --> ' + str(round(Tinflow_Q85, 4)))
    print('Q95 --> ' + str(round(Tinflow_Q95, 4)))
    print('Min --> ' + str(round(Tinflow_min, 4)))
    print('Max --> ' + str(round(Tinflow_max, 4)))


print('Total release from use and EoL in one period, including stocked from previous years')
# equals to immediate release from consumption and stock release from consumption"

Total_PMCrel = simulator.getLoggedCategoryOutflowSum('PMC')

Total_PMCrel_Mean = np.mean(Total_PMCrel[:,Speriod])
Total_PMCrel_SD = np.std(Total_PMCrel[:,Speriod])
Total_PMCrel_Q50 = np.percentile(Total_PMCrel[:,Speriod], 50)
Total_PMCrel_Q5 = np.percentile(Total_PMCrel[:,Speriod], 5)
Total_PMCrel_Q15 = np.percentile(Total_PMCrel[:,Speriod], 15)
Total_PMCrel_Q25 = np.percentile(Total_PMCrel[:,Speriod], 25)
Total_PMCrel_Q75 = np.percentile(Total_PMCrel[:,Speriod], 75)
Total_PMCrel_Q85 = np.percentile(Total_PMCrel[:,Speriod], 85)
Total_PMCrel_Q95 = np.percentile(Total_PMCrel[:,Speriod], 95)
Total_PMCrel_min = min(Total_PMCrel[:,Speriod])
Total_PMCrel_max = max(Total_PMCrel[:,Speriod])

print('From use and EoL immediate total release')
print('Total_PMCrel_Mean --> ' + str(round(Total_PMCrel_Mean, RoundV)))
print('Total_PMCrel_SD --> ' + str(round(Total_PMCrel_SD, RoundV)))
print('Total_PMCrel_Q50 --> ' + str(round(Total_PMCrel_Q50, RoundV)))
print('Total_PMCrel_Q5 --> ' + str(round(Total_PMCrel_Q5, RoundV)))
print('Total_PMCrel_Q15 --> ' + str(round(Total_PMCrel_Q15, RoundV)))
print('Total_PMCrel_Q25 --> ' + str(round(Total_PMCrel_Q25, RoundV)))
print('Total_PMCrel_Q75 --> ' + str(round(Total_PMCrel_Q75, RoundV)))
print('Total_PMCrel_Q85 --> ' + str(round(Total_PMCrel_Q85, RoundV)))
print('Total_PMCrel_Q95 --> ' + str(round(Total_PMCrel_Q95, RoundV)))
print('Total_PMCrel_min --> ' + str(round(Total_PMCrel_min, RoundV)))
print('Total_PMCrel_max --> ' + str(round(Total_PMCrel_max, RoundV)))


print('Current cumulative PMC stock')
Total_PMCstock = simulator.getLoggedCategoryStock('PMC')
Total_PMCstock_Mean = np.mean(Total_PMCstock[:,Speriod])
Total_PMCstock_SD = np.std(Total_PMCstock[:,Speriod])
Total_PMCstock_Q50 = np.percentile(Total_PMCstock[:,Speriod], 50)
Total_PMCstock_Q5 = np.percentile(Total_PMCstock[:,Speriod], 5)
Total_PMCstock_Q15 = np.percentile(Total_PMCstock[:,Speriod], 15)
Total_PMCstock_Q25 = np.percentile(Total_PMCstock[:,Speriod], 25)
Total_PMCstock_Q75 = np.percentile(Total_PMCstock[:,Speriod], 75)
Total_PMCstock_Q85 = np.percentile(Total_PMCstock[:,Speriod], 85)
Total_PMCstock_Q95 = np.percentile(Total_PMCstock[:,Speriod], 95)
Total_PMCstock_Min = min(Total_PMCstock[:,Speriod])
Total_PMCstock_Max = max(Total_PMCstock[:,Speriod])

print('Total_PMCstock_Mean --> ' + str(round(Total_PMCstock_Mean, RoundV)))
print('Total_PMCstock_SD --> ' + str(round(Total_PMCstock_SD, RoundV)))
print('Total_PMCstock_Q50 --> ' + str(round(Total_PMCstock_Q50, RoundV)))
print('Total_PMCstock_Q5 --> ' + str(round(Total_PMCstock_Q5, RoundV)))
print('Total_PMCstock_Q15 --> ' + str(round(Total_PMCstock_Q15, RoundV)))
print('Total_PMCstock_Q25 --> ' + str(round(Total_PMCstock_Q25, RoundV)))
print('Total_PMCstock_Q75 --> ' + str(round(Total_PMCstock_Q75, RoundV)))
print('Total_PMCstock_Q85 --> ' + str(round(Total_PMCstock_Q85, RoundV)))
print('Total_PMCstock_Q95 --> ' + str(round(Total_PMCstock_Q95, RoundV)))
print('Total_PMCstock_Min --> ' + str(round(Total_PMCstock_Min, RoundV)))
print('Total_PMCstock_Max --> ' + str(round(Total_PMCstock_Max, RoundV)))

PMC_Stock_DF = pd.DataFrame({
    'Specific year': Speriod + startYear,
    'PMC_Stock_Mean': round(Total_PMCstock_Mean, RoundV),
    'PMC_Stock_SD': round(Total_PMCstock_SD, RoundV),
    'PMC_Stock_Q50': round(Total_PMCstock_Q50, RoundV),
    'PMC_Stock_Q5': round(Total_PMCstock_Q5, RoundV),
    'PMC_Stock_Q15': round(Total_PMCstock_Q15, RoundV),
    'PMC_Stock_Q25': round(Total_PMCstock_Q25, RoundV),
    'PMC_Stock_Q75': round(Total_PMCstock_Q75, RoundV),
    'PMC_Stock_Q85': round(Total_PMCstock_Q85, RoundV),
    'PMC_Stock_Q95': round(Total_PMCstock_Q95, RoundV),
    'PMC_Stock_MIN': round(Total_PMCstock_Min, RoundV),
    'PMC_Stock_MAX': round(Total_PMCstock_Max, RoundV)},
    index=[0]
)

PMC_Stock_DF_Transp = pd.DataFrame.transpose(PMC_Stock_DF)
toCSV_under_directory(PMC_Stock_DF_Transp, 'RESULTS_PMC_CUMULATIVE_STOCK - S1', True, ',')

    
# # SOFT FALL SURFACES STOCK
# Total_SFSrel = simulator.getLoggedCategoryOutflowSum('SFS_stock')

# Total_SFSrel_Mean = np.mean(Total_SFSrel[:,Speriod])
# Total_SFSrel_SD = np.std(Total_SFSrel[:,Speriod])
# Total_SFSrel_Q50 = np.percentile(Total_SFSrel[:,Speriod], 50)
# Total_SFSrel_Q25 = np.percentile(Total_SFSrel[:,Speriod], 25)
# Total_SFSrel_Q75 = np.percentile(Total_SFSrel[:,Speriod], 75)
# Total_SFSrel_min = min(Total_SFSrel[:,Speriod])
# Total_SFSrel_max = max(Total_SFSrel[:,Speriod])

# print('From SoftFallSurfaces immediate total release')
# print('Total_SFSrel_Mean --> ' + str(round(Total_SFSrel_Mean, RoundV)))
# print('Total_SFSrel_SD --> ' + str(round(Total_SFSrel_SD, RoundV)))
# print('Total_SFSrel_Q50 --> ' + str(round(Total_SFSrel_Q50, RoundV)))
# print('Total_SFSrel_Q25 --> ' + str(round(Total_SFSrel_Q25, RoundV)))
# print('Total_SFSrel_Q75 --> ' + str(round(Total_SFSrel_Q75, RoundV)))
# print('Total_SFSrel_min --> ' + str(round(Total_SFSrel_min, RoundV)))
# print('Total_SFSrel_max --> ' + str(round(Total_SFSrel_max, RoundV)))

# SFS_Stock_Rel_DF = pd.DataFrame({
#     'Specific year': Speriod + startYear,
#     'SoftFallSurfaces_Stock_Release_Mean': round(Total_SFSrel_Mean, RoundV),
#     'SoftFallSurfaces_Stock_Release_SD': round(Total_SFSrel_SD, RoundV),
#     'SoftFallSurfaces_Stock_Release_Q50': round(Total_SFSrel_Q50, RoundV),
#     'SoftFallSurfaces_Stock_Release_Q25': round(Total_SFSrel_Q25, RoundV),
#     'SoftFallSurfaces_Stock_Release_Q75': round(Total_SFSrel_Q75, RoundV),
#     'SoftFallSurfaces_Stock_Release_MIN': round(Total_SFSrel_min, RoundV),
#     'SoftFallSurfaces_Stock_Release_MAX': round(Total_SFSrel_max, RoundV)},
#     index=[0]
# )

# SFS_Stock_Rel_DF_Transp = pd.DataFrame.transpose(SFS_Stock_Rel_DF)
# toCSV_under_directory(SFS_Stock_Rel_DF_Transp, 'RESULTS_SFS_CUMULATIVE_RELEASE - Tires', True, ',')

# print('Current cumulative SoftFallSurfaces stock')
# Total_SFSstock = simulator.getLoggedCategoryStock('SFS_stock')
# Total_SFSstock_Mean = np.mean(Total_SFSstock[:,Speriod])
# Total_SFSstock_SD = np.std(Total_SFSstock[:,Speriod])
# Total_SFSstock_Q50 = np.percentile(Total_SFSstock[:,Speriod], 50)
# Total_SFSstock_Q25 = np.percentile(Total_SFSstock[:,Speriod], 25)
# Total_SFSstock_Q75 = np.percentile(Total_SFSstock[:,Speriod], 75)
# Total_SFSstock_Min = min(Total_SFSstock[:,Speriod])
# Total_SFSstock_Max = max(Total_SFSstock[:,Speriod])

# print('Total_SFSstock_Mean --> ' + str(round(Total_SFSstock_Mean, RoundV)))
# print('Total_SFSstock_SD --> ' + str(round(Total_SFSstock_SD, RoundV)))
# print('Total_SFSstock_Q50 --> ' + str(round(Total_SFSstock_Q50, RoundV)))
# print('Total_SFSstock_Q25 --> ' + str(round(Total_SFSstock_Q25, RoundV)))
# print('Total_SFSstock_Q75 --> ' + str(round(Total_SFSstock_Q75, RoundV)))
# print('Total_SFSstock_Min --> ' + str(round(Total_SFSstock_Min, RoundV)))
# print('Total_SFSstock_Max --> ' + str(round(Total_SFSstock_Max, RoundV)))

# SFS_Stock_DF = pd.DataFrame({
#     'Specific year': Speriod + startYear,
#     'SoftFallSurfaces_Stock_Mean': round(Total_SFSstock_Mean, RoundV),
#     'SoftFallSurfaces_Stock_SD': round(Total_SFSstock_SD, RoundV),
#     'SoftFallSurfaces_Stock_Q50': round(Total_SFSstock_Q50, RoundV),
#     'SoftFallSurfaces_Stock_Q25': round(Total_SFSstock_Q25, RoundV),
#     'SoftFallSurfaces_Stock_Q75': round(Total_SFSstock_Q75, RoundV),
#     'SoftFallSurfaces_Stock_MIN': round(Total_SFSstock_Min, RoundV),
#     'SoftFallSurfaces_Stock_MAX': round(Total_SFSstock_Max, RoundV)},
#     index=[0]
# )

# SFS_Stock_DF_Transp = pd.DataFrame.transpose(SFS_Stock_DF)
# toCSV_under_directory(SFS_Stock_DF_Transp, 'RESULTS_SFS_CUMULATIVE_STOCK - Tires', True, ',')


print('Immediate total release from (consumption) the inflow in current year (to annual release)')
PMC2TE_IM = simulator.getCategoryImmediateFlowFromStockSum('PMC')
PMC2TE_IM_Mean = np.mean(PMC2TE_IM[:,Speriod])
PMC2TE_IM_SD = np.std(PMC2TE_IM[:,Speriod])
PMC2TE_IM_Q50 = np.percentile(PMC2TE_IM[:,Speriod], 50)
PMC2TE_IM_Q5 = np.percentile(PMC2TE_IM[:,Speriod], 5)
PMC2TE_IM_Q15 = np.percentile(PMC2TE_IM[:,Speriod], 15)
PMC2TE_IM_Q25 = np.percentile(PMC2TE_IM[:,Speriod], 25)
PMC2TE_IM_Q75 = np.percentile(PMC2TE_IM[:,Speriod], 75)
PMC2TE_IM_Q85 = np.percentile(PMC2TE_IM[:,Speriod], 85)
PMC2TE_IM_Q95 = np.percentile(PMC2TE_IM[:,Speriod], 95)
PMC2TE_IM_Min = min(PMC2TE_IM[:,Speriod])
PMC2TE_IM_Max = max(PMC2TE_IM[:,Speriod])

print('PMC2TE_IM_Mean --> ' + str(round(PMC2TE_IM_Mean, RoundV)))
print('PMC2TE_IM_SD --> ' + str(round(PMC2TE_IM_SD, RoundV)))
print('PMC2TE_IM_Q50 --> ' + str(round(PMC2TE_IM_Q50, RoundV)))
print('PMC2TE_IM_Q5 --> ' + str(round(PMC2TE_IM_Q5, RoundV)))
print('PMC2TE_IM_Q15 --> ' + str(round(PMC2TE_IM_Q15, RoundV)))
print('PMC2TE_IM_Q25 --> ' + str(round(PMC2TE_IM_Q25, RoundV)))
print('PMC2TE_IM_Q75 --> ' + str(round(PMC2TE_IM_Q75, RoundV)))
print('PMC2TE_IM_Q85 --> ' + str(round(PMC2TE_IM_Q85, RoundV)))
print('PMC2TE_IM_Q95 --> ' + str(round(PMC2TE_IM_Q95, RoundV)))
print('PMC2TE_IM_Min --> ' + str(round(PMC2TE_IM_Min, RoundV)))
print('PMC2TE_IM_Max --> ' + str(round(PMC2TE_IM_Max, RoundV)))



print('From In-use stock to (current) annual release')
print('From PMC immediate total release')
Stock2CR = Total_PMCrel[:,Speriod] - PMC2TE_IM[:,Speriod]
Stock2CR_Mean = np.mean(Stock2CR)
Stock2CR_SD = np.std(Stock2CR)
Stock2CR_Q50 = np.percentile(Stock2CR, 50)
Stock2CR_Q5 = np.percentile(Stock2CR, 5)
Stock2CR_Q15 = np.percentile(Stock2CR, 15)
Stock2CR_Q25 = np.percentile(Stock2CR, 25)
Stock2CR_Q75 = np.percentile(Stock2CR, 75)
Stock2CR_Q85 = np.percentile(Stock2CR, 85)
Stock2CR_Q95 = np.percentile(Stock2CR, 95)
Stock2CR_Min = min(Stock2CR)
Stock2CR_Max = max(Stock2CR)

print('Stock2CR_Mean --> ' + str(round(Stock2CR_Mean, RoundV)))
print('Stock2CR_SD --> ' + str(round(Stock2CR_SD, RoundV)))
print('Stock2CR_Q50 --> ' + str(round(Stock2CR_Q50, RoundV)))
print('Stock2CR_Q5 --> ' + str(round(Stock2CR_Q5, RoundV)))
print('Stock2CR_Q15 --> ' + str(round(Stock2CR_Q15, RoundV)))
print('Stock2CR_Q25 --> ' + str(round(Stock2CR_Q25, RoundV)))
print('Stock2CR_Q75 --> ' + str(round(Stock2CR_Q75, RoundV)))
print('Stock2CR_Q85 --> ' + str(round(Stock2CR_Q85, RoundV)))
print('Stock2CR_Q95 --> ' + str(round(Stock2CR_Q95, RoundV)))
print('Stock2CR_Min --> ' + str(round(Stock2CR_Min, RoundV)))
print('Stock2CR_Max --> ' + str(round(Stock2CR_Max, RoundV)))


print('From Manufacture to Consumption')
Manu2Con = simulator.getLoggedCategoryInflows('consumption')
Manu2Con_Mean = np.mean(Manu2Con[:,Speriod])
Manu2Con_SD = np.std(Manu2Con[:,Speriod])
Manu2Con_Q50 = np.percentile(Manu2Con[:,Speriod], 50)
Manu2Con_Q5 = np.percentile(Manu2Con[:,Speriod], 5)
Manu2Con_Q15 = np.percentile(Manu2Con[:,Speriod], 15)
Manu2Con_Q25 = np.percentile(Manu2Con[:,Speriod], 25)
Manu2Con_Q75 = np.percentile(Manu2Con[:,Speriod], 75)
Manu2Con_Q85 = np.percentile(Manu2Con[:,Speriod], 85)
Manu2Con_Q95 = np.percentile(Manu2Con[:,Speriod], 95)
Manu2Con_Min = min(Manu2Con[:,Speriod])
Manu2Con_Max = max(Manu2Con[:,Speriod])

print('Manu2Con_Mean --> ' + str(round(Manu2Con_Mean, RoundV)))
print('Manu2Con_SD --> ' + str(round(Manu2Con_SD, RoundV)))
print('Manu2Con_Q50 --> ' + str(round(Manu2Con_Q50, RoundV)))
print('Manu2Con_Q5 --> ' + str(round(Manu2Con_Q5, RoundV)))
print('Manu2Con_Q15 --> ' + str(round(Manu2Con_Q15, RoundV)))
print('Manu2Con_Q25 --> ' + str(round(Manu2Con_Q25, RoundV)))
print('Manu2Con_Q75 --> ' + str(round(Manu2Con_Q75, RoundV)))
print('Manu2Con_Q85 --> ' + str(round(Manu2Con_Q85, RoundV)))
print('Manu2Con_Q95 --> ' + str(round(Manu2Con_Q95, RoundV)))
print('Manu2Con_Min --> ' + str(round(Manu2Con_Min, RoundV)))
print('Manu2Con_Max --> ' + str(round(Manu2Con_Max, RoundV)))


print('From Consumption to In-use Stock')
Cons2Stock = Manu2Con[:,Speriod] - PMC2TE_IM[:,Speriod]
Cons2Stock_Mean = np.mean(Cons2Stock)
Cons2Stock_SD = np.std(Cons2Stock)
Cons2Stock_Q50 = np.percentile(Cons2Stock, 50)
Cons2Stock_Q5 = np.percentile(Cons2Stock, 5)
Cons2Stock_Q15 = np.percentile(Cons2Stock, 15)
Cons2Stock_Q25 = np.percentile(Cons2Stock, 25)
Cons2Stock_Q75 = np.percentile(Cons2Stock, 75)
Cons2Stock_Q85 = np.percentile(Cons2Stock, 85)
Cons2Stock_Q95 = np.percentile(Cons2Stock, 95)
Cons2Stock_Min = min(Cons2Stock)
Cons2Stock_Max = max(Cons2Stock)

print('Cons2Stock_Mean --> ' + str(round(Cons2Stock_Mean, RoundV)))
print('Cons2Stock_SD --> ' + str(round(Cons2Stock_SD, RoundV)))
print('Cons2Stock_Q50 --> ' + str(round(Cons2Stock_Q50, RoundV)))
print('Cons2Stock_Q5 --> ' + str(round(Cons2Stock_Q5, RoundV)))
print('Cons2Stock_Q15 --> ' + str(round(Cons2Stock_Q15, RoundV)))
print('Cons2Stock_Q25 --> ' + str(round(Cons2Stock_Q25, RoundV)))
print('Cons2Stock_Q75 --> ' + str(round(Cons2Stock_Q75, RoundV)))
print('Cons2Stock_Q85 --> ' + str(round(Cons2Stock_Q85, RoundV)))
print('Cons2Stock_Q95 --> ' + str(round(Cons2Stock_Q95, RoundV)))
print('Cons2Stock_Min --> ' + str(round(Cons2Stock_Min, RoundV)))
print('Cons2Stock_Max --> ' + str(round(Cons2Stock_Max, RoundV)))



# In[22]:
# Sum up all inflows to sinks
## CHANGE: battW_RecMat_GBM is a FlowCompartment allocated to EoW_alternatives. Thus, only 11 SINK compartments.
# battW_RecMat_GBM_sink = Comp_with_loggedInflows['Batteries Rec_Recovered_GBM']

EoW_alternatives_sink = Comp_with_loggedInflows['End-of-Waste Alternatives']
Reuse_sink = Comp_with_loggedInflows['Reuse']
Landfill_sink = Comp_with_loggedInflows['Landfill']
Elimination_sink = Comp_with_loggedInflows['Elimination']
Export_sink = Comp_with_loggedInflows['Export']
Unknown_sink = Comp_with_loggedInflows['Unknown']
Roadside_sink = Comp_with_loggedInflows['Roadside']
STsoil_sink = Comp_with_loggedInflows['Sludge treated soil']
NUsoil_sink = Comp_with_loggedInflows['Natural&Urban soil']
Surfacewater_sink = Comp_with_loggedInflows['Surfacewater']
Subsurface_sink = Comp_with_loggedInflows['Subsurface']

sink_matrices = {
    "End-of-Waste Alternatives": EoW_alternatives_sink,
    "Reuse": Reuse_sink,
    "Landfill": Landfill_sink,
    "Elimination": Elimination_sink,
    "Export": Export_sink,
    "Unknown": Unknown_sink,
    "Roadside": Roadside_sink,
    "Sludge treated soil": STsoil_sink,
    "Natural&Urban soil": NUsoil_sink,
    "Surfacewater": Surfacewater_sink,
    "Subsurface": Subsurface_sink,
}

# Cumulative over time from the beginning of the model, 2014–2050.
export_cumulative_sink_results(
    sink_matrices=sink_matrices,
    scenario_label=f"{SCENARIO_LABEL}_{VARIANT_LABEL}",
    output_dir=dirsname,
    start_year=startYear,
    t=t,
    cumulative_start_year=None
)

# Cumulative over time for the prospective/policy period, 2025–2050.
export_cumulative_sink_results(
    sink_matrices=sink_matrices,
    scenario_label=f"{SCENARIO_LABEL}_{VARIANT_LABEL}",
    output_dir=dirsname,
    start_year=startYear,
    t=t,
    cumulative_start_year=2025
)

# Master summary for full temporal boundary, 2014–2050.
export_cumulative_sinks_master_summary(
    sink_matrices=sink_matrices,
    scenario_label=f"{SCENARIO_LABEL}_{VARIANT_LABEL}",
    output_dir=dirsname,
    start_year=startYear,
    cumulative_from_year=2014,
    cumulative_to_year=2050,
    t=t,
    save_run_level_vectors=True
)

# Master summary for prospective/policy period, 2025–2050.
export_cumulative_sinks_master_summary(
    sink_matrices=sink_matrices,
    scenario_label=f"{SCENARIO_LABEL}_{VARIANT_LABEL}",
    output_dir=dirsname,
    start_year=startYear,
    cumulative_from_year=2025,
    cumulative_to_year=2050,
    t=t,
    save_run_level_vectors=True
)

# In[23]:
# Extract means to sinks for each year
## CHANGE to 13 Sink compartments
Means_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    Means_To_Sinks[0,i]=np.mean(EoW_alternatives_sink[:,i])
    Means_To_Sinks[1,i]=np.mean(Reuse_sink[:,i])
    Means_To_Sinks[2,i]=np.mean(Landfill_sink[:,i])
    Means_To_Sinks[3,i]=np.mean(Elimination_sink[:,i])
    Means_To_Sinks[4,i]=np.mean(Export_sink[:,i])
    Means_To_Sinks[5,i]=np.mean(Unknown_sink[:,i])
    Means_To_Sinks[6,i]=np.mean(Roadside_sink [:,i])  
    Means_To_Sinks[7,i]=np.mean(STsoil_sink[:,i])
    Means_To_Sinks[8,i]=np.mean(NUsoil_sink[:,i])
    Means_To_Sinks[9,i]=np.mean(Surfacewater_sink[:,i])
    Means_To_Sinks[10,i]=np.mean(Subsurface_sink[:,i])
  
Means_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Means_To_Sinks[0,:], 
                                  'Reuse':Means_To_Sinks[1,:],
                                  'Landfill':Means_To_Sinks[2,:],
                                  'Elimination':Means_To_Sinks[3,:],
                                  'Export':Means_To_Sinks[4,:],
                                  'Unknown':Means_To_Sinks[5,:], 
                                  'Roadside':Means_To_Sinks[6,:],
                                  'Sludge treated soil':Means_To_Sinks[7,:],
                                  'Natural&Urban soil':Means_To_Sinks[8,:],
                                  'Surfacewater':Means_To_Sinks[9,:],
                                  'Subsurface':Means_To_Sinks[10,:],
                                 })

Means_To_Sinks_DF_Transp = pd.DataFrame.transpose(Means_To_Sinks_DF)
toCSV_under_directory(Means_To_Sinks_DF_Transp, 'RESULTS_Means to Sinks - S1', True, ',')


SDs_To_Sinks = np.zeros(shape=(11,Tperiods))
  
for i in np.arange(0,Tperiods):
    SDs_To_Sinks[0,i]=np.std(EoW_alternatives_sink[:,i])
    SDs_To_Sinks[1,i]=np.std(Reuse_sink[:,i])
    SDs_To_Sinks[2,i]=np.std(Landfill_sink[:,i])
    SDs_To_Sinks[3,i]=np.std(Elimination_sink[:,i])
    SDs_To_Sinks[4,i]=np.std(Export_sink[:,i])
    SDs_To_Sinks[5,i]=np.std(Unknown_sink[:,i])
    SDs_To_Sinks[6,i]=np.std(Roadside_sink [:,i])  
    SDs_To_Sinks[7,i]=np.std(STsoil_sink[:,i])
    SDs_To_Sinks[8,i]=np.std(NUsoil_sink[:,i])
    SDs_To_Sinks[9,i]=np.std(Surfacewater_sink[:,i])
    SDs_To_Sinks[10,i]=np.std(Subsurface_sink[:,i])
  
SDs_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':SDs_To_Sinks[0,:], 
                                  'Reuse':SDs_To_Sinks[1,:],
                                  'Landfill':SDs_To_Sinks[2,:],
                                  'Elimination':SDs_To_Sinks[3,:],
                                  'Export':SDs_To_Sinks[4,:],
                                  'Unknown':SDs_To_Sinks[5,:], 
                                  'Roadside':SDs_To_Sinks[6,:],
                                  'Sludge treated soil':SDs_To_Sinks[7,:],
                                  'Natural&Urban soil':SDs_To_Sinks[8,:],
                                  'Surfacewater':SDs_To_Sinks[9,:],
                                  'Subsurface':SDs_To_Sinks[10,:],
                                 })

SDs_To_Sinks_DF_Transp = pd.DataFrame.transpose(SDs_To_Sinks_DF)
toCSV_under_directory(SDs_To_Sinks_DF_Transp, 'RESULTS_SDs to Sinks - S1', True, ',')

# In[24]:
# Extract minimums to sinks for each year

Min_To_Sinks = np.zeros(shape=(11,Tperiods))
  
for i in np.arange(0,Tperiods):
    Min_To_Sinks[0,i]=min(EoW_alternatives_sink[:,i])
    Min_To_Sinks[1,i]=min(Reuse_sink[:,i])
    Min_To_Sinks[2,i]=min(Landfill_sink[:,i])
    Min_To_Sinks[3,i]=min(Elimination_sink[:,i])
    Min_To_Sinks[4,i]=min(Export_sink[:,i])
    Min_To_Sinks[5,i]=min(Unknown_sink[:,i])
    Min_To_Sinks[6,i]=min(Roadside_sink [:,i])  
    Min_To_Sinks[7,i]=min(STsoil_sink[:,i])
    Min_To_Sinks[8,i]=min(NUsoil_sink[:,i])
    Min_To_Sinks[9,i]=min(Surfacewater_sink[:,i])
    Min_To_Sinks[10,i]=min(Subsurface_sink[:,i])
  
Min_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Min_To_Sinks[0,:], 
                                  'Reuse':Min_To_Sinks[1,:],
                                  'Landfill':Min_To_Sinks[2,:],
                                  'Elimination':Min_To_Sinks[3,:],
                                  'Export':Min_To_Sinks[4,:],
                                  'Unknown':Min_To_Sinks[5,:], 
                                  'Roadside':Min_To_Sinks[6,:],
                                  'Sludge treated soil':Min_To_Sinks[7,:],
                                  'Natural&Urban soil':Min_To_Sinks[8,:],
                                  'Surfacewater':Min_To_Sinks[9,:],
                                  'Subsurface':Min_To_Sinks[10,:],
                                 })


Min_To_Sinks_DF_Transp = pd.DataFrame.transpose(Min_To_Sinks_DF)
toCSV_under_directory(Min_To_Sinks_DF_Transp, 'RESULTS_Min to Sinks - S1', True, ',')

# In[25]:
# Extract minimums to sinks for each year

Max_To_Sinks = np.zeros(shape=(11,Tperiods))
  
for i in np.arange(0,Tperiods):
    Max_To_Sinks[0,i]=max(EoW_alternatives_sink[:,i])
    Max_To_Sinks[1,i]=max(Reuse_sink[:,i])
    Max_To_Sinks[2,i]=max(Landfill_sink[:,i])
    Max_To_Sinks[3,i]=max(Elimination_sink[:,i])
    Max_To_Sinks[4,i]=max(Export_sink[:,i])
    Max_To_Sinks[5,i]=max(Unknown_sink[:,i])
    Max_To_Sinks[6,i]=max(Roadside_sink [:,i])  
    Max_To_Sinks[7,i]=max(STsoil_sink[:,i])
    Max_To_Sinks[8,i]=max(NUsoil_sink[:,i])
    Max_To_Sinks[9,i]=max(Surfacewater_sink[:,i])
    Max_To_Sinks[10,i]=max(Subsurface_sink[:,i])
  
Max_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Max_To_Sinks[0,:], 
                                  'Reuse':Max_To_Sinks[1,:],
                                  'Landfill':Max_To_Sinks[2,:],
                                  'Elimination':Max_To_Sinks[3,:],
                                  'Export':Max_To_Sinks[4,:],
                                  'Unknown':Max_To_Sinks[5,:], 
                                  'Roadside':Max_To_Sinks[6,:],
                                  'Sludge treated soil':Max_To_Sinks[7,:],
                                  'Natural&Urban soil':Max_To_Sinks[8,:],
                                  'Surfacewater':Max_To_Sinks[9,:],
                                  'Subsurface':Max_To_Sinks[10,:],
                                 })


Max_To_Sinks_DF_Transp = pd.DataFrame.transpose(Max_To_Sinks_DF)
toCSV_under_directory(Max_To_Sinks_DF_Transp, 'RESULTS_Max to Sinks - S1', True, ',')


# In[26]:
# Extract medians to sinks for each year
Medians_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    #Medians_To_Sinks[0,i]=np.percentile(reprocessing_sink[:,i],50)
    Medians_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],50)
    Medians_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],50)
    Medians_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],50)
    Medians_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],50)
    Medians_To_Sinks[4,i]= np.percentile(Export_sink[:,i],50)
    Medians_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],50)
    Medians_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],50)
    Medians_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],50)
    Medians_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],50)
    Medians_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],50)
    Medians_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],50)
  
Medians_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Medians_To_Sinks[0,:], 
                                  'Reuse':Medians_To_Sinks[1,:],
                                  'Landfill':Medians_To_Sinks[2,:],
                                  'Elimination':Medians_To_Sinks[3,:],
                                  'Export':Medians_To_Sinks[4,:],
                                  'Unknown':Medians_To_Sinks[5,:], 
                                  'Roadside':Medians_To_Sinks[6,:],
                                  'Sludge treated soil':Medians_To_Sinks[7,:],
                                  'Natural&Urban soil':Medians_To_Sinks[8,:],
                                  'Surfacewater':Medians_To_Sinks[9,:],
                                  'Subsurface':Medians_To_Sinks[10,:],
                                 })


Medians_To_Sinks_DF_Transp = pd.DataFrame.transpose(Medians_To_Sinks_DF)
toCSV_under_directory(Medians_To_Sinks_DF_Transp, 'RESULTS_Medians to Sinks - S1', True, ',')


# In[27]:
# Extract Q5 to sinks for each year
Q5_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    #Q5_To_Sinks[0,i]=np.percentile(reprocessing_sink[:,i],5)
    Q5_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],5)
    Q5_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],5)
    Q5_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],5)
    Q5_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],5)
    Q5_To_Sinks[4,i]= np.percentile(Export_sink[:,i],5)
    Q5_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],5)
    Q5_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],5)
    Q5_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],5)
    Q5_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],5)
    Q5_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],5)
    Q5_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],5)
  
Q5_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q5_To_Sinks[0,:], 
                                  'Reuse':Q5_To_Sinks[1,:],
                                  'Landfill':Q5_To_Sinks[2,:],
                                  'Elimination':Q5_To_Sinks[3,:],
                                  'Export':Q5_To_Sinks[4,:],
                                  'Unknown':Q5_To_Sinks[5,:], 
                                  'Roadside':Q5_To_Sinks[6,:],
                                  'Sludge treated soil':Q5_To_Sinks[7,:],
                                  'Natural&Urban soil':Q5_To_Sinks[8,:],
                                  'Surfacewater':Q5_To_Sinks[9,:],
                                  'Subsurface':Q5_To_Sinks[10,:],
                                 })


Q5_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q5_To_Sinks_DF)
toCSV_under_directory(Q5_To_Sinks_DF_Transp, 'RESULTS_Q5 to Sinks - S1', True, ',')


# In[28]:
# Extract Q15 to sinks for each year
Q15_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    Q15_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],15)
    Q15_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],15)
    Q15_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],15)
    Q15_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],15)
    Q15_To_Sinks[4,i]= np.percentile(Export_sink[:,i],15)
    Q15_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],15)
    Q15_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],15)
    Q15_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],15)
    Q15_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],15)
    Q15_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],15)
    Q15_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],15)
  
Q15_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q15_To_Sinks[0,:], 
                                  'Reuse':Q15_To_Sinks[1,:],
                                  'Landfill':Q15_To_Sinks[2,:],
                                  'Elimination':Q15_To_Sinks[3,:],
                                  'Export':Q15_To_Sinks[4,:],
                                  'Unknown':Q15_To_Sinks[5,:], 
                                  'Roadside':Q15_To_Sinks[6,:],
                                  'Sludge treated soil':Q15_To_Sinks[7,:],
                                  'Natural&Urban soil':Q15_To_Sinks[8,:],
                                  'Surfacewater':Q15_To_Sinks[9,:],
                                  'Subsurface':Q15_To_Sinks[10,:],
                                 })


Q15_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q15_To_Sinks_DF)
toCSV_under_directory(Q15_To_Sinks_DF_Transp, 'RESULTS_Q15 to Sinks - S1', True, ',')
    
    
# In[29]:    
# Extract Q25 to sinks for each year
Q25_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    Q25_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],25)
    Q25_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],25)
    Q25_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],25)
    Q25_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],25)
    Q25_To_Sinks[4,i]= np.percentile(Export_sink[:,i],25)
    Q25_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],25)
    Q25_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],25)
    Q25_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],25)
    Q25_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],25)
    Q25_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],25)
    Q25_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],25)
  
Q25_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q25_To_Sinks[0,:], 
                                  'Reuse':Q25_To_Sinks[1,:],
                                  'Landfill':Q25_To_Sinks[2,:],
                                  'Elimination':Q25_To_Sinks[3,:],
                                  'Export':Q25_To_Sinks[4,:],
                                  'Unknown':Q25_To_Sinks[5,:], 
                                  'Roadside':Q25_To_Sinks[6,:],
                                  'Sludge treated soil':Q25_To_Sinks[7,:],
                                  'Natural&Urban soil':Q25_To_Sinks[8,:],
                                  'Surfacewater':Q25_To_Sinks[9,:],
                                  'Subsurface':Q25_To_Sinks[10,:],
                                 })


Q25_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q25_To_Sinks_DF)
toCSV_under_directory(Q25_To_Sinks_DF_Transp, 'RESULTS_Q25 to Sinks - S1', True, ',')

# In[30]:
# Extract Q75 to sinks for each year
Q75_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    Q75_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],75)
    Q75_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],75)
    Q75_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],75)
    Q75_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],75)
    Q75_To_Sinks[4,i]= np.percentile(Export_sink[:,i],75)
    Q75_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],75)
    Q75_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],75)
    Q75_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],75)
    Q75_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],75)
    Q75_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],75)
    Q75_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],75)
  
Q75_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q75_To_Sinks[0,:], 
                                  'Reuse':Q75_To_Sinks[1,:],
                                  'Landfill':Q75_To_Sinks[2,:],
                                  'Elimination':Q75_To_Sinks[3,:],
                                  'Export':Q75_To_Sinks[4,:],
                                  'Unknown':Q75_To_Sinks[5,:], 
                                  'Roadside':Q75_To_Sinks[6,:],
                                  'Sludge treated soil':Q75_To_Sinks[7,:],
                                  'Natural&Urban soil':Q75_To_Sinks[8,:],
                                  'Surfacewater':Q75_To_Sinks[9,:],
                                  'Subsurface':Q75_To_Sinks[10,:],
                                 })


Q75_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q75_To_Sinks_DF)
toCSV_under_directory(Q75_To_Sinks_DF_Transp, 'RESULTS_Q75 to Sinks - S1', True, ',')


# In[31]:
# Extract Q85 to sinks for each year
Q85_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
  Q85_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],85)
  Q85_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],85)
  Q85_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],85)
  Q85_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],85)
  Q85_To_Sinks[4,i]= np.percentile(Export_sink[:,i],85)
  Q85_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],85)
  Q85_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],85)
  Q85_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],85)
  Q85_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],85)
  Q85_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],85)
  Q85_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],85)
  
Q85_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q85_To_Sinks[0,:], 
                                  'Reuse':Q85_To_Sinks[1,:],
                                  'Landfill':Q85_To_Sinks[2,:],
                                  'Elimination':Q85_To_Sinks[3,:],
                                  'Export':Q85_To_Sinks[4,:],
                                  'Unknown':Q85_To_Sinks[5,:], 
                                  'Roadside':Q85_To_Sinks[6,:],
                                  'Sludge treated soil':Q85_To_Sinks[7,:],
                                  'Natural&Urban soil':Q85_To_Sinks[8,:],
                                  'Surfacewater':Q85_To_Sinks[9,:],
                                  'Subsurface':Q85_To_Sinks[10,:],
                                 })


Q85_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q85_To_Sinks_DF)
toCSV_under_directory(Q85_To_Sinks_DF_Transp, 'RESULTS_Q85 to Sinks - S1', True, ',')


# In[32]:
# Extract Q95 to sinks for each year
Q95_To_Sinks = np.zeros(shape=(11,Tperiods))

for i in np.arange(0,Tperiods):
    Q95_To_Sinks[0,i]= np.percentile(EoW_alternatives_sink[:,i],95)
    Q95_To_Sinks[1,i]= np.percentile(Reuse_sink[:,i],95)
    Q95_To_Sinks[2,i]= np.percentile(Landfill_sink[:,i],95)
    Q95_To_Sinks[3,i]= np.percentile(Elimination_sink[:,i],95)
    Q95_To_Sinks[4,i]= np.percentile(Export_sink[:,i],95)
    Q95_To_Sinks[5,i]= np.percentile(Unknown_sink[:,i],95)
    Q95_To_Sinks[6,i]= np.percentile(Roadside_sink [:,i],95)
    Q95_To_Sinks[7,i]= np.percentile(STsoil_sink[:,i],95)
    Q95_To_Sinks[8,i]= np.percentile(NUsoil_sink[:,i],95)
    Q95_To_Sinks[9,i]= np.percentile(Surfacewater_sink[:,i],95)
    Q95_To_Sinks[10,i]= np.percentile(Subsurface_sink[:,i],95)
  
Q95_To_Sinks_DF = pd.DataFrame({
                                  'End-of-Waste Alternatives':Q95_To_Sinks[0,:], 
                                  'Reuse':Q95_To_Sinks[1,:],
                                  'Landfill':Q95_To_Sinks[2,:],
                                  'Elimination':Q95_To_Sinks[3,:],
                                  'Export':Q95_To_Sinks[4,:],
                                  'Unknown':Q95_To_Sinks[5,:], 
                                  'Roadside':Q95_To_Sinks[6,:],
                                  'Sludge treated soil':Q95_To_Sinks[7,:],
                                  'Natural&Urban soil':Q95_To_Sinks[8,:],
                                  'Surfacewater':Q95_To_Sinks[9,:],
                                  'Subsurface':Q95_To_Sinks[10,:],
                                 })

Q95_To_Sinks_DF_Transp = pd.DataFrame.transpose(Q95_To_Sinks_DF)
toCSV_under_directory(Q95_To_Sinks_DF_Transp, 'RESULTS_Q95 to Sinks - S1', True, ',')


### Sink dynamnics of SW, STsoil and NUsoil
SW_sink_currentyear_data = pd.DataFrame(Surfacewater_sink[:,i])
toCSV_under_directory(SW_sink_currentyear_data, 'RESULTS_Distribution to Surfacewater - S1', True, ',')

STs_sink_currentyear_data = pd.DataFrame(STsoil_sink[:,i])
toCSV_under_directory(STs_sink_currentyear_data, 'RESULTS_Distribution to STsoil - S1', True, ',')

NUs_sink_currentyear_data = pd.DataFrame(NUsoil_sink[:,i])
toCSV_under_directory(NUs_sink_currentyear_data, 'RESULTS_Distribution to NUsoil - S1', True, ',')

###
dynamic_EoW_alternatives_sink = pd.DataFrame(EoW_alternatives_sink)
dynamic_Reuse_sink = pd.DataFrame(Reuse_sink)
dynamic_Landfill_sink = pd.DataFrame(Landfill_sink)
dynamic_Elimination_sink = pd.DataFrame(Elimination_sink)
dynamic_Export_sink = pd.DataFrame(Export_sink)
dynamic_Unknown_sink = pd.DataFrame(Unknown_sink)
dynamic_Roadside_sink = pd.DataFrame(Roadside_sink)
dynamic_STsoil_sink = pd.DataFrame(STsoil_sink)
dynamic_NUsoil_sink = pd.DataFrame(NUsoil_sink)
dynamic_Surfacewater_sink = pd.DataFrame(Surfacewater_sink)
dynamic_Subsurface_sink = pd.DataFrame(Subsurface_sink)


toCSV_under_directory(dynamic_EoW_alternatives_sink, 'Results_SinkDynamic_EoW_alternatives-S1', True, ',')
toCSV_under_directory(dynamic_Reuse_sink , 'Results_SinkDynamic_Reuse-S1', True, ',')
toCSV_under_directory(dynamic_Landfill_sink, 'Results_SinkDynamic_Landfill-S1', True, ',')
toCSV_under_directory(dynamic_Elimination_sink, 'Results_SinkDynamic_Elimination-S1', True, ',')
toCSV_under_directory(dynamic_Export_sink, 'Results_SinkDynamic_Export-S1', True, ',')
toCSV_under_directory(dynamic_Unknown_sink, 'Results_SinkDynamic_Unknown-S1', True, ',')
toCSV_under_directory(dynamic_Roadside_sink, 'Results_SinkDynamic_Roadside-S1', True, ',')
toCSV_under_directory(dynamic_STsoil_sink, 'Results_SinkDynamic_STsoil-S1', True, ',')
toCSV_under_directory(dynamic_NUsoil_sink, 'Results_SinkDynamic_NUsoil-S1', True, ',')
toCSV_under_directory(dynamic_Surfacewater_sink, 'Results_SinkDynamic_Surfacewater-S1', True, ',')
toCSV_under_directory(dynamic_Subsurface_sink, 'Results_SinkDynamic_Subsurface-S1', True, ',')


### Probability distribution: product category allocation 2050 (Inflows to PC)
batteries_PC = Comp_with_loggedInflows['batteries']
tires_PC = Comp_with_loggedInflows['tires']
construction_PC = Comp_with_loggedInflows['construction']

dynamic_PC_batteries = pd.DataFrame(batteries_PC)
toCSV_under_directory(dynamic_PC_batteries, 'Results_dynamic_PC_batteries - S1', True, ',')

dynamic_PC_tires = pd.DataFrame(tires_PC)
toCSV_under_directory(dynamic_PC_tires, 'Results_dynamic_PC_tires - S1', True, ',')

dynamic_PC_construction = pd.DataFrame(construction_PC)
toCSV_under_directory(dynamic_PC_construction, 'Results_dynamic_PC_construction - S1', True, ',')

###  Time series: In-use Stock (line plot): In use stock (Sum of Inflow-outflow of each PC)
inUseStockSize_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    inUseStockSize_matrix[:,i] = Total_PMCstock[:,i]

inUseStockSize_DF = pd.DataFrame(inUseStockSize_matrix)
toCSV_under_directory(inUseStockSize_DF, 'RESULTS_inUseStockSize- S1', True, ',')

### Consumption vs Stock vs annual release: mean of Release immediate release (Outflow) vs stock (Inflow-outflow) of each year
consumption_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    consumption_matrix[:,i] = Manu2Con[:,i]

consumption_DF = pd.DataFrame(consumption_matrix)
toCSV_under_directory(consumption_DF, 'RESULTS_consumptionOfYear - S1', True, ',')


#for stock, inUseStockSize_2014_2050 which created before should be used

annualRelease_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    annualRelease_matrix[:,i] = Total_PMCrel[:,i]

annualRelease_DF = pd.DataFrame(annualRelease_matrix)
toCSV_under_directory(annualRelease_DF, 'RESULTS_annualReleaseOfYear - S1', True, ',')



