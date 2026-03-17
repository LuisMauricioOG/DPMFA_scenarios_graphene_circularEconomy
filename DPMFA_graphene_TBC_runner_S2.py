# -*- coding: utf-8 -*-
#######################################
# @author: Luis Mauricio Ortiz-Galvez
# Last modification on October 16, 2025
 
startYear = 2014 #from 2014 to 2050
Tperiods = 37     # total time range considered, also for np.arange, important for development of  Accumulation, release, Concentration
                    # if Tperiods is changed, the import settings of the production data in the model must also be changed (if Tperiods get higher)
RUNs = 10000 # 1000
Speriod = 36 # specific period for concentrations and Mass-flows
########################################



# In[2]:
# from datetime import datetime
import time
import os
import DPMFA_graphene_TBC_model_S2_20251016_final as ScII_Model # here space is not allowed in model name
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
dirsname = os.path.join('results_runner', t + "_" + tYear)

print('dirsname' + dirsname)
if not os.path.exists(dirsname):
    os.makedirs(dirsname)

## prepare function to export csv at the right location
def withOpen_under_directory(filename, mode, delimiter_setting):
    with open(dirsname + '\\' + filename + t + '.csv', mode) as RM : #w:write mode, 'a':append mode
        a = csv.writer(RM, delimiter = delimiter_setting)
        data = np.asarray(Flows)
        a.writerows(data)
# use example: withOpen_under_directory('RESULTS_AllFlows_UC2_mean_2025_2050_', 'w', ' ')

## prepare function to export csv at the right location
def toCSV_under_directory(dataname, filename, index_value, seperator):
    # index_value = True or False
    # seperator = ',' or ' '
    dataname.to_csv(dirsname + '\\' + filename + t + '.csv', index = index_value, sep = seperator)
# use example: toCSV_under_directory(Min_To_Sinks_DF_Transp, 'RESULTS_Min to Sinks - UC2', True, ',')


# In[4]:
print('Results here in console for period: ' + str(Speriod))
plt.close() #
model = ScII_Model.model
model.checkModelValidity()

# In[5]:
simulator = sc.Simulator(RUNs, Tperiods, 2250, True, True) # 2250 is just a seed
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
        np.savetxt(f'C:\\Users\\ortl\\Documents\\Empa\\Articles\\pDPMFA_AdMa_EU_CE\\ProspectiveMFA_scenarios\\results_runner\\2025_10_16_Results\\Flows_txt\\RawData_S2\\{Comp.name}TO{str(Target_name)}.txt',raw,delimiter = ',')
                   
        print('Mean--> ' + str(round(np.mean(value[:,Speriod]),RoundV)))
        print('SD--> ' + str(round(np.std(value[:,Speriod]),RoundV)))
        print('Q5--> ' + str(round(np.percentile(value[:,Speriod],5),RoundV)))
        print('Q15--> ' + str(round(np.percentile(value[:,Speriod],15),RoundV)))
        print('Q25--> ' + str(round(np.percentile(value[:,Speriod],25),RoundV)))
        print('Q50--> ' + str(round(np.percentile(value[:,Speriod],50),RoundV)))
        print('Q5--> ' + str(round(np.percentile(value[:,Speriod],75),RoundV)))
        print('Q5--> ' + str(round(np.percentile(value[:,Speriod],85),RoundV)))
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

    withOpen_under_directory('RESULTS_AllFlows_S2_mean_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_std_2014_2050', 'w', ' ')


# In[12]:

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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q5_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q15_2014_2050', 'w', ' ')    
    
    
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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q25_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q75_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q85_2014_2050', 'w', ' ')
    
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
    withOpen_under_directory('RESULTS_AllFlows_S2_Q95_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_min_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_max_2014_2050', 'w', ' ')

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
    withOpen_under_directory('RESULTS_AllFlows_S2_medians_2014_2050', 'w', ' ')


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
toCSV_under_directory(PMC_Stock_DF_Transp, 'RESULTS_PMC_CUMULATIVE_STOCK - S2', True, ',')

    
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
toCSV_under_directory(Means_To_Sinks_DF_Transp, 'RESULTS_Means to Sinks - S2', True, ',')


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
toCSV_under_directory(SDs_To_Sinks_DF_Transp, 'RESULTS_SDs to Sinks - S2', True, ',')

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
toCSV_under_directory(Min_To_Sinks_DF_Transp, 'RESULTS_Min to Sinks - S2', True, ',')

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
toCSV_under_directory(Max_To_Sinks_DF_Transp, 'RESULTS_Max to Sinks - S2', True, ',')


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
    Medians_To_Sinks[9,i]= np.percentile(NUsoil_sink[:,i],50)
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
toCSV_under_directory(Medians_To_Sinks_DF_Transp, 'RESULTS_Medians to Sinks - S2', True, ',')


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
toCSV_under_directory(Q5_To_Sinks_DF_Transp, 'RESULTS_Q5 to Sinks - S2', True, ',')


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
toCSV_under_directory(Q15_To_Sinks_DF_Transp, 'RESULTS_Q15 to Sinks - S2', True, ',')
    
    
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
toCSV_under_directory(Q25_To_Sinks_DF_Transp, 'RESULTS_Q25 to Sinks - S2', True, ',')

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
toCSV_under_directory(Q75_To_Sinks_DF_Transp, 'RESULTS_Q75 to Sinks - S2', True, ',')


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
toCSV_under_directory(Q85_To_Sinks_DF_Transp, 'RESULTS_Q85 to Sinks - S2', True, ',')


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
toCSV_under_directory(Q95_To_Sinks_DF_Transp, 'RESULTS_Q95 to Sinks - S2', True, ',')


### Sink dynamnics of SW, STsoil and NUsoil
SW_sink_currentyear_data = pd.DataFrame(Surfacewater_sink[:,i])
toCSV_under_directory(SW_sink_currentyear_data, 'RESULTS_Distribution to Surfacewater - S2', True, ',')

STs_sink_currentyear_data = pd.DataFrame(STsoil_sink[:,i])
toCSV_under_directory(STs_sink_currentyear_data, 'RESULTS_Distribution to STsoil - S2', True, ',')

NUs_sink_currentyear_data = pd.DataFrame(NUsoil_sink[:,i])
toCSV_under_directory(NUs_sink_currentyear_data, 'RESULTS_Distribution to NUsoil - S2', True, ',')

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


toCSV_under_directory(dynamic_EoW_alternatives_sink, 'Results_SinkDynamic_EoW_alternatives-S2', True, ',')
toCSV_under_directory(dynamic_Reuse_sink , 'Results_SinkDynamic_Reuse-S2', True, ',')
toCSV_under_directory(dynamic_Landfill_sink, 'Results_SinkDynamic_Landfill-S2', True, ',')
toCSV_under_directory(dynamic_Elimination_sink, 'Results_SinkDynamic_Elimination-S2', True, ',')
toCSV_under_directory(dynamic_Export_sink, 'Results_SinkDynamic_Export-S2', True, ',')
toCSV_under_directory(dynamic_Unknown_sink, 'Results_SinkDynamic_Unknown-S2', True, ',')
toCSV_under_directory(dynamic_Roadside_sink, 'Results_SinkDynamic_Roadside-S2', True, ',')
toCSV_under_directory(dynamic_STsoil_sink, 'Results_SinkDynamic_STsoil-S2', True, ',')
toCSV_under_directory(dynamic_NUsoil_sink, 'Results_SinkDynamic_NUsoil-S2', True, ',')
toCSV_under_directory(dynamic_Surfacewater_sink, 'Results_SinkDynamic_Surfacewater-S2', True, ',')
toCSV_under_directory(dynamic_Subsurface_sink, 'Results_SinkDynamic_Subsurface-S2', True, ',')


### Probability distribution: product category allocation 2050 (Inflows to PC)
batteries_PC = Comp_with_loggedInflows['batteries']
tires_PC = Comp_with_loggedInflows['tires']
construction_PC = Comp_with_loggedInflows['construction']

dynamic_PC_batteries = pd.DataFrame(batteries_PC)
toCSV_under_directory(dynamic_PC_batteries, 'Results_dynamic_PC_batteries - S2', True, ',')

dynamic_PC_tires = pd.DataFrame(tires_PC)
toCSV_under_directory(dynamic_PC_tires, 'Results_dynamic_PC_tires - S2', True, ',')

dynamic_PC_construction = pd.DataFrame(construction_PC)
toCSV_under_directory(dynamic_PC_construction, 'Results_dynamic_PC_construction - S2', True, ',')

###  Time series: In-use Stock (line plot): In use stock (Sum of Inflow-outflow of each PC)
inUseStockSize_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    inUseStockSize_matrix[:,i] = Total_PMCstock[:,i]

inUseStockSize_DF = pd.DataFrame(inUseStockSize_matrix)
toCSV_under_directory(inUseStockSize_DF, 'RESULTS_inUseStockSize- S2', True, ',')

### Consumption vs Stock vs annual release: mean of Release immediate release (Outflow) vs stock (Inflow-outflow) of each year
consumption_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    consumption_matrix[:,i] = Manu2Con[:,i]

consumption_DF = pd.DataFrame(consumption_matrix)
toCSV_under_directory(consumption_DF, 'RESULTS_consumptionOfYear - S2', True, ',')


#for stock, inUseStockSize_2014_2050 which created before should be used

annualRelease_matrix = np.zeros(shape=(RUNs,Tperiods))
for i in np.arange(0,Tperiods):
    annualRelease_matrix[:,i] = Total_PMCrel[:,i]

annualRelease_DF = pd.DataFrame(annualRelease_matrix)
toCSV_under_directory(annualRelease_DF, 'RESULTS_annualReleaseOfYear - S2', True, ',')



