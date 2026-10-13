# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:54:28 2025

@author: Luis Mauricio Ortiz-Galvez, Lora Dameska (MSc. thesis)

Description: This DPMFA code works to further investigate the mass flows 
entering the recycling streams. However, there is specific information needed 
to adapt it to a specific advanced material in a product, and then into the 
corresponding waste category. Therefore, the practitioner would need to 
check the decision tree and recycling streams diagrams developed to get 
guidance on potential flows on their corresponding case study.

System boundaries: EU27(+CH+NO+UK), focus on recycling streams from 2025 to 2050.
Based on the graphene model (H. Hong et al., Environ. Sci. Technol. 2022, 56, 19, 13798–13809.
DOI: 10.1021/acs.est.2c04002)

Note: retrospective mass flows from the market entry year (i.e., batteries in 2014) to 2025 are considered. 
"""
import sys
import os

project_path = r'C:\Users\ortl\Documents\Empa\Articles\pDPMFA_AdMa_EU_CE\ProspectiveMFA_scenarios'

if project_path not in sys.path:
    sys.path.insert(0, project_path)
    
import components as cp
import numpy.random as nr
import numpy as np
import Model
import pandas as pd
from scipy.stats import norm
from matplotlib import pyplot as plt
import TruncatingFunctions as tf
import TriangularBalanced as tb


###define model
model = Model.Model('Extended MFA of graphene in Tires-Batteries-Construction from 2014 to 2050, Scenario 1')

###PMC
# split of the total inflow to different compartments
totalInflow = cp.FlowCompartment('Total Inflow', logInflows=True, logOutflows=True)

# Production of the nano-material and Manufacturing of the products it is used in
production = cp.Stock('production',  logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['production'])
manufacture = cp.Stock('manufacture', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['manufacture'])
consumption = cp.FlowCompartment('consumption', logInflows=True, logOutflows=True, categories=['consumption'])

### Definition of product categories (compartments)
batteries = cp.FlowCompartment('batteries', logInflows=True, logOutflows=True, categories=['Products'])
batteries_Use = cp.Stock('batteries_Use', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['PMC','Use','PMCrel'])
batteries_EoL = cp.Stock('batteries_EoL', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['PMC','EoL'])

tires = cp.FlowCompartment('tires', logInflows=True,logOutflows=True,categories=['Products'])
tires_Use = cp.Stock('tires_Use', logInflows=True,logOutflows=True,logImmediateFlows =True,categories=['PMC','Use','PMCrel'])
tires_EoL = cp.Stock('tires_EoL', logInflows=True,logOutflows=True,logImmediateFlows =True,categories=['PMC','EoL'])

construction = cp.FlowCompartment('construction', logInflows=True, logOutflows=True, categories=['Products'])
construction_Use =cp.Stock('construction_Use', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['PMC','Use','PMCrel'])
construction_EoL =cp.Stock('construction_EoL', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['PMC','EoL'])
#construction refers to the application of graphene in concrete


#### Compartments for solid waste flows
### Definition of solid waste categories
#CDW stands for Construction and Demolition waste
CDW = cp.FlowCompartment('CDW', logInflows=True,logOutflows=True, categories =['WasteStream', 'SolidWaste', 'PMCrel'])
tiresW = cp.FlowCompartment('tiresW', logInflows=True,logOutflows=True, categories =['WasteStream','SolidWaste', 'PMCrel'])
battW = cp.FlowCompartment('battW', logInflows=True,logOutflows=True, categories =['WasteStream','SolidWaste', 'PMCrel'])


### Definition of Batteries EoL sorting and Recycling
#From EoL, the undefined streams are equal to 1 - batteries_collection (TimeDep), corresponding to unknown, disposal, and export. ADD EXPLANATION (based on portable batteries and accumulators eurostat's statistics and MACRAME UC2 - case study on BTMS of EVs)
# CHANGE: Name of battW_Undefined to battW_Unidentified
battW_Unidentified = cp.FlowCompartment('Batteries Unidentified Streams', logInflows=True, logOutflows=True, categories=['Batteries Unidentified Streams'])
battW_Dismantling = cp.FlowCompartment('Batteries Dismantling', logInflows=True, logOutflows=True, categories=['Batteries Dismantling'])
# Reuse of batteries can be added in the future
battW_Recycling = cp.FlowCompartment('Batteries Recycling', logInflows=True, logOutflows=True, categories=['Batteries Recycling'])

battW_SizeReduction_like = cp.FlowCompartment('Batteries Size Reduction', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Air'])

# Definition of pyrometallurgical recycling process' compartment and subcompartments
battW_Pyro = cp.FlowCompartment('Batteries Pyro', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Pyro'])
battW_Pyro_like = cp.FlowCompartment('Batteries Pyro to Elimination', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Pyro', 'Elimination'])
battW_Pyro_losses = cp.FlowCompartment('Batteries Pyro_losses', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Pyro', 'Batteries Rec_Residue'])
battW_Pyro_output = cp.FlowCompartment('Batteries Pyro_output', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Pyro', 'Batteries Rec_RecoveredMat'])

# Definition of hydrometallurgical (hydroleaching) recycling process' compartment and subcompartments
battW_HydroLeaching = cp.FlowCompartment('Batteries Hydro(Leaching)', logInflows=True, logOutflows=True, categories=['Batteries Recycling'])
battW_HydroLeaching_like = cp.FlowCompartment('Batteries Hydro(Leaching) to Elimination', logInflows=True, logOutflows=True, categories=['Batteries Recycling','Batteries Hydro(Leaching)', 'Elimination'])
battW_HydroLeaching_losses = cp.FlowCompartment('Batteries Hydro(Leaching)_losses', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Hydro(Leaching)', 'Batteries Rec_Residue'])
battW_HydroLeaching_output = cp.FlowCompartment('Batteries Hydro(Leaching)_output', logInflows=True, logOutflows=True, categories=['Batteries Recycling', 'Batteries Hydro(Leaching)', 'Batteries Rec_RecoveredMat'])

# Definition of direct recycling process' compartment (for scenario 2 and 3)
#battW_DirectRecycling = cp.FlowCompartment('Batteries Direct Recycling', logInflows=True, logOutflows=True, categories=['Batteries Recycling'])

# Definition of battery recycling residues' and recovered material/electrode's compartments
## CHANGE: allocated to Unknown (now as FlowCompartment)
battW_Recycling_Residue = cp.FlowCompartment('Batteries Recycling_Residue', logInflows=True, logOutflows=True, categories=['Batteries Rec_Residue', 'Reuse'])

## CHANGE: disaggregated: before was Recovered Material(RecMat), now it is Secondary Material with GBM impurities (RecMat_impurities) + Recovered GBM (RecMat_GBM)
## change: to FlowCompartment (allocated to EoW_alternatives)
# battW_RecMat_GBM = cp.FlowCompartment('Batteries Rec_Recovered_GBM', logInflows=True, logOutflows=True, categories=['Batteries Rec_RecoveredMat','Batteries Rec_Recovered_GBM'])
## CHANGE: battW_RecMat_impurities as FlowCompartment 
battW_RecMat_impurities = cp.FlowCompartment('Batteries Rec_RecoveredMat_impurities', logInflows=True, logOutflows=True, categories=['Batteries Rec_RecoveredMat', 'Batteries Rec_RecoveredMat_impurities'])

### Definition of Construction EoL sorting and Recycling (mainly concrete/mineral waste)
sorting_CDW_Mineral = cp.FlowCompartment('sorting_CDW_Mineral', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting1'])

## Definition of CDW - Mineral waste subcompartments (energyRec cover WIP and energy recovery)
## CHANGE: backfilling of mineral waste is allocated to End-of-Waste alternatives (FlowCompartment)
CDW_Mineral_Backf = cp.FlowCompartment('CDW_Mineral waste to Backfilling', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2'])
# CDW_Mineral_disposal = cp.FlowCompartment('CDW_Mineral_Disposal', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2'])
CDW_Mineral_CementKilns = cp.FlowCompartment('CDW_Mineral to Cement Kilns', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2'])
recycling_CDW_Mineral = cp.FlowCompartment('Mineral waste to Recycling', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2', 'Recycling CDW'])

reprocessing_CDW_Mineral = cp.FlowCompartment('Mineral waste to reprocess',logInflows=True, logOutflows=True, categories=['Recycling CDW', 'Reprocessed Minerals'])

## CHANGE: Add this compartment
Backfilling = cp.FlowCompartment('Backfilling', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2', 'Backfilling'])

## ADD (Recycled) Aggregates. Allocated to cementFill due to use in mixed concrete (FlowCompartment)
Aggregates = cp.FlowCompartment('Recycled Aggregates', logInflows=True, logOutflows=True, categories=['Recycling CDW', 'Reprocessed Minerals', 'Recycled Aggregates'])

## CHANGE: cementFill is allocated to EoW_alternatives (FlowCompartment)
cementFill = cp.FlowCompartment('Cement Filler', logInflows=True, logOutflows=True, categories=['Cement Filler'])


## Definition of Tires waste EoL and Recycling streams
#Tires_reuse cover direct reuse, retreading, and regrooving
#Tires_unknown covers improper disposal, undetermined, and storage
tiresW_preRecycling = cp.FlowCompartment('Tires preRecycling', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting1'])

tiresW_EnergyRec = cp.FlowCompartment('Tires Energy Recovery', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2', 'Energy Recovery'])
cementKilns = cp.FlowCompartment('Cement Kilns', logInflows=True, logOutflows=True, categories=['Cement Kilns'])
## CHANGE: energyRecovery To FlowCompartment, and allocated to WIP (eliminating the Likelihood approach, and summing up with WIP)
energyRecovery = cp.FlowCompartment('Energy recovery', logInflows=True, logOutflows=True, categories=['Energy Recovery'])

tiresW_MatRecovery = cp.FlowCompartment('Tires Material Recovery', logInflows=True, logOutflows=True, categories=['Sorting', 'Sorting2', 'Material Recovery'])
## CHANGE: CivilEng To FlowCompartment
CivilEng = cp.FlowCompartment('Tires in Civil_Eng', logInflows=True, logOutflows=True, categories=['Material Recovery', 'Civil Eng'])
tiresW_Recycling = cp.FlowCompartment('Tires Recycling', logInflows=True, logOutflows=True, categories=['Material Recovery', 'Tires Recycling'])


Pyrolysis = cp.FlowCompartment('Tires Pyrolysis', logInflows=True, logOutflows=True, categories=['Tires Recycling', 'Pyrolysis'])
## CHANGE: CarbonBlack To FlowCompartment
CarbonBlack = cp.FlowCompartment('Tires Carbon_Black', logInflows=True, logOutflows=True, categories=['Carbon Black', 'Pyrolysis', 'Material Recovery'])


MechanicalRecycling_sorting  = cp.FlowCompartment('Tires Mechanical_Recycling', logInflows=True, logOutflows=True, categories=['Tires Recycling', 'Tires Granulation'])
tiresW_MechanicalRecycling_like  =  cp.FlowCompartment('Tires Mechanical_Recycling to Elimination', logInflows=True, logOutflows=True, categories=['Tires Recycling', 'Tires Granulation', 'Elimination'])
        
## CHANGE: DockFenders To FlowCompartment
DockFenders = cp.FlowCompartment('Tires in Dock Fenders', logInflows=True, logOutflows=True, categories=['Tires Recycling', 'Tires in Dock Fenders'])


Reuse = cp.Sink('Reuse', logInflows=True, categories=['Reuse'])
# In the future, it should be a stock for New Tires, and other products
#  Reuse = cp.Stock ('Reuse', logInflows=True, logOutflows=True, logImmediateFlows=True, categories=['Reuse'])

## CHANGE: OtherUses To FlowCompartment, allocated to Reuse
OtherUses = cp.FlowCompartment('Tires Other_Uses', logInflows=True, logOutflows=True, categories=['Tires Granulation', 'Tires Other_Uses'])

## CHANGE: MouldedProducts To FlowCompartment
MouldedProducts = cp.FlowCompartment('Tires Moulded_Products', logInflows=True, logOutflows=True, categories=['Tires Granulation', 'Tires MouldedProducts'])

## CHANGE: SoftFallSurfaces To FlowCompartment
SoftFallSurfaces = cp.FlowCompartment('Tires SoftFall_Surfaces', logInflows=True, logOutflows=True, categories=['Tires Granulation', 'Tires SoftFall_Surfaces'])
#SoftFallSurfaces covers artificial turfs, and sports and children playgrounds

## CHANGE: Pavings To FlowCompartment
Pavings = cp.FlowCompartment('Pavings', logInflows=True, logOutflows=True, categories=['Tires Granulation','Pavings'])

## ADD: End-of-Waste alternatives, referring to the concept of the WFD. 
EoW_alternatives = cp.Sink('End-of-Waste Alternatives', logInflows=True, categories=['End-of-Waste Alternatives'])

# Other waste management and environmental compartments
Landfill = cp.Sink('Landfill', logInflows=True, categories=['Landfill'])
Elimination = cp.Sink('Elimination', logInflows=True, categories=['Elimination'])
Export = cp.Sink('Export', logInflows=True, categories=['Export'])
Unknown = cp.Sink('Unknown', logInflows=True, categories=['Unknown'])

NUsoil = cp.Sink('Natural&Urban soil', logInflows=True,categories=['Soil', 'NU Soil'])
STsoil = cp.Sink('Sludge treated soil', logInflows=True, categories=['Soil', 'ST Soil'])

Air = cp.FlowCompartment('Air', logInflows=True, logOutflows=True, categories=['Air'])
Surfacewater = cp.Sink('Surfacewater', logInflows=True, categories=['Surfacewater'])
Subsurface = cp.Sink('Subsurface', logInflows=True, categories=['Subsurface'])
Roadside = cp.Sink('Roadside', logInflows=True, categories=['Roadside'])

# WIP, STP and WWTP
WIP = cp.FlowCompartment('Waste Incineration Plant', logInflows=True, logOutflows=True, categories=['WIP'])
Burning = cp.FlowCompartment('Burning', logInflows=True, logOutflows=True, categories=['WIP'])
Remains_burning = cp.FlowCompartment('Remains_burning', logInflows=True, logOutflows=True, categories=['WIP'])
Bottomash = cp.FlowCompartment('Bottom ash', logInflows=True, logOutflows=True, categories=['WIP'])
Flyash = cp.FlowCompartment('Flyash', logInflows=True, logOutflows=True, categories=['WIP'])
Filterash = cp.FlowCompartment('Filterash', logInflows=True, logOutflows=True, categories=['WIP'])

# WWTP and onsite treatments
WWTP = cp.FlowCompartment('WWTP', logInflows=True, logOutflows=True, categories=['WWTP'])
TreatIOnsite = cp.FlowCompartment('TreatIOnsite', logInflows=True, logOutflows=True, categories=['Treat', 'TreatIOnsite'])
TreatI = cp.FlowCompartment('TreatI', logInflows=True, logOutflows=True, categories=['STP', 'Treat', 'TreatI'])
TreatII = cp.FlowCompartment('TreatII', logInflows=True, logOutflows=True, categories=['STP', 'Treat', 'TreatII'])
TreatIII = cp.FlowCompartment('TreatIII', logInflows=True, logOutflows=True, categories=['STP', 'Treat', 'TreatIII'])


## CHANGE: OnsiteSludge allocated to STPSludge
OnsiteSludge = cp.FlowCompartment('OnsiteSludge', logInflows=True, logOutflows=True, categories=['onsiteTreat', 'OnsiteSludge', 'STPSludge'])
Wastewater = cp.FlowCompartment('Wastewater', logInflows=True, logOutflows=True, categories=['AllWastewater'])
NoSewageSystem = cp.FlowCompartment('NoSewageSystem', logInflows=True, logOutflows=True, categories=['NoSewageSystem'])
SewageSystem = cp.FlowCompartment('SewageSystem', logInflows=True, logOutflows=True, categories=['SewageSystem'])
onsiteTreat = cp.FlowCompartment('onsiteTreat', logInflows=True, logOutflows=True, categories=['onsiteTreat'])
STPoverflow = cp.FlowCompartment('STPoverflow', logInflows=True, logOutflows=True, categories=['STP', 'Overflow', 'Surfacewater'])
STPeffluent = cp.FlowCompartment('STPeffluent', logInflows=True, logOutflows=True, categories=['STP', 'STPeffluent', 'Surfacewater'])
STPsludge = cp.FlowCompartment('STPsludge', logInflows=True, logOutflows=True, categories=['STP', 'STPSludge'])

## CHANGE: reprocessing To FlowCompartment, and allocated to Reuse
reprocessing = cp.FlowCompartment('Reprocessing from ashes', logInflows=True, logOutflows=True, categories=['Reprocessing Ashes'])

#ADD ALL
compartmentList = [totalInflow, production, manufacture, consumption,
                   batteries, batteries_Use, batteries_EoL,
                   tires, tires_Use, tires_EoL,
                   construction, construction_Use, construction_EoL,
                   
                   battW, battW_Unidentified, battW_Dismantling,
                   
                   battW_Recycling, battW_SizeReduction_like,
                   #battW_DirectRecycling,
                   battW_HydroLeaching, battW_HydroLeaching_like, battW_HydroLeaching_losses, battW_HydroLeaching_output, 
                   battW_Pyro, battW_Pyro_like, battW_Pyro_losses, battW_Pyro_output, 
                   #battW_RecMat_GBM, 
                   battW_RecMat_impurities, battW_Recycling_Residue, 
            
                   tiresW, tiresW_preRecycling, 
                   
                   tiresW_EnergyRec, cementKilns,
                   
                   tiresW_MatRecovery, tiresW_Recycling, CivilEng, 
                   tiresW_MechanicalRecycling_like, MechanicalRecycling_sorting,
                   Pyrolysis, CarbonBlack,
                   
                   DockFenders, OtherUses, MouldedProducts, SoftFallSurfaces,
                   
                   CDW, sorting_CDW_Mineral,
                   CDW_Mineral_Backf, CDW_Mineral_CementKilns, 
                   recycling_CDW_Mineral, reprocessing_CDW_Mineral, 
                   Backfilling, Aggregates, cementFill,
                   
                   Landfill, Elimination, Export, Unknown, Reuse, EoW_alternatives,

                   Roadside, Pavings, energyRecovery,

                   NUsoil, STsoil,

                   Wastewater, Air, Surfacewater, 
                   
                   WIP, Burning, Remains_burning, Bottomash, Flyash, Filterash,
                   WWTP, TreatI, TreatII, TreatIII, TreatIOnsite, OnsiteSludge,
                   NoSewageSystem, SewageSystem, onsiteTreat,
                   STPoverflow, STPeffluent, STPsludge,
                   reprocessing, Subsurface
                   ]



model.setCompartments(compartmentList)

# Load production data
productionData_TBC = pd.read_csv('C:\\Users\\ortl\\Documents\\Empa\\Articles\\pDPMFA_AdMa_EU_CE\\ProspectiveMFA_scenarios\\inputData\\for_py_graphene_TBC.csv', header=0) 
# 2014 -> 2050
# From production data, set as inflow
productionData_TBC_EU = productionData_TBC.values

Productionvolume = []
for i in np.arange(0,37):
    Productionvolume.append((productionData_TBC_EU[:, i])) 
periodRange = np.arange(0,37)
model.addInflow(cp.ExternalListInflow(totalInflow, [cp.RandomChoiceInflow(Productionvolume[x]) for x in periodRange]))

# productionData_EU = productionData.values #in earlier code as " Productiondata = Productiondata.as_matrix(Productiondata)

# Productionvolume=[]
# for i in np.arange(0,6):
#     Productionvolume.append(productionData_EU[:,i])
# periodRange = np.arange(0,6)
# model.addInflow(cp.ExternalListInflow(totalInflow, [cp.RandomChoiceInflow(Productionvolume[x]) for x in periodRange]))

## CHANGE (CONTINUE HERE)
# PMC transfer coefficients from previous studies (Hong et al., 2022)
totalInflow.transfers = [cp.ConstTransfer(1, production)]
## Production and Manufacture
production.localRelease = cp.ListRelease([1.0])

manufacture.localRelease = cp.ListRelease([1.0])

production.transfers = [cp.StochasticTransfer(nr.triangular, [0.371*0.005, 0.75*0.005, 1*0.005], Wastewater, priority=2),
                        cp.StochasticTransfer(nr.triangular, [0.125*0.005, 0.25*0.005, 0.375*0.005], Air, priority=2),
                        cp.ConstTransfer(1, manufacture, priority=1)]

manufacture.transfers = [cp.StochasticTransfer(nr.triangular, [0.165*(0.005/0.995), 0.33*(0.005/0.995), 0.495*(0.005/0.995)], Wastewater, priority=2),
                         cp.StochasticTransfer(nr.triangular, [0.165*(0.005/0.995), 0.33*(0.005/0.995), 0.495*(0.005/0.995)], WIP, priority=2),
                         cp.StochasticTransfer(nr.triangular, [0.17*(0.005/0.995), 0.34*(0.005/0.995), 0.51*(0.005/0.995)], Air, priority=2),
                         cp.ConstTransfer(1, consumption, priority=1)]


################## Sampling of the removal data from different sources ########
s = 1000 #sampling size used for sampling in uncertainty ranges of several sources/values

############# Dynamic product categories ##########
# Batteries
# Start in 2014 (Hong et al., 2022) but adapted to only cover the Tires-Batteries-Construction system (equal to 1 per year)
batteries_y0 = cp.TransferConstant(1) #2014
batteries_y1 = cp.TransferDistribution(tf.TriangTruncDet, [0.9797*0.5, 0.9797, 0.9797*1.5, 1, 0, 1])
batteries_y2 = cp.TransferDistribution(tf.TriangTruncDet, [0.9731*0.5, 0.9731, 0.9731*1.5, 1, 0, 1])
batteries_y3 = cp.TransferDistribution(tf.TriangTruncDet, [0.9698*0.5, 0.9698, 0.9698*1.5, 1, 0, 1])
batteries_y4 = cp.TransferDistribution(tf.TriangTruncDet, [0.9556*0.5, 0.9556, 0.9556*1.5, 1, 0, 1])
batteries_y5 = cp.TransferDistribution(tf.TriangTruncDet, [0.9447*0.5, 0.9447, 0.9447*1.5, 1, 0, 1])
batteries_y6 = cp.TransferDistribution(tf.TriangTruncDet, [0.9135*0.5, 0.9135, 0.9135*1.5, 1, 0, 1])
batteries_y7 = cp.TransferDistribution(tf.TriangTruncDet, [0.8824*0.5, 0.8824, 0.8824*1.5, 1, 0, 1])
batteries_y8 = cp.TransferDistribution(tf.TriangTruncDet, [0.8525*0.5, 0.8525, 0.8525*1.5, 1, 0, 1])
batteries_y9 = cp.TransferDistribution(tf.TriangTruncDet, [0.8186*0.5, 0.8186, 0.8186*1.5, 1, 0, 1])
batteries_y10 = cp.TransferDistribution(tf.TriangTruncDet, [0.7856*0.5, 0.7856, 0.7856*1.5, 1, 0, 1])
batteries_y11 = cp.TransferDistribution(tf.TriangTruncDet, [0.7672*0.5, 0.7672, 0.7672*1.5, 1, 0, 1])
batteries_y12 = cp.TransferDistribution(tf.TriangTruncDet, [0.7502*0.5, 0.7502, 0.7502*1.5, 1, 0, 1])
batteries_y13 = cp.TransferDistribution(tf.TriangTruncDet, [0.7398*0.5, 0.7398, 0.7398*1.5, 1, 0, 1])
batteries_y14 = cp.TransferDistribution(tf.TriangTruncDet, [0.7462*0.5, 0.7462, 0.7462*1.5, 1, 0, 1])
batteries_y15 = cp.TransferDistribution(tf.TriangTruncDet, [0.7151*0.5, 0.7151, 0.7151*1.5, 1, 0, 1])
batteries_y16 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
# From 2030 to 2050, the TC is the same as the one from 2030 (y16)
batteries_y17 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y18 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y19 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y20 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y21 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y22 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y23 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y24 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y25 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y26 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y27 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y28 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y29 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y30 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y31 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y32 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y33 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y34 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y35 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])
batteries_y36 = cp.TransferDistribution(tf.TriangTruncDet, [0.7259*0.5, 0.7259, 0.7259*1.5, 1, 0, 1])

# Tires
# Start in 2015 (Hong et al., 2022)
tires_y0 = cp.TransferConstant(0)
tires_y1 = cp.TransferDistribution(tf.TriangTruncDet, [0.0203*0.5, 0.0203, 0.0203*1.5, 1, 0, 1])
tires_y2 = cp.TransferDistribution(tf.TriangTruncDet, [0.0269*0.5, 0.0269, 0.0269*1.5, 1, 0, 1])
tires_y3 = cp.TransferDistribution(tf.TriangTruncDet, [0.0302*0.5, 0.0302, 0.0302*1.5, 1, 0, 1])
tires_y4 = cp.TransferDistribution(tf.TriangTruncDet, [0.0317*0.5, 0.0317, 0.0317*1.5, 1, 0, 1])
tires_y5 = cp.TransferDistribution(tf.TriangTruncDet, [0.0418*0.5, 0.0418, 0.0418*1.5, 1, 0, 1])
tires_y6 = cp.TransferDistribution(tf.TriangTruncDet, [0.0649*0.5, 0.0649, 0.0649*1.5, 1, 0, 1])
tires_y7 = cp.TransferDistribution(tf.TriangTruncDet, [0.0891*0.5, 0.0891, 0.0891*1.5, 1, 0, 1])
tires_y8 = cp.TransferDistribution(tf.TriangTruncDet, [0.1033*0.5, 0.1033, 0.1033*1.5, 1, 0, 1])
tires_y9 = cp.TransferDistribution(tf.TriangTruncDet, [0.1134*0.5, 0.1134, 0.1134*1.5, 1, 0, 1])
tires_y10 = cp.TransferDistribution(tf.TriangTruncDet, [0.1367*0.5, 0.1367, 0.1367*1.5, 1, 0, 1])
tires_y11 = cp.TransferDistribution(tf.TriangTruncDet, [0.1397*0.5, 0.1397, 0.1397*1.5, 1, 0, 1])
tires_y12 = cp.TransferDistribution(tf.TriangTruncDet, [0.1353*0.5, 0.1353, 0.1353*1.5, 1, 0, 1])
tires_y13 = cp.TransferDistribution(tf.TriangTruncDet, [0.1260*0.5, 0.1260, 0.1260*1.5, 1, 0, 1])
tires_y14 = cp.TransferDistribution(tf.TriangTruncDet, [0.1106*0.5, 0.1106, 0.1106*1.5, 1, 0, 1])
tires_y15 = cp.TransferDistribution(tf.TriangTruncDet, [0.1041*0.5, 0.1041, 0.1041*1.5, 1, 0, 1])
tires_y16 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
# From 2030 to 2050, the TC is the same as the one from 2030 (y16)
tires_y17 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y18 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y19 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y20 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y21 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y22 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y23 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y24 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y25 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y26 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y27 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y28 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y29 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y30 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y31 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y32 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y33 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y34 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y35 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])
tires_y36 = cp.TransferDistribution(tf.TriangTruncDet, [0.1000*0.5, 0.1000, 0.1000*1.5, 1, 0, 1])


# Construction
# Start in 2018 (Hong et al., 2022)
construction_y0 = cp.TransferConstant(0)
construction_y1 = cp.TransferConstant(0)
construction_y2 = cp.TransferConstant(0)
construction_y3 = cp.TransferConstant(0)
construction_y4 = cp.TransferDistribution(tf.TriangTruncDet, [0.0127*0.5, 0.0127, 0.0127*1.5, 1, 0, 1])
construction_y5 = cp.TransferDistribution(tf.TriangTruncDet, [0.0136*0.5, 0.0136, 0.0136*1.5, 1, 0, 1])
construction_y6 = cp.TransferDistribution(tf.TriangTruncDet, [0.0216*0.5, 0.02164, 0.0216*1.5, 1, 0, 1])
construction_y7 = cp.TransferDistribution(tf.TriangTruncDet, [0.0285*0.5, 0.0285, 0.0285*1.5, 1, 0, 1])
construction_y8 = cp.TransferDistribution(tf.TriangTruncDet, [0.0443*0.5, 0.0443, 0.0443*1.5, 1, 0, 1])
construction_y9 = cp.TransferDistribution(tf.TriangTruncDet, [0.0680*0.5, 0.0680, 0.0680*1.5, 1, 0, 1])
construction_y10 = cp.TransferDistribution(tf.TriangTruncDet, [0.0777*0.5, 0.0777, 0.0777*1.5, 1, 0, 1])
construction_y11 = cp.TransferDistribution(tf.TriangTruncDet, [0.0931*0.5, 0.0931, 0.0931*1.5, 1, 0, 1])
construction_y12 = cp.TransferDistribution(tf.TriangTruncDet, [0.1145*0.5, 0.1145, 0.1145*1.5, 1, 0, 1])
construction_y13 = cp.TransferDistribution(tf.TriangTruncDet, [0.1341*0.5, 0.1341, 0.1341*1.5, 1, 0, 1])
construction_y14 = cp.TransferDistribution(tf.TriangTruncDet, [0.1432*0.5, 0.1432, 0.1432*1.5, 1, 0, 1])
construction_y15 = cp.TransferDistribution(tf.TriangTruncDet, [0.1808*0.5, 0.1808, 0.1808*1.5, 1, 0, 1])
construction_y16 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
# From 2030 to 2050, the TC is the same as the one from 2030 (y16)
construction_y17 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y18 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y19 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y20 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y21 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y22 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y23 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y24 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y25 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y26 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y27 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y28 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y29 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y30 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y31 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y32 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y33 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y34 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y35 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])
construction_y36 = cp.TransferDistribution(tf.TriangTruncDet, [0.1741*0.5, 0.1741, 0.1741*1.5, 1, 0, 1])



## Definition of consumption transfers based on the values above
consumption.transfers = [cp.TimeDependendDistributionTransfer([batteries_y0,batteries_y1,batteries_y2,batteries_y3,batteries_y4,batteries_y5,batteries_y6,batteries_y7,batteries_y8,batteries_y9,batteries_y10,batteries_y11,batteries_y12,batteries_y13,batteries_y14,batteries_y15,batteries_y16,batteries_y17,batteries_y18,batteries_y19,batteries_y20,batteries_y21,batteries_y22,batteries_y23,batteries_y24,batteries_y25,batteries_y26,batteries_y27,batteries_y28,batteries_y29,batteries_y30,
                                                             batteries_y31,batteries_y32,batteries_y33,batteries_y34,batteries_y35,batteries_y36],
                                                             batteries, consumption, priority=2),
                        cp.TimeDependendDistributionTransfer([tires_y0,tires_y1,tires_y2,tires_y3,tires_y4,tires_y5,tires_y6,tires_y7,tires_y8,tires_y9,tires_y10,tires_y11,tires_y12,tires_y13,tires_y14,tires_y15,tires_y16,tires_y17,tires_y18,tires_y19,tires_y20,tires_y21,tires_y22,tires_y23,tires_y24,tires_y25,tires_y26,tires_y27,tires_y28,tires_y29,tires_y30,tires_y31,tires_y32,tires_y33,tires_y34,tires_y35,tires_y36],
                                                             tires, consumption, priority=2),
                        cp.TimeDependendDistributionTransfer([construction_y0,construction_y1,construction_y2,construction_y3,construction_y4,construction_y5,construction_y6,construction_y7,construction_y8,construction_y9,construction_y10,construction_y11,construction_y12,construction_y13,construction_y14,construction_y15,construction_y16,construction_y17,construction_y18,construction_y19,construction_y20,construction_y21,construction_y22,construction_y23,construction_y24,construction_y25,construction_y26,construction_y27,construction_y28,construction_y29,construction_y30,construction_y31,construction_y32,construction_y33,construction_y34,construction_y35,construction_y36],
                                                             construction, consumption, priority=2),
                        ]



### Products categories to their USE and EoL compartments ###
# Batteries
batteries.transfers = [cp.ConstTransfer(1, batteries_EoL,priority=1)]


# Tires
tires.transfers = [cp.StochasticTransfer(nr.triangular, [0.105, 0.2112, 0.315], tires_Use, priority=2),
                   cp.ConstTransfer(1, tires_EoL, priority=1)]

# Construction
construction.transfers = [cp.StochasticTransfer(nr.triangular, [0.005, 0.01, 0.015], construction_Use, priority=2),
                          cp.ConstTransfer(1, construction_EoL, priority=1)]


### Release schedules for Use and EoL
# Batteries release schedules
batteries_Use.localRelease = cp.ListRelease([0.0])
batteries_EoL.localRelease = cp.ListRelease([0,0,0,0,0,0,0,0,0,0.003466974,0.01439745,0.04894278,0.1172529,
                                             0.1980285,0.2358228,0.1980285,0.1172529,0.04894278,0.01439745,
                                             0.003466974])

# Tires release schedules
tires_Use.localRelease = cp.ListRelease([0.25,0.25,0.25, 0.25])
tires_EoL.localRelease = cp.ListRelease([0.25,0.25,0.25, 0.25])

# Construction release schedules
construction_Use.localRelease = cp.ListRelease([0.9,0.01666667,0.01666667,0.01666667,0.01666667,
                                          0.01666667,0.01666667,0.01666667])

construction_EoL.localRelease = cp.ListRelease([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.001722281,
                                          0.001037804,0.001572363,0.00232936,0.003374171,0.004779078,
                                          0.006618624,0.00896268,0.01186737,0.01536449,0.01945038,
                                          0.02407602,0.0291399,0.03448561,0.03990567,0.04515209,
                                          0.04995371,0.05403864,0.0571594,0.05911772,0.05978529,
                                          0.05911772,0.0571594,0.05403864,0.04995371,0.04515209,
                                          0.03990567,0.03448561,0.0291399,0.02407602,0.01945038,
                                          0.01536449,0.01186737,0.00896268,0.006618624,0.004779078,
                                          0.003374171,0.00232936,0.001572363,0.001037804,0.001722281])



### Release from Use and EoL and allocation
# Batteries 
batteries_Use.transfers = [cp.ConstTransfer(1, NUsoil, priority=1)] # anyway nothing should go out to NUsoil as inflow to use is 0
batteries_EoL.transfers = [cp.ConstTransfer(1, battW, priority =1)]
# batteries_EoL.transfers = [cp.StochasticTransfer(nr.triangular, [0.5 * 0.4186, 0.4186, 1.5 * 0.4186], battW, priority=1),
#                            cp.StochasticTransfer(nr.triangular, [0.5 * 0.2450, 0.2450, 1.5 * 0.2450], Export, priority=2),
#                            cp.StochasticTransfer(nr.triangular, [0.5 * 0.4550, 0.4550, 1.5 * 0.4550], batt_unknown, priority=2),
#                            cp.StochasticTransfer(nr.triangular, [0.5 * 0.17, 0.17, 1.5 * 0.17], batt_improperDisposal, priority=2)]

# Tires
tires_Use.transfers   = [cp.StochasticTransfer(nr.triangular, [0.0781, 0.1562, 0.2343], Surfacewater, priority=1),
                         cp.StochasticTransfer(nr.triangular, [0.28135, 0.5627, 0.84405], Roadside, priority=1),
                         cp.StochasticTransfer(nr.triangular, [0.1333, 0.2666, 0.3999], Elimination, priority=1),
                         cp.StochasticTransfer(nr.triangular, [0.00725, 0.0145, 0.02175], NUsoil, priority=1)] 
tires_EoL.transfers = [cp.ConstTransfer(1, tiresW, priority=1)]


# Construction
construction_Use.transfers = [cp.ConstTransfer(1, Wastewater, priority=1)]
construction_EoL.transfers = [cp.ConstTransfer(1, CDW, priority =1)]

### Release schedules for reprocessing and reuse and corresponding allocation
#Assumption: those materials would be reused ONE year after they entered EoL and are recycled/reprocessed
# reuse_CDW_Mineral.localRelease = cp.ListRelease([0.25,0.25,0.25, 0.25])

# reuse_CDW_Mineral.transfers = [cp.ConstTransfer(1, Reuse, priority=1)]


###########################  COLLECTION RATE AND RECYCLING STREAMS  #################################
# Values are from CNT data from Rajkovic et al. (2020) in Hong et al. (2022)
# Battery waste
""" Done  """
## 2004 was Hong et al's model: y4, see below
# BW_S_y0 = cp.TransferDistribution (nr.triangular, [0.0485, 0.097, 0.1455])
# BW_S_y1 = cp.TransferDistribution (nr.triangular, [0.0599625, 0.117, 0.1740375])
# BW_S_y2 = cp.TransferDistribution (nr.triangular, [0.095625, 0.17, 0.244375])
# BW_S_y3 = cp.TransferDistribution (nr.triangular, [0.0917375, 0.179, 0.2662625])
## 2004
# BW_S_y4 = cp.TransferDistribution (nr.triangular, [0.094, 0.188, 0.282])
# BW_S_y5 = cp.TransferDistribution (nr.triangular, [0.0985, 0.197, 0.2955])
# BW_S_y6 = cp.TransferDistribution (nr.triangular, [0.105575, 0.206, 0.306425])
# BW_S_y7 = cp.TransferDistribution (nr.triangular, [0.1209375, 0.215, 0.3090625])
# BW_S_y8 = cp.TransferDistribution (nr.triangular, [0.138053833333333, 0.227, 0.315946166666667])
# BW_S_y9 = cp.TransferDistribution (nr.triangular, [0.167512333333333, 0.238, 0.308487666666667])
# BW_S_y10 = cp.TransferDistribution (nr.triangular, [0.212375, 0.25, 0.287625])
# BW_S_y11 = cp.TransferDistribution (nr.triangular, [0.280335, 0.33, 0.379665])
# BW_S_y12 = cp.TransferDistribution (nr.triangular, [0.30582, 0.36, 0.41418])
# BW_S_y13 = cp.TransferDistribution (nr.triangular, [0.32281, 0.38, 0.43719])
## 2014: 
# BW_S_y14 = cp.TransferDistribution (nr.triangular, [0.3398, 0.4, 0.4602])
# BW_S_y15 = cp.TransferDistribution (nr.triangular, [0.348295, 0.41, 0.471705])
# BW_S_y16 = cp.TransferDistribution (nr.triangular, [0.3410732, 0.428, 0.5149268])
# BW_S_y17 = cp.TransferDistribution (nr.triangular, [0.3319578, 0.446, 0.5600422])
# BW_S_y18 = cp.TransferDistribution (nr.triangular, [0.3441488, 0.464, 0.5838512])
# BW_S_y19 = cp.TransferDistribution (nr.triangular, [0.3803462, 0.482, 0.5836538])
# BW_S_y20 = cp.TransferDistribution (nr.triangular, [0.41825, 0.5, 0.58175])

# BW_S_y21 = cp.TransferConstant ( 0.52)
# BW_S_y22 = cp.TransferConstant ( 0.54)
# BW_S_y23 = cp.TransferConstant ( 0.56)
# BW_S_y24 = cp.TransferConstant ( 0.58)
# BW_S_y25 = cp.TransferConstant ( 0.6)
# BW_S_y26 = cp.TransferConstant ( 0.62)
# BW_S_y27 = cp.TransferConstant ( 0.64)
# BW_S_y28 = cp.TransferConstant ( 0.66)
# BW_S_y29 = cp.TransferConstant ( 0.68)
# BW_S_y30 = cp.TransferConstant ( 0.7)

BW_S_y0 = cp.TransferDistribution (nr.triangular, [0.3398, 0.4, 0.4602])
BW_S_y1 = cp.TransferDistribution (nr.triangular, [0.348295, 0.41, 0.471705])
BW_S_y2 = cp.TransferDistribution (nr.triangular, [0.3410732, 0.428, 0.5149268])
BW_S_y3 = cp.TransferDistribution (nr.triangular, [0.3319578, 0.446, 0.5600422])
BW_S_y4 = cp.TransferDistribution (nr.triangular, [0.3441488, 0.464, 0.5838512])
BW_S_y5 = cp.TransferDistribution (nr.triangular, [0.3803462, 0.482, 0.5836538])
BW_S_y6 = cp.TransferDistribution (nr.triangular, [0.41825, 0.5, 0.58175])
# #Change to triangular distributions
# BW_S_y7 = cp.TransferConstant ( 0.52)
# BW_S_y8 = cp.TransferConstant ( 0.54)
# BW_S_y9 = cp.TransferConstant ( 0.56)
# BW_S_y10 = cp.TransferConstant ( 0.58)
# BW_S_y11 = cp.TransferConstant ( 0.6)
# BW_S_y12 = cp.TransferConstant ( 0.62)
# BW_S_y13 = cp.TransferConstant ( 0.64)
# BW_S_y14 = cp.TransferConstant ( 0.66)
# BW_S_y15 = cp.TransferConstant ( 0.68)
# BW_S_y16 = cp.TransferConstant ( 0.7)

#Change to triangular distributions
#Example filtrationSeparation_y11 = cp.TransferDistribution(tf.TriangTruncDet, [0.000399136*0.5, 0.000399136, 0.000399136*1.5, 1, 0, 1])

## CHANGE: TC (mean) = 0.5, and adding uncertainty of +/-50% TC in triangular distributions
BW_S_y7 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
## CHANGE: constant after 2020. 
BW_S_y8 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y9 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y10 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
#2025
BW_S_y11 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y12 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y13 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y14 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y15 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y16 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
# Same value as the one from 2030 (y16) from now on
BW_S_y17 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y18 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y19 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y20 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y21 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y22 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y23 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y24 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y25 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y26 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y27 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y28 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y29 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y30 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y31 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y32 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y33 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y34 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y35 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])
BW_S_y36 = cp.TransferDistribution(nr.triangular, [0.25, 0.5, 0.75])



# CDW
""" Done """
# CDW_LF_y0 = cp.TransferDistribution (nr.triangular, [0.402781, 0.646, 0.889219])
# CDW_LF_y1 = cp.TransferDistribution (nr.triangular, [0.4155495, 0.617, 0.8184505])
# CDW_LF_y2 = cp.TransferDistribution (nr.triangular, [0.425418, 0.588, 0.750582])
# CDW_LF_y3 = cp.TransferDistribution (nr.triangular, [0.4323865, 0.559, 0.6856135])
# CDW_LF_y4 = cp.TransferDistribution (nr.triangular, [0.44469, 0.54, 0.63531])
# CDW_LF_y5 = cp.TransferDistribution (nr.triangular, [0.37091275, 0.507, 0.64308725])
# CDW_LF_y6 = cp.TransferDistribution (nr.triangular, [0.303202, 0.474, 0.644798])
# CDW_LF_y7 = cp.TransferDistribution (nr.triangular, [0.24155775, 0.441, 0.64044225])
# CDW_LF_y8 = cp.TransferDistribution (nr.triangular, [0.22678, 0.408, 0.58922])
# CDW_LF_y9 = cp.TransferDistribution (nr.triangular, [0.21146875, 0.375, 0.53853125])
# CDW_LF_y10 = cp.TransferDistribution (nr.triangular, [0.195624, 0.342, 0.488376])
# CDW_LF_y11 = cp.TransferDistribution (nr.triangular, [0.164952, 0.316, 0.467048])
# CDW_LF_y12 = cp.TransferDistribution (nr.triangular, [0.166452, 0.291, 0.415548])
# CDW_LF_y13 = cp.TransferDistribution (nr.triangular, [0.14355, 0.275, 0.40645])

CDW_LF_y0 = cp.TransferDistribution (nr.triangular, [0.148148, 0.259, 0.369852])
CDW_LF_y1 = cp.TransferDistribution (nr.triangular, [0.10962, 0.21, 0.31038])
CDW_LF_y2 = cp.TransferDistribution (nr.triangular, [0.0905, 0.181, 0.2715])
CDW_LF_y3 = cp.TransferDistribution (nr.triangular, [0.076, 0.152, 0.228])
CDW_LF_y4 = cp.TransferDistribution (nr.triangular, [0.0615, 0.123, 0.1845])
CDW_LF_y5 = cp.TransferDistribution (nr.triangular, [0.0465, 0.093, 0.1395])
CDW_LF_y6 = cp.TransferDistribution (nr.triangular, [0.032, 0.0639999999999999, 0.0959999999999999])

#Change to triangular distributions
CDW_LF_y7 = cp.TransferDistribution (nr.triangular, [0.061*0.5, 0.061, 0.061*1.5])
CDW_LF_y8 = cp.TransferDistribution (nr.triangular, [0.058*0.5, 0.058, 0.058*1.5])
CDW_LF_y9 = cp.TransferDistribution (nr.triangular, [0.054*0.5, 0.054, 0.054*1.5])
CDW_LF_y10 = cp.TransferDistribution (nr.triangular, [0.051*0.5, 0.051, 0.051*1.5])
CDW_LF_y11 = cp.TransferDistribution (nr.triangular, [0.047*0.5, 0.047, 0.047*1.5])
CDW_LF_y12 = cp.TransferDistribution (nr.triangular, [0.044*0.5, 0.044, 0.044*1.5])
CDW_LF_y13 = cp.TransferDistribution (nr.triangular, [0.04*0.5, 0.04, 0.04*1.5])
CDW_LF_y14 = cp.TransferDistribution (nr.triangular,[0.037*0.5, 0.037, 0.037*1.5])
CDW_LF_y15 = cp.TransferDistribution (nr.triangular, [0.033*0.5, 0.033, 0.033*1.5])
CDW_LF_y16 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
#Same as y16
CDW_LF_y17 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y18 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y19 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y20 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y21 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y22 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y23 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y24 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y25 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y26 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y27 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y28 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y29 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y30 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y31 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y32 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y33 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y34 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y35 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])
CDW_LF_y36 = cp.TransferDistribution (nr.triangular, [0.03*0.5, 0.03, 0.03*1.5])


# TiresW
#Only data regarding tires collection for 2016, based on ETRma (2018). Therefore, it is used as a constant transfer coefficient in the corresponding flows, below.



################## Transfers from waste streams

### NOTE: For hypothesis of release or transformation in the EoL and Recycling processes, it is recommended to have multiple expert opinions and concatenate them to create a probability distribultion of the Likelihoods and pathways:
  # Something similar as creating a probability distribution from different sources, in this case each person's opinion (likelihood) would be in a uniform distribution. Example from unburnt material in WIP below:
  # Unburnt = np.concatenate([ s*[0], [0], nr.uniform(1e-9, 0.02, int(s*0.454)), nr.uniform(1e-9, 0.25, int(s*0.454)), nr.triangular(0.25, 0.25, 0.3, int(s*0.092))])

  # Burning.transfers =[cp.RandomChoiceTransfer(Unburnt, Remains_burning, priority=2),
  #                     cp.ConstTransfer(1,Elimination, priority=1)]

  # Remains_burning.transfers =[cp.StochasticTransfer(nr.triangular,[0.10,0.19,0.29], Flyash, priority=2),
  #                             cp.ConstTransfer(1,Bottomash, priority=1)]

  # Flyash.transfers =[cp.StochasticTransfer(nr.triangular, [0.000025, 0.00005, 0.000075], Air, priority=2), # 1- 0.99995 (filter efficiency)
  #                    cp.StochasticTransfer(nr.triangular, [1e-9, 0.0001, 0.00015], Wastewater, priority=2),
  #                    cp.ConstTransfer(1,Filterash,priority=1)]



## Likelihood of releases (or transformation), using the Decision Tree (considering the questions related to the material and its application, see below) and Recycling Streams diagrams from the corresponding post-consumer waste category (MACRAME D4.2)
    # Logarithmic (factor of TEN):
        # Low = 0 - 0.01
        # Medium = 0.01 - 0.1
        # High =  0.1 - 0.99
        # Very High = 0.99 - 1

#TO DO: the values could be modified to avoid overlapping of a value between categories. For instance, define Medium as >=0.01, 0.1 and High as >=0.1, 1. 
#Moreover, the likelihood approach could be implemented as a component and only call the corresponding category based on the hypothesis.

    # Questions: 
        # 1. Is the AdMa coating or embedded into the product?
        # 2. What temperatures can the AdMa withstand?
        # 3. Which pH can the AdMa withstand?
        # 4. Is there anything else that can eliminate (transformed) the AdMa?

## CDW flows
CDW.transfers = [cp.TimeDependendDistributionTransfer([CDW_LF_y0,CDW_LF_y1,CDW_LF_y2,CDW_LF_y3,CDW_LF_y4,CDW_LF_y5,CDW_LF_y6,CDW_LF_y7,CDW_LF_y8,CDW_LF_y9,CDW_LF_y10,
                                                       CDW_LF_y11,CDW_LF_y12,CDW_LF_y13,CDW_LF_y14,CDW_LF_y15,CDW_LF_y16,CDW_LF_y17,CDW_LF_y18,CDW_LF_y19,CDW_LF_y20,
                                                       CDW_LF_y21,CDW_LF_y22,CDW_LF_y23,CDW_LF_y24,CDW_LF_y25,CDW_LF_y26,CDW_LF_y27,CDW_LF_y28,CDW_LF_y29,CDW_LF_y30,
                                                       CDW_LF_y31,CDW_LF_y32,CDW_LF_y33,CDW_LF_y34,CDW_LF_y35,CDW_LF_y36],
                                                       Landfill, CDW, priority=2),
                 cp.ConstTransfer(1, sorting_CDW_Mineral, priority=1)
                 ]

# Since previously CDW was divided into different waste types and here we focus on mineral waste (concrete), we did not use same as the previous models. But we used a constant transfer coefficient, based on the fact that the management of mineral waste has not been mantained similar over the time boundary (2014-2020) based on EC:JRC (2023)
sorting_CDW_Mineral.transfers =  [cp.ConstTransfer(1, recycling_CDW_Mineral, priority=1), 
                                  cp.StochasticTransfer(nr.triangular,[0.0679, 0.1160, 0.1641], Landfill, priority=2),
                                  cp.StochasticTransfer(nr.triangular,[0.0289, 0.0495, 0.07004], CDW_Mineral_Backf, priority=2),
                                  cp.StochasticTransfer(nr.triangular,[0.0039, 0.0067, 0.0095], CDW_Mineral_CementKilns, priority=2)
                                 ]
 
CDW_Mineral_CementKilns.transfers =[cp.ConstTransfer(1, cementFill, priority=1)]  

# CHANGED to priority 1 (+/-50% uncertainty)
recycling_CDW_Mineral.transfers = [cp.StochasticTransfer(nr.triangular, [0.0816, 0.1631, 0.2447], Air, priority=1),
                         cp.StochasticTransfer(nr.triangular, [0.1537, 0.3074, 0.4611], reprocessing_CDW_Mineral, priority=1),
                         cp.StochasticTransfer(nr.triangular, [0.2648, 0.5295, 0.7943], Landfill, priority=1), 
                        ]

## CHANGE: Backfilling is not included as a Compartment, but the flow between CDW_Mineral_Backf to 
CDW_Mineral_Backf.transfers = [cp.ConstTransfer(1, Backfilling, priority=1)]

Backfilling.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]

reprocessing_CDW_Mineral.transfers = [cp.ConstTransfer(1, Aggregates, priority=1)]
## NOTE (from Qie Wu, 2019, MSc thesis): Aggregates are defined as the final products of mineral recycling plants, serving as the recycled aggregates for ready-mixed concrete plant.

## CHANGE: ADDING this flow
Aggregates.transfers = [cp.ConstTransfer(1, cementFill, priority=1)]


## TiresW flows (NEW)
tiresW.transfers = [cp.StochasticTransfer(nr.triangular, [0.0554, 0.0947, 0.1340], Reuse, priority=2), #Reuse + Retreading
                    cp.StochasticTransfer(nr.triangular, [0.0405, 0.0693, 0.0981], Export, priority=2),
                    cp.ConstTransfer(1, tiresW_preRecycling, priority=1) ## ELT Arising (ETRma). CHANGE: to ConstTransfer
                   ]

## CHANGE: before as Improper disposal, here as Unknown (based on ETRma (2018), here it represents the sum of Landfill/Unknown + Stocks)
tiresW_preRecycling.transfers = [cp.StochasticTransfer(nr.triangular,[0.3430, 0.5864, 0.8298], tiresW_MatRecovery, priority=1),
                                 cp.StochasticTransfer(nr.triangular,[0.2006, 0.3429,0.4852], tiresW_EnergyRec, priority=2),
                                 cp.StochasticTransfer(nr.triangular,[0.0414, 0.0708, 0.1002], Unknown, priority=2), #Before as ImproperDisposal
                                ]

# From Figure of ETRma's ELTs Report (2016), but not much explanation from where the values come from. Therefore, adding extra conversion factor (+/-30% unc.).
tiresW_EnergyRec.transfers = [cp.ConstTransfer(1, cementKilns, priority=1),
                              cp.StochasticTransfer(nr.triangular,[0.0506, 0.09, 0.1294], energyRecovery, priority=2),
                             ]


# Hypothesis: (very) high likelihood of release (elimination/transformation). In this case depending on the temperatures, but assuming working temperatures conditions of around 1450°C, 
#graphene might be transformed = eliminated. Assuming that in the Energy Industry (incineration) compartment, the same hypothesis applies. But more variation in temperature; thus, allocated to WIP (simplification of the model).
   
cementKilns_dist = np.concatenate([nr.uniform(0.99, 1, s)])
cementKilns.transfers = [cp.RandomChoiceTransfer(cementKilns_dist, Elimination, priority=2),
                         cp.StochasticTransfer(tf.TriangTruncDet,[0.375, 0.75, 1.125, 1, 0, 1], energyRecovery, priority=1), # CHANGE to StochasticTransfer, previously Const
                         cp.StochasticTransfer(nr.triangular,[0.125, 0.25, 0.375], tiresW_MatRecovery, priority=1),
                        ]

## CHANGE: assumption = energyRecovery would have similar path as WIP (simplification of the model).
energyRecovery.transfers = [cp.ConstTransfer(1, WIP, priority=1)]


tiresW_MatRecovery.transfers = [cp.StochasticTransfer(nr.triangular, [0.0383, 0.0654, 0.0925], CivilEng, priority=2),
                                cp.ConstTransfer(1, tiresW_Recycling, priority=1),
                               ]

## CHANGE: allocated to EoW_alternatives
CivilEng.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
                                                 
tiresW_Recycling.transfers = [cp.ConstTransfer(1, tiresW_MechanicalRecycling_like, priority=1),
                              cp.StochasticTransfer(nr.triangular,[0.0988, 0.1756, 0.2524], cementFill, priority=2), #Process not further investigated
                              cp.StochasticTransfer(nr.triangular,[0.0037, 0.0065, 0.0093], Pyrolysis, priority=2),
                              cp.StochasticTransfer(nr.triangular,[0.0122, 0.0217, 0.0312],  DockFenders, priority=2) #Process not further investigated
                             ]

# Hypothesis: high likelihood of release (elimination). In this case, the majority is expected to remain in the Carbon Black content (releases to Air are neglected but further research is needed)                                            

Pyrolisis_dist = np.concatenate([nr.uniform(0.1, 0.99, s)])
Pyrolysis.transfers = [cp.RandomChoiceTransfer(Pyrolisis_dist, Elimination, priority=2),
                       cp.ConstTransfer(1, CarbonBlack, priority=1)                                 
                      ]                     


## CHANGE: allocation to EoW_alternatives
cementFill.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
DockFenders.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
CarbonBlack.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]


# Hypothesis: (very) low likelihood of release (elimination) In this case, it might be emitted to air (elimination is neglected but further research is needed) and remain in rubber                                             
tiresW_MechanicalRecycling_dist = np.concatenate([nr.uniform(0, 0.01, s)])
tiresW_MechanicalRecycling_like.transfers = [cp.RandomChoiceTransfer(tiresW_MechanicalRecycling_dist, Air, priority=2),
                                             cp.ConstTransfer(1, MechanicalRecycling_sorting, priority=1)                                      
                                             ]
                                        
MechanicalRecycling_sorting.transfers = [cp.StochasticTransfer(nr.triangular,[0.3094, 0.55, 0.7906], SoftFallSurfaces, priority=1),
                                        cp.StochasticTransfer(nr.triangular,[0.1238, 0.22 ,0.3163], MouldedProducts, priority=1),
                                        cp.StochasticTransfer(nr.triangular,[0.0956,0.17,0.2444], Unknown, priority=1), #Undetermined (Export, Trade) in EuRIC (2022) Fact sheet. https://euric.org/resource-hub/reports-studies/mechanical-tyre-recycling-factsheet
                                        cp.StochasticTransfer(nr.triangular,[0.0113, 0.02, 0.0288], Pavings, priority=1),
                                        cp.StochasticTransfer(nr.triangular,[0.0225, 0.04, 0.0575], OtherUses, priority=1)
                                        ]


## CHANGE: allocation to EoW_alternatives
SoftFallSurfaces.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
MouldedProducts.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
Pavings.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]
OtherUses.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]


# tires_EoL.transfers   = [cp.StochasticTransfer(nr.triangular, [0.580808081*0.5, 0.580808081, 0.580808081*1.5], tires_preRecycling, priority=2),
#                          cp.StochasticTransfer(nr.triangular, [0.01540404*0.5, 0.01540404, 0.01540404*1.5], tires_improperDisposal, priority=2),
#                          cp.StochasticTransfer(nr.triangular, [0.373737374*0.5, 0.373737374, 0.373737374*1.5], tires_Export, priority=2),
#                          cp.StochasticTransfer(nr.triangular, [1e-9*0.5, 1e-9, 1e-9*1.5], tires_unknown, priority=1)]



## BattW flows
battW.transfers = [cp.TimeDependendDistributionTransfer([BW_S_y0,BW_S_y1,BW_S_y2,BW_S_y3,BW_S_y4,BW_S_y5,BW_S_y6,BW_S_y7,BW_S_y8,BW_S_y9,BW_S_y10,
                                                          BW_S_y11,BW_S_y12,BW_S_y13,BW_S_y14,BW_S_y15,BW_S_y16,BW_S_y17,BW_S_y18,BW_S_y19,BW_S_y20,
                                                          BW_S_y21,BW_S_y22,BW_S_y23,BW_S_y24,BW_S_y25,BW_S_y26,BW_S_y27,BW_S_y28,BW_S_y29,BW_S_y30,
                                                          BW_S_y31,BW_S_y32,BW_S_y33,BW_S_y34,BW_S_y35, BW_S_y36
                                                         ], battW_Dismantling, battW, priority=2),
                    cp.ConstTransfer(1, battW_Unidentified, priority=1)
                   ]

# It could be changed to the surplus that is not possible to treat in existing plants as [REF] mentioned
# CHANGE: Export is not included. According to Leon et al. (2020) portable Li-ion batteries:
# Using the relationship: 42% collection-to-recycling rate, 17% disposal rate (allocated to Landfill in this model), and 41% gap (allocated to Unknown in this model). Thus, Unidentified = 0.29 Landfill + 0.71 Unknown
# CV of 34.5% to each value (scientific article + weighting) due to lack of data, no possibility to add the dynamics; assuming it has remained constant during the temporal boundary
battW_Unidentified.transfers = [cp.StochasticTransfer(nr.triangular,[0.1900, 0.29, 0.3901], Landfill, priority=1), #CV=0.10005
                                cp.ConstTransfer(1, Unknown, priority=1) 
                               ]

battW_Dismantling.transfers = [#batteries_Reuse, batteries_Recycling, priority=2),
                               cp.ConstTransfer(1, battW_SizeReduction_like, priority=1)
                              ]


# Hypothesis = low likelihood to eliminate/release material. In this case, releases to air. However, further research is needed
battW_SizeReduction_dist = np.concatenate([nr.uniform(0, 0.01, s)])

battW_SizeReduction_like.transfers = [cp.RandomChoiceTransfer(battW_SizeReduction_dist, Air, priority=2),
                                      cp.ConstTransfer(1, battW_Recycling, priority=1)
                                     ]

#batteries_Reuse.transfers = [cp.ConstTransfer(1, consumption_batteries)]
#batteries_Reuse.localRelease = cp.ListRelease([1.0], 1)

# According to Bruno & Fiore (2023). battW_Pyro stands for Pyrometallurgical recycling process for batteries. Assumption: general MFA apply to graphite (in this case, similarly to GRMs)

battW_Recycling.transfers = [cp.StochasticTransfer(nr.triangular, [0.105, 0.21, 0.315], battW_HydroLeaching_like, priority=2),
                             cp.StochasticTransfer(tf.TriangTruncDet, [0.395, 0.79, 1.185, 1, 0, 1], battW_Pyro_like, priority=1)
                            ]


# Hypothesis = high likelihood. In this case, eliminated
battW_Pyro_dist = np.concatenate([nr.uniform(0.1, 0.99, s)])
battW_Pyro_like.transfers = [cp.RandomChoiceTransfer(battW_Pyro_dist, Elimination, priority=2),
                             cp.ConstTransfer(1, battW_Pyro, priority=1)
                            ]

battW_Pyro.transfers = [cp.StochasticTransfer(nr.triangular,[0.145,0.29,0.435], battW_HydroLeaching, priority=2),
                        cp.StochasticTransfer(nr.triangular,[0.035,0.07,0.105], battW_Pyro_losses, priority= 2),
                        cp.ConstTransfer(1, battW_Pyro_output, priority= 1)                 
                       ]



# Hypothesis = medium, to elimination (or leached out) but further research is needed.
battW_HydroLeaching_dist = np.concatenate([nr.uniform(0.01, 0.1, s)])
battW_HydroLeaching_like.transfers = [cp.RandomChoiceTransfer(battW_HydroLeaching_dist, Elimination, priority=2),
                                      cp.ConstTransfer(1, battW_HydroLeaching, priority=1)
                                     ]

battW_HydroLeaching.transfers = [cp.StochasticTransfer(nr.triangular, [0.045, 0.09, 0.135], battW_HydroLeaching_losses, priority=2),
                                 cp.ConstTransfer(1, battW_HydroLeaching_output, priority=1)
                                ]

#Assuming that graphene flows within the residue and losses material (losses), and to recovered material (output)
battW_Pyro_losses.transfers = [cp.ConstTransfer(1, battW_Recycling_Residue, priority=1)]

battW_HydroLeaching_losses.transfers = [cp.ConstTransfer(1, battW_Recycling_Residue, priority=1)]

## CHANGE: battW_Recycling_Residue allocated to Unknown
battW_Recycling_Residue.transfers = [cp.ConstTransfer(1, Unknown, priority=1)]


battW_HydroLeaching_output.transfers = [cp.ConstTransfer(1, battW_RecMat_impurities, priority=1)]

battW_Pyro_output.transfers = [cp.ConstTransfer(1, battW_RecMat_impurities, priority=1)]

## CHANGE: battW_RecMat_impurities NOT SINK compartment
battW_RecMat_impurities.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]


#battW_DirectRecycling.transfers = [cp.ConstTransfer(1, battW_Recycling_Residue, priority=2), # (the flow to battW_Recycling_Residue was deleted because we want to compare the maximal mass recovered [assumption] to the rest)
                                #cp.StochasticTransfer(tf.TrapezTrunc(0.775, 0.9, 0.5, 0.5, s, 0, 1), NewBatt, priority=1)]
# The compartment could be named reprocessing_battW_RecMat OR directly NewBatt. In addition, CV needs to be modified according to the references. Now, it is considered from expert opinion ( = 0.5)      



############################ Wastewater, WWTP, STP systems ############################
### Dynamic values from Hong et al. (2022), but they built the model up based on the CNT Model for EU28+CH+NO (author:Sun.T.Y, Rajkovic S.). 
# Now taking y0 as 2014, in Sun et al.'s study case is y14. They started from the year 2000.
#################################   WWTP   ####################################
## Waste water
""" DONE """
## UNCHANGED
#WW_NS_y0 = cp.TransferDistribution (nr.triangular, [0.0815, 0.163, 0.2445])
#WW_NS_y1 = cp.TransferDistribution (nr.triangular, [0.075, 0.15, 0.225])
#WW_NS_y2 = cp.TransferDistribution (nr.triangular, [0.0775, 0.155, 0.2325])
#WW_NS_y3 = cp.TransferDistribution (nr.triangular, [0.0765, 0.153, 0.2295])
#WW_NS_y4 = cp.TransferDistribution (nr.triangular, [0.076, 0.152, 0.228])
#WW_NS_y5 = cp.TransferDistribution (nr.triangular, [0.0785, 0.157, 0.2355])
#WW_NS_y6 = cp.TransferDistribution (nr.triangular, [0.067, 0.134, 0.201])
#WW_NS_y7 = cp.TransferDistribution (nr.triangular, [0.0715, 0.143, 0.2145])
#WW_NS_y8 = cp.TransferDistribution (nr.triangular, [0.071, 0.142, 0.213])
#WW_NS_y9 = cp.TransferDistribution (nr.triangular, [0.072, 0.144, 0.216])
#WW_NS_y10 = cp.TransferDistribution (nr.triangular, [0.0715, 0.143, 0.2145])
#WW_NS_y11 = cp.TransferDistribution (nr.triangular, [0.074, 0.148, 0.222])
#WW_NS_y12 = cp.TransferDistribution (nr.triangular, [0.0725, 0.145, 0.2175])
#WW_NS_y13 = cp.TransferDistribution (nr.triangular, [0.072, 0.144, 0.216])
#WW_NS_y14 = cp.TransferDistribution(nr.triangular, [0.071, 0.142, 0.213])
WW_NS_y0 = cp.TransferDistribution(nr.triangular, [0.071, 0.142, 0.213]) #2014
WW_NS_y1 = cp.TransferDistribution(nr.triangular, [0.0715, 0.143, 0.2145])
WW_NS_y2 = cp.TransferDistribution(nr.triangular, [0.0715, 0.143, 0.2145])
WW_NS_y3 = cp.TransferDistribution(nr.triangular, [0.0715, 0.143, 0.2145])
WW_NS_y4 = cp.TransferDistribution(nr.triangular, [0.0715, 0.143, 0.2145])
WW_NS_y5 = cp.TransferDistribution(nr.triangular, [0.0715, 0.143, 0.2145])
WW_NS_y6 = cp.TransferDistribution(nr.triangular, [0.072, 0.144, 0.216])

WW_NS_y7 = cp.TransferConstant(0.128)
WW_NS_y8 = cp.TransferConstant(0.121)
WW_NS_y9 = cp.TransferConstant(0.113)
WW_NS_y10 = cp.TransferConstant(0.106)
WW_NS_y11 = cp.TransferConstant(0.098)
WW_NS_y12 = cp.TransferConstant(0.09)
WW_NS_y13 = cp.TransferConstant(0.083)
WW_NS_y14 = cp.TransferConstant(0.075)
WW_NS_y15 = cp.TransferConstant(0.068)
WW_NS_y16 = cp.TransferConstant(0.06) #2030
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
WW_NS_y17 = cp.TransferConstant(0.06)
WW_NS_y18 = cp.TransferConstant(0.06)
WW_NS_y19 = cp.TransferConstant(0.06)
WW_NS_y20 = cp.TransferConstant(0.06)
WW_NS_y21 = cp.TransferConstant(0.06)
WW_NS_y22 = cp.TransferConstant(0.06)
WW_NS_y23 = cp.TransferConstant(0.06)
WW_NS_y24 = cp.TransferConstant(0.06)
WW_NS_y25 = cp.TransferConstant(0.06)
WW_NS_y26 = cp.TransferConstant(0.06)
WW_NS_y27 = cp.TransferConstant(0.06)
WW_NS_y28 = cp.TransferConstant(0.06)
WW_NS_y29 = cp.TransferConstant(0.06)
WW_NS_y30 = cp.TransferConstant(0.06)
WW_NS_y31 = cp.TransferConstant(0.06)
WW_NS_y32 = cp.TransferConstant(0.06)
WW_NS_y33 = cp.TransferConstant(0.06)
WW_NS_y34 = cp.TransferConstant(0.06)
WW_NS_y35 = cp.TransferConstant(0.06)
WW_NS_y36 = cp.TransferConstant(0.06)



## Sewage system
# SS_NT_y0 = cp.TransferDistribution(nr.triangular, [0.0325, 0.065, 0.0975])
# SS_NT_y1 = cp.TransferDistribution(nr.triangular, [0.034, 0.068, 0.102])
# SS_NT_y2 = cp.TransferDistribution(nr.triangular, [0.028, 0.056, 0.084])
# SS_NT_y3 = cp.TransferDistribution(nr.triangular, [0.0245, 0.049, 0.0735])
# SS_NT_y4 = cp.TransferDistribution(nr.triangular, [0.0225, 0.045, 0.0675])
# SS_NT_y5 = cp.TransferDistribution(nr.triangular, [0.018, 0.036, 0.054])
# SS_NT_y6 = cp.TransferDistribution(nr.triangular, [0.0285, 0.057, 0.0855])
# SS_NT_y7 = cp.TransferDistribution(nr.triangular, [0.024, 0.048, 0.072])
# SS_NT_y8 = cp.TransferDistribution(nr.triangular, [0.0235, 0.047, 0.0705])
# SS_NT_y9 = cp.TransferDistribution(nr.triangular, [0.0205, 0.041, 0.0615])
# SS_NT_y10 = cp.TransferDistribution(nr.triangular, [0.0185, 0.037, 0.0555])
# SS_NT_y11 = cp.TransferDistribution(nr.triangular, [0.016, 0.032, 0.048])
# SS_NT_y12 = cp.TransferDistribution(nr.triangular, [0.014, 0.028, 0.042])
# SS_NT_y13 = cp.TransferDistribution(nr.triangular, [0.0125, 0.025, 0.0375])
# SS_NT_y14 = cp.TransferDistribution(nr.triangular, [0.0105, 0.021, 0.0315])
SS_NT_y0 = cp.TransferDistribution(nr.triangular, [0.0105, 0.021, 0.0315])
SS_NT_y1 = cp.TransferDistribution(nr.triangular, [0.0085, 0.017, 0.0255])
SS_NT_y2 = cp.TransferDistribution(nr.triangular, [0.007, 0.014, 0.021])
SS_NT_y3 = cp.TransferDistribution(nr.triangular, [0.005, 0.01, 0.015])
SS_NT_y4 = cp.TransferDistribution(nr.triangular, [0.0035, 0.007, 0.0105])
SS_NT_y5 = cp.TransferDistribution(nr.triangular, [0.00152318181818182, 0.003, 0.00447681818181818])
SS_NT_y6 = cp.TransferDistribution(nr.triangular, [0, 0, 0.0001])  # changed by hand to be different from 0 by HH

SS_NT_y7 = cp.TransferConstant(0)
SS_NT_y8 = cp.TransferConstant(0)
SS_NT_y9 = cp.TransferConstant(0)
SS_NT_y10 = cp.TransferConstant(0)
SS_NT_y11 = cp.TransferConstant(0)
SS_NT_y12 = cp.TransferConstant(0)
SS_NT_y13 = cp.TransferConstant(0)
SS_NT_y14 = cp.TransferConstant(0)
SS_NT_y15 = cp.TransferConstant(0)
SS_NT_y16 = cp.TransferConstant(0)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
SS_NT_y17 = cp.TransferConstant(0)
SS_NT_y18 = cp.TransferConstant(0)
SS_NT_y19 = cp.TransferConstant(0)
SS_NT_y20 = cp.TransferConstant(0)
SS_NT_y21 = cp.TransferConstant(0)
SS_NT_y22 = cp.TransferConstant(0)
SS_NT_y23 = cp.TransferConstant(0)
SS_NT_y24 = cp.TransferConstant(0)
SS_NT_y25 = cp.TransferConstant(0)
SS_NT_y26 = cp.TransferConstant(0)
SS_NT_y27 = cp.TransferConstant(0)
SS_NT_y28 = cp.TransferConstant(0)
SS_NT_y29 = cp.TransferConstant(0)
SS_NT_y30 = cp.TransferConstant(0)
SS_NT_y31 = cp.TransferConstant(0)
SS_NT_y32 = cp.TransferConstant(0)
SS_NT_y33 = cp.TransferConstant(0)
SS_NT_y34 = cp.TransferConstant(0)
SS_NT_y35 = cp.TransferConstant(0)
SS_NT_y36 = cp.TransferConstant(0)

### WWTP
""" Done """
# WWTP_TI_y0 = cp.TransferDistribution(nr.triangular, [0.0295, 0.059, 0.0885])
# WWTP_TI_y1 = cp.TransferDistribution(nr.triangular, [0.0255, 0.051, 0.0765])
# WWTP_TI_y2 = cp.TransferDistribution(nr.triangular, [0.021, 0.042, 0.063])
# WWTP_TI_y3 = cp.TransferDistribution(nr.triangular, [0.0175, 0.035, 0.0525])
# WWTP_TI_y4 = cp.TransferDistribution(nr.triangular, [0.014, 0.028, 0.042])
# WWTP_TI_y5 = cp.TransferDistribution(nr.triangular, [0.011, 0.022, 0.033])
# WWTP_TI_y6 = cp.TransferDistribution(nr.triangular, [0.012, 0.024, 0.036])
# WWTP_TI_y7 = cp.TransferDistribution(nr.triangular, [0.0105, 0.021, 0.0315])
# WWTP_TI_y8 = cp.TransferDistribution(nr.triangular, [0.011, 0.022, 0.033])
# WWTP_TI_y9 = cp.TransferDistribution(nr.triangular, [0.008, 0.016, 0.024])
# WWTP_TI_y10 = cp.TransferDistribution(nr.triangular, [0.0075, 0.015, 0.0225])
# WWTP_TI_y11 = cp.TransferDistribution(nr.triangular, [0.0065, 0.013, 0.0195])
# WWTP_TI_y12 = cp.TransferDistribution(nr.triangular, [0.0055, 0.011, 0.0165])
# WWTP_TI_y13 = cp.TransferDistribution(nr.triangular, [0.005, 0.01, 0.015])
# WWTP_TI_y14 = cp.TransferDistribution(nr.triangular, [0.004, 0.008, 0.012])
WWTP_TI_y0 = cp.TransferDistribution(nr.triangular, [0.004, 0.008, 0.012])
WWTP_TI_y1 = cp.TransferDistribution(nr.triangular, [0.0035, 0.007, 0.0105])
WWTP_TI_y2 = cp.TransferDistribution(nr.triangular, [0.0025, 0.005, 0.0075])
WWTP_TI_y3 = cp.TransferDistribution(nr.triangular, [0.002, 0.004, 0.006])
WWTP_TI_y4 = cp.TransferDistribution(nr.triangular, [0.0015, 0.003, 0.0045])
WWTP_TI_y5 = cp.TransferDistribution(nr.triangular, [0.000507727272727272, 0.001, 0.00149227272727273])
WWTP_TI_y6 = cp.TransferDistribution(nr.triangular, [0, 0, 0.0001])  # adapted by hand to be not equal 0

WWTP_TI_y7 = cp.TransferConstant(0)
WWTP_TI_y8 = cp.TransferConstant(0)
WWTP_TI_y9 = cp.TransferConstant(0)
WWTP_TI_y10 = cp.TransferConstant(0)
WWTP_TI_y11 = cp.TransferConstant(0)
WWTP_TI_y12 = cp.TransferConstant(0)
WWTP_TI_y13 = cp.TransferConstant(0)
WWTP_TI_y14 = cp.TransferConstant(0)
WWTP_TI_y15 = cp.TransferConstant(0)
WWTP_TI_y16 = cp.TransferConstant(0)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
WWTP_TI_y17 = cp.TransferConstant(0)
WWTP_TI_y18 = cp.TransferConstant(0)
WWTP_TI_y19 = cp.TransferConstant(0)
WWTP_TI_y20 = cp.TransferConstant(0)
WWTP_TI_y21 = cp.TransferConstant(0)
WWTP_TI_y22 = cp.TransferConstant(0)
WWTP_TI_y23 = cp.TransferConstant(0)
WWTP_TI_y24 = cp.TransferConstant(0)
WWTP_TI_y25 = cp.TransferConstant(0)
WWTP_TI_y26 = cp.TransferConstant(0)
WWTP_TI_y27 = cp.TransferConstant(0)
WWTP_TI_y28 = cp.TransferConstant(0)
WWTP_TI_y29 = cp.TransferConstant(0)
WWTP_TI_y30 = cp.TransferConstant(0)
WWTP_TI_y31 = cp.TransferConstant(0)
WWTP_TI_y32 = cp.TransferConstant(0)
WWTP_TI_y33 = cp.TransferConstant(0)
WWTP_TI_y34 = cp.TransferConstant(0)
WWTP_TI_y35 = cp.TransferConstant(0)
WWTP_TI_y36 = cp.TransferConstant(0)

# WWTP_TII_y0 = cp.TransferDistribution(nr.triangular, [0.1725, 0.345, 0.5175])
# WWTP_TII_y1 = cp.TransferDistribution(nr.triangular, [0.1615, 0.323, 0.4845])
# WWTP_TII_y2 = cp.TransferDistribution(nr.triangular, [0.161, 0.322, 0.483])
# WWTP_TII_y3 = cp.TransferDistribution(nr.triangular, [0.1535, 0.307, 0.4605])
# WWTP_TII_y4 = cp.TransferDistribution(nr.triangular, [0.145, 0.29, 0.435])
# WWTP_TII_y5 = cp.TransferDistribution(nr.triangular, [0.123, 0.246, 0.369])
# WWTP_TII_y6 = cp.TransferDistribution(nr.triangular, [0.1305, 0.261, 0.3915])
# WWTP_TII_y7 = cp.TransferDistribution(nr.triangular, [0.124, 0.248, 0.372])
# WWTP_TII_y8 = cp.TransferDistribution(nr.triangular, [0.118, 0.236, 0.354])
# WWTP_TII_y9 = cp.TransferDistribution(nr.triangular, [0.1125, 0.225, 0.3375])
# WWTP_TII_y10 = cp.TransferDistribution(nr.triangular, [0.108, 0.216, 0.324])
# WWTP_TII_y11 = cp.TransferDistribution(nr.triangular, [0.105, 0.21, 0.315])
# WWTP_TII_y12 = cp.TransferDistribution(nr.triangular, [0.1025, 0.205, 0.3075])
# WWTP_TII_y13 = cp.TransferDistribution(nr.triangular, [0.0975, 0.195, 0.2925])
# WWTP_TII_y14 = cp.TransferDistribution(nr.triangular, [0.0925, 0.185, 0.2775])
WWTP_TII_y0 = cp.TransferDistribution(nr.triangular, [0.0925, 0.185, 0.2775])
WWTP_TII_y1 = cp.TransferDistribution(nr.triangular, [0.0875, 0.175, 0.2625])
WWTP_TII_y2 = cp.TransferDistribution(nr.triangular, [0.083, 0.166, 0.249])
WWTP_TII_y3 = cp.TransferDistribution(nr.triangular, [0.0785, 0.157, 0.2355])
WWTP_TII_y4 = cp.TransferDistribution(nr.triangular, [0.074, 0.148, 0.222])
WWTP_TII_y5 = cp.TransferDistribution(nr.triangular, [0.0705740909090909, 0.139, 0.207425909090909])
WWTP_TII_y6 = cp.TransferDistribution(nr.triangular, [0.060931, 0.109, 0.157069])

WWTP_TII_y7 = cp.TransferConstant(0.117)
WWTP_TII_y8 = cp.TransferConstant(0.113)
WWTP_TII_y9 = cp.TransferConstant(0.11)
WWTP_TII_y10 = cp.TransferConstant(0.107)
WWTP_TII_y11 = cp.TransferConstant(0.103)
WWTP_TII_y12 = cp.TransferConstant(0.1)
WWTP_TII_y13 = cp.TransferConstant(0.096)
WWTP_TII_y14 = cp.TransferConstant(0.093)
WWTP_TII_y15 = cp.TransferConstant(0.09)
WWTP_TII_y16 = cp.TransferConstant(0.086)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
WWTP_TII_y17 = cp.TransferConstant(0.086)
WWTP_TII_y18 = cp.TransferConstant(0.086)
WWTP_TII_y19 = cp.TransferConstant(0.086)
WWTP_TII_y20 = cp.TransferConstant(0.086)
WWTP_TII_y21 = cp.TransferConstant(0.086)
WWTP_TII_y22 = cp.TransferConstant(0.086)
WWTP_TII_y23 = cp.TransferConstant(0.086)
WWTP_TII_y24 = cp.TransferConstant(0.086)
WWTP_TII_y25 = cp.TransferConstant(0.086)
WWTP_TII_y26 = cp.TransferConstant(0.086)
WWTP_TII_y27 = cp.TransferConstant(0.086)
WWTP_TII_y28 = cp.TransferConstant(0.086)
WWTP_TII_y29 = cp.TransferConstant(0.086)
WWTP_TII_y30 = cp.TransferConstant(0.086)
WWTP_TII_y31 = cp.TransferConstant(0.086)
WWTP_TII_y32 = cp.TransferConstant(0.086)
WWTP_TII_y33 = cp.TransferConstant(0.086)
WWTP_TII_y34 = cp.TransferConstant(0.086)
WWTP_TII_y35 = cp.TransferConstant(0.086)
WWTP_TII_y36 = cp.TransferConstant(0.086)

### STPsludge
""" Done """
# SG_WIIP_y0 = cp.TransferDistribution(nr.triangular, [0.188727, 0.209, 0.229273])
# SG_WIIP_y1 = cp.TransferDistribution(nr.triangular, [0.1822428, 0.213, 0.2437572])
# SG_WIIP_y2 = cp.TransferDistribution(nr.triangular, [0.1753794, 0.217, 0.2586206])
# SG_WIIP_y3 = cp.TransferDistribution(nr.triangular, [0.1791868, 0.221, 0.2628132])
# SG_WIIP_y4 = cp.TransferDistribution(nr.triangular, [0.194265, 0.225, 0.255735])
# SG_WIIP_y5 = cp.TransferDistribution(nr.triangular, [0.209764, 0.229, 0.248236])
# SG_WIIP_y6 = cp.TransferDistribution(nr.triangular, [0.21696, 0.256, 0.29504])
# SG_WIIP_y7 = cp.TransferDistribution(nr.triangular, [0.222794, 0.286, 0.349206])
# SG_WIIP_y8 = cp.TransferDistribution(nr.triangular, [0.2410785, 0.317, 0.3929215])
# SG_WIIP_y9 = cp.TransferDistribution(nr.triangular, [0.252648, 0.319, 0.385352])
# SG_WIIP_y10 = cp.TransferDistribution(nr.triangular, [0.2643435, 0.321, 0.3776565])
# SG_WIIP_y11 = cp.TransferDistribution(nr.triangular, [0.22309875, 0.315, 0.40690125])
# SG_WIIP_y12 = cp.TransferDistribution(nr.triangular, [0.186795, 0.315, 0.443205])
# SG_WIIP_y13 = cp.TransferDistribution(nr.triangular, [0.177947, 0.308, 0.438053])
# SG_WIIP_y14 = cp.TransferDistribution(nr.triangular, [0.1693125, 0.301, 0.4326875])
SG_WIIP_y0 = cp.TransferDistribution(nr.triangular, [0.1693125, 0.301, 0.4326875])
SG_WIIP_y1 = cp.TransferDistribution(nr.triangular, [0.1614375, 0.315, 0.4685625])
SG_WIIP_y2 = cp.TransferDistribution(nr.triangular, [0.1645, 0.329, 0.4935])
SG_WIIP_y3 = cp.TransferDistribution(nr.triangular, [0.1715, 0.343, 0.5145])
SG_WIIP_y4 = cp.TransferDistribution(nr.triangular, [0.1785, 0.357, 0.5355])
SG_WIIP_y5 = cp.TransferDistribution(nr.triangular, [0.1901375, 0.371, 0.5518625])
SG_WIIP_y6 = cp.TransferDistribution(nr.triangular, [0.217125, 0.386, 0.554875])

SG_WIIP_y7 = cp.TransferConstant(0.376)
SG_WIIP_y8 = cp.TransferConstant(0.359)
SG_WIIP_y9 = cp.TransferConstant(0.343)
SG_WIIP_y10 = cp.TransferConstant(0.326)
SG_WIIP_y11 = cp.TransferConstant(0.309)
SG_WIIP_y12 = cp.TransferConstant(0.292)
SG_WIIP_y13 = cp.TransferConstant(0.275)
SG_WIIP_y14 = cp.TransferConstant(0.259)
SG_WIIP_y15 = cp.TransferConstant(0.242)
SG_WIIP_y16 = cp.TransferConstant(0.225)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
SG_WIIP_y17 = cp.TransferConstant(0.225)
SG_WIIP_y18 = cp.TransferConstant(0.225)
SG_WIIP_y19 = cp.TransferConstant(0.225)
SG_WIIP_y20 = cp.TransferConstant(0.225)
SG_WIIP_y21 = cp.TransferConstant(0.225)
SG_WIIP_y22 = cp.TransferConstant(0.225)
SG_WIIP_y23 = cp.TransferConstant(0.225)
SG_WIIP_y24 = cp.TransferConstant(0.225)
SG_WIIP_y25 = cp.TransferConstant(0.225)
SG_WIIP_y26 = cp.TransferConstant(0.225)
SG_WIIP_y27 = cp.TransferConstant(0.225)
SG_WIIP_y28 = cp.TransferConstant(0.225)
SG_WIIP_y29 = cp.TransferConstant(0.225)
SG_WIIP_y30 =  cp.TransferConstant(0.225)
SG_WIIP_y31 =  cp.TransferConstant(0.225)
SG_WIIP_y32 =  cp.TransferConstant(0.225)
SG_WIIP_y33 =  cp.TransferConstant(0.225)
SG_WIIP_y34 =  cp.TransferConstant(0.225)
SG_WIIP_y35 =  cp.TransferConstant(0.225)
SG_WIIP_y36 = cp.TransferConstant(0.225)


# SG_LF_y0 = cp.TransferDistribution(nr.triangular, [0.188727, 0.209, 0.229273])
# SG_LF_y1 = cp.TransferDistribution(nr.triangular, [0.1779648, 0.208, 0.2380352])
# SG_LF_y2 = cp.TransferDistribution(nr.triangular, [0.1681056, 0.208, 0.2478944])
# SG_LF_y3 = cp.TransferDistribution(nr.triangular, [0.1678356, 0.207, 0.2461644])
# SG_LF_y4 = cp.TransferDistribution(nr.triangular, [0.1778604, 0.206, 0.2341396])
# SG_LF_y5 = cp.TransferDistribution(nr.triangular, [0.18778, 0.205, 0.22222])
# SG_LF_y6 = cp.TransferDistribution(nr.triangular, [0.17289, 0.204, 0.23511])
# SG_LF_y7 = cp.TransferDistribution(nr.triangular, [0.158137, 0.203, 0.247863])
# SG_LF_y8 = cp.TransferDistribution(nr.triangular, [0.1543815, 0.203, 0.2516185])
# SG_LF_y9 = cp.TransferDistribution(nr.triangular, [0.145728, 0.184, 0.222272])
# SG_LF_y10 = cp.TransferDistribution(nr.triangular, [0.1375245, 0.167, 0.1964755])
# SG_LF_y11 = cp.TransferDistribution(nr.triangular, [0.10694575, 0.151, 0.19505425])
# SG_LF_y12 = cp.TransferDistribution(nr.triangular, [0.081834, 0.138, 0.194166])
# SG_LF_y13 = cp.TransferDistribution(nr.triangular, [0.06990775, 0.121, 0.17209225])
# SG_LF_y14 = cp.TransferDistribution(nr.triangular, [0.0590625, 0.105, 0.1509375])
SG_LF_y0 = cp.TransferDistribution(nr.triangular, [0.0590625, 0.105, 0.1509375])
SG_LF_y1 = cp.TransferDistribution(nr.triangular, [0.056712, 0.102, 0.147288])
SG_LF_y2 = cp.TransferDistribution(nr.triangular, [0.053851, 0.098, 0.142149])
SG_LF_y3 = cp.TransferDistribution(nr.triangular, [0.051585, 0.095, 0.138415])
SG_LF_y4 = cp.TransferDistribution(nr.triangular, [0.0579215, 0.091, 0.1240785])
SG_LF_y5 = cp.TransferDistribution(nr.triangular, [0.06424, 0.088, 0.11176])
SG_LF_y6 = cp.TransferDistribution(nr.triangular, [0.069174, 0.084, 0.098826])

SG_LF_y7 = cp.TransferConstant(0.078)
SG_LF_y8 = cp.TransferConstant(0.072)
SG_LF_y9 = cp.TransferConstant(0.066)
SG_LF_y10 = cp.TransferConstant(0.06)
SG_LF_y11 = cp.TransferConstant(0.054)
SG_LF_y12 = cp.TransferConstant(0.048)
SG_LF_y13 = cp.TransferConstant(0.043)
SG_LF_y14 = cp.TransferConstant(0.037)
SG_LF_y15 = cp.TransferConstant(0.031)
SG_LF_y16 = cp.TransferConstant(0.025)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
SG_LF_y17 = cp.TransferConstant(0.025)
SG_LF_y18 = cp.TransferConstant(0.025)
SG_LF_y19 = cp.TransferConstant(0.025)
SG_LF_y20 = cp.TransferConstant(0.025)
SG_LF_y21 = cp.TransferConstant(0.025)
SG_LF_y22 = cp.TransferConstant(0.025)
SG_LF_y23 = cp.TransferConstant(0.025)
SG_LF_y24 = cp.TransferConstant(0.025)
SG_LF_y25 = cp.TransferConstant(0.025)
SG_LF_y26 = cp.TransferConstant(0.025)
SG_LF_y27 = cp.TransferConstant(0.025)
SG_LF_y28 = cp.TransferConstant(0.025)
SG_LF_y29 = cp.TransferConstant(0.025)
SG_LF_y30 = cp.TransferConstant(0.025)
SG_LF_y31 = cp.TransferConstant(0.025)
SG_LF_y32 = cp.TransferConstant(0.025)
SG_LF_y33 = cp.TransferConstant(0.025)
SG_LF_y34 = cp.TransferConstant(0.025)
SG_LF_y35 = cp.TransferConstant(0.025)
SG_LF_y36 = cp.TransferConstant(0.025)


### Flyash
# FA_LF_y0 = cp.TransferConstant(1)
# FA_LF_y1 = cp.TransferConstant(1)
# FA_LF_y2 = cp.TransferConstant(1)
# FA_LF_y3 = cp.TransferConstant(1)
# FA_LF_y4 = cp.TransferConstant(1)
# FA_LF_y5 = cp.TransferConstant(1)
# FA_LF_y6 = cp.TransferConstant(1)
# FA_LF_y7 = cp.TransferConstant(1)
# FA_LF_y8 = cp.TransferConstant(1)
# FA_LF_y9 = cp.TransferConstant(1)
# FA_LF_y10 = cp.TransferConstant(1)
# FA_LF_y11 = cp.TransferConstant(1)
# FA_LF_y12 = cp.TransferConstant(1)
# FA_LF_y13 = cp.TransferConstant(1)
# FA_LF_y14 = cp.TransferConstant(1)
FA_LF_y0 = cp.TransferConstant(1) #2014
FA_LF_y1 = cp.TransferConstant(1)
FA_LF_y2 = cp.TransferConstant(1)
FA_LF_y3 = cp.TransferConstant(1)
FA_LF_y4 = cp.TransferConstant(1)
FA_LF_y5 = cp.TransferConstant(1)
FA_LF_y6 = cp.TransferConstant(1)
FA_LF_y7 = cp.TransferConstant(0.96)
FA_LF_y8 = cp.TransferConstant(0.92)
FA_LF_y9 = cp.TransferConstant(0.88)
FA_LF_y10 = cp.TransferConstant(0.84)
FA_LF_y11 = cp.TransferConstant(0.8)
FA_LF_y12 = cp.TransferConstant(0.76)
FA_LF_y13 = cp.TransferConstant(0.72)
FA_LF_y14 = cp.TransferConstant(0.68)
FA_LF_y15 = cp.TransferConstant(0.64)
FA_LF_y16 = cp.TransferConstant(0.6)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
FA_LF_y17 = cp.TransferConstant(0.6)
FA_LF_y18 = cp.TransferConstant(0.6)
FA_LF_y19 = cp.TransferConstant(0.6)
FA_LF_y20 = cp.TransferConstant(0.6)
FA_LF_y21 = cp.TransferConstant(0.6)
FA_LF_y22 = cp.TransferConstant(0.6)
FA_LF_y23 = cp.TransferConstant(0.6)
FA_LF_y24 = cp.TransferConstant(0.6)
FA_LF_y25 = cp.TransferConstant(0.6)
FA_LF_y26 = cp.TransferConstant(0.6)
FA_LF_y27 = cp.TransferConstant(0.6)
FA_LF_y28 = cp.TransferConstant(0.6)
FA_LF_y29 = cp.TransferConstant(0.6)
FA_LF_y30 = cp.TransferConstant(0.6)
FA_LF_y31 = cp.TransferConstant(0.6)
FA_LF_y32 = cp.TransferConstant(0.6)
FA_LF_y33 = cp.TransferConstant(0.6)
FA_LF_y34 = cp.TransferConstant(0.6)
FA_LF_y35 = cp.TransferConstant(0.6)
FA_LF_y36 = cp.TransferConstant(0.6)



### Bottomash
# BA_R_y0 = cp.TransferDistribution(nr.triangular, [0.2445, 0.489, 0.7335])
# BA_R_y1 = cp.TransferDistribution(nr.triangular, [0.243, 0.486, 0.729])
# BA_R_y2 = cp.TransferDistribution(nr.triangular, [0.247025, 0.482, 0.716975])
# BA_R_y3 = cp.TransferDistribution(nr.triangular, [0.253125, 0.45, 0.646875])
# BA_R_y4 = cp.TransferDistribution(nr.triangular, [0.253125, 0.45, 0.646875])
# BA_R_y5 = cp.TransferDistribution(nr.triangular, [0.236478333333333, 0.46, 0.683521666666667])
# BA_R_y6 = cp.TransferDistribution(nr.triangular, [0.235, 0.47, 0.705])
# BA_R_y7 = cp.TransferDistribution(nr.triangular, [0.24, 0.48, 0.72])
# BA_R_y8 = cp.TransferDistribution(nr.triangular, [0.245, 0.49, 0.735])
# BA_R_y9 = cp.TransferDistribution(nr.triangular, [0.260208333333333, 0.5, 0.739791666666667])
# BA_R_y10 = cp.TransferDistribution(nr.triangular, [0.29172, 0.51, 0.72828])
# BA_R_y11 = cp.TransferDistribution(nr.triangular, [0.24012, 0.46, 0.67988])
# BA_R_y12 = cp.TransferDistribution(nr.triangular, [0.23452, 0.41, 0.58548])
# BA_R_y13 = cp.TransferDistribution(nr.triangular, [0.215586, 0.413, 0.610414])
# BA_R_y14 = cp.TransferDistribution(nr.triangular, [0.237952, 0.416, 0.594048])
BA_R_y0 = cp.TransferDistribution(nr.triangular, [0.237952, 0.416, 0.594048])
BA_R_y1 = cp.TransferDistribution(nr.triangular, [0.227592, 0.436, 0.644408])
BA_R_y2 = cp.TransferDistribution(nr.triangular, [0.2165, 0.433, 0.6495])
BA_R_y3 = cp.TransferDistribution(nr.triangular, [0.2145, 0.429, 0.6435])
BA_R_y4 = cp.TransferDistribution(nr.triangular, [0.2125, 0.425, 0.6375])
BA_R_y5 = cp.TransferDistribution(nr.triangular, [0.211, 0.422, 0.633])
BA_R_y6 = cp.TransferDistribution(nr.triangular, [0.209, 0.418, 0.627])

BA_R_y7 = cp.TransferConstant(0.535)
BA_R_y8 = cp.TransferConstant(0.565)
BA_R_y9 = cp.TransferConstant(0.594)
BA_R_y10 = cp.TransferConstant(0.623)
BA_R_y11 = cp.TransferConstant(0.653)
BA_R_y12 = cp.TransferConstant(0.682)
BA_R_y13 = cp.TransferConstant(0.712)
BA_R_y14 = cp.TransferConstant(0.741)
BA_R_y15 = cp.TransferConstant(0.771)
BA_R_y16 = cp.TransferConstant(0.8)
# From 2030 forward, the dynamic TC remain constant using the last value (y16 in this code) until 2050 (y36)
BA_R_y17 = cp.TransferConstant(0.8)
BA_R_y18 = cp.TransferConstant(0.8)
BA_R_y19 = cp.TransferConstant(0.8)
BA_R_y20 = cp.TransferConstant(0.8)
BA_R_y21 = cp.TransferConstant(0.8)
BA_R_y22 = cp.TransferConstant(0.8)
BA_R_y23 = cp.TransferConstant(0.8)
BA_R_y24 = cp.TransferConstant(0.8)
BA_R_y25 = cp.TransferConstant(0.8)
BA_R_y26 = cp.TransferConstant(0.8)
BA_R_y27 = cp.TransferConstant(0.8)
BA_R_y28 = cp.TransferConstant(0.8)
BA_R_y29 = cp.TransferConstant(0.8)
BA_R_y30 = cp.TransferConstant(0.8)
BA_R_y31 = cp.TransferConstant(0.8)
BA_R_y32 = cp.TransferConstant(0.8)
BA_R_y33 = cp.TransferConstant(0.8)
BA_R_y34 = cp.TransferConstant(0.8)
BA_R_y35 = cp.TransferConstant(0.8)
BA_R_y36 = cp.TransferConstant(0.8)


Wastewater.transfers = [
    cp.TimeDependendDistributionTransfer([WW_NS_y0, WW_NS_y1, WW_NS_y2, WW_NS_y3, WW_NS_y4, WW_NS_y5, WW_NS_y6, WW_NS_y7, WW_NS_y8, WW_NS_y9, WW_NS_y10,
                                          WW_NS_y11, WW_NS_y12, WW_NS_y13, WW_NS_y14, WW_NS_y15, WW_NS_y16, 
                                          WW_NS_y17, WW_NS_y18, WW_NS_y19, WW_NS_y20, WW_NS_y21, WW_NS_y22, WW_NS_y23, WW_NS_y24, WW_NS_y25, WW_NS_y26, WW_NS_y27,
                                          WW_NS_y28, WW_NS_y29, WW_NS_y30, WW_NS_y31, WW_NS_y32, WW_NS_y33, WW_NS_y34, WW_NS_y35, WW_NS_y36],
                                         NoSewageSystem, Wastewater, priority=2),
    cp.ConstTransfer(1, SewageSystem, priority=1)]

NoSewageSystem.transfers = [cp.ConstTransfer(0.03, Surfacewater, priority=1),
                            cp.ConstTransfer(0.11, onsiteTreat, priority=1)]

onsiteTreat.transfers = [cp.ConstTransfer(1, TreatIOnsite, priority=1)]  # STPoverflow  => TreatIOnsite

SewageSystem.transfers = [cp.StochasticTransfer(nr.triangular, [0.016, 0.032, 0.048], STPoverflow, priority=2),
                          # overflow rate
                          cp.StochasticTransfer(nr.triangular, [0.0325, 0.065, 0.0975], Subsurface, priority=2),
                          # Exfiltration
                          cp.TimeDependendDistributionTransfer(
                              [SS_NT_y0, SS_NT_y1, SS_NT_y2, SS_NT_y3, SS_NT_y4, SS_NT_y5, SS_NT_y6, SS_NT_y7, SS_NT_y8, SS_NT_y9, SS_NT_y10,
                               SS_NT_y11, SS_NT_y12, SS_NT_y13, SS_NT_y14, SS_NT_y15, SS_NT_y16, 
                               SS_NT_y17, SS_NT_y18, SS_NT_y19, SS_NT_y20, SS_NT_y21, SS_NT_y22, SS_NT_y23, SS_NT_y24, SS_NT_y25, SS_NT_y26, SS_NT_y27, SS_NT_y28,
                               SS_NT_y29, SS_NT_y30, SS_NT_y31, SS_NT_y32, SS_NT_y33, SS_NT_y34, SS_NT_y35, SS_NT_y36],
                              Surfacewater, SewageSystem, priority=2),
                          # discharges of sewer to nature (when no wwtp is installed) #value from Sana's study year 19
                          cp.ConstTransfer(1, WWTP, priority=1)]

WWTP.transfers = [cp.TimeDependendDistributionTransfer([
                                                        WWTP_TI_y0, WWTP_TI_y1, WWTP_TI_y2, WWTP_TI_y3, WWTP_TI_y4, WWTP_TI_y5, WWTP_TI_y6, WWTP_TI_y7, WWTP_TI_y8, WWTP_TI_y9, WWTP_TI_y10,
                                                        WWTP_TI_y11, WWTP_TI_y12, WWTP_TI_y13, WWTP_TI_y14, WWTP_TI_y15, WWTP_TI_y16, 
                                                        WWTP_TI_y17, WWTP_TI_y18, WWTP_TI_y19, WWTP_TI_y20, WWTP_TI_y21, WWTP_TI_y22, WWTP_TI_y23, WWTP_TI_y24, WWTP_TI_y25, WWTP_TI_y26, WWTP_TI_y27, WWTP_TI_y28,
                                                        WWTP_TI_y29, WWTP_TI_y30, WWTP_TI_y31, WWTP_TI_y32, WWTP_TI_y33, WWTP_TI_y34, WWTP_TI_y35, WWTP_TI_y36
                                                        ], TreatI, WWTP, priority=2),
                  cp.TimeDependendDistributionTransfer([
                                                        WWTP_TII_y0, WWTP_TII_y1, WWTP_TII_y2, WWTP_TII_y3, WWTP_TII_y4, WWTP_TII_y5, WWTP_TII_y6, WWTP_TII_y7, WWTP_TII_y8, WWTP_TII_y9, WWTP_TII_y10,
                                                        WWTP_TII_y11, WWTP_TII_y12, WWTP_TII_y13, WWTP_TII_y14, WWTP_TII_y15, WWTP_TII_y16, 
                                                        WWTP_TII_y17, WWTP_TII_y18, WWTP_TII_y19, WWTP_TII_y20, WWTP_TII_y21, WWTP_TII_y22, WWTP_TII_y23, WWTP_TII_y24, WWTP_TII_y25, WWTP_TII_y26, WWTP_TII_y27,
                                                        WWTP_TII_y28, WWTP_TII_y29, WWTP_TII_y30, WWTP_TII_y31, WWTP_TII_y32, WWTP_TII_y33, WWTP_TII_y34, WWTP_TII_y35, WWTP_TII_y36
                                                        ], TreatII, WWTP, priority=2),
                  cp.ConstTransfer(1, TreatIII, priority=1)]

# NonRemovalData - represents the inversed removal data means what actually stays in the water
# only one source for removal data was found and thus applied to all treatment levels
NonRemovalData_AllTreat = np.concatenate(
    [[0], nr.uniform(0, 0.1, int(s * 0.91)), nr.triangular(0.1, 0.1, 0.12, int(s * 0.09))])

TreatI.transfers = [cp.RandomChoiceTransfer(NonRemovalData_AllTreat, STPeffluent, priority=2),
                    cp.ConstTransfer(1, STPsludge, priority=1)]

TreatIOnsite.transfers = [cp.RandomChoiceTransfer(NonRemovalData_AllTreat, Subsurface, priority=2),
                          cp.ConstTransfer(1, OnsiteSludge, priority=1)]


## CHANGE: allocating OnsiteSludge to STPsludge flows (assuming this follows the same waste management streams)
OnsiteSludge.transfers = [cp.ConstTransfer(1, STPsludge, priority=1)]

TreatII.transfers = [cp.RandomChoiceTransfer(NonRemovalData_AllTreat, STPeffluent, priority=2),
                     cp.ConstTransfer(1, STPsludge, priority=1)]

TreatIII.transfers = [cp.RandomChoiceTransfer(NonRemovalData_AllTreat, STPeffluent, priority=2),
                      cp.ConstTransfer(1, STPsludge, priority=1)]
# STP
STPsludge.transfers = [cp.TimeDependendDistributionTransfer(
                                                            [SG_WIIP_y0, SG_WIIP_y1, SG_WIIP_y2, SG_WIIP_y3, SG_WIIP_y4, SG_WIIP_y5, SG_WIIP_y6, SG_WIIP_y7, SG_WIIP_y8, SG_WIIP_y9, SG_WIIP_y10,
                                                             SG_WIIP_y11, SG_WIIP_y12, SG_WIIP_y13, SG_WIIP_y14, SG_WIIP_y15, SG_WIIP_y16, 
                                                             SG_WIIP_y17, SG_WIIP_y18, SG_WIIP_y19, SG_WIIP_y20, SG_WIIP_y21, SG_WIIP_y22, SG_WIIP_y23, SG_WIIP_y24, SG_WIIP_y25, SG_WIIP_y26, SG_WIIP_y27, SG_WIIP_y28,
                                                             SG_WIIP_y29, SG_WIIP_y30, SG_WIIP_y31, SG_WIIP_y32, SG_WIIP_y33, SG_WIIP_y34, SG_WIIP_y35, SG_WIIP_y36],
                                                            WIP, STPsludge, priority=2),
                       cp.TimeDependendDistributionTransfer(
                                                           [SG_LF_y0, SG_LF_y1, SG_LF_y2, SG_LF_y3, SG_LF_y4, SG_LF_y5, SG_LF_y6, SG_LF_y7, SG_LF_y8, SG_LF_y9, SG_LF_y10,
                                                            SG_LF_y11, SG_LF_y12, SG_LF_y13, SG_LF_y14, SG_LF_y15, SG_LF_y16, 
                                                            SG_LF_y17, SG_LF_y18, SG_LF_y19, SG_LF_y20, SG_LF_y21, SG_LF_y22, SG_LF_y23, SG_LF_y24, SG_LF_y25, SG_LF_y26, SG_LF_y27, SG_LF_y28,
                                                            SG_LF_y29, SG_LF_y30, SG_LF_y31, SG_LF_y32, SG_LF_y33, SG_LF_y34, SG_LF_y35, SG_LF_y36],
                                                            Landfill, STPsludge, priority=2),
                       cp.ConstTransfer(1, STsoil, priority=1)]

STPoverflow.transfers = [cp.ConstTransfer(1, Surfacewater, priority=1)]

STPeffluent.transfers = [cp.ConstTransfer(1, Surfacewater, priority=1)]

# WIP (Elimination while burning)
WIP.transfers = [cp.ConstTransfer(1, Burning, priority=1)]

# Elimination while burning
# special case where müller  gave min, mean and max - thus two uniform distributions around that 'mean'
Unburnt = np.concatenate([s * [0], [0], nr.uniform(1e-9, 0.02, int(s * 0.454)), nr.uniform(1e-9, 0.25, int(s * 0.454)),
                          nr.triangular(0.25, 0.25, 0.3, int(s * 0.092))])

Burning.transfers = [cp.RandomChoiceTransfer(Unburnt, Remains_burning, priority=2),
                     cp.ConstTransfer(1, Elimination, priority=1)]

Remains_burning.transfers = [cp.StochasticTransfer(nr.triangular, [0.10, 0.19, 0.29], Flyash, priority=2),
                             cp.ConstTransfer(1, Bottomash, priority=1)]

Flyash.transfers = [cp.StochasticTransfer(nr.triangular, [0.000025, 0.00005, 0.000075], Air, priority=2),
                    # 1- 0.99995 (filter efficiency)
                    cp.StochasticTransfer(nr.triangular, [1e-9, 0.0001, 0.00015], Wastewater, priority=2),
                    cp.ConstTransfer(1, Filterash, priority=1)]

Filterash.transfers = [
    cp.TimeDependendDistributionTransfer([FA_LF_y0, FA_LF_y1, FA_LF_y2, FA_LF_y3, FA_LF_y4, FA_LF_y5, FA_LF_y6, FA_LF_y7, FA_LF_y8, FA_LF_y9, FA_LF_y10,
                                          FA_LF_y11, FA_LF_y12, FA_LF_y13, FA_LF_y14, FA_LF_y15, FA_LF_y16, 
                                          FA_LF_y17, FA_LF_y18, FA_LF_y19, FA_LF_y20, FA_LF_y21, FA_LF_y22, FA_LF_y23, FA_LF_y24, FA_LF_y25, FA_LF_y26, FA_LF_y27,
                                          FA_LF_y28, FA_LF_y29, FA_LF_y30, FA_LF_y31, FA_LF_y32, FA_LF_y33, FA_LF_y34, FA_LF_y35, FA_LF_y36],
                                         Landfill, Filterash, priority=2),
    cp.ConstTransfer(1, reprocessing, priority=1)]

Bottomash.transfers = [
    cp.TimeDependendDistributionTransfer([BA_R_y0, BA_R_y1, BA_R_y2, BA_R_y3, BA_R_y4, BA_R_y5, BA_R_y6, BA_R_y7, BA_R_y8, BA_R_y9, BA_R_y10,
                                          BA_R_y11, BA_R_y12, BA_R_y13, BA_R_y14, BA_R_y15, BA_R_y16, 
                                          BA_R_y17, BA_R_y18, BA_R_y19, BA_R_y20, BA_R_y21, BA_R_y22, BA_R_y23, BA_R_y24, BA_R_y25, BA_R_y26, BA_R_y27,
                                          BA_R_y28, BA_R_y29, BA_R_y30, BA_R_y31, BA_R_y32, BA_R_y33, BA_R_y34, BA_R_y35, BA_R_y36],
                                         reprocessing, Bottomash, priority=2),
    cp.ConstTransfer(1, Landfill, priority=1)]

## CHANGE: reprocessing (ashes) to EoW_alternatives (in this case, after reprocessing, ashes could be used for different sectors, e.g., construction)
reprocessing.transfers = [cp.ConstTransfer(1, EoW_alternatives, priority=1)]

Air.transfers = [cp.StochasticTransfer(nr.triangular, [0.047 * 0.5, 0.047, 0.047 * 1.5], Surfacewater, priority=2),
                 cp.ConstTransfer(1, NUsoil, priority=1)]

model.checkModelValidity()
