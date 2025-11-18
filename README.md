# Extended_DPMFA_scenarios_graphene
These simulated scenarios are developed as part of the following publication:
Title: 

Authors:

Journal:

These prospective scenarios generate probability distributions associated with the flows of graphene within the lifecycle of three product categories: batteries, tires, and construction, and potential releases to the environment, at the European level (EU27, Norway, the United Kingdom, and Switzerland) from 2014 to 2050 with yearly intervals.

The tool is written using the DPMFA Python package (github.com/empa-tsl/dpmfa), consisting of:
* Components.py: the different classes of compartments and flows are defined;
* Model.py: the required methods to link the compartments together with the flows are defined;
* Simulator.py: includes the methods required to run the simulation and extract results;
* [NAME OF EACH MODEL AND RUN]

Each model.py takes as input a csv file containing the probability distributions of the produced amounts of graphene in each year in the corresponding product categories: "[NAME].csv"

For further details, the user is referred to the above publication.
