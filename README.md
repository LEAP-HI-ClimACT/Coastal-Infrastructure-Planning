# Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change
This repository includes the codes used in the paper: **Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change**.

**TLDR**: Climate change adaptation is often addressed with static policies that ignore future climate uncertainties. This work develops a climate change adaptation framework using stochastic-control methodologies, offering optimal solutions that address model uncertainty, social costs of carbon, and a wide variety of adaptation actions.

This repository contains the steps along with the necessary codes for generating components of the MDP and POMDP models for flood management and for generating the input
for dynamic programming solvers like FRTDP.  The workflow for generating the MDP input for the coastal city setting application with the two floodwall options, as described in the paper, using the provided code is outlined below.

# Steps to generate adaptation policies for coastal applications
## Datasets used for sea-level-rise simulations
- The IPCC projection dataset (for Battery tide gauge in New York) used to simulate SLR models is publicly available at: [IPCC slr projections](https://zenodo.org/records/6382554)
- The dataset for SLR trends (for Battery tide gauge in New York) used to fit the noise present in the past SLR observations is taken from: [NOAA slr trends](https://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?id=8518750)
## Sea-level-rise simulations
The SLR trajectories are simulated using the **slr_simulations.m** file. This file requires the SLR projections to be stored in the format as is provided in **SLR_245.mat** file. This file provides sea-level rise (SLR) projections across various years (**years.mat**) and quantile values (**quantiles.mat**) for the SSP245 scenario at the New York Battery tidal gauge (downloaded from the links above). The generated simulations are used to inform SLR state transitions.

## Datasets used for storm surge models
The dataset of annual extremes (for Battery tide gauge in New York) used to fit a GEV model for storm surges is taken from: [NOAA Extremes](https://tidesandcurrents.noaa.gov/est/est_station.shtml?stnid=8518750)

## Generate transition models
To generate discrete-state transition probabilities, the SLR and storm surge levels are suitably discretized into 77 and 72 discrete states, respectively. The discrete-state transitions for Sea Level Rise (SLR) and storm surge are generated using the codes, as:
1. **slr transitions.m** for SLR transitions
2. **surge_transitions.m** for storm surge transitions

The SLR transition probability matrices derived and used in this paper can be found at **t_slr_245.mat** file. The storm surge transitions derived from the GEV model is uploaded here as **t_surge.mat** file.

The full state space of the problem needs to account for non-stationary sea level rise trends and different system configurations possible with different sets of actions. The complete state space of the problem is derived in the paper. 
After generating the SLR and surge transitions, the full transition model over the complete MDP state space is generated considering state-space augmentation with time and systems, using **transition_model_MDP.m** file.

## Generate reward models
Generate the rewards for states in the different systems for all the actions in the application setting, along with social cost of carbon considerations, using the **rewards_with_scc_systems.m** file.
After generating the rewards for the different systems, the rewards for the full MDP state space are generated using the **rewards_model_MDP.m** file. This file calls **rewards_fp.m** file to shape the rewards in such a way so as to use the fast parser for the input file in the solver used.

## Generate the MDP input file for dynamic programming solvers
After generating the transitions and rewards components of the MDP model in the above steps, the MDP input file is generated using **inp_mdp_file_generator.m** file. The generated input is solved using the FRTDP solver. The solver is publicly available at: [FRTDP](https://github.com/trey0/zmdp). The details related to the solver and the input file formats can be found in
this [GitHub](https://github.com/trey0/zmdp).

## Citation

If you use this repository in your research or projects, please cite our paper as follows:


