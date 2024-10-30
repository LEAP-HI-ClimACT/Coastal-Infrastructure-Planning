# Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change
This repository includes the codes used in the paper: **Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change**.

**TLDR**:![Coastal Infrastructure Adaptation under Climate Change](https://github.com/ashmita695/Coastal_Infrastructure_Planning/blob/main/coastal%20application%20framework.pdf)

Coastal infrastructure adaptation against climate change effects is often addressed with static policies that ignore future climate uncertainties. This work develops a climate change adaptation framework using stochastic-control methodologies, offering optimal solutions that address model uncertainty, social costs of carbon, and a wide variety of adaptation actions.

This repository contains the steps along with the necessary codes for generating components of the MDP and POMDP models for flood management and for generating the input
for dynamic programming solvers like FRTDP.  The workflow for generating the MDP input for the coastal city setting application with the two floodwall options for the SSP:245 scenario, as described in the paper, using the provided code, is outlined below.

# Steps to generate adaptation policies for coastal applications
## Datasets used for sea-level-rise simulations
- The IPCC projection dataset (for Battery tide gauge in New York) used to simulate SLR models is available at: [IPCC slr projections](https://zenodo.org/records/6382554). These datasets, if used directly, are to be cited properly as documented in [IPCC SLR projections](https://zenodo.org/records/6382554).
- The dataset for SLR trends (for Battery tide gauge in New York) used to fit the noise present in the past SLR observations is taken from: [NOAA slr trends](https://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?id=8518750)
## Sea-level-rise simulations
The SLR trajectories are simulated using the **slr_simulations.m** file. This file requires the SLR projections to be stored in the format as is provided in **SLR_245.mat** file. This file provides sea-level rise (SLR) projections across various years (**years.mat**) and quantile values (**quantiles.mat**) for the SSP245 scenario at the New York Battery tidal gauge (downloaded from [IPCC slr projections](https://zenodo.org/records/6382554)). The SLR projections for other SSP scenarios can also be found at [IPCC slr projections](https://zenodo.org/records/6382554). The generated simulations are used to inform SLR state transitions.

## Datasets used for storm surge models
The dataset of annual extremes (for Battery tide gauge in New York) used to fit a GEV model for storm surges is taken from: [NOAA Extremes](https://tidesandcurrents.noaa.gov/est/est_station.shtml?stnid=8518750)

## Generate transition models
To generate discrete-state transition probabilities, the SLR and storm surge levels are suitably discretized into 77 and 72 discrete states, respectively. The discrete-state transitions for Sea Level Rise (SLR) and storm surge are generated using the codes, as:
1. **slr_transitions.m** for SLR transitions
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

## Data Sources for SLR projections
- Fox-Kemper, B., H.T. Hewitt, C. Xiao, G. Aðalgeirsdóttir, S.S. Drijfhout, T.L. Edwards, N.R. Golledge, M. Hemer, R.E. Kopp, G. Krinner, A. Mix, D. Notz, S. Nowicki, I.S. Nurhati, L. Ruiz, J.-B. Sallée, A.B.A. Slangen, and Y. Yu, 2021: Ocean, Cryosphere and Sea Level Change. In Climate Change 2021: The Physical Science Basis. Contribution of Working Group I to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change [Masson-Delmotte, V., P. Zhai, A. Pirani, S.L. Connors, C. Péan, S. Berger, N. Caud, Y. Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. Lonnoy, J.B.R. Matthews, T.K. Maycock, T. Waterfield, O. Yelekçi, R. Yu, and B. Zhou (eds.)]. Cambridge University Press, Cambridge, United Kingdom and New York, NY, USA, pp. 1211–1362, doi:10.1017/9781009157896.011.
- Kopp, R. E., Garner, G. G., Hermans, T. H. J., Jha, S., Kumar, P., Reedy, A., Slangen, A. B. A., Turilli, M., Edwards, T. L., Gregory, J. M., Koubbe, G., Levermann, A., Merzky, A., Nowicki, S., Palmer, M. D., & Smith, C. (2023). The Framework for Assessing Changes To Sea-Level (FACTS) v1.0: A platform for characterizing parametric and structural uncertainty in future global, relative, and extreme sea-level change. Geoscientific Model Development, 16, 7461–7489. https://doi.org/10.5194/gmd-16-7461-2023
- Garner, G. G., T. Hermans, R. E. Kopp, A. B. A. Slangen, T. L. Edwards, A. Levermann, S. Nowikci, M. D. Palmer, C. Smith, B. Fox-Kemper, H. T. Hewitt, C. Xiao, G. Aðalgeirsdóttir, S. S. Drijfhout, T. L. Edwards, N. R. Golledge, M. Hemer, G. Krinner, A. Mix, D. Notz, S. Nowicki, I. S. Nurhati, L. Ruiz, J-B. Sallée, Y. Yu, L. Hua, T. Palmer, B. Pearson, 2021. IPCC AR6 Sea Level Projections. Version 20210809. Dataset accessed [2024-10-27] at https://doi.org/10.5281/zenodo.5914709.
  
## Citation

If you use this repository in your research or projects, please cite our paper as follows:


