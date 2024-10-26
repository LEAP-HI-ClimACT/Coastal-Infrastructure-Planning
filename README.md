# Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change
This repository includes the codes used in the paper: **Optimal Life-Cycle Adaptation of Coastal Infrastructure under Climate Change**.

**TLDR**: Climate change adaptation is often addressed with static policies that ignore future climate uncertainties. This work develops a climate change adaptation framework using stochastic-control methodologies, offering optimal solutions that address model uncertainty, social costs of carbon, and a wide variety of adaptation actions.

This repository contains the steps along with the necessary codes for generating components of the MDP and POMDP models for flood management and for generating the input
for dynamic programming solvers like FRTDP.  The workflow for generating the MDP input for the coastal city setting application with the two floodwall options, as described in the paper, using the provided code is outlined below.

# Steps to generate adaptation policies for coastal applications
## Datasets used for sea-level-rise simulations
- The IPCC projection dataset (for Battery tide gauge in New York) used to simulate SLR models is publicly available at: [IPCC slr projections](https://zenodo.org/records/6382554)
- The dataset for SLR trends (for Battery tide gauge in New York) used to fit the noise present in the past SLR observations is taken from: [NOAA slr trends](https://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?id=8518750)
## Datasets used for storm surge simulations
The dataset of annual extremes (for Battery tide gauge in New York) used to fit storm surges is taken from: [NOAA Extremes](https://tidesandcurrents.noaa.gov/est/est_station.shtml?stnid=8518750)

## Generate transition models

## Generate reward models

