This folder contains scripts and outputs for the greenocean group's Global Carbon Budget submission for the year 2023.

Contents:

- `./submissionData2023`: contains data (not available on github), we are using model version RIV12
- `GCB2023_ocean_model_protocol_final.pdf`: protocol as received from J. Hauck
- `./scripts` contains scripts that produce output, with the batch submission file `runGCB.bsub`, that can be run using `sbatch < runGCB.bsub` 
    - Note: currently `runGCB.bsub` requires some python packages that are contained in an environment called swamp2 (notably arrow), if you are running this in the future and it's not working, check imports.
    - scripts used to make each variable are documented in SPREADSHEET (see References)
- `./notebooks` contains notebooks for visualizing outputs
- `./ModelCodeAndRunFiles` contains the archive of the model Fortran code and associated run files 


References:
- The following notion page documents the code, forcing, run id's, etc: [NOTION](https://nettle-pajama-b85.notion.site/GCB-2023-run-specifications-and-monitor-setup-81add128b34a4d87a7630e03fc542ee3?pvs=4)
- The following google spreadsheet, editable by TJŠJ and viewable by all, contains a list of all variables required by the GCB, what the output files are, where they are made, etc: [SPREADSHEET](https://docs.google.com/spreadsheets/d/186mFWSIaPWu7X_RldiGPrPzIuqyNupUR7XwoCBkMaIk/edit?usp=sharing)
- The following notebook contains a visualization of the output of this year's model runs and a comparison with last year's, where applicable:
[NOTEBOOK](https://nbviewer.org/github/tjarnikova/GCB2023/blob/main/notebooks/visualiseAllGCBOutputs.ipynb)
- The following notebook contains a visualization of the submitted model air-sea CO2 flux (version CAL12) and a comparison with the TOM12 version and the 2022 version:
[NOTEBOOK](https://nbviewer.org/github/tjarnikova/GCB2023/blob/main/notebooks/visualiseModelFlux.ipynb)
