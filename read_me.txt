This repository is made for:
'Potential influence of the marine carbon cycle on the multiequilibria window of the Atlantic Meridional Overturning Circulation'

Authors: Boot, A., von der Heydt, A.S., and Dijkstra, H.A.
Submitted to Geophysical Research Letters (06/2023)

Contact person: Amber Boot (she/they; d.boot@uu.nl)

This repository contains several subfolders:
CMIP6 fit:
- Citations for the used ESM models to construct the CMIP6 fits 
- Mask used in the analysis representing box t on a 1x1 degree ESM grid
- Script to integrate ESM data over box t region
- Processed data sets

* Unprocessed data sets can be downloaded from https://esgf-node.llnl.gov/search/cmip6/

Continuation sequences:
- .txt files with the continuation sequences used in AUTO to get to results in the main text
- Using AUTO I used the following methodology:
- Load in data: a = load('coupled_model_v1')
- Time integrate into steady state: r1 = run(a)
- Continue in Ea using the sequences provided in the .txt file

Data:
- All data files from all the experiments used in the paper

Observations:
- Observational data used to check whether the used fits are in a reasonable range for Es.

Plotting scripts:
- Scripts used to make Figures 2, 3 S1, S2 and S3