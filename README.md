# H5N1_viral_kinetics

This GitHub repository provides the code and data to reproduce the analysis in [Viral kinetics of H5N1 infections in dairy cattle (to be uploaded to biorxiv)](https://www.biorxiv.org).

## Running the code
All statistical models can be run using the scripts 'ct_value_model.R' and 'pcr_to_inf_virus.R'. These scripts will fit the stan models to the data (see below) and output the fitted stan model to the 'fit_stan_models' subdirectory. The 'ct_values_figure.R' script can then be run to produce the figures included in the manuscript. The numerical values quoted in the manuscript are also printed at the end of the script. Please note that: the CT value models are fit with 4,000 iterations and a burn in of 1,000 iterations, with 4 chains; and the log-titre models are fit with 10,000 iterations and a burn in of 2,000 iterations, with 4 chains. This means that the figure script can take a moderate time to run. The models will provide reasonable posterior estimates (for most of the models) with only 2000 iterations (and a burn in period of 500 iterations). This will speed up the run time if just checking that the code works.



## Data
The data required to run the analysis is available in the 'data' subdirectory. The data has been obtained from [Halwe et.al 2024](https://www.nature.com/articles/s41586-024-08063-y), [Caserta et.al. 2024](https://www.nature.com/articles/s41586-024-07849-4), [Baker et.al. 2024](https://www.nature.com/articles/s41586-024-08166-6), and [Facciuolo et al. 2025](https://www.nature.com/articles/s41564-025-01998-6)
