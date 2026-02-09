# turbulence_AUSCAT

Project undertaken by Gen Tolhurst (gen.tolhurst@unimelb.edu.au) under the supervision of Prof Todd Lane. 
Research is undertaken at The School of Geography, Earth and Atmospheric Sciences (SGEAS) at The University of Melbourne.
Project duration from August 2025 to February 2026.

<img width="2013" height="3420" alt="flow_chart_pic" src="https://github.com/user-attachments/assets/e336accd-4125-4981-b406-658ce61d008a" />


## Module for calculating indices
cat_indices.py

## Module for evaluation and trend tables
cat_evaluation.py

## Module for plotting
auscat_plots.py  
continuous_colormaps_rgb_0-1 is for colourmaps used in asucat_plots.py

## Calculate and summarise the data
streamlined_turbulence_code.ipynb - this prepares BARRA-R data for evaluation, uses cat_indices.py to calculate the requested index, calculates the 99th percentile from BARRA-R baseline 1990-2009 period, calculates the monthly percentiles for historical evalaution and summarises all experiments in terms of frequency of exceeding the 99th percentile as defined by BARRA-R. Data is produced for subsequent analysis. This process is quite time/ memory intensive. Developed using 28 CPUs and 126 GB of memory. Takes several hours to compute.

## Pre-prepare data
preprocess_BARRA_R.ipynb - obsolete

## Calculate indices
CAT_indices.ipynb - obsolete. Used to develop index calculations
contains formulas used for calculating indices

## Model validation and selection
CAT_evaluation.ipynb - obsolete. replaced by cat_evaluation.py module. Used to develop and explore evaluation for models. Lots of different exploratory evaluations.


## Supplementary info
(BARPA_BARRA_map_extent.ipynb takes a sample of data and visualises the extents of the BARRA-R and BARPA-R domains on the world map)
(list_of_model_variables.txt)


## Small sample of BARPA data to test index calculations
(sample_BARPA-R_data.ipynb for a small time slice of many layers to test indices calculations)
We compare our index calculations in compare_CAT_indices_using_BARRA-R_sample.ipynb
