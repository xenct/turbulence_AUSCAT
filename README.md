# turbulence_AUSCAT

Project undertaken by Gen Tolhurst (gen.tolhurst@unimelb.edu.au or gentolhurst@gmail.com) under the supervision of Prof Todd Lane. 
Research is undertaken at The School of Geography, Earth and Atmospheric Sciences (SGEAS) at The University of Melbourne.
Project duration from August 2025 to February 2026.

## Flow chart describes analysis process
<img width="1000" height="1400" alt="flow_chart_pic" src="https://github.com/user-attachments/assets/87a4fcca-b77d-440c-9ab9-d9337ab83e82" />

## Module for calculating indices
[cat_indices.py](https://github.com/xenct/turbulence_AUSCAT/blob/main/cat_indices.py)

## Module for evaluation and trend tables
[cat_evaluation.py](https://github.com/xenct/turbulence_AUSCAT/blob/main/cat_evaluation.py)

https://github.com/xenct/turbulence_AUSCAT/blob/main/MOG_plots.ipynb

## Module for plotting
[auscat_plots.py](https://github.com/xenct/turbulence_AUSCAT/blob/main/auscat_plots.py)

[continuous_colormaps_rgb_0-1](https://github.com/xenct/turbulence_AUSCAT/blob/main/continuous_colormaps_rgb_0-1) is for colourmaps used in [asucat_plots.py](https://github.com/xenct/turbulence_AUSCAT/blob/main/asucat_plots.py)

(For context, some functions are adapted from this ACS repo I worked on https://github.com/AusClimateService/plotting_maps )

## Calculate and summarise the data
[streamlined_turbulence_code.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/streamlined_turbulence_code.ipynb) - this prepares BARRA-R data for evaluation, uses cat_indices.py to calculate the requested index, calculates the 99th percentile from BARRA-R baseline 1990-2009 period, calculates the monthly percentiles for historical evalaution and summarises all experiments in terms of frequency of exceeding the 99th percentile as defined by BARRA-R. Data is produced for subsequent analysis. This process is quite time/ memory intensive. Developed using 28 CPUs and 126 GB of memory. Takes several hours to compute.

## Plot the figures and summarise in tables
[MOG_plots.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/MOG_plots.ipynb)

## Pre-prepare data
[preprocess_BARRA_R.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/preprocess_BARRA_R.ipynb) - obsolete

## Calculate indices
[CAT_indices.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/CAT_indices.ipynb) - obsolete. Used to develop index calculations
contains formulas used for calculating indices

## Model validation and selection
[CAT_evaluation.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/CAT_evaluation.ipynb) - obsolete. replaced by cat_evaluation.py module. Used to develop and explore evaluation for models. Lots of different exploratory evaluations.

## Supplementary info
([BARPA_BARRA_map_extent.ipynb](https://github.com/xenct/turbulence_AUSCAT/blob/main/BARPA_BARRA_map_extent.ipynb) takes a sample of data and visualises the extents of the BARRA-R and BARPA-R domains on the world map)
([list_of_model_variables.txt](https://github.com/xenct/turbulence_AUSCAT/blob/main/list_of_model_variables.txt) for all available variable in BARPA-R)

## Small sample of BARPA data to test index calculations
(sample_BARPA-R_data.ipynb for a small time slice of many layers to test indices calculations)
We compare our index calculations in compare_CAT_indices_using_BARRA-R_sample.ipynb

# Ideas To Do
can be found under [issues](https://github.com/xenct/turbulence_AUSCAT/issues/1)
