#!/usr/bin/env python
# coding: utf-8

# after you have run "streamlined_turbulence_code.ipynb", run this to perform the model evaluation
# you should have these files: glob.glob(f"/scratch/v46/gt3409/{turbulence_index}/{P}hPa/freq-above-p99/{turbulence_index}-{P}hPa-monthly-freq-above-p99_AUS-15_*_BOM_BARPA-R_v1-r1_6hr.nc") 
# and glob.glob(f"/scratch/v46/gt3409/{turbulence_index}/{P}hPa/percentiles/{turbulence_index}-{P}hPa-monthly-percentiles_AUS-15_*_BOM_BARPA-R_v1-r1_6hr.nc")

import xarray as xr
import glob
import numpy as np
import seaborn as sns
from scipy import stats
import calendar
import pandas as pd
from turbulence_AUSCAT.cat_indices import calc_turbulence_indices, windspeed, VWS
from xarray.groupers import SeasonResampler
import matplotlib.pyplot as plt

path= "/scratch/v46/gt3409/turbulence_AUSCAT/"

mid_lat_slice = slice(-50,-25)
lon_slice = slice(90,195)

baseline_time_range = np.arange(1990,2009+1)
baseline_time_slice = slice("1990", "2009")

P0 = 1000
step_size=  0.1545


list_evaluation = ['evaluation_BARRA-R_r1i1p1f1',
                   # 'evaluation_ERA5_r1i1p1f1',
                  ]

list_historical = ['historical_ACCESS-CM2_r4i1p1f1', 
                   'historical_ACCESS-ESM1-5_r6i1p1f1',
                   'historical_CESM2_r11i1p1f1', 
                   'historical_CMCC-ESM2_r1i1p1f1',
                   'historical_EC-Earth3_r1i1p1f1',
                   'historical_MPI-ESM1-2-HR_r1i1p1f1',
                   'historical_NorESM2-MM_r1i1p1f1',
                  ]

list_ssp126 = [
                 'ssp126_ACCESS-CM2_r4i1p1f1', 
                 'ssp126_ACCESS-ESM1-5_r6i1p1f1',
                 'ssp126_CESM2_r11i1p1f1',
                 'ssp126_CMCC-ESM2_r1i1p1f1',
                 'ssp126_EC-Earth3_r1i1p1f1',
                 'ssp126_MPI-ESM1-2-HR_r1i1p1f1',
                 'ssp126_NorESM2-MM_r1i1p1f1',
              ]

list_ssp370 = ['ssp370_ACCESS-CM2_r4i1p1f1',
                 'ssp370_ACCESS-ESM1-5_r6i1p1f1',
                 'ssp370_CESM2_r11i1p1f1',
                 'ssp370_CMCC-ESM2_r1i1p1f1',
                 'ssp370_EC-Earth3_r1i1p1f1',
                 'ssp370_MPI-ESM1-2-HR_r1i1p1f1',
                 'ssp370_NorESM2-MM_r1i1p1f1',
              ]

list_ssp585 = ['ssp585_ACCESS-CM2_r4i1p1f1',
                 'ssp585_EC-Earth3_r1i1p1f1']

list_future = list_ssp126 + list_ssp370 + list_ssp585


# # Make tables for evaluation:

# ## Kovmogolov-Smirnov two sample test for model evaluation

# evaluation
evaluation = ['evaluation_BARRA-R_r1i1p1f1']

def plot_kstest(ds_eval, ds_hist, turbulence_index, time_selection, evaluation, run, pvalue, significance, ax=None, show_years=False):
    if ax is None:
        fig, ax = plt.subplots()
    ax.set_title(f"{run} in {time_selection}")        

    if show_years:
        # show each single year
        for year in baseline_time_range:
            sns.lineplot(data = ds_eval.sel({"time": (ds_eval.time.dt.year==year)}).mean("time"),
                         x=turbulence_index,
                         y="quantile",
                         color="r",
                         alpha=0.1,
                         ax=ax)
        for year in baseline_time_range:
            sns.lineplot(data = ds_hist.sel({"time": (ds_hist.time.dt.year==year)}).mean("time"),
                         x=turbulence_index,
                         y="quantile",
                         color="k",
                         alpha=0.1,
                         ax=ax)
            
    # empirical CDF for ks test        
    sns.ecdfplot(data = ds_eval.mean("time")[turbulence_index], color="k", label=evaluation, ax=ax,)
    sns.ecdfplot(data = ds_hist.mean("time")[turbulence_index], color="r", label=run, ax=ax,)
    
    # add pvalue to plot
    ax.text(x=0.98, y=0.02, s= f"pvalue = {pvalue:.3f}{significance}", transform = ax.transAxes, va="bottom", ha="right",)
    ax.set_xlabel(turbulence_index)
    ax.legend(loc="upper left", fontsize=10)
    return 

def calc_kstest_table(ds = None, 
                      turbulence_index=None,
                      pressure=None,
                      evaluation = evaluation, 
                      list_historical=list_historical, 
                      outfile=None,
                      show_plot = False):
    """"""
    kstest_table = [["sample1", "sample2", "time_selection", "mon/seas/yr", "pvalue", "significance"]]
    for run in list_historical:
        # annual
        ds_eval = ds.sel({"run":evaluation})
        ds_hist = ds.sel({"run":run})
    
        # calculate pvalue
        ks_statistic, pvalue = stats.ks_2samp(ds_eval.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),
                                              ds_hist.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),)
        
        if pvalue < 0.01:
            significance = "**"
        elif pvalue < 0.05:
            significance = "*"
        else:
            significance = ""
    
        # update ks test stats
        kstest_table.append([evaluation, run, "annual", "year", f"{pvalue:.3f}", significance])
    
        if show_plot:
            plot_kstest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)

        # cool v warm season
        for time_selection in ["MJJASO", "NDJFMA"]:
            ds_hist = ds.sel({"run":run,}).resample({"time":SeasonResampler([time_selection], drop_incomplete=True)}).mean()
            ds_eval = ds.sel({"run":evaluation,}).resample({"time":SeasonResampler([time_selection], drop_incomplete=True)}).mean()
    
            # calculate pvalue
            ks_statistic, pvalue = stats.ks_2samp(ds_eval.resample({"time":"YE-OCT"}).mean()[turbulence_index].values.flatten(),
                                                  ds_hist.resample({"time":"YE-OCT"}).mean()[turbulence_index].values.flatten(),)
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
    
            # update ks test stats
            kstest_table.append([evaluation, run, time_selection, "6M", f"{pvalue:.3f}", significance])
    
            if show_plot:
                plot_kstest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
    
        # seasonal
        for time_selection in ["DJF", "MAM", "JJA", "SON"]:
            ds_hist = ds.sel({"run":run, "time": (ds.time.dt.season==time_selection)}).resample({"time":"QS-DEC"}).mean()
            ds_eval = ds.sel({"run":evaluation, "time": (ds.time.dt.season==time_selection)}).resample({"time":"QS-DEC"}).mean()
    
            # calculate pvalue
            ks_statistic, pvalue = stats.ks_2samp(ds_eval.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),
                                                  ds_hist.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),)
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
    
            # update ks test stats
            kstest_table.append([evaluation, run, time_selection, "season", f"{pvalue:.3f}", significance])
    
            if show_plot:
                plot_kstest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
    
        # monthly
        for time_selection in np.arange(1,12+1):
            ds_hist = ds.sel({"run":run, "time": (ds.time.dt.month==time_selection)})
            ds_eval = ds.sel({"run":evaluation, "time": (ds.time.dt.month==time_selection)})
    
            # calculate pvalue
            ks_statistic, pvalue = stats.ks_2samp(ds_eval.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),
                                                  ds_hist.resample({"time":"YE"}).mean()[turbulence_index].values.flatten(),)
            # print(f"{run}: pvalue {pvalue:.3f}, ks_stat = {ks_statistic}")
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
    
            kstest_table.append([evaluation, run, calendar.month_name[time_selection], "month", f"{pvalue:.3f}", significance])     
    
            if show_plot:
                plot_kstest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
    
    kstest_df = pd.DataFrame(columns = kstest_table[0], data = kstest_table[1:])
    if outfile is None:
        outfile = f"{path}/evaluation_kstest_table_{turbulence_index}-{pressure}hPa.csv"
    kstest_df.to_csv(outfile)
    return kstest_df


def calc_ttest_table(ds=None, q=0.33, turbulence_index=None, P = None, evaluation = evaluation, outfile=None,):
    """"""
    # q ttest value in the mid lat box to evaluation spatial conherence
    ttest_table = [["sample1", "sample2", "time_selection", "mon/seas/yr", f"p{int(q*100)}_pvalue", "significance"]]
    for run in list_historical:
        # annual
        ds_eval = ds.sel({"run":evaluation}).squeeze()
        ds_hist = ds.sel({"run":run}).squeeze()
    
        # calculate pvalue
        stat, pvalue = stats.ttest_ind(ds_eval.resample({"time":"YE"}).mean()[turbulence_index],
                                       ds_hist.resample({"time":"YE"}).mean()[turbulence_index])
    
        pvalue = xr.Dataset(data_vars={"pval" : (["lat", "lon"],  pvalue),}, 
                            coords= {"lat": ds_hist.lat, "lon": ds_hist.lon},
                           )["pval"].quantile(q, dim=["lat", "lon"])
        
        # mark significance
        if pvalue < 0.01:
            significance = "**"
        elif pvalue < 0.05:
            significance = "*"
        else:
            significance = ""
    
        # update t test stats
        ttest_table.append([evaluation, run, "annual", "year", f"{pvalue:.3f}", significance])
    
        if False:
            plot_kde_test(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
        
        # seasonal
        ds_hist_seas = ds.sel({"run":run,}).resample({"time":"QS-DEC"}).mean().sel({"time":baseline_time_slice})
        ds_eval_seas = ds.sel({"run":evaluation,}).resample({"time":"QS-DEC"}).mean().sel({"time":baseline_time_slice})
        
        for time_selection in ["DJF", "MAM", "JJA", "SON"]:
            ds_hist = ds_hist_seas.sel({"time": (ds_hist_seas.time.dt.season==time_selection)}).squeeze()
            ds_eval = ds_eval_seas.sel({"time": (ds_eval_seas.time.dt.season==time_selection)}).squeeze()
    
            # calculate pvalue
            stat, pvalue = stats.ttest_ind(ds_eval.resample({"time":"YE"}).mean()[turbulence_index],
                                       ds_hist.resample({"time":"YE"}).mean()[turbulence_index])
    
            pvalue = xr.Dataset(data_vars={"pval" : (["lat", "lon"],  pvalue),}, 
                            coords= {"lat": ds_hist.lat, "lon": ds_hist.lon},
                           )["pval"].quantile(q, dim=["lat", "lon"])
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
    
            # update t test stats
            ttest_table.append([evaluation, run, time_selection, "season", f"{pvalue:.3f}", significance])
    
            if False:
                plot_kde_ttest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)

        # cool v warm season
        for time_selection in ["MJJASO", "NDJFMA"]:
            ds_hist = ds.sel({"run":run,}).resample({"time":SeasonResampler([time_selection], drop_incomplete=True)}).mean()
            ds_eval = ds.sel({"run":evaluation,}).resample({"time":SeasonResampler([time_selection], drop_incomplete=True)}).mean()
    

            # calculate pvalue
            stat, pvalue = stats.ttest_ind(ds_eval.resample({"time":"YE-OCT"}).mean()[turbulence_index],
                                       ds_hist.resample({"time":"YE-OCT"}).mean()[turbulence_index])
    
            pvalue = xr.Dataset(data_vars={"pval" : (["lat", "lon"],  pvalue),}, 
                            coords= {"lat": ds_hist.lat, "lon": ds_hist.lon},
                           )["pval"].quantile(q, dim=["lat", "lon"])
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
                
            # update t test stats
            ttest_table.append([evaluation, run, time_selection, "6M", f"{pvalue:.3f}", significance])
    
            if False:
                plot_kde_ttest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
    
        # monthly
        for time_selection in np.arange(1,12+1):
            ds_hist = ds.sel({"run":run, "time": (ds.time.dt.month==time_selection)}).squeeze()
            ds_eval = ds.sel({"run":evaluation, "time": (ds.time.dt.month==time_selection)}).squeeze()
    
            # calculate pvalue
            stat, pvalue = stats.ttest_ind(ds_eval.resample({"time":"YE"}).mean()[turbulence_index],
                                       ds_hist.resample({"time":"YE"}).mean()[turbulence_index])
    
            pvalue = xr.Dataset(data_vars={"pval" : (["lat", "lon"],  pvalue),}, 
                            coords= {"lat": ds_hist.lat, "lon": ds_hist.lon},
                           )["pval"].quantile(q, dim=["lat", "lon"])
            
            
            if pvalue < 0.01:
                significance = "**"
            elif pvalue < 0.05:
                significance = "*"
            else:
                significance = ""
    
            ttest_table.append([evaluation, run, calendar.month_name[time_selection], "month", f"{pvalue:.3f}", significance])     
    
            if False:
                plot_kde_ttest(ds_eval, ds_hist, time_selection, evaluation, run, pvalue)
    
    
    ttest_df = pd.DataFrame(columns = ttest_table[0], data = ttest_table[1:])
    if outfile is None:
        outfile = f"{path}/evaluation_ttest_{turbulence_index}-{P}hPa_p{int(q*100)}_table.csv"
    ttest_df.to_csv(outfile)
    return ttest_df

def combined_significance_table(turbulence_index, P, path=None):
    ds = xr.open_mfdataset(glob.glob(f"/scratch/v46/gt3409/{turbulence_index}/{P}hPa/percentiles/{turbulence_index}-{P}hPa-monthly-percentiles_AUS-15_*_BOM_BARPA-R_v1-r1_6hr.nc"),
                      combine="nested",
                      concat_dim="run",
                          compat='no_conflicts')
    ktest_df = calc_kstest_table(ds = ds, 
                          turbulence_index=turbulence_index,
                          pressure=P,
                          evaluation = evaluation[0], 
                          list_historical=list_historical, 
                          outfile=None)
    
    
    ds = xr.open_mfdataset(glob.glob(f"/scratch/v46/gt3409/{turbulence_index}/{P}hPa/freq-above-p99/{turbulence_index}-{P}hPa-monthly-freq-above-p99_AUS-15_*_BOM_BARPA-R_v1-r1_6hr.nc"),
                      combine="nested",
                      concat_dim="run",
                           coords="different",
                          join="outer",
                          compat='no_conflicts').sel({"time":baseline_time_slice, "lon":lon_slice })
    q=0.33
    ttest_df = calc_ttest_table(ds = ds,
                                 q=q, 
                                 turbulence_index=turbulence_index,
                                 P = P, 
                                 evaluation = evaluation[0], )

    # read the test p values from the saved files
    evaluation_kstest = pd.read_csv(f"{path}/evaluation_kstest_table_{turbulence_index}-{P}hPa.csv", index_col=0,).replace(np.nan, "")
    evaluation_ttest = pd.read_csv(f"{path}/evaluation_ttest_{turbulence_index}-{P}hPa_p{int(q*100)}_table.csv", index_col=0,)[["p33_pvalue", "significance"]].replace(np.nan, "")

    evaluation_combined = pd.concat([evaluation_kstest, evaluation_ttest], axis=1, join="inner")
    evaluation_combined["combined_significance"] = evaluation_kstest["significance"] + evaluation_ttest["significance"]
    evaluation_combined.to_csv(f"{path}/evaluation_combined_tests_table_{turbulence_index}-{P}hPa.csv")
    return evaluation_combined





