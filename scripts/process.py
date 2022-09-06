# Script to parse IBLP output.
# The script may produce auxiliary folder/files in run directory.
# Some os function requires Python 3.5+

import sys
import os
import collections
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parseOutput(outputDir, versions, scenarios, writeCSV=True, filename="summary.csv"):
    """
    The function parse the output file in the given directory.
    Assume the subfolders hierarchy: outputDir/param_scenario_name/testset_name/file.out.
    The result will also be written to a .csv file if not specified.
    Input:
        outputDir: string, a path to the parent output directory
        writeCSV: boolean, whether to save the results in structured format
    Return:
        a pandas dataframe containing results
    """
    # may move those up as input in the future...
    # then need to match keywords and fields
    keywords = {
        'obj' : 'Objective value',
        'status' : 'Result',
        'gap' : 'Gap',
        'nodes' : 'Enumerated nodes',
        'cpu' : 'CPU seconds',
        'wall' : 'Wallclock seconds'
    }

    status = {
        'OPT': 'Optimal solution found',
        'TIME' : 'Stopped on time limit'
    }

    results = collections.defaultdict(list)
    etol = np.finfo(float).eps

    # iterate over versions, scenarios, datasets, and files in each folder to read results
    for v in versions:
        for s in scenarios:
            resultDir = os.path.join(outputDir, v, s)
            if os.path.isdir(resultDir) == False:
                continue
            # iterate over different datasets available
            with os.scandir(resultDir) as dataset_it:
                for d_entry in dataset_it:
                    if d_entry.name not in dataSets:
                        continue
                    # iterate over files in the folder
                    with os.scandir(d_entry.path) as output_it:
                        for o_entry in output_it:
                            if o_entry.name.endswith(".out"):
                                # start to write result to the dictionary
                                results["dataset"].append(d_entry.name)
                                results["scenario"].append(scenarios[s])
                                results["version"].append(v)
                                results["instance"].append(
                                    os.path.splitext(o_entry.name)[0]
                                )
                                results['obj'].append('None')
                                results['status'].append('UNK')
                                results['gap'].append(100000)
                                results['nodes'].append(-1)
                                results['cpu'].append(100000)
                                results['wall'].append(-1)
                                
                                with open(o_entry.path, "r") as file:
                                    for line in file.read().splitlines():
                                        if len(line) == 0:
                                            continue
                                        if keywords["obj"] in line:
                                            results["obj"][-1] = float(line.split()[2])
                                        elif keywords["status"] in line:
                                            for stat in status:
                                                if status[stat] in line:
                                                    results["status"][-1] = stat
                                        elif (keywords["gap"] in line and
                                              "inf" not in line):
                                            results["gap"][-1] = float(line.split()[1])
                                        elif keywords["nodes"] in line:
                                            results["nodes"][-1] = line.split()[2]
                                        elif keywords["cpu"] in line and "Total" not in line:
                                            results["cpu"][-1] = float(line.split()[3])
                                        elif keywords["wall"] in line and "Total" not in line:
                                            results["wall"][-1] = float(line.split()[3])
                                if results["status"][-1] == 'OPT':
                                    results["gap"][-1] = 0.0
                                elif results["status"][-1] == 'TIME':
                                    results["cpu"][-1] = 3600
                                    results["wall"][-1] = 3600
                                                    


    #for k in results:
    #   print (k)
    #   print(len(results[k]))
    df_result = pd.DataFrame(results)

    print (df_result)
    
    # write results to .csv file
    if writeCSV:
        # df_result.to_csv(filename, mode='a', header=False, index=False) # append results only
        df_result.to_csv(filename, index=False)

    return df_result

def processTable(df, displayCols, writeLTX=False, filename="ltx_tb.txt"):
    """
    Print a summary table for required columns.
    Input:
        df: a dataframe with all info from parseOutput
        displayCol: columns to print
    """

    # separate instance to different tables
    # convert each instance related data into a dictionary
    # each data field can print to a table
    # or print a summary table where instance by row

    # obtain the list of instances
    instList = list(df.instance.unique())
    scnList = list(df.scenario.unique())
    versionList = list(df.version.unique())
    # print(instList)

    # collect required info into dict
    rsltDict = {}
    for inst in instList:
        rsltDict[inst] = {}
        if 'nw04' in inst:
            continue
        for scn in scnList:
            for v in versionList:
                cond = (
                    (df["scenario"] == scn)
                    & (df["instance"] == inst)
                    & (df["version"] == v)
                )
                df_temp = df[cond]
                if len(df_temp["dataset"].values) > 0:
                    ds = df_temp["dataset"].values[0]
                    rsltDict[inst].update(
                        {(scn, v, ds, col): df_temp[col].values[0] for col in displayCols}
                )

    # convert dict to structured df: change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient="index")
    df_forprint.columns.names = ["scn", "v", "datasets", "fields"]
    df_forprint = df_forprint.sort_index()

    # OPTION 1: print results to a single table: suggest to use when display col number < 2
    # with open('ltx_tb1.txt', 'w') as file:
    #     file.write(df_forprint.to_latex())

    # OPTION 2: for each displayCol, print a table; using slicer indexing
    if writeLTX:
        with open(filename, "w") as file:
            for col in displayCols:
                for scn in scnList:
                    file.write(df_forprint.loc[:, (scn, slice(None), col)].to_latex())

    # OPTION 3: just process table, do not print latex table to file
    # pass

    return df_forprint

def dropFilter(df, scenarios, ds):
    """
    Prepare data for plotting performance profile; running time only.
    Input:
        df: pandas dataframe output from processTable
        plotCol: columns to make single plots
        scenarios: scenarios on one plot
    """
    df = df[scenarios.values()]
    # apply index filter on solution time
    df_time = df.xs(
        (ds, "cpu"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    for c in df_time.columns:
        df_time[c] = pd.to_numeric(df_time[c],
                                   errors='coerce').replace(np.nan, 36000)
    df_solved = df.xs(
        (ds, "status"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    # filter out cases where time is < 5'' or > 3600'' for all methods
    col_list = df_time.columns.values.tolist()
    drop_easy = df_time[(df_time[col_list] < 5).all(axis=1)].index.tolist()
    drop_unsolved = df_solved[(df_solved[col_list] != 'OPT').all(axis=1)].index.tolist()
    drop_list_time = list(set(drop_easy) | set(drop_unsolved))
    #drop_list_time.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    df_solved = df.drop(drop_list_time)

    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    for c in df_gap.columns:
        df_gap[c] = pd.to_numeric(df_gap[c],
                                  errors='coerce').replace(np.nan, 100000)
    drop_no_gap = df_gap[(df_gap[col_list] > 1000).all(axis=1)].index.tolist()
    drop_list_gap = list(drop_no_gap)
    #drop_list_gap.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    df_has_soln = df.drop(drop_list_gap)

    return df_solved, df_has_soln

def plotPerfProf(
        df, plotname="perf_profile", plottitle="Performance Profile",
        xmin=0.0, xmax=None, legendnames={}, versionlegend=False
):
    """
    Generate a performance profile plot for the given dataframe.
    Assume data given are in number types.
    x-axis label: multiple of virtual best;
    y-axis label: franction of instances.
    Input:
        df: instances as index, field-to-plot as columns
        plotname: name of the plot
        fixmin: the base value used to compute ratio; using df min if not given
        xmin: the smallest x-ticker to display; set by xlim
        xmax: the largest x-ticker to display; set by xlim
        displaynames: a dictionary contains legend name; using df col name if not given
    """

    fig, ax = plt.subplots(1, 1)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()
    df["virtual_best"] = df[col_list].min(axis=1)

    for col in col_list:
        print(col)
        # for each col, compute ratio
        ratios = df[col] / df["virtual_best"]
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)

        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            np.append(uniq_ratios, xmax)  # append array at the boundary point
            np.append(cum_frac, cum_frac[-1])

        # add turning points and form series to plot
        x_val = []
        y_val = []
        x_val.append(1.0)
        y_val.append(0.0)
        if uniq_ratios[0] > 1:
            x_val.append(uniq_ratios[0])
            y_val.append(0)
        x_val.append(uniq_ratios[0])
        y_val.append(cum_frac[0])
        for j, r in enumerate(uniq_ratios[1:]):
            x_val.extend([r, r])
            y_val.extend([cum_frac[j], cum_frac[j + 1]])
        if cum_frac[-1] == 1.0:
            x_val.append(xmax)
            y_val.append(1.0)

        if legendnames:
            # , color=colors[i])
            plt.plot(x_val, y_val, label=legendnames[col])
        elif versionlegend:
            plt.plot(x_val, y_val, label=col)  # , color=colors[i])
        else:
            plt.plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.02, 1.05)
    ax.tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax.set_title(plottitle)
    ax.set_xlabel("Multiple of virtual best")
    ax.set_ylabel("Fraction of instances")
    ax.legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.9),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)


def plotCumProf(df, plotname="cum_profile", plottitle = "Cumulative Profile",
                legendnames={}, versionlegend=False):

    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, wspace=0)
    ax = gs.subplots(sharey=True)

    df_time = df.xs(
        (ds, "cpu"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    for c in df_time.columns:
        df_time[c] = pd.to_numeric(df_time[c],
                                   errors='coerce').replace(np.nan, 100000)
    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    for c in df_gap.columns:
        df_gap[c] = pd.to_numeric(df_gap[c],
                                  errors='coerce').replace(np.nan, 100000)
    col_list = df_time.columns.values.tolist()
    time_buckets = range(0, 3600)

    for col in col_list:
        print(col)
        times = df_time[col]
        print(times)
        cum_cnt = np.sum(np.array([times <= t for t in time_buckets]), axis=1)
        print("##############")
        print("############## printing cum_cnt for CumProf")
        print("##############")
        print(cum_cnt)
        cum_frac = cum_cnt / len(df)
        x_val = []
        if legendnames:
            ax[0].plot(time_buckets, cum_frac, label=legendnames[col])
        elif versionlegend:
            ax[0].plot(time_buckets, cum_frac, label=col)  # , color=colors[i])
        else:
            ax[0].plot(time_buckets, cum_frac, label=col[0])

    ax[0].set_xlim(0, 3599)
    ax[0].set_ylim(0.0, 1)
    ax[0].tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Fraction of instances")

    gap_buckets = np.linspace(0, 100, 1000)

    for col in col_list:
        print(col)
        gaps = df_gap[col]
        print(gaps)
        cum_cnt = np.sum(np.array([gaps <= g for g in gap_buckets]), axis=1)
        cum_frac = cum_cnt / len(df_gap)
        #print(cum_frac)
        x_val = []
        if legendnames:
            ax[1].plot(gap_buckets, cum_frac, label=legendnames[col])
        elif versionlegend:
            ax[1].plot(gap_buckets, cum_frac, label=col)  # , color=colors[i])
        else:
            ax[1].plot(gap_buckets, cum_frac, label=col[0])

    ax[1].set_xlim(0.0, 20)
    ax[1].tick_params(axis="both", direction="in", right=True)
    ax[1].legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.95),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    # set other figure elements
    ax[1].set_xlabel("Gap")
    ax[1].label_outer()

    fig.suptitle(plottitle)
    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)
    # fig.savefig("./performance/barchart/"+plotname+'.eps', format='eps', dpi=600)

def plotBaselineProf(
        df, baseline, plotname="base_profile", plottitle="Baseline Profile",
        xmin=0.0, xmax=None, legendnames={}, versionlegend=False
):
    """
    Generate a performance profile plot for the given dataframe.
    Assume data given are in number types.
    x-axis label: multiple of virtual best;
    y-axis label: franction of instances.
    Input:
        df: instances as index, field-to-plot as columns
        plotname: name of the plot
        fixmin: the base value used to compute ratio; using df min if not given
        xmin: the smallest x-ticker to display; set by xlim
        xmax: the largest x-ticker to display; set by xlim
        displaynames: a dictionary contains legend name; using df col name if not given
    """

    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, wspace=0)
    ax = gs.subplots(sharey=True)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()

    for col in col_list:
        if col == baseline or col[0] == "virtual_best":
            continue
        print(col)
        # for each col, compute ratio
        ratios = df[col] / df[baseline]
        print(df[col])
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place
        print(uniq_ratios)
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)

        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            uniq_ratios = np.append(uniq_ratios, xmax)  # append array at the boundary point
            cum_frac = np.append(cum_frac, cum_frac[-1])

        print(uniq_ratios)
        #print(cum_frac)

        # Values less than one are scaled differently
        if uniq_ratios[0] < 1:
            x_val = []
            y_val = []
            x_val.append(0.0)
            y_val.append(0.0)
            x_val.append(uniq_ratios[0])
            y_val.append(0)
            x_val.append(uniq_ratios[0])
            y_val.append(cum_frac[0])
            for j, r in enumerate(uniq_ratios[1:]):
                if r > 1:
                    x_val.append(r)
                    y_val.append(cum_frac[j])
                    break
                x_val.extend([r, r])
                # j is indexed starting at zero, not one!
                y_val.extend([cum_frac[j], cum_frac[j + 1]])

            if legendnames:
                # , color=colors[i])
                ax[0].plot(x_val, y_val, label=legendnames[col])
            elif versionlegend:
                ax[0].plot(x_val, y_val, label=col)  # , color=colors[i])
            else:
                ax[0].plot(x_val, y_val, label=col[0])  # , color=colors[i])

        # add turning points and form series to plot
        x_val = []
        y_val = []
        if uniq_ratios[0] >= 1:
            x_val.append(1.0)
            y_val.append(0.0)
            j = 0
        if uniq_ratios[0] > 1:
            x_val.append(uniq_ratios[0])
            y_val.append(0)
        x_val.append(uniq_ratios[j])
        y_val.append(cum_frac[j])
        
        for k, r in enumerate(uniq_ratios[j+1:]):
            x_val.extend([r, r])
            y_val.extend([cum_frac[k+j], cum_frac[k+j+1]])

        if legendnames:
            # , color=colors[i])
            ax[1].plot(x_val, y_val, label=legendnames[col])
        elif versionlegend:
            ax[1].plot(x_val, y_val, label=col)  # , color=colors[i])
        else:
            ax[1].plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(-0.02, 1.05)
    ax[0].tick_params(axis="both", direction="in", right=True)

    ax[1].set_xlim(1, xmax)
    ax[1].label_outer()
    ax[1].tick_params(axis="both", direction="in", right=True)
    ax[1].legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.05),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    fig.supxlabel("Ratio of baseline")
    fig.supylabel("Fraction of instances")
    fig.suptitle(plottitle)
    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)

if __name__ == "__main__":

    dataSets = [
        'MIPLIB2017'
    ]

    versions = [
        "2.10",
        "before",
        "after-old-style",
        "after",
        "f462681a"
    ]
    
    # Output parent path
    outputDir = "/mnt/c/Users/tkral/Documents/Projects/Cbc/output"

    scenarios = {
        'default' : 'Default Settings',
    #    'p4': 'Four Cores'
    }
    ################# Process & Save | Load from CSV ###################
    # specify summary file name
    file_csv = "summary_"+dataSets[0]+".csv"

    # if len(args) == 0:
    if 1:
        df_r = parseOutput(
            outputDir, versions, scenarios, writeCSV=True, filename=file_csv
        )
    else:
        try:
            df_r = pd.read_csv(file_csv)
            set_cond = (df_r["scenario"].isin(scenarios.values())) | (
                df_r["dataset"].isin(dataSets)
            )
            df_r = df_r[set_cond]
        except FileNotFoundError:
            print("{} does not exist in current directory.".format(file_csv))
        else:
            print("Reading from", file_csv)

    # columns to process and print
    displayCols = {
        "cpu": "CPU Search Time",
        "nodes": "Number of Processed Nodes",
        "gap": "Final Gap",
        "status": "Status"
    }

    df_proc = processTable(df_r, displayCols)

    ################### Make Performance Profile ####################
    # columns to compare in the plot
    plotCols = {
        "cpu": ["CPU Time", 25],
        "nodes": ["Nodes Processed", 50],
    }
    # plotCols = {}

    # manual input example:
    # for k in scenarios:
    #     if '01' in k:
    #         scenarios[k] = 'linkingBranchStrategy'
    #     else:
    #         scenarios[k] = 'fractionalBranchStrategy'

    #baseline = None
    baseline = ('Default Settings', '2.10')
    if len(versions) > 1:
        versionlegend = True
    else:
        versionlegend = False
        
    for ds in dataSets:
        df_solved, df_has_soln = dropFilter(df_proc, scenarios, ds)
        for col in plotCols:
            df_sub = df_solved.xs(
                (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()
            for c in df_sub.columns:
                df_sub[c] = pd.to_numeric(df_sub[c],
                                          errors='coerce').replace(np.nan, 100000)
            print("")
            print("Creating performance profile for "+col)
            print("")
            plotPerfProf(
                df_sub, plotname="perf_" + col + "_" + ds,
                plottitle = "Performance Profile: "+plotCols[col][0]+" ("+ds+")",
                xmin = 0.0, xmax=plotCols[col][1],
                versionlegend = versionlegend
            )
            if baseline is not None: 
                print("")
                print("Creating baseline profile for "+col)
                print("")
                plotBaselineProf(
                    df_sub, baseline = baseline,
                    plotname="base_"+baseline[0]+"_"+col+"_"+ds,
                    plottitle = "Baseline Profile: "+plotCols[col][0]+" vs "+baseline[1],
                    xmax=plotCols[col][1],
                    versionlegend = versionlegend
                )
            plotCumProf(df_has_soln, plotname="cum_" + col + "_" + ds,
                        plottitle="Cumulative Profile: "+plotCols[col][0]+" ("+ds+")",
                        versionlegend = versionlegend
)
        if baseline is not None: 
            print("")
            print("Creating baseline profile for gap")
            print("")
            df_gap = df_has_soln.xs(
                (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()
            df_baseline_has_gap = df_gap.drop(df_gap[df_gap[baseline] == 0].index.to_list())
            # plotBaselineProf(
            #     df_baseline_has_gap, baseline = baseline,
            #     plotname="base_" + baseline[0] + "_" + "gap_" + ds,
            #     plottitle = "Baseline Profile: Gap ("+ds+")",
            #     xmax=25, versionlegend = versionlegend
            # )
    
