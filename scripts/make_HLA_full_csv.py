import os
import re

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

recs = [0,1,2]
errors = [0,3,5]

def parse_logs(directory):
    times = {}
    memory = {}
    cells_explored = {}
    for gene in genes:
        for rec in recs:
            for error in errors:
                dir_path = f"{directory}/{gene}/{rec}/reads_{error}_split/"
                count = 0
                count_gaf = 0
                for filename in os.listdir(dir_path):
                    if filename.endswith(".log"):
                        count += 1
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f:
                                if "Done in" in line:
                                    match = re.search(r'Done in (\d+\.?\d*)([a-z]+)', line)
                                    if match:
                                        value = float(match.group(1))
                                        unit = match.group(2)
                                        if unit == "s":
                                            times[(gene, rec, error, count)] = value
                                        elif unit == "ms":
                                            times[(gene, rec, error, count)] = value / 1000
                                        elif unit == "µs":
                                            times[(gene, rec, error, count)] = value / 1000000
                                    else:
                                        match = re.search(r'Done in (\d+)', line)
                                        if match:
                                            times[(gene, rec, error, count)] = int(match.group(1))
                                if "Maximum resident set size (kbytes)" in line:
                                    memory[(gene, rec, error, count)] = int(line.split()[-1])
                    if filename.endswith(".gaf"):
                        count_gaf += 1
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f.readlines():
                                if line.startswith("@CO"):
                                    cells_explored[(gene, rec, error, count_gaf)] = float(line.split("\t")[-1][:-2])
                
    return times, memory, cells_explored

times, memory, cells_explored = parse_logs("output/HLA_full/new_recgraph")

# make csv with mean values for time, memory and cells explored grouped by gene, rec and error

import pandas as pd
import numpy as np

def calculate_means(times, memory, cells_explored):
    # Convert dictionaries to DataFrame
    times_df = pd.DataFrame(list(times.items()), columns=["key", "time"])
    memory_df = pd.DataFrame(list(memory.items()), columns=["key", "memory"])
    cells_explored_df = pd.DataFrame(list(cells_explored.items()), columns=["key", "cells_explored"])
    
    # Split key tuple into separate columns
    times_df[['gene', 'rec', 'error', 'count']] = pd.DataFrame(times_df['key'].tolist(), index=times_df.index)
    memory_df[['gene', 'rec', 'error', 'count']] = pd.DataFrame(memory_df['key'].tolist(), index=memory_df.index)
    cells_explored_df[['gene', 'rec', 'error', 'count']] = pd.DataFrame(cells_explored_df['key'].tolist(), index=cells_explored_df.index)
    
    # Drop the 'key' and 'count' columns as they are not needed for grouping
    times_df = times_df.drop(columns=['key', 'count'])
    memory_df = memory_df.drop(columns=['key', 'count'])
    cells_explored_df = cells_explored_df.drop(columns=['key', 'count'])
    
    # Group by gene, rec, and error and calculate the mean
    times_mean = times_df.groupby(['gene', 'rec', 'error']).mean().reset_index()
    memory_mean = memory_df.groupby(['gene', 'rec', 'error']).mean().reset_index()
    cells_explored_mean = cells_explored_df.groupby(['gene', 'rec', 'error']).mean().reset_index()
    
    # Merge the three dataframes
    merged_df = pd.merge(times_mean, memory_mean, on=['gene', 'rec', 'error'], how='outer')
    merged_df = pd.merge(merged_df, cells_explored_mean, on=['gene', 'rec', 'error'], how='outer')
    
    return merged_df

# Parse logs and calculate means
times, memory, cells_explored = parse_logs("output/HLA_full/new_recgraph")
mean_values_df = calculate_means(times, memory, cells_explored)

mean_values_df = mean_values_df.applymap(lambda x: f"{x:.2e}" if isinstance(x, float) else x)

# Save the dataframe to a CSV file
mean_values_df.to_csv("mean_values.csv", index=False)
