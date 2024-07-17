import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

genes = []
with open("genes_HLA_short.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

modes = ["0", "3", "5", "10", "rec"]

def parse_logs(directory):
    times = {}
    memory = {}
    edit_distance = {}
    for gene in genes:
        for mode in modes:
            dir_path = f"{directory}/{gene}/reads_{mode}_split/"
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
                                        times[(gene, mode, count)] = value
                                    elif unit == "ms":
                                        times[(gene, mode, count)] = value / 1000
                                    elif unit == "µs":
                                        times[(gene, mode, count)] = value / 1000000
                                else:
                                    match = re.search(r'Done in (\d+)', line)
                                    if match:
                                        times[(gene, mode, count)] = int(match.group(1))
                            if "Maximum resident set size (kbytes)" in line:
                                memory[(gene, mode, count)] = int(line.split()[-1])
                if filename.endswith(".gaf"):
                    count_gaf += 1
                    if "old" in directory:
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f.readlines():
                                cigar = line.strip().split("\t")[12].split(",")[0]
                                edit_score = edit_distance_from_cigar(cigar)
                                if "recombination" in line:
                                    edit_score += 4
                            edit_distance[(gene, mode, count_gaf)] = edit_score
                    else:
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f.readlines():
                                if line.startswith("@CO"):
                                    edit_distance[(gene, mode, count_gaf)] = int(line.split("\t")[1])
                    
    return times, memory, edit_distance

def parse_cigar(cigar_string):
    
    operations = []
    current_length = ""
    for char in cigar_string:
        if char.isdigit():
            current_length += char
        else:
            operations.append((int(current_length), char))
            current_length = ""
    return operations

def edit_distance_from_cigar(cigar_string):
    
    operations = parse_cigar(cigar_string)
    edit_distance = 0
    for length, op in operations:
        if op in 'IDX':
            edit_distance += length
    return edit_distance

old_times, old_memory, old_edit = parse_logs("output/HLA/old_recgraph")
new_times, new_memory, new_edit = parse_logs("output/HLA/new_recgraph")

data_time = []
data_memory = []

for gene in genes:
    for mode in modes:
        for count in range(1, max(len(old_times), len(new_times)) + 1):
            if (gene, mode, count) in old_times and (gene, mode, count) in new_times:
                data_time.append((gene, mode, old_times[(gene, mode, count)], "old"))
                data_time.append((gene, mode, new_times[(gene, mode, count)], "new"))

            if (gene, mode, count) in old_memory and (gene, mode, count) in new_memory:
                data_memory.append((gene, mode, old_memory[(gene, mode, count)], "old"))
                data_memory.append((gene, mode, new_memory[(gene, mode, count)], "new"))



data_edit = []

for gene, mode, count in old_edit:
    data_edit.append((gene, mode, old_edit[(gene, mode, count)], "old"))
    data_edit.append((gene, mode, new_edit[(gene, mode, count)], "new"))

print(data_edit)
df_edit = pd.DataFrame(data_edit, columns=["Gene", "Mode", "Edit", "Version"])

                    
fig, axes = plt.subplots(nrows=len(genes), ncols=1, figsize=(10, 5 * len(genes)))

for i, gene in enumerate(genes):
    df_edit_gene = df_edit[df_edit["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Edit", hue="Version", data=df_edit_gene)
    axes[i].set_title(f"Edit distance comparison for Gene: {gene}")
    axes[i].set_yscale("symlog")
    axes[i].set_ylim(-1, 10000)
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Edit distance")
    axes[i].legend(title='Version')
    
plt.tight_layout()
plt.show()

df_time = pd.DataFrame(data_time, columns=["Gene", "Mode", "Time", "Version"])

fig, axes = plt.subplots(nrows=len(genes), ncols=1, figsize=(10, 5 * len(genes)))

for i, gene in enumerate(genes):
    df_time_gene = df_time[df_time["Gene"] == gene]

    sns.boxplot(ax=axes[i], x="Mode", y="Time", hue="Version", data=df_time_gene)
    axes[i].set_title(f"Execution time comparison for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Time (s)")
    axes[i].legend(title='Version')

plt.tight_layout()
plt.show()

df_memory = pd.DataFrame(data_memory, columns=["Gene", "Mode", "Memory", "Version"])
fig, axes = plt.subplots(nrows=len(genes), ncols=1, figsize=(10, 5 * len(genes)))

for i, gene in enumerate(genes):
    df_memory_gene = df_memory[df_memory["Gene"] == gene]

    sns.boxplot(ax=axes[i], x="Mode", y="Memory", hue="Version", data=df_memory_gene)
    axes[i].set_title(f"Memory usage comparison for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Memory used (KB)")
    axes[i].legend(title='Version')

plt.tight_layout()
plt.show()