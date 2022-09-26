#!/usr/bin/env python3

import collections
import random
import sys

random.seed(int(220923))

def edges_from_path(path: str):
    vertices = tuple([-2] + [ int(v.rstrip("+-")) for v in path.split(",") ] + [-1])
    return (list(zip(vertices[:-1], vertices[1:])), vertices)

paths = {}
path_set = {}
edge_counts = collections.Counter()

vlens = { -1: 0, -2: 0 }

dup_paths = {}

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('S'):
        prefix, vid, seq = line.split(sep=None)
        vlens[int(vid)] = len(seq)

    if not line.startswith('P'):
        print(line)
        continue
    
    prefix, pid, path, rest = line.split(sep=None, maxsplit=3)

    edges, vertices = edges_from_path(path)
    if vertices in path_set:
        print("# Path", pid, "is a duplicate of path", path_set[vertices], file=sys.stderr)
        dup_paths[pid] = ((prefix, pid, path, rest), edges)
        continue
    path_set[vertices] = pid

    paths[pid] = ((prefix, pid, path, rest), edges)
    edge_counts.update(edges)

# Sample paths
to_delete = []
pids = list(paths.keys())
for pid in random.sample(pids, len(pids)):
    _, edges = paths[pid]
    # Check if counters are all gt 1
    if all((edge_counts[edge]>1 for edge in edges)):
        to_delete.append(pid)
        edge_counts.subtract(edges)

to_delete = set(to_delete)
to_keep = paths.keys() - to_delete

print("# Input paths:", len(paths) + len(dup_paths), file=sys.stderr)
print("# Duplicated paths:", len(dup_paths), file=sys.stderr)
print("# Output paths: ", len(paths) - len(to_delete), file=sys.stderr)

i = 0
path_id = {}
edge_to_paths = collections.defaultdict(set)
for pid in to_keep:
    (line, edges) = paths[pid]
    for edge in edges:
        edge_to_paths[edge].update([pid])
    i += 1
    path_id[pid] = i
    print("# PATH", i, " -> ", pid, file=sys.stderr)
    print("\t".join(line))

def ok_path(pid, edges, cedge, curr, path, curr_dist, max_recomb, min_dist):
    if cedge >= len(edges):
        return path
    edge = edges[cedge]
    candidates = edge_to_paths[edge]
    if curr in candidates:
        result = ok_path(pid, edges, cedge+1, curr, path + [path_id[curr]], curr_dist + vlens[edge[1]], max_recomb, min_dist)
        if result is not None:
            return result
    # We have to introduce a recombination
    if curr_dist < min_dist or max_recomb == 0:
        return None
    for candidate in candidates:
        if candidate == pid or candidate == curr:
            continue
        result = ok_path(pid, edges, cedge+1, candidate, path + [path_id[candidate]], vlens[edge[1]], max_recomb-1, min_dist)
        if result is not None:
            return result
    return None

min_dist = 100
for pid in to_delete:
    (_, edges) = paths[pid]
    result = None
    for max_recomb in range(1, 10):
        for candidate in edge_to_paths[edges[0]]:
            if candidate == pid:
                continue
            result = ok_path(pid, edges, 1, candidate, [], 0, max_recomb, min_dist)
            if result is not None:
                print(f"# Path {pid} is a mosaic with {max_recomb} recombinations of {result}", file=sys.stderr)
                break
        if result is not None:
            break
    if result is None:
        print(f"# Path {pid} is a NOT a mosaic of the path cover with minimum recombination distance of {min_dist}", file=sys.stderr)

