#!/usr/bin/env python3

import collections
import random
import sys

random.seed(int(220923))

def edges_from_path(path: str):
    vertices = [ int(v.rstrip("+-")) for v in path.split(",") ]
    return (list(zip(vertices[:-1], vertices[1:])), vertices)

paths = []
edge_counts = collections.Counter()

discarded_single_vertex_paths = 0
single_vertex_paths = 0
single_vertices = set()

for line in sys.stdin:
    line = line.rstrip()
    if not line.startswith('P'):
        print(line)
        continue
    
    prefix, pid, path, rest = line.split(sep=None, maxsplit=3)

    edges, vertices = edges_from_path(path)
    if len(vertices) == 1:
        single_vertex_paths += 1

        if vertices[0] not in single_vertices:
            print(line)
            single_vertices.add(vertices[0])
        else:
            discarded_single_vertex_paths += 1
    else:
        paths.append(((prefix, pid, path, rest), edges))
        edge_counts.update(edges)

# Sample paths
to_delete = []
for i in random.sample(range(len(paths)), len(paths)):
    _, edges = paths[i]
    # Check if counters are all gt 1
    if all((edge_counts[edge]>1 for edge in edges)):
        to_delete.append(i)
        edge_counts.subtract(edges)

to_delete = set(to_delete)

for i, (line, _) in enumerate(paths):
    if i in to_delete:
        continue
    print("\t".join(line))

print("# Input paths:", single_vertex_paths+len(paths), file=sys.stderr)
print("#    --> single-vertex:", single_vertex_paths, file=sys.stderr)
print("#    --> multiple-vertex:", len(paths), file=sys.stderr)
print("# Output paths: ", len(paths) +single_vertex_paths - len(to_delete) - discarded_single_vertex_paths, file=sys.stderr)
print("#    --> single-vertex:", single_vertex_paths - discarded_single_vertex_paths, file=sys.stderr)
print("#    --> multiple-vertex:", len(paths) - len(to_delete), file=sys.stderr)

