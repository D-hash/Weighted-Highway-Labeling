# Weighted Highway Labeling

Source code accompanying the paper **“Provably Minimal and Reconfigurable Highway Cover Labeling on Weighted Digraphs,”** submitted to the *Data Mining and Knowledge Discovery* journal.

This repository provides implementations for constructing minimal **Highway Labelings (HLs)** on:

- Unweighted, undirected graphs
- Weighted, directed graphs

It also includes algorithms for efficiently updating the labeling index after inserting or deleting landmarks.

## Dependencies

The project uses the [NetworKit](https://github.com/networkit/networkit) library for graph processing.

## Constructing a Highway Labeling

After loading a graph using NetworKit, create an HL object with:

```cpp
HighwayLabelling(
    graph,
    numberOfLandmarks,
    landmarkSelectionStrategy,
    numberOfLandmarkChanges,
    dynamicExperimentType
);
```

The constructor arguments are:

1. `graph` — The loaded NetworKit graph.
2. `numberOfLandmarks` — The number of landmarks to select.
3. `landmarkSelectionStrategy` — An integer from `1` to `4`:
   - `1` — Degree ordering
   - `2` — Approximate betweenness
   - `3` — Distance-\(\log(n)\)-bounded dominating set
   - `4` — Approximate closeness
4. `numberOfLandmarkChanges` — The number of landmark updates (which can be also 0).
5. `dynamicExperimentType` — The type of dynamic experiment:
   - `1` — Incremental
   - `2` — Decremental
   - `3` — Mixed

## Building the Index

Use the following methods to construct the labeling index:

- `FarhanConstruction()` — Constructs an HL for unweighted, undirected graphs.
- `ConstructDirWeighHL()` — Constructs an HL for weighted, directed graphs.

`FarhanConstruction()` implements the construction algorithm presented in:

> Muhammad Farhan, Qing Wang, Yu Lin, and Brendan D. McKay.  
> “A Highly Scalable Labelling Approach for Exact Distance Queries in Complex Networks.”  
> *EDBT 2019*, pp. 13–24.

## Querying Distances

After constructing the index, query distances using:

- `QueryDistance()` — For unweighted, undirected graphs.
- `DirectedQueryDistance()` — For weighted, directed graphs.

## Updating Landmarks

### Adding Landmarks

- `AddLandmarkUnweighted()` — Promotes a vertex to a landmark in an unweighted, undirected graph.
- `AddLandmarkDirected()` — Promotes a vertex to a landmark in a weighted, directed graph.

### Removing Landmarks

- `RemoveLandmarkUnweighted()` — Demotes a landmark to a standard vertex in an unweighted, undirected graph.
- `RemoveLandmarkDirected()` — Demotes a landmark to a standard vertex in a weighted, directed graph.
