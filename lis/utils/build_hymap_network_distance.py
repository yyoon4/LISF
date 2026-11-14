#!/usr/bin/env python3
"""Add explicit HYMAP routing distances to a COPY of a localization map.

This is an offline utility, not part of the LIS time-step loop. It preserves
the input weight/sites variables. Requires numpy and netCDF4. No model or
private research data are included in this utility.
"""
from pathlib import Path
import argparse
import json
import shutil

import numpy as np
from netCDF4 import Dataset

FILL = -9999.0
SPATIAL_DIMS = ("north_south", "east_west")
WEIGHT_DIMS = ("count",) + SPATIAL_DIMS


def read_grid(dataset, name):
    var = dataset.variables[name]
    if var.dimensions != SPATIAL_DIMS:
        raise ValueError(f"{name}: expected dimensions {SPATIAL_DIMS}, got {var.dimensions}")
    if "scale_factor" in var.ncattrs() or "add_offset" in var.ncattrs():
        raise ValueError(f"{name}: packed variables are not supported")
    return np.asarray(np.ma.filled(var[:].astype(np.float64), np.nan))


def routing_forest(next_x, next_y, edge_m, active):
    """Index a directed single-downstream network. Edge u->next(u) costs edge_m[u].

    Coordinates are one-based, as in HYMAP. Nonpositive x is an outlet;
    missing coordinates/inactive cells do not supply network connections.
    Depth is path length to the tree root, NOT straight-line distance.
    """
    shape = next_x.shape
    if any(a.shape != shape for a in (next_y, edge_m, active)):
        raise ValueError("Routing array shapes differ")
    ny, nx = shape
    n = nx * ny
    x, y, edge = next_x.ravel(), next_y.ravel(), edge_m.ravel()
    valid = active.ravel().astype(bool) & np.isfinite(x) & (x != FILL)
    parent = np.full(n, -1, dtype=np.int64)
    links = np.flatnonzero(valid & (x > 0))
    if np.any(~np.isfinite(y[links])) or np.any(x[links] != np.floor(x[links])) or np.any(y[links] != np.floor(y[links])):
        raise ValueError("Noninteger or missing downstream coordinates")
    if np.any((x[links] > nx) | (y[links] < 1) | (y[links] > ny)):
        raise ValueError("Downstream coordinates outside the static grid")
    destination = (y[links].astype(np.int64) - 1) * nx + x[links].astype(np.int64) - 1
    links, destination = links[valid[destination]], destination[valid[destination]]
    if np.any(~np.isfinite(edge[links])) or np.any(edge[links] <= 0):
        raise ValueError("Connected routing edges require finite positive HYMAP_grid_distance in metres")
    parent[links] = destination
    children = [[] for _ in range(n)]
    for child in links:
        children[parent[child]].append(int(child))
    tin = np.full(n, -1, dtype=np.int64)
    tout = np.full(n, -1, dtype=np.int64)
    depth_m = np.full(n, np.nan)
    clock = 0
    for root in np.flatnonzero(valid & (parent < 0)):
        depth_m[root] = 0.0
        stack = [(int(root), False)]
        while stack:
            node, closing = stack.pop()
            if closing:
                tout[node] = clock
                continue
            tin[node] = clock
            clock += 1
            stack.append((node, True))
            for child in reversed(children[node]):
                depth_m[child] = depth_m[node] + edge[child]
                stack.append((child, False))
    if np.any(valid & (tin < 0)):
        raise ValueError("Directed routing cycle detected; cannot define a single river-path distance")
    return shape, valid, tin, tout, depth_m


def distance_from_source(forest, source):
    shape, valid, tin, tout, depth_m = forest
    if source < 0 or source >= valid.size or not valid[source]:
        raise ValueError("Localization source is outside the active routing network")
    # Ancestor OR descendant on the same directed chain. Sharing a downstream
    # confluence does NOT connect two separate tributaries to one another.
    upstream = (tin >= tin[source]) & (tin < tout[source])
    downstream = (tin <= tin[source]) & (tout > tin[source])
    connected = valid & (upstream | downstream)
    distance = np.full(valid.size, FILL, dtype=np.float32)
    distance[connected] = (np.abs(depth_m[connected] - depth_m[source]) / 1000.0).astype(np.float32)
    return distance.reshape(shape)


def source_cells(sites, count):
    flat = sites.ravel()
    positive = np.isfinite(flat) & (flat > 0)
    values = flat[positive]
    if np.any(values != np.floor(values)) or np.any(values > count):
        raise ValueError("sites must contain one-based integer weight-slice IDs")
    result = np.full(count, -1, dtype=np.int64)
    for cell, value in zip(np.flatnonzero(positive), values.astype(np.int64)):
        if result[value - 1] >= 0:
            raise ValueError(f"Ambiguous source: site ID {value} occurs at multiple cells")
        result[value - 1] = cell
    if np.any(result < 0):
        raise ValueError("Some weight-slice IDs have no source coordinate in sites")
    return result


def build(weight_map, static_file, output):
    weight_map, static_file, output = map(Path, (weight_map, static_file, output))
    if output.exists() or output.resolve() in (weight_map.resolve(), static_file.resolve()):
        raise ValueError("Output must be a new file; source files are never modified")
    with Dataset(weight_map) as src, Dataset(static_file) as static:
        weight = src.variables["weight"]
        if weight.dimensions != WEIGHT_DIMS or weight.dtype.kind != "f":
            raise ValueError(f"weight must be floating point with dimensions {WEIGHT_DIMS}")
        if "scale_factor" in weight.ncattrs() or "add_offset" in weight.ncattrs():
            raise ValueError("Packed weights are not supported")
        for attribute in ("_FillValue", "missing_value"):
            if attribute in weight.ncattrs():
                missing = np.asarray(weight.getncattr(attribute), dtype=np.float64)
                if missing.size != 1 or not np.all(np.isfinite(missing)) or np.any(missing != FILL):
                    raise ValueError(f"weight {attribute} must be -9999 for HYMAP")
        # Match Fortran nf90_get_var: validate the raw stored values. Automatic
        # valid_range/missing-value masking must not hide values HYMAP will read.
        weight.set_auto_maskandscale(False)
        if "network_distance_km" in src.variables:
            raise ValueError("Source already contains network_distance_km; use the original weight map")
        sites = read_grid(src, "sites")
        next_x = read_grid(static, "HYMAP_flow_direction_x")
        next_y = read_grid(static, "HYMAP_flow_direction_y")
        edge_m = read_grid(static, "HYMAP_grid_distance")
        basin_mask = read_grid(static, "HYMAP_basin_mask")
        if sites.shape != next_x.shape or weight.shape[1:] != sites.shape:
            raise ValueError("Localization map and HYMAP static grid shapes differ")
        units = str(getattr(static.variables["HYMAP_grid_distance"], "units", "m")).strip().lower()
        if units not in ("m", "meter", "meters", "metre", "metres"):
            raise ValueError(f"HYMAP_grid_distance is interpreted in metres, not {units!r}")
        # Where both files carry explicit spatial coordinates, require agreement.
        for name in ("latitude", "longitude", "lat", "lon", *SPATIAL_DIMS):
            if name in src.variables and name in static.variables:
                a = np.asarray(np.ma.filled(src.variables[name][:].astype(np.float64), np.nan))
                b = np.asarray(np.ma.filled(static.variables[name][:].astype(np.float64), np.nan))
                if a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-7, equal_nan=True):
                    raise ValueError(f"Grid coordinate mismatch: {name}")
        forest = routing_forest(next_x, next_y, edge_m, basin_mask > 0)
        cells = source_cells(sites, weight.shape[0])
        if np.any(~forest[1][cells]):
            raise ValueError("A localization source lies outside the HYMAP active network")
        # The source is copied byte-for-byte before appending the new variable.
        # Work is one source plane at a time; no source*global-grid distance cube.
        try:
            shutil.copyfile(weight_map, output)
            with Dataset(output, "a") as dst:
                kwargs = {"fill_value": FILL}
                if dst.data_model.startswith("NETCDF4"):
                    kwargs.update(zlib=True, complevel=4, chunksizes=(1, min(sites.shape[0], 32), min(sites.shape[1], 32)))
                distances = dst.createVariable("network_distance_km", "f4", WEIGHT_DIMS, **kwargs)
                distances.units = "km"
                distances.long_name = "Directed HYMAP routing path distance between weight source and target"
                distances.edge_convention = "u -> next(u), edge length = HYMAP_grid_distance[u] in metres"
                distances.connectivity = "same directed chain only; existing weight support retained"
                distances.source_index = "one-based sites value indexes count"
                dst.network_distance_builder = Path(__file__).name
                dst.network_distance_static_file = static_file.name
                kept = missing_path = 0
                for site, cell in enumerate(cells):
                    plane = np.asarray(weight[site], dtype=np.float64)
                    support = np.isfinite(plane) & (plane != FILL)
                    if np.any(support & ((plane < 0) | (plane > 1))):
                        raise ValueError(f"Site {site + 1}: area-distance-max requires weights in [0,1]")
                    if np.any(~np.isfinite(plane)):
                        raise ValueError(f"Site {site + 1}: nonfinite weights")
                    distance = distance_from_source(forest, int(cell))
                    missing_path += int(np.count_nonzero(support & (distance == FILL)))
                    distance[~support] = FILL
                    kept += int(np.count_nonzero(distance != FILL))
                    distances[site] = distance
                    if (site + 1) % 100 == 0 or site + 1 == cells.size:
                        print(f"Network distances: {site + 1}/{cells.size}", flush=True)
        except BaseException:
            output.unlink(missing_ok=True)
            raise
    return {"status": "COMPLETE", "sites": int(cells.size), "connected_pairs": kept,
            "supported_pairs_without_direct_path": missing_path, "output": str(output),
            "source_map_modified": False, "model_run": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight-map", type=Path, required=True)
    parser.add_argument("--static-file", type=Path, required=True,
                        help="The same HYMAP grid used by the weight map; HYMAP_grid_distance is in metres")
    parser.add_argument("--output", type=Path, required=True, help="New map containing unchanged weights plus network_distance_km")
    args = parser.parse_args()
    print(json.dumps(build(args.weight_map, args.static_file, args.output), indent=2))


if __name__ == "__main__":
    main()
