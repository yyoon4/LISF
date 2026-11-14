# Optional river-distance weighting for HYMAP SWOT updates

This option changes how an existing localization map weights a source's analysis
increment. It does not change SWOT QC/BC, EnKF gains, source increment caps,
restart files, routing physics, or MPI communication.

## Methods

| Method | Effective weight | Combination of sources |
| --- | --- | --- |
| `area-max` (default) | Existing map weight | Select the greatest weight; retain the first source on ties |
| `area-distance-max` (opt-in) | Existing map weight times a cosine river-distance taper | Select the greatest effective weight; retain the first source on ties |

The names describe the intended area-based input map. The code **does not
recalculate or invert an area ratio**; it uses the supplied `weight` values.
The maximum is over weights, not over absolute increments. For a recipient
and ensemble member, the selected increment is `source_dH * effective_weight`.
There is no averaging or summation of multiple sources.

For network distance `d` and cutoff length `L`, both in km:

```text
T(d/L) = (1 + cos(pi*d/L))/2    for 0 <= d < L
         0                     for d >= L
effective_weight = existing_weight * T(d/L)
```

This is a cosine taper, not linear attenuation. At distances 0, L/4, L/2,
3L/4, and L, the factors are 1, about 0.854, 0.5, about 0.146, and 0.
Each source's weight is reduced or unchanged. The *selected* source can change,
so the final increment, added storage, and discharge need not decrease.

## Configuration

An old configuration needs no edits. Omitting the new method selects `area-max`
and does not read any distance variable or apply a taper.

To explicitly select the existing behavior:

```text
HYMAP3 localization weighting method: "area-max"
```

To select the new behavior, retain the existing localization enable/window
settings and point to a map with the additional distance variable:

```text
HYMAP3 use localization update in DA: 1
HYMAP3 localization update window size: 10
HYMAP3 localization weight map: ./input/localization_with_network_distance.nc
HYMAP3 localization weighting method: "area-distance-max"
HYMAP3 localization distance cutoff (km): 100.0
```

`100.0` is an illustration, **not a recommended or validated length**. A finite
positive cutoff is required only for active `area-distance-max` nests. Use the
`HYMAP2` prefix for HYMAP2. Supply one value per nest, as with the existing
localization settings. Localization-off nests do not read or taper a map.

The existing pixel window still applies: DX=5, 10, or 15 searches a square
of up to 11x11, 21x21, or 31x31 cells, respectively. `L` is a separate physical
distance along routing links; DX is not a count of river links and is not
converted to km. The effective support is the intersection of the existing
map, pixel window, current MPI-local availability, and `d < L`. Increasing L
does not restore sources missing at MPI boundaries. This patch does not fix
that separate cross-rank exchange limitation.

## Distance map

The new method requires this additional variable in the **same NetCDF file**
as `weight` and `sites`:

```text
float network_distance_km(count, north_south, east_west)
    units = "km"
    _FillValue = -9999.f
```

It uses the same source index and target grid as `weight`. A source's `sites`
value is a one-based index on `count`, not a routing index or station ID.
Unknown or disconnected paths use -9999. No Euclidean/pixel distance fallback
or inferred interpretation of a variable named `distance` is used. Source-to-
self distance is zero. Missing original weights stay missing. Weights in the
new mode must be finite values in [0,1], apart from the existing missing value.
Invalid weights are not silently clamped. The default method retains existing
weight semantics, including its original behavior for unusual finite values.

During initialization, the taper is applied once in place to the loaded
weights, using one temporary 2-D distance plane. No second permanent distance
cube or per-time-step network search is added. HYMAP3 retains rank hyperslab
reads; HYMAP2 retains its existing global weight allocation. There is no change
to the maximum-weight update routine.

## Offline map preparation

`build_hymap_network_distance.py` writes a new copy of a weight map, preserving
all original weight/sites values, and appends explicit routing distances.
It requires Python, NumPy, and netCDF4:

```bash
python -B lis/utils/build_hymap_network_distance.py \
  --weight-map input/localization_weights.nc \
  --static-file input/lis_input.nc \
  --output input/localization_with_network_distance.nc
```

Use the **same grid and routing topology** that generated the localization
map. The static variables are `HYMAP_flow_direction_x`,
`HYMAP_flow_direction_y`, `HYMAP_grid_distance`, and `HYMAP_basin_mask`, on
`(north_south,east_west)`. Next-cell coordinates are one-based. For each
directed edge `u -> next(u)`, the length is `HYMAP_grid_distance[u]` in metres.
The distance is the sum of these lengths, divided by 1000. This is the model's
routing distance, not an independently surveyed channel-centerline length.

Only ancestors or descendants on the same directed chain are connected.
Two tributaries that merely share a downstream confluence are not joined.
Existing map support determines which of the upstream/downstream connections
can be used; the utility does not create new weight support. An outlet's
unused edge is not added. Masked cells break connectivity. Duplicate source
IDs, routing cycles, and invalid connected edges must be reconciled explicitly.
The utility does not guess a source location or repair the static grid.

The output file must not already exist. Neither the source map nor the static
file is modified. The output includes one additional compressed distance
cube; it is prepared once and can be reused for different cutoff lengths.

## Multiple-source mixing is a separate future option

Positive normalized weights produce an average within the source increment
range, not a sum that grows with source count. They can reduce independent
errors, but shared bias, error correlation, and different hydrologic phases
remain. Area/distance weights measure spatial representativeness, not
observation accuracy. Opposite increment signs are not proof of a bad input.

In particular, `sum(w*T*dH)/sum(w*T)` cancels the taper when there is only
one source, and cancels any common attenuation across sources. A future
mixing method must specify both the mixing proportions and the attenuation
of the total update. It must not accidentally undo distance damping or count
correlated observations as independent. No mixing option is implemented here.

## Validation scope

The synthetic tests cover the distance kernel and initialization of the actual
HYMAP2/HYMAP3 modules with minimal LIS/ESMF/NetCDF API stubs, including legacy
defaults, rank offsets, and malformed option inputs. The standalone builder
uses real NetCDF files in its synthetic tests and checks the original map is
unchanged, directed-chain membership, edge indexing, units, and rejection of
ambiguous source IDs/cycles. The existing actual update-routine regressions
remain applicable because those routines are unchanged.

These checks are not a full LIS/ESMF/MPI build or a hydrologic skill evaluation.
Before using the new option operationally, compile the full model and compare
both methods from the same restart, input observations, and MPI partition.
No length scale is selected by this patch.
