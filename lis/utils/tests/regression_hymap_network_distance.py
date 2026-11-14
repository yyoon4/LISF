#!/usr/bin/env python3
"""Synthetic directed-network and real NetCDF tests; no research data needed."""
from pathlib import Path
import hashlib
import importlib.util
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset

SCRIPT = Path(__file__).resolve().parents[1] / "build_hymap_network_distance.py"
spec = importlib.util.spec_from_file_location("network_distance", SCRIPT)
nd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nd)


def grid():
    # A->C(10), B->C(40), C->D(20), D->E(30); F->G(7); H inactive.
    parents = [2, 2, 3, 4, -1, 6, -1, -1]
    x = np.array([(p % 4 + 1) if p >= 0 else 0 for p in parents], float).reshape(2, 4)
    y = np.array([(p // 4 + 1) if p >= 0 else 0 for p in parents], float).reshape(2, 4)
    edge = np.array([10, 40, 20, 30, 999, 7, 777, 0.], float).reshape(2, 4)
    active = np.array([1] * 7 + [0], bool).reshape(2, 4)
    return x, y, edge, active


class ForestTests(unittest.TestCase):
    def test_edge_source_and_branches(self):
        forest = nd.routing_forest(*grid())
        expected = {0: {0: 0, 2: 10, 3: 30, 4: 60},
                    2: {0: 10, 1: 40, 2: 0, 3: 20, 4: 50},
                    3: {0: 30, 1: 60, 2: 20, 3: 0, 4: 30},
                    4: {0: 60, 1: 90, 2: 50, 3: 30, 4: 0},
                    5: {5: 0, 6: 7}}
        for source, values in expected.items():
            distance = nd.distance_from_source(forest, source).ravel()
            self.assertEqual(set(np.flatnonzero(distance != nd.FILL)), set(values))
            for target, metres in values.items():
                self.assertAlmostEqual(float(distance[target]), metres / 1000., places=7)

    def test_cycle_and_self_loop(self):
        for destination in [0, 2]:
            x, y, edge, active = grid()
            x.flat[2] = destination % 4 + 1
            y.flat[2] = destination // 4 + 1
            with self.assertRaisesRegex(ValueError, "cycle"):
                nd.routing_forest(x, y, edge, active)

    def test_bad_connected_edge(self):
        for bad in [0, -1, np.nan, np.inf]:
            x, y, edge, active = grid()
            edge.flat[0] = bad
            with self.assertRaisesRegex(ValueError, "positive"):
                nd.routing_forest(x, y, edge, active)

    def test_bad_coordinates(self):
        for bad_x, bad_y in [(5, 1), (1.5, 1), (1, 0), (1, np.nan)]:
            x, y, edge, active = grid()
            x.flat[0], y.flat[0] = bad_x, bad_y
            with self.assertRaises(ValueError):
                nd.routing_forest(x, y, edge, active)

    def test_inactive_barrier(self):
        x, y, edge, active = grid()
        active.flat[2] = False
        forest = nd.routing_forest(x, y, edge, active)
        self.assertEqual(set(np.flatnonzero(nd.distance_from_source(forest, 0).ravel() != nd.FILL)), {0})
        with self.assertRaises(ValueError):
            nd.distance_from_source(forest, 2)

    def test_source_ids(self):
        np.testing.assert_array_equal(nd.source_cells(np.array([[2., -9999., 1.]]), 2), [2, 0])
        for sites, count in [([[1, 1]], 1), ([[1]], 2), ([[1.5]], 1), ([[2]], 1)]:
            with self.assertRaises(ValueError):
                nd.source_cells(np.asarray(sites, float), count)


class FileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="hymap-network-distance-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.weight = self.root / "weight.nc"
        self.static = self.root / "static.nc"
        self.output = self.root / "with_distance.nc"
        with Dataset(self.weight, "w") as ds:
            ds.createDimension("count", 2)
            ds.createDimension("north_south", 2)
            ds.createDimension("east_west", 4)
            w = ds.createVariable("weight", "f4", nd.WEIGHT_DIMS, fill_value=nd.FILL, zlib=True)
            values = np.full((2, 2, 4), .8, np.float32)
            values[0].flat[1] = nd.FILL  # source C must not revive this supported ancestor.
            values[0].flat[3] = 0.0
            w[:] = values
            s = ds.createVariable("sites", "i4", nd.SPATIAL_DIMS, fill_value=int(nd.FILL))
            values = np.full((2, 4), int(nd.FILL), np.int32)
            values.flat[2], values.flat[0] = 1, 2
            s[:] = values
            ds.description = "synthetic fixture"
        with Dataset(self.static, "w") as ds:
            ds.createDimension("north_south", 2)
            ds.createDimension("east_west", 4)
            for name, values in zip(["HYMAP_flow_direction_x", "HYMAP_flow_direction_y", "HYMAP_grid_distance", "HYMAP_basin_mask"], grid()):
                v = ds.createVariable(name, "f8", nd.SPATIAL_DIMS)
                v[:] = values
            ds.variables["HYMAP_grid_distance"].units = "m"

    def test_real_netcdf_and_unchanged_original(self):
        before = hashlib.sha256(self.weight.read_bytes()).hexdigest()
        result = nd.build(self.weight, self.static, self.output)
        self.assertEqual(result["status"], "COMPLETE")
        self.assertGreater(result["supported_pairs_without_direct_path"], 0)
        self.assertEqual(before, hashlib.sha256(self.weight.read_bytes()).hexdigest())
        with Dataset(self.weight) as old, Dataset(self.output) as new:
            self.assertEqual(new.description, old.description)
            for name in ["sites", "weight"]:
                old.variables[name].set_auto_mask(False)
                new.variables[name].set_auto_mask(False)
                np.testing.assert_array_equal(old.variables[name][:], new.variables[name][:])
                self.assertEqual(old.variables[name].dtype, new.variables[name].dtype)
            d = new.variables["network_distance_km"]
            d.set_auto_mask(False)
            self.assertEqual(d.units, "km")
            self.assertEqual(d.dimensions, nd.WEIGHT_DIMS)
            self.assertEqual(d[0].flat[1], nd.FILL)
            self.assertEqual(d[1].flat[1], nd.FILL)  # sibling branch
            self.assertAlmostEqual(float(d[0].flat[0]), .01, places=7)
            self.assertAlmostEqual(float(d[0].flat[3]), .02, places=7)
            self.assertEqual(float(new.variables["weight"][0].flat[3]), 0.)

    def test_output_collision(self):
        with self.assertRaises(ValueError):
            nd.build(self.weight, self.static, self.weight)

    def test_wrong_distance_units(self):
        with Dataset(self.static, "a") as ds:
            ds.variables["HYMAP_grid_distance"].units = "km"
        with self.assertRaisesRegex(ValueError, "metres"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_bad_weight_cleans_output(self):
        with Dataset(self.weight, "a") as ds:
            ds.variables["weight"][0, 0, 0] = 2.
        with self.assertRaisesRegex(ValueError, "weights"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_duplicate_source_rejected_before_copy(self):
        with Dataset(self.weight, "a") as ds:
            ds.variables["sites"][0, 1] = 1
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_other_missing_value_rejected(self):
        with Dataset(self.weight, "a") as ds:
            ds.variables["weight"].missing_value = -1.
        with self.assertRaisesRegex(ValueError, "-9999"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_valid_range_cannot_hide_invalid_weight(self):
        with Dataset(self.weight, "a") as ds:
            ds.variables["weight"][0, 0, 0] = 2.
            ds.variables["weight"].valid_range = np.asarray([0., 1.], np.float32)
        with self.assertRaisesRegex(ValueError, "weights"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_raw_nonfinite_weight_rejected(self):
        with Dataset(self.weight, "a") as ds:
            ds.variables["weight"][0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "Nonfinite|nonfinite"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_coordinate_mask_mismatch(self):
        for file, masked in [(self.weight, False), (self.static, True)]:
            with Dataset(file, "a") as ds:
                v = ds.createVariable("latitude", "f4", nd.SPATIAL_DIMS, fill_value=nd.FILL)
                v[:] = 40.
                if masked:
                    v[0, 0] = nd.FILL
        with self.assertRaisesRegex(ValueError, "coordinate mismatch"):
            nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())

    def test_partial_copy_cleanup(self):
        def failed_copy(source, destination):
            Path(destination).write_bytes(b"partial copy")
            raise OSError("synthetic copy failure")
        with patch.object(nd.shutil, "copyfile", failed_copy):
            with self.assertRaisesRegex(OSError, "synthetic copy failure"):
                nd.build(self.weight, self.static, self.output)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
