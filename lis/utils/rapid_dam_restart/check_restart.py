#!/usr/bin/env python3
"""Compare continuous/split-run RAPID restart states; never execute LIS.

Requires numpy and netCDF4. --self-test checks the comparison utility with
temporary synthetic files, not the Fortran model or a reservoir simulation.
"""

import argparse
from pathlib import Path
import shutil
import tempfile

import netCDF4
import numpy as np


IDS = ("rapid_dam_id", "rapid_dam_river_id")
STATES = ("rapid_dam_parameters", "rapid_dam_storage", "rapid_dam_inflow_next",
          "rapid_dam_outflow_prev", "rapid_dam_river_qout")
ATTRS = ("rapid_dam_state_version", "rapid_dam_tau", "rapid_dam_dt")


def compare(reference, candidate, rtol=1.e-12, atol=1.e-8):
    """Raise ValueError for missing, invalid or different restart state."""
    with netCDF4.Dataset(reference) as left, netCDF4.Dataset(candidate) as right:
        for ds in (left, right):
            for name in ("rapid_dam", "rapid_dam_member", "rivid", "rapid_dam_parameter"):
                if name not in ds.dimensions or len(ds.dimensions[name]) <= 0:
                    raise ValueError(f"invalid dimension: {name}")
            if len(ds.dimensions["rapid_dam_parameter"]) != 4:
                raise ValueError("dam parameter dimension must have length four")
            for name in ATTRS:
                if name not in ds.ncattrs():
                    raise ValueError(f"missing attribute: {name}")
                if not np.isfinite(ds.getncattr(name)):
                    raise ValueError(f"invalid attribute: {name}")
            if ds.getncattr("rapid_dam_state_version") != 1:
                raise ValueError("unsupported dam restart schema")
        for name in ATTRS:
            if left.getncattr(name) != right.getncattr(name):
                raise ValueError(f"changed attribute: {name}")
        for name in (*IDS, *STATES, "Qout"):
            if name not in left.variables or name not in right.variables:
                raise ValueError(f"missing variable: {name}")
            a, b = left[name], right[name]
            if a.dimensions != b.dimensions or a.shape != b.shape:
                raise ValueError(f"changed dimensions: {name}")
            if name in STATES and any(v.dtype.kind != "f" or v.dtype.itemsize != 8 for v in (a, b)):
                raise ValueError(f"dam state must use double precision: {name}")
            av, bv = a[:], b[:]
            if np.any(np.ma.getmaskarray(av)) or np.any(np.ma.getmaskarray(bv)):
                raise ValueError(f"unwritten/missing values: {name}")
            av, bv = np.asarray(av), np.asarray(bv)
            if not np.all(np.isfinite(av)) or not np.all(np.isfinite(bv)):
                raise ValueError(f"nonfinite values: {name}")
            if name in IDS or name == "rapid_dam_parameters":
                same = np.array_equal(av, bv)
            else:
                same = np.allclose(av, bv, rtol=rtol, atol=atol)
            if not same:
                raise ValueError(f"state differs: {name}; maximum absolute error={np.max(np.abs(av-bv)):g}")


def fixture(path):
    """A small state-file fixture for testing this utility's failure detection."""
    with netCDF4.Dataset(path, "w") as ds:
        for name, size in (("rivid", 3), ("rapid_dam", 1), ("rapid_dam_member", 2),
                           ("rapid_dam_parameter", 4)):
            ds.createDimension(name, size)
        ds.rapid_dam_state_version = 1
        ds.rapid_dam_tau = 10800.0
        ds.rapid_dam_dt = 900.0
        ds.createVariable(IDS[0], "i4", ("rapid_dam",))[:] = [20]
        ds.createVariable(IDS[1], "i4", ("rivid",))[:] = [10, 20, 30]
        ds.createVariable(STATES[0], "f8", ("rapid_dam_parameter", "rapid_dam"))[:] = [[100], [0], [1.5], [.01]]
        for name in STATES[1:4]:
            ds.createVariable(name, "f8", ("rapid_dam_member", "rapid_dam"))[:] = [[12], [15]]
        ds.createVariable(STATES[4], "f8", ("rapid_dam_member", "rivid"))[:] = [[1, 2, 3], [4, 5, 6]]
        ds.createVariable("Qout", "f4", ("rapid_dam_member", "rivid"))[:] = [[1, 2, 3], [4, 5, 6]]


def self_test():
    rejected = 0
    with tempfile.TemporaryDirectory(prefix="rapid-restart-check-") as tmp:
        reference, candidate = Path(tmp)/"reference.nc", Path(tmp)/"candidate.nc"
        fixture(reference)
        shutil.copyfile(reference, candidate)
        compare(reference, candidate)
        for name in (*IDS, *STATES, "Qout"):
            shutil.copyfile(reference, candidate)
            with netCDF4.Dataset(candidate, "r+") as ds:
                values = ds[name][:]
                values.flat[0] += 1
                ds[name][:] = values
            try:
                compare(reference, candidate)
            except ValueError:
                rejected += 1
            else:
                raise AssertionError(f"missed corruption: {name}")
        for attr, value in (("rapid_dam_state_version", 99), ("rapid_dam_tau", 1.),
                            ("rapid_dam_dt", 1.)):
            shutil.copyfile(reference, candidate)
            with netCDF4.Dataset(candidate, "r+") as ds:
                ds.setncattr(attr, value)
            try:
                compare(reference, candidate)
            except ValueError:
                rejected += 1
            else:
                raise AssertionError(f"missed changed metadata: {attr}")
        shutil.copyfile(reference, candidate)
        with netCDF4.Dataset(candidate, "r+") as ds:
            ds["rapid_dam_storage"][0, 0] = np.nan
        try:
            compare(reference, candidate)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("missed NaN storage")
    print(f"PASS: comparator accepted identical state and rejected {rejected} corruptions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", nargs="?")
    parser.add_argument("candidate", nargs="?")
    parser.add_argument("--rtol", type=float, default=1.e-12)
    parser.add_argument("--atol", type=float, default=1.e-8)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.reference and args.candidate:
        compare(args.reference, args.candidate, args.rtol, args.atol)
        print("PASS: full RAPID dam restart states agree")
    else:
        parser.error("provide reference and candidate restart files, or --self-test")


if __name__ == "__main__":
    main()
