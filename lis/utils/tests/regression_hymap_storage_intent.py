#!/usr/bin/env python3
"""Check storage argument contracts and execute actual HYMAP2/3 updates.

Synthetic single-cell budgets only; no LIS inputs or experiment are used.
The legacy OUT declaration is a contract error. A numeric failure with any
particular compiler is not assumed: --allow-legacy-out can execute that code
for comparison, without asserting that its behavior is Fortran-conforming.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


DRIVER = r"""
program storage_regression
  use storage_test_module
  implicit none
  real :: river,flood,total_before
  integer :: version,j,checks
  checks=0
  do version=2,3
    ! No-flow calls must preserve both input reservoirs.
    river=100.; flood=40.
    call route(version,10.,0.,0.,0.,0.,river,flood)
    call pair(river,flood,100.,40.,'quiescent routing')

    ! Net external hydrologic flux: +30 river, +20 flood.
    river=100.; flood=40.
    total_before=river+flood
    call route(version,10.,2.,1.,5.,3.,river,flood)
    call pair(river,flood,130.,60.,'routing preserves input plus flux')
    call close(river+flood-total_before,50.,'routing total budget')

    ! Signed edge fluxes are retained, without triggering zero clipping.
    river=100.; flood=40.
    call route(version,10.,-2.,-1.,-1.,-2.,river,flood)
    call pair(river,flood,110.,30.,'signed routing flux')

    ! Withdraw exactly the available donor storage.
    river=100.; flood=40.
    call route(version,10.,10.,4.,0.,0.,river,flood)
    call pair(river,flood,0.,0.,'empty donor at supported flux')

    ! Preserve the existing clipping behavior for unsupported excess flux.
    river=3.; flood=2.
    call route(version,10.,10.,10.,0.,0.,river,flood)
    call pair(river,flood,0.,0.,'routing existing nonnegative clip')

    river=100.; flood=40.
    call runoff(version,10.,.25,0.,river,flood)
    call pair(river,flood,100.,40.,'zero runoff preserves storage')

    river=100.; flood=40.
    total_before=river+flood
    call runoff(version,10.,.25,8.,river,flood)
    call pair(river,flood,160.,60.,'positive runoff partition')
    call close(river+flood-total_before,80.,'positive runoff total budget')

    river=100.; flood=40.
    call runoff(version,10.,.25,-4.,river,flood)
    call pair(river,flood,70.,30.,'negative runoff within both reservoirs')

    ! Flood storage is insufficient: the remaining sink comes from river.
    river=100.; flood=2.
    total_before=river+flood
    call runoff(version,10.,.5,-4.,river,flood)
    call pair(river,flood,62.,0.,'negative runoff redistribution')
    call close(river+flood-total_before,-40.,'negative runoff total budget')

    ! Requested loss exceeds total storage: retain the existing lower bound.
    river=10.; flood=2.
    call runoff(version,10.,.5,-4.,river,flood)
    call pair(river,flood,0.,0.,'sink limited by available storage')

    river=100.; flood=40.
    call runoff(version,10.,0.,8.,river,flood)
    call pair(river,flood,180.,40.,'runoff zero flooded fraction')
    river=100.; flood=40.
    call runoff(version,10.,1.,8.,river,flood)
    call pair(river,flood,100.,120.,'runoff full flooded fraction')

    river=100.; flood=40.
    total_before=river+flood
    call route(version,10.,2.,1.,5.,3.,river,flood)
    call runoff(version,10.,.25,8.,river,flood)
    call pair(river,flood,190.,80.,'routing then runoff')
    call close(river+flood-total_before,130.,'combined total budget')

    river=100.; flood=40.
    do j=1,100
      call route(version,1.,1.,1.,4.,2.,river,flood)
      call runoff(version,1.,.25,2.,river,flood)
    enddo
    call pair(river,flood,550.,190.,'repeated updates retain previous state')
    call close(river+flood,740.,'repeated total budget')
  enddo
  print *, 'PASS actual HYMAP2/3 storage updates; checks=',checks

contains
  subroutine route(version,dt,rout,fout,rin,fin,rs,fs)
    integer, intent(in) :: version
    real, intent(in) :: dt,rout,fout,rin,fin
    real, intent(inout) :: rs,fs
    if(version==2) then
      call HYMAP2_calc_stonxt(dt,rout,fout,rin,fin,rs,fs)
    else
      call HYMAP3_calc_stonxt(dt,rout,fout,rin,fin,rs,fs)
    endif
  end subroutine

  subroutine runoff(version,dt,fraction,q,rs,fs)
    integer, intent(in) :: version
    real, intent(in) :: dt,fraction,q
    real, intent(inout) :: rs,fs
    if(version==2) then
      call HYMAP2_calc_runoff(dt,fraction,q,rs,fs)
    else
      call HYMAP3_calc_runoff(dt,fraction,q,rs,fs)
    endif
  end subroutine

  subroutine pair(rs,fs,expected_rs,expected_fs,label)
    real, intent(in) :: rs,fs,expected_rs,expected_fs
    character(*), intent(in) :: label
    call close(rs,expected_rs,label//' river')
    call close(fs,expected_fs,label//' flood')
  end subroutine

  subroutine close(actual,expected,label)
    real, intent(in) :: actual,expected
    character(*), intent(in) :: label
    if(.not.(abs(actual-expected)<=1.e-5*max(1.,abs(expected)))) then
      print *, 'FAIL ',label,' H',version,' actual=',actual,' expected=',expected
      stop 1
    endif
    checks=checks+1
  end subroutine
end program
"""


def extract_subroutine(source, name):
    pattern = rf"(?ims)^\s*subroutine\s+{re.escape(name)}\s*\(.*?^\s*end\s+subroutine\s+{re.escape(name)}\s*$"
    matches = re.findall(pattern, source)
    if len(matches) != 1:
        raise ValueError(f"Expected one {name}, found {len(matches)}")
    return matches[0]


def locate_root(explicit):
    if explicit:
        return Path(explicit).resolve()
    for candidate in (Path.cwd(), *Path(__file__).resolve().parents):
        if (candidate / "lis/routing/HYMAP2_router/HYMAP2_modelMod.F90").is_file():
            return candidate
    raise ValueError("Cannot locate LIS source; pass --source-root REPOSITORY")


def compiler_command(explicit):
    if explicit:
        command = shlex.split(explicit)
    elif os.environ.get("FC"):
        command = shlex.split(os.environ["FC"])
    else:
        command = next(([name] for name in ("gfortran", "gfortran-15", "gfortran-14", "gfortran-13")
                        if shutil.which(name)), [])
    if not command or not shutil.which(command[0]):
        return None
    return command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", help="LIS repository root")
    parser.add_argument("--compiler", help="GNU Fortran command (defaults to FC or gfortran)")
    parser.add_argument("--check-extraction", action="store_true",
                        help="Check extraction/contracts only; no Fortran execution")
    parser.add_argument("--allow-legacy-out", action="store_true",
                        help="Permit legacy OUT declarations for numeric comparison only")
    args = parser.parse_args()
    root = locate_root(args.source_root)
    routines = []
    for version in (2, 3):
        source = (root / f"lis/routing/HYMAP{version}_router/HYMAP{version}_modelMod.F90").read_text()
        for operation in ("runoff", "stonxt"):
            name = f"HYMAP{version}_calc_{operation}"
            routine = extract_subroutine(source, name)
            if not args.allow_legacy_out:
                for argument in ("rivsto", "fldsto"):
                    if not re.search(rf"(?im)^\s*real\s*,\s*intent\(inout\)\s*::\s*{argument}\b", routine):
                        raise ValueError(f"{name}: input storage {argument} must have INTENT(INOUT)")
            routines.append(routine)
    if args.check_extraction:
        print("Extracted four actual storage subroutines; Fortran execution NOT performed")
        return 0
    compiler = compiler_command(args.compiler)
    if compiler is None:
        print("SKIP: GNU Fortran unavailable; Fortran regression NOT executed", file=sys.stderr)
        return 77
    wrapper = ("module LIS_coreMod\n  implicit none\nend module\n"
               "module storage_test_module\n  implicit none\ncontains\n" +
               "\n".join(routines) + "\nend module\n")
    with tempfile.TemporaryDirectory(prefix="hymap-storage-regression-") as directory:
        work = Path(directory)
        (work / "storage.F90").write_text(wrapper)
        (work / "driver.F90").write_text(DRIVER)
        for optimization in ("-O0", "-O2"):
            command = compiler + [optimization, "-g", "-fcheck=all", "-fbacktrace", "-Wall", "-Wextra",
                                  "-ffree-line-length-none", "storage.F90", "driver.F90", "-o", "regression"]
            subprocess.run(command, cwd=work, check=True)
            print(f"Storage regression {optimization}", flush=True)
            subprocess.run([str(work / "regression")], cwd=work, check=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
