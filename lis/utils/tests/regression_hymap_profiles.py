#!/usr/bin/env python3
"""Compile and test the actual HYMAP2/3 volume/elevation profile subroutines.

This isolated regression uses synthetic geometry. It does not build LIS, read
model inputs, or execute an experiment. GNU Fortran is required; unavailable
compilers produce exit 77 rather than a successful test result.
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
program profile_regression
  use profile_test_module
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer, parameter :: dp=real64, nz=3
  real(dp), parameter :: bed=90._dp, bank=100._dp
  real(dp), parameter :: length=100._dp, width=10._dp
  real(dp), parameter :: ar=length*width, vmax=ar*(bank-bed)
  real(dp), parameter :: dh=1._dp/1024._dp, tol=1.e-10_dp
  real(dp), parameter :: heights(nz)=[1._dp,2._dp,5._dp]
  real(dp), parameter :: areas(5)=[100000._dp,2000._dp,1001._dp,1000._dp,500._dp]
  real(dp) :: area, bins(nz), width_increment, previous_height, previous_volume
  real(dp) :: h, v, h1, h2, v1, v2, repeat_volume
  real(dp) :: stage_points(5+nz), volume_points(3+nz)
  integer :: version, geometry, j, side, repeat, checks

  checks=0
  do version=2,3
    do geometry=1,size(areas)
      area=areas(geometry)
      width_increment=(area-ar)/length/nz
      previous_height=0._dp
      previous_volume=vmax
      do j=1,nz
        bins(j)=previous_volume+length*(width+width_increment*(real(j,dp)-.5_dp))* &
                (heights(j)-previous_height)
        previous_volume=bins(j)
        previous_height=heights(j)
      enddo

      ! Below bed, dry bed, within channel, bankfull, every table boundary,
      ! highest boundary, and well above the table. Binary offsets avoid
      ! conflating the pre-existing REAL() channel conversion with this fix.
      stage_points=[bed-1._dp,bed,bed+.25_dp,bed+5._dp,bank,bank+heights]
      do j=1,size(stage_points)
        do side=-1,1
          h=stage_points(j)+real(side,dp)*dh
          call volume(version,h,v)
          call close(v,expected_volume(h),'forward profile',version,geometry)
          call elevation(version,v,h1)
          call close(h1,max(h,bed),'elevation round trip',version,geometry)
          call volume(version,h1,v1)
          call close(v1,v,'volume round trip from elevation',version,geometry)
        enddo
      enddo
      ! Interior of each floodplain segment and above the table.
      previous_height=0._dp
      do j=1,nz
        h=bank+(previous_height+heights(j))*.5_dp
        call check_elevation(h,version,geometry)
        previous_height=heights(j)
      enddo
      call check_elevation(bank+heights(nz)+2._dp,version,geometry)
      call check_elevation(bank+heights(nz)+100._dp,version,geometry)

      if(area>ar) then
        volume_points=[0._dp,.25_dp*vmax,vmax,bins]
      else
        volume_points=[0._dp,.25_dp*vmax,vmax,vmax+ar*heights]
      endif
      do j=1,size(volume_points)
        do side=-1,1
          v=volume_points(j)+real(side,dp)*ar*dh
          if(v<0._dp) cycle
          call elevation(version,v,h)
          call volume(version,h,v1)
          call close(v1,v,'volume round trip',version,geometry)
        enddo
      enddo
      if(area>ar) then
        v=bins(nz)+2._dp*area
      else
        v=ar*(bank+heights(nz)+2._dp-bed)
      endif
      repeat_volume=v
      do repeat=1,100
        call elevation(version,repeat_volume,h)
        call volume(version,h,repeat_volume)
      enddo
      call close(repeat_volume,v,'repeated no-change round trip',version,geometry)

      ! Negative elevation/volume behavior is intentionally unchanged:
      ! below-bed WSE clips to zero; a negative input volume maps below bed,
      ! and the subsequent forward mapping clips it to zero.
      call volume(version,bed-1._dp,v)
      call close(v,0._dp,'negative depth clips to zero',version,geometry)
      call elevation(version,-ar,h)
      call close(h,bed-1._dp,'negative volume existing behavior',version,geometry)
      call volume(version,h,v)
      call close(v,0._dp,'negative volume round trip clips',version,geometry)

      ! Continuity and monotonicity at bankfull and every stage boundary.
      do j=5,size(stage_points)
        h=stage_points(j)
        call volume(version,h-dh,v1)
        call volume(version,h+dh,v2)
        if(v2<v1.or.v2-v1>2._dp*dh*max(area,ar)*(1._dp+tol)) then
          print *, 'FAIL forward continuity/monotonicity', version,geometry,j,v1,v2
          stop 1
        endif
        checks=checks+1
        call volume(version,h,v)
        call elevation(version,v-ar*dh,h1)
        call elevation(version,v+ar*dh,h2)
        if(h2<h1.or.h2-h1>2._dp*dh*(1._dp+tol)) then
          print *, 'FAIL inverse continuity/monotonicity', version,geometry,j,h1,h2
          stop 1
        endif
        checks=checks+1
      enddo
    enddo
  enddo
  print *, 'PASS actual HYMAP2/3 profile subroutines; checks=',checks

contains
  subroutine volume(version,h,v)
    integer, intent(in) :: version
    real(dp), intent(in) :: h
    real(dp), intent(out) :: v
    if(version==2) then
      call HYMAP2_get_volume_profile(nz,bank,heights,bins,area,vmax,bed,length,width,h,v)
    else
      call HYMAP3_get_volume_profile(nz,bank,heights,bins,area,vmax,bed,length,width,h,v)
    endif
  end subroutine

  subroutine elevation(version,v,h)
    integer, intent(in) :: version
    real(dp), intent(in) :: v
    real(dp), intent(out) :: h
    if(version==2) then
      call HYMAP2_get_elevation_profile(nz,bank,heights,bins,area,vmax,bed,length,width,h,v)
    else
      call HYMAP3_get_elevation_profile(nz,bank,heights,bins,area,vmax,bed,length,width,h,v)
    endif
  end subroutine

  function expected_volume(h) result(v)
    real(dp), intent(in) :: h
    real(dp) :: v,lo_h,lo_v
    integer :: k
    if(h<=bank.or.area<=ar) then
      v=max(h-bed,0._dp)*ar
    elseif(h>bank+heights(nz)) then
      v=bins(nz)+(h-bank-heights(nz))*area
    else
      lo_h=bank
      lo_v=vmax
      do k=1,nz
        if(h<=bank+heights(k)) then
          v=lo_v+(bins(k)-lo_v)*(h-lo_h)/(bank+heights(k)-lo_h)
          return
        endif
        lo_h=bank+heights(k)
        lo_v=bins(k)
      enddo
      stop 2
    endif
  end function

  subroutine close(actual,expected,label,version,geometry)
    real(dp), intent(in) :: actual,expected
    character(*), intent(in) :: label
    integer, intent(in) :: version,geometry
    if(.not.(abs(actual-expected)<=tol*max(1._dp,abs(expected)))) then
      print *, 'FAIL ',label,' H',version,' geometry',geometry
      print *, 'actual=',actual,' expected=',expected
      stop 1
    endif
    checks=checks+1
  end subroutine

  subroutine check_elevation(h,version,geometry)
    real(dp), intent(in) :: h
    integer, intent(in) :: version,geometry
    real(dp) :: v,round_h
    call volume(version,h,v)
    call close(v,expected_volume(h),'segment forward',version,geometry)
    call elevation(version,v,round_h)
    call close(round_h,h,'segment round trip',version,geometry)
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
                        help="Check extraction only; this does not validate Fortran behavior")
    args = parser.parse_args()
    root = locate_root(args.source_root)
    routines = []
    for version in (2, 3):
        path = root / f"lis/routing/HYMAP{version}_router/HYMAP{version}_modelMod.F90"
        source = path.read_text()
        for operation in ("volume", "elevation"):
            routines.append(extract_subroutine(source, f"HYMAP{version}_get_{operation}_profile"))
    if args.check_extraction:
        print("Extracted four actual Fortran subroutines; Fortran execution NOT performed")
        return 0
    compiler = compiler_command(args.compiler)
    if compiler is None:
        print("SKIP: GNU Fortran unavailable; Fortran regression NOT executed", file=sys.stderr)
        return 77
    wrapper = ("module LIS_logMod\n  implicit none\n  integer :: LIS_logunit=6\nend module\n"
               "module profile_test_module\n  implicit none\ncontains\n" + "\n".join(routines) +
               "\nend module\nsubroutine LIS_endrun()\n  stop 99\nend subroutine\n")
    with tempfile.TemporaryDirectory(prefix="hymap-profile-regression-") as directory:
        work = Path(directory)
        (work / "profiles.F90").write_text(wrapper)
        (work / "driver.F90").write_text(DRIVER)
        command = compiler + ["-O0", "-g", "-fcheck=all", "-fbacktrace", "-Wall", "-Wextra",
                              "-ffree-line-length-none", "profiles.F90", "driver.F90", "-o", "regression"]
        subprocess.run(command, cwd=work, check=True)
        subprocess.run([str(work / "regression")], cwd=work, check=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
