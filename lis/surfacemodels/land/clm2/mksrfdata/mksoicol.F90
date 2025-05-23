!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "preproc.h"
#include "LIS_misc.h"

subroutine mksoicol (fsoicol, ndiag, pctgla_o, soil_color_o)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Make soil color classes for LSM grid from BATS T42 data
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: mksoicol.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------

  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_precisionMod
  use clm2_varpar        !lsm parameters
  use clm2_varsur        !lsm surface variables
  use fileutils, only : getfil
  use clm2_areaMod          !area averaging routines 
  use clm2_shr_sys_mod, only : shr_sys_flush 
  implicit none

! ------------------------ arguments ------------------------------
  character(len=*), intent(in) :: fsoicol              !input soicol dataset file name
  integer , intent(in) :: ndiag                        !unit number for diagnostic output
  real(r8), intent(in) :: pctgla_o(lsmlon,lsmlat)      !output grid: percent glacier
  integer , intent(out):: soil_color_o(lsmlon,lsmlat)  !LSM soil color classes on output grid
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  character(len=LIS_CONST_PATH_LEN) locfn     !local dataset file name

  integer :: nlon_i                           !input grid : longitude points (read in)
  integer :: nlat_i                           !input grid : latitude  points (read in)
  integer :: ncid,dimid,varid                 !input netCDF id's
  integer :: ier                              !error status   

  integer :: ii                               !longitude index for BATS grid
  integer :: io                               !longitude index for LSM grid
  integer :: ji                               !latitude  index for BATS grid
  integer :: jo                               !latitude  index for LSM grid
  integer :: k                                !temporary BATS or LSM soil color
  integer :: miss = 99999                     !missing data indicator
  integer :: n                                !loop index
  integer, parameter :: num=2                 !get 1st and 2nd largest areas of overlap
  integer :: wsti(num)                        !index to 1st and 2nd largest values in wst

  integer, parameter :: nsoicol=8             !number of LSM color classes
  real(r8) :: wst(0:nsoicol)                  !overlap weights, by surface type
  real(r8) :: gast_i(0:nsoicol)               !input grid : global area, by surface type
  real(r8) :: gast_o(0:nsoicol)               !output grid: global area, by surface type
  character(len=35) :: col(0:nsoicol)         !name of each color

  real(r8) ::edge_i(4)                        !input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: latixy_i(:,:)      !input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)      !input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)        !input grid: number longitude points by lat
  real(r8), allocatable :: lon_i(:,:)         !input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lon_i_offset(:,:)  !input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)           !input grid: latitude, south edge (degrees)
  real(r8), allocatable :: area_i(:,:)        !input grid: cell area
  real(r8), allocatable :: mask_i(:,:)        !input grid: mask (0, 1)
  integer , allocatable :: soil_color_i(:,:)  !input grid: BATS soil color

  real(r8) :: mask_o                          !output grid: mask (0, 1)
  integer  :: novr_i2o                        !number of overlapping input cells
  integer  :: iovr_i2o(maxovr)                !lon index of overlap input cell
  integer  :: jovr_i2o(maxovr)                !lat index of overlap input cell
  real(r8) :: wovr_i2o(maxovr)                !weight    of overlap input cell
  real(r8) :: offset                          !used to shift x-grid 360 degrees

  real(r8) :: fld_o(lsmlon,lsmlat)            !output grid: dummy field 
  real(r8) :: fld_i                           !input grid: dummy field 
  real(r8) :: sum_fldo                        !global sum of dummy output field
  real(r8) :: sum_fldi                        !global sum of dummy input field
  real(r8) :: relerr = 0.00001                !max error: sum overlap weights ne 1
! -----------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call clm2_shr_sys_flush(6)

! -----------------------------------------------------------------
! Define the LSM color classes: 0 to nsoicol
! -----------------------------------------------------------------

  col(0) = 'no soil                            '
  col(1) = 'class 1: light                     '
  col(2) = 'class 2:                           '
  col(3) = 'class 3:                           '
  col(4) = 'class 4:                           '
  col(5) = 'class 5:                           '
  col(6) = 'class 6:                           '
  col(7) = 'class 7:                           '
  col(8) = 'class 8: dark                      '
! col(9) = 'class 9: very light North Africa   '

! -----------------------------------------------------------------
! Read input soil colors
! -----------------------------------------------------------------

! Obtain input grid info

  call getfil (fsoicol, locfn, 0)
  call wrap_open(locfn, 0, ncid)

  call wrap_inq_dimid  (ncid, 'lon', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlon_i)

  call wrap_inq_dimid  (ncid, 'lat', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlat_i)

  allocate (latixy_i(nlon_i,nlat_i), stat=ier)   
  if (ier/=0) call endrun
  allocate (longxy_i(nlon_i,nlat_i), stat=ier)   
  if (ier/=0) call endrun
  allocate (numlon_i(nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (lon_i(nlon_i+1,nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (lon_i_offset(nlon_i+1,nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (lat_i(nlat_i+1), stat=ier)        
  if (ier/=0) call endrun
  allocate (area_i(nlon_i,nlat_i), stat=ier)  
  if (ier/=0) call endrun
  allocate (mask_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun
  allocate (soil_color_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun

  call wrap_inq_varid (ncid, 'LATIXY', varid)
  call wrap_get_var_realx (ncid, varid, latixy_i)

  call wrap_inq_varid (ncid, 'LONGXY', varid)
  call wrap_get_var_realx (ncid, varid, longxy_i)

  call wrap_inq_varid (ncid, 'EDGEN', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(1))

  call wrap_inq_varid (ncid, 'EDGEE', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(2))

  call wrap_inq_varid (ncid, 'EDGES', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(3))

  call wrap_inq_varid (ncid, 'EDGEW', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(4))

! Obtain input data

  call wrap_inq_varid (ncid, 'SOIL_COLOR', varid)
  call wrap_get_var_int (ncid, varid, soil_color_i)

  call wrap_close(ncid)

! -----------------------------------------------------------------
! Map data from input grid to land model grid. Get:
! -----------------------------------------------------------------

! Determine input grid cell and cell areas

  numlon_i(:) = nlon_i

  call celledge (nlat_i    , nlon_i    , numlon_i  , longxy_i  ,  &
                 latixy_i  , edge_i(1) , edge_i(2) , edge_i(3) ,  &
                 edge_i(4) , lat_i     , lon_i     )

  call cellarea (nlat_i    , nlon_i    , numlon_i  , lat_i     ,  &
                 lon_i     , edge_i(1) , edge_i(2) , edge_i(3) ,  &
                 edge_i(4) , area_i    )

  do ji = 1, nlat_i
     do ii = 1, numlon_i(ji)
        mask_i(ii,ji) = 1.
     end do
  end do

! Shift x-grid to locate periodic grid intersections. This
! assumes that all lon_i(1,j) have the same value for all
! latitudes j and that the same holds for lon_o(1,j)

  if (lon_i(1,1) < lonw(1,1)) then
     offset = 360.0
  else
     offset = -360.0
  end if
  
  do ji = 1, nlat_i
     do ii = 1, numlon_i(ji) + 1
        lon_i_offset(ii,ji) = lon_i(ii,ji) + offset
     end do
  end do

! Process each cell on land model grid
! novr_i2o - number of input grid cells that overlap each land grid cell
! iovr_i2o - longitude index of overlapping input grid cell
! jovr_i2o - latitude  index of overlapping input grid cell
! wovr_i2o - fraction of land grid cell overlapped by input grid cell

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,k,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i, &
!$OMP & wst, wsti)
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i   , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats     , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o , &
                           wovr_i2o)                             

! Sum overlap weights by color class

        do k = 0, nsoicol
           wst(k) = 0.
        end do
        do n = 1, novr_i2o         !overlap cell index
           ii = iovr_i2o(n)        !lon index (input grid) of overlap cell
           ji = jovr_i2o(n)        !lat index (input grid) of overlap cell
           k = soil_color_i(ii,ji) !color class (input grid)
           wst(k) = wst(k) + wovr_i2o(n)
        end do

! Rank non-zero weights by color type. wsti(1) is the most extensive
! color type. wsti(2) is the second most extensive color type

        call mkrank (nsoicol, wst, miss, wsti, num)
        soil_color_o(io,jo) = wsti(1)

! If land but no color, set color to 4

        if (landmask(io,jo)==1 .and. soil_color_o(io,jo)==0) soil_color_o(io,jo) = 4

! Set ocean colors to zero

        if (landmask(io,jo) == 0) soil_color_o(io,jo) = 0

! Set color for grid cells that are 100% glacier to zero. Otherwise,
! must have a soil color for the non-glacier portion of grid cell.

        if (abs(pctgla_o(io,jo)-100.)<1.e-06) soil_color_o(io,jo)=0

! Error checks

        if (soil_color_o(io,jo) < 0 .or. soil_color_o(io,jo) > nsoicol) then
           write (6,*) 'MKSOICOL error: land model soil color = ', &
                soil_color_o(io,jo),' is not valid for lon,lat = ',io,jo
           call endrun
        end if

! Global sum of output field -- must multiply by fraction of
! output grid that is land as determined by input grid

        fld_o(io,jo) = 0.
        do n = 1, novr_i2o
           ii = iovr_i2o(n)
           ji = jovr_i2o(n)
           fld_i = ((ji-1)*nlon_i + ii) 
           fld_o(io,jo) = fld_o(io,jo) + wovr_i2o(n) * fld_i 
        end do

     end do  !end of output longitude loop
  end do     !end of output latitude  loop
!$OMP END PARALLEL DO

! -----------------------------------------------------------------
! Error check1
! Compare global sum fld_o to global sum fld_i. 
! -----------------------------------------------------------------

! This check is true only if both grids span the same domain. 
! To obtain global sum of input field must multiply by 
! fraction of input grid that is land as determined by input grid

  sum_fldo = 0.
  do jo = 1,lsmlat
     do io = 1,numlon(jo)
        sum_fldo = sum_fldo + area(io,jo) * fld_o(io,jo) 
     end do
  end do

  sum_fldi = 0.
  do ji = 1, nlat_i      
     do ii = 1, numlon_i(ji)
        fld_i = ((ji-1)*nlon_i + ii) 
        sum_fldi = sum_fldi + area_i(ii,ji) * fld_i
     end do
  end do

  if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
     write (6,*) 'MKGLACIER error: input field not conserved'
     write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
     write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
     call endrun
  end if

! -----------------------------------------------------------------
! Error check2
! Compare global area of each soil color on input and output grids
! -----------------------------------------------------------------

! input grid

  gast_i(:) = 0.
  do ji = 1, nlat_i
     do ii = 1, nlon_i
        k = soil_color_i(ii,ji)
        gast_i(k) = gast_i(k) + area_i(ii,ji)
     end do
  end do

! output grid

  gast_o(:) = 0.
  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        k = soil_color_o(io,jo)
        gast_o(k) = gast_o(k) + area(io,jo)
     end do
  end do

! area comparison

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Soil Color Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001)
1001 format (1x,'soil color type',20x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)

  do k = 0, nsoicol
     write (ndiag,1002) col(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002 format (1x,a35,f16.3,f17.3)
  end do

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
     
  write (6,*) 'Successfully made soil color classes'
  write (6,*)
  call clm2_shr_sys_flush(6)


! Deallocate dynamic memory

  deallocate (latixy_i)
  deallocate (longxy_i)
  deallocate (numlon_i)
  deallocate (lon_i)
  deallocate (lon_i_offset)
  deallocate (lat_i)
  deallocate (area_i)
  deallocate (mask_i)
  deallocate (soil_color_i)

  return
end subroutine mksoicol


