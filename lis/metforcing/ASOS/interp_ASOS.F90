!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: interp_ASOS.F90
! \label{interp_ASOS}
!
! !REVISION HISTORY:
!  3 Jul 2013: Soni Yatheendradas: Modified to add findex argument towards
!              spatial interpolation choice
!  6 Sep 2016: Minwook Kim: based on new TRMM 3B42V7 code that has changes
!
! !INTERFACE: 

subroutine interp_ASOS(n, nASOS, f, lb, lis_gds, nc, nr, varfield, &
                             findex) 
! !USES:
  use ASOS_forcingMod, only : ASOS_struc
  use LIS_coreMod,           only : LIS_rc , LIS_domain
  use LIS_logMod,            only : LIS_logunit, LIS_endrun 

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer                :: nc
  integer                :: nr
  integer                :: nASOS
  real                   :: lis_gds(50)
  real                   :: f(nASOS)
  logical*1              :: lb(nASOS)
  real, dimension(nc,nr) :: varfield
  integer, intent(in) :: findex ! SY

!
! !DESCRIPTION:
!   This subroutine interpolates a given ASOS field
!   to the LIS grid.
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nASOS]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description}
!

!
!  The routines invoked are:
!  \begin{description}
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP


  integer :: iret
  integer :: mo
  integer :: count,i,j

  real, dimension(nc*nr) :: lis1d

  logical*1 :: geogmask(nc,nr)
  logical*1 :: lo(nc*nr)

  !=== End variable declarations
  !--------------------------------------------------------------------
  ! Setting interpolation options (ip=0,bilinear)
  ! (km=1, one parameter, ibi=1,use undefined bitmap
  ! (needed for soil moisture and temperature only)
  ! Use budget bilinear (ip=3) for precip forcing fields
  !--------------------------------------------------------------------
  mo = nc*nr

  !--------------------------------------------------------------------  
  ! Initialize output bitmap. Important for soil moisture and temp.
  !--------------------------------------------------------------------  
  lb = .true.
  lo = .true.


  ASOS_struc(n)%mi = nASOS

  if(LIS_rc%met_interp(findex).eq."budget-bilinear") then ! SY
   call conserv_interp(lis_gds,lb,f,lo,lis1d,ASOS_struc(n)%mi,mo,&
        LIS_domain(n)%lat, LIS_domain(n)%lon, &
        ASOS_struc(n)%w112,&
        ASOS_struc(n)%w122,ASOS_struc(n)%w212,&
        ASOS_struc(n)%w222,&
        ASOS_struc(n)%n112,ASOS_struc(n)%n122,&
        ASOS_struc(n)%n212,&
        ASOS_struc(n)%n222,-9999.0,iret)
  elseif(LIS_rc%met_interp(findex).eq."neighbor") then ! SY
  ! SY: begin for nearest neighbor
   call neighbor_interp(lis_gds,lb,f,lo,lis1d,ASOS_struc(n)%mi,mo,&
        LIS_domain(n)%lat, LIS_domain(n)%lon, &
        ASOS_struc(n)%n112, LIS_rc%udef,iret)
  ! SY: end for nearest neighbor
  else
   write(LIS_logunit,*) 'This interpolation not defined for ASOS data'
   write(LIS_logunit,*) 'Program stopping ... '
   call LIS_endrun()
  endif

  !--------------------------------------------------------------------    
  ! Create 2D array for main program. Also define a "soil" mask
  ! due to different geography between GDAS & LDAS. For LDAS land 
  ! points not included in GDAS geography dataset only.
  !--------------------------------------------------------------------    
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count)
        geogmask(i,j) = lo(i+count)
     enddo
     count = count + nc
  enddo
end subroutine interp_ASOS

